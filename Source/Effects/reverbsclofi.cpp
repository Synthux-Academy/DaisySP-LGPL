/*
Copyright (c) 2023 Electrosmith, Corp, Sean Costello, Istvan Varga, Paul Batchelor

Use of this source code is governed by the LGPL V2.1
license that can be found in the LICENSE file or at
https://opensource.org/license/lgpl-2-1/
*/

// Lower-CPU variant of daisysp::ReverbSc - see reverbsclofi.h for the list of changes.

#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "reverbsclofi.h"

#define REVSC_OK 0
#define REVSC_NOT_OK 1

#define DEFAULT_SRATE 48000.0
#define MIN_SRATE 5000.0
#define MAX_SRATE 1000000.0
#define MAX_PITCHMOD 20.0
#define DELAYPOS_SHIFT 28
#define DELAYPOS_SCALE 0x10000000
#define DELAYPOS_MASK 0x0FFFFFFF

#ifndef M_PI
#define M_PI 3.14159265358979323846 /* pi */
#endif

using namespace daisysp;

/* kReverbParams[n][0] = delay time (in seconds)                     */
/* kReverbParams[n][1] = random variation in delay time (in seconds) */
/* kReverbParams[n][2] = random variation frequency (in 1/sec)       */
/* kReverbParams[n][3] = random seed (0 - 32767)                     */

/* First kReverbScLoFiLines rows of the stock 8-line parameter table. */
static const float kReverbParams[kReverbScLoFiLines][4]
    = {{(2473.0 / DEFAULT_SRATE), 0.0010, 3.100, 1966.0},
       {(2767.0 / DEFAULT_SRATE), 0.0011, 3.500, 29491.0},
       {(3217.0 / DEFAULT_SRATE), 0.0017, 1.110, 22937.0},
       {(3557.0 / DEFAULT_SRATE), 0.0006, 3.973, 9830.0},
       {(3907.0 / DEFAULT_SRATE), 0.0010, 2.341, 20643.0},
       {(4127.0 / DEFAULT_SRATE), 0.0011, 1.897, 22937.0}};

static int DelayLineMaxSamples(float sr, float i_pitch_mod, int n);
static const float kOutputGain = 0.35;
// Junction scattering 2/N keeps the line network lossless, so the decay is
// set by feedback_ and the damping filter alone (stock: 0.25 for 8 lines).
static const float kJpScale    = 2.0f / (float)kReverbScLoFiLines;

/* Lookup table replacing cosf() in the damping filter coefficient calc
   below; that calc only runs when lpfreq_ changes, but avoiding cosf keeps
   the cost bounded and predictable. */
#define REVSC_COS_TABLE_SIZE 1024
static float kCosTable[REVSC_COS_TABLE_SIZE];
static bool  kCosTableInit = false;

static void InitCosTable()
{
    if(kCosTableInit)
        return;
    for(int i = 0; i < REVSC_COS_TABLE_SIZE; i++)
    {
        kCosTable[i]
            = cosf((2.0f * (float)M_PI * (float)i) / (float)REVSC_COS_TABLE_SIZE);
    }
    kCosTableInit = true;
}

static float LookupCos(float x)
{
    const float two_pi = 2.0f * (float)M_PI;
    x                  = fmodf(x, two_pi);
    if(x < 0.0f)
        x += two_pi;

    float findex = x * ((float)REVSC_COS_TABLE_SIZE / two_pi);
    int   idx0   = (int)findex;
    float frac   = findex - (float)idx0;
    int   idx1   = (idx0 + 1) & (REVSC_COS_TABLE_SIZE - 1);
    idx0 &= (REVSC_COS_TABLE_SIZE - 1);

    return kCosTable[idx0] + frac * (kCosTable[idx1] - kCosTable[idx0]);
}

int ReverbScLoFi::Init(float sr, float *buf)
{
    InitCosTable();
    aux_           = buf;
    i_sample_rate_ = sr;
    sample_rate_   = sr;
    feedback_      = 0.97;
    lpfreq_        = 10000;
    i_pitch_mod_   = 1;
    i_skip_init_   = 0;
    damp_fact_     = 1.0;
    prv_lpfreq_    = 0.0;
    init_done_     = 1;
    // NB: upstream advanced `aux_` by a *byte* count while treating it as a
    // float index, so its aux_ was ~4x oversized. Here the lines are packed
    // tightly by sample count and DSY_REVERBSCLOFI_MAX_SIZE is the real float
    // capacity of aux_, in samples.
    int n_samps = 0;
    for(int i = 0; i < kReverbScLoFiLines; i++)
    {
        const int len = DelayLineMaxSamples(sr, 1, i);
        if(n_samps + len > DSY_REVERBSCLOFI_MAX_SIZE)
            return 1;
        delay_lines_[i].buf = aux_ + n_samps;
        InitDelayLine(&delay_lines_[i], i);
        n_samps += len;
    }
    return 0;
}

static int DelayLineMaxSamples(float sr, float i_pitch_mod, int n)
{
    float max_del;

    max_del = kReverbParams[n][0];
    max_del += (kReverbParams[n][1] * (float)i_pitch_mod * 1.125);
    return (int)(max_del * sr + 16.5);
}

void ReverbScLoFi::NextRandomLineseg(ReverbScLoFiDl *lp, int n)
{
    float prv_del, nxt_del, phs_inc_val;

    /* update random seed */
    if(lp->seed_val < 0)
        lp->seed_val += 0x10000;
    lp->seed_val = (lp->seed_val * 15625 + 1) & 0xFFFF;
    if(lp->seed_val >= 0x8000)
        lp->seed_val -= 0x10000;
    /* length of next segment in samples */
    lp->rand_line_cnt = (int)((sample_rate_ / kReverbParams[n][2]) + 0.5);
    prv_del           = (float)lp->write_pos;
    prv_del -= ((float)lp->read_pos
                + ((float)lp->read_pos_frac / (float)DELAYPOS_SCALE));
    while(prv_del < 0.0)
        prv_del += lp->buffer_size;
    prv_del = prv_del / sample_rate_; /* previous delay time in seconds */
    nxt_del = (float)lp->seed_val * kReverbParams[n][1] / 32768.0;
    /* next delay time in seconds */
    nxt_del = kReverbParams[n][0] + (nxt_del * (float)i_pitch_mod_);
    /* calculate phase increment per sample */
    phs_inc_val           = (prv_del - nxt_del) / (float)lp->rand_line_cnt;
    phs_inc_val           = phs_inc_val * sample_rate_ + 1.0;
    lp->read_pos_frac_inc = (int)(phs_inc_val * DELAYPOS_SCALE + 0.5);
}

int ReverbScLoFi::InitDelayLine(ReverbScLoFiDl *lp, int n)
{
    float read_pos;
    /* int     i; */

    /* calculate length of delay line */
    lp->buffer_size = DelayLineMaxSamples(sample_rate_, 1, n);
    lp->dummy       = 0;
    lp->write_pos   = 0;
    /* set random seed */
    lp->seed_val = (int)(kReverbParams[n][3] + 0.5);
    /* set initial delay time */
    read_pos     = (float)lp->seed_val * kReverbParams[n][1] / 32768;
    read_pos     = kReverbParams[n][0] + (read_pos * (float)i_pitch_mod_);
    read_pos     = (float)lp->buffer_size - (read_pos * sample_rate_);
    lp->read_pos = (int)read_pos;
    read_pos     = (read_pos - (float)lp->read_pos) * (float)DELAYPOS_SCALE;
    lp->read_pos_frac = (int)(read_pos + 0.5);
    /* initialise first random line segment */
    NextRandomLineseg(lp, n);
    /* clear delay line to zero */
    lp->filter_state = 0.0;
    for(int i = 0; i < lp->buffer_size; i++)
    {
        lp->buf[i] = 0;
    }
    return REVSC_OK;
}

int ReverbScLoFi::Process(const float &in1,
                      const float &in2,
                      float *      out1,
                      float *      out2)
{
    float       a_in_l, a_in_r, a_out_l, a_out_r;
    float       v0, v1, frac;
    ReverbScLoFiDl *lp;
    int         read_pos;
    uint32_t    n;
    int         buffer_size; /* Local copy */
    float       damp_fact = damp_fact_;

    if(init_done_ <= 0)
        return REVSC_NOT_OK;

    /* calculate tone filter coefficient if frequency changed */
    if(lpfreq_ != prv_lpfreq_)
    {
        prv_lpfreq_ = lpfreq_;
        damp_fact
            = 2.0f - LookupCos(prv_lpfreq_ * (2.0f * (float)M_PI) / sample_rate_);
        damp_fact = damp_fact_
            = damp_fact - sqrtf(damp_fact * damp_fact - 1.0f);
    }

    /* calculate "resultant junction pressure" and mix to input signals */

    a_in_l = a_out_l = a_out_r = 0.0;
    for(n = 0; n < (uint32_t)kReverbScLoFiLines; n++)
    {
        a_in_l += delay_lines_[n].filter_state;
    }
    a_in_l *= kJpScale;
    a_in_r = a_in_l + in2;
    a_in_l = a_in_l + in1;

    /* loop through all delay lines */

    for(n = 0; n < (uint32_t)kReverbScLoFiLines; n++)
    {
        lp          = &delay_lines_[n];
        buffer_size = lp->buffer_size;

        /* send input signal and feedback to delay line */

        lp->buf[lp->write_pos]
            = (float)((n & 1 ? a_in_r : a_in_l) - lp->filter_state);
        if(++lp->write_pos >= buffer_size)
        {
            lp->write_pos -= buffer_size;
        }

        /* read from delay line with linear interpolation */

        if(lp->read_pos_frac >= DELAYPOS_SCALE)
        {
            lp->read_pos += (lp->read_pos_frac >> DELAYPOS_SHIFT);
            lp->read_pos_frac &= DELAYPOS_MASK;
        }
        if(lp->read_pos >= buffer_size)
            lp->read_pos -= buffer_size;
        read_pos = lp->read_pos;
        frac     = (float)lp->read_pos_frac * (1.0f / (float)DELAYPOS_SCALE);

        /* read two adjacent samples and interpolate */

        v0 = lp->buf[read_pos];
        if(++read_pos >= buffer_size)
            read_pos -= buffer_size;
        v1 = lp->buf[read_pos];
        v0 = v0 + frac * (v1 - v0);

        /* update buffer read position */

        lp->read_pos_frac += lp->read_pos_frac_inc;

        /* apply feedback gain and lowpass filter */

        v0 *= (float)feedback_;
        v0               = (lp->filter_state - v0) * damp_fact + v0;
        lp->filter_state = v0;

        /* mix to output */

        if(n & 1)
        {
            a_out_r += v0;
        }
        else
        {
            a_out_l += v0;
        }

        /* start next random line segment if current one has reached endpoint */

        if(--(lp->rand_line_cnt) <= 0)
        {
            NextRandomLineseg(lp, n);
        }
    }
    /* someday, use a_out_r for multimono out */

    *out1 = a_out_l * kOutputGain;
    *out2 = a_out_r * kOutputGain;
    return REVSC_OK;
}
