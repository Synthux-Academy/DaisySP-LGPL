/*
Copyright (c) 2023 Electrosmith, Corp, Sean Costello, Istvan Varga, Paul Batchelor

Use of this source code is governed by the LGPL V2.1
license that can be found in the LICENSE file or at
https://opensource.org/license/lgpl-2-1/
*/

// ReverbScLoFi: a lighter variant of daisysp::ReverbSc (see reverbsc.h).
// Differences from ReverbSc:
//   - 6 internal delay lines instead of 8 (lower CPU)
//   - linear delay-line interpolation instead of 4-point cubic (lower CPU)
//   - delay storage sized for those 6 lines at 48kHz (20544 floats, ~80KB)
//     instead of the stock 98936 floats. Sample rates above ~48kHz don't
//     fit and Init() returns 1.
//   - delay storage is a buffer the caller provides to Init() (so it can be
//     placed in a chosen RAM region, e.g. DTCM/SRAM or SDRAM) instead of an
//     array inside the object; the object itself is small. The buffer must
//     hold kReverbScLoFiBufSamples floats and outlive the reverb.
//   - junction scaling 2/kReverbScLoFiLines (upstream hard-codes 0.25 = 2/8;
//     kept at 0.25 with 6 lines the network loses energy every pass and
//     the tail is much shorter than the feedback setting implies)
//   - delay lines are packed by sample count (upstream's offset math mixed
//     bytes and floats, oversizing its buffer ~4x)
//   - damping coefficient uses a cosine lookup table instead of cosf()
// Everything else (per-line random modulation, damping filter, parameter
// ranges, Process() API) is unchanged, so it can be swapped in for ReverbSc
// with only the Init() call changing.
//
// Usage:
//   static float rev_buf[daisysp::kReverbScLoFiBufSamples];
//   daisysp::ReverbScLoFi rev;
//   rev.Init(sample_rate, rev_buf);

#pragma once
#ifndef DSYSP_REVERBSCLOFI_H
#define DSYSP_REVERBSCLOFI_H

namespace daisysp
{
/** Number of parallel comb/allpass delay lines. */
static constexpr int kReverbScLoFiLines = 6;

/** Delay storage in floats. The 6 delay lines below are packed tightly and
 *  sum to 20494 samples at 48kHz; this leaves a few samples for float
 *  rounding in DelayLineMaxSamples(). Init() returns 1 if they don't fit. */
#define DSY_REVERBSCLOFI_MAX_SIZE 20544
static constexpr int kReverbScLoFiBufSamples = DSY_REVERBSCLOFI_MAX_SIZE;

/** Delay line for internal reverb use */
typedef struct
{
    int    write_pos;         /**< write position */
    int    buffer_size;       /**< buffer size */
    int    read_pos;          /**< read position */
    int    read_pos_frac;     /**< fractional component of read pos */
    int    read_pos_frac_inc; /**< increment for fractional */
    int    dummy;             /**<  dummy var */
    int    seed_val;          /**< randseed */
    int    rand_line_cnt;     /**< number of random lines */
    float  filter_state;      /**< state of filter */
    float  *buf;              /**< buffer ptr */
} ReverbScLoFiDl;

/** Stereo Reverb */
class ReverbScLoFi
{
  public:
    ReverbScLoFi() {}
    ~ReverbScLoFi() {}
    /** Initializes the reverb module, and sets the sample_rate at which the Process function will be called.
        buf holds kReverbScLoFiBufSamples samples of delay storage; Init() clears it.
        Returns 0 if all good, or 1 if it runs out of delay times exceed maximum allowed.
    */
    int Init(float sample_rate, float *buf);

    /** Process the input through the reverb, and updates values of out1, and out2 with the new processed signal.
    */
    int Process(const float &in1, const float &in2, float *out1, float *out2);

    /** controls the reverb time. reverb tail becomes infinite when set to 1.0
        \param fb - sets reverb time. range: 0.0 to 1.0
    */
    inline void SetFeedback(const float &fb) { feedback_ = fb; }
    /** controls the internal dampening filter's cutoff frequency.
        \param freq - low pass frequency. range: 0.0 to sample_rate / 2
    */
    inline void SetLpFreq(const float &freq) { lpfreq_ = freq; }

  private:
    void       NextRandomLineseg(ReverbScLoFiDl *lp, int n);
    int        InitDelayLine(ReverbScLoFiDl *lp, int n);
    float      feedback_, lpfreq_;
    float      i_sample_rate_, i_pitch_mod_, i_skip_init_;
    float      sample_rate_;
    float      damp_fact_;
    float      prv_lpfreq_;
    int        init_done_;
    ReverbScLoFiDl delay_lines_[kReverbScLoFiLines];
    float     *aux_;
};


} // namespace daisysp
#endif
