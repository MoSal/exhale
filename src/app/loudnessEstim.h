/* loudnessEstim.h - header file for class with ITU-R BS.1770-4 loudness level estimation
 * written by C. R. Helmrich, last modified in 2020 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under a Modified BSD-Style License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2020 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#ifndef _LOUDNESS_ESTIM_H_
#define _LOUDNESS_ESTIM_H_

#include "exhaleAppPch.h"

// constants, experimental macros
#define LE_THRESH_ABS   (15.0f / 268435456.0f) // absolute threshold for -70 LUFS
#define LE_THRESH_REL       0.1f // second stage, relative threshold 10dB below L
#define LE_LUFS_OFFSET  2.53125f // to get -3.01 LUFS for mono 997-Hz 0-dBFS sine

// ITU-R loudness estimator class
class LoudnessEstimator
{
private:

  // member variables
  int32_t  m_filterMemoryI[8]; // channel-wise preceding K-weighting filter input
  int32_t  m_filterMemoryO[8]; // channel-wise previous K-weighting filter output
  uint64_t m_powerValue[4][8]; // channel-wise power in each gating block quarter
  float    m_gbNormFactor; // 64-sample normalization factor, 1/(4*m_gbHopSize64)
  uint8_t  m_filterFactor; // sampling rate dependent K-weighting filter constant
  uint8_t  m_gbHopLength64;  // number of 64-sample units in gating block quarter
  uint8_t  m_gbHopSize64;  // hop-size between gating blocks, 25% of block length
  uint8_t  m_inputChannels;
  uint32_t m_inputMaxValue;
  uint32_t m_inputPeakValue;
  int32_t* m_inputPcmData;
  std::vector <uint32_t> m_gbRmsValues; // sqrt of power average per gating block

public:

  // constructor
  LoudnessEstimator (int32_t* const inputPcmData,       const unsigned bitDepth = 24,
                     const unsigned sampleRate = 44100, const unsigned numChannels = 2);
  // destructor
  ~LoudnessEstimator () { reset (); }
  // public functions
  uint32_t addNewPcmData (const unsigned samplesPerChannel);
  uint32_t getStatistics (const bool includeWarmUp = false);
  void     reset () { m_gbHopLength64 = m_inputPeakValue = 0; m_gbRmsValues.clear (); memset (m_powerValue, 0, sizeof (m_powerValue)); }

}; // LoudnessEstimator

#endif // _LOUDNESS_ESTIM_H_