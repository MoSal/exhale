/* specAnalysis.h - header file for class providing spectral analysis of MCLT signals
 * written by C. R. Helmrich, last modified in 2020 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under the exhale Copyright License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2020 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#ifndef _SPEC_ANALYSIS_H_
#define _SPEC_ANALYSIS_H_

#include "exhaleLibPch.h"
#include "linearPrediction.h"

// constants, experimental macros
#define SA_BW_SHIFT             5
#define SA_BW  (1 << SA_BW_SHIFT)
#define SA_EPS               1024
#define SA_EXACT_COMPLEX_ABS    0
#define SA_IMPROVED_SFM_ESTIM   1
#define SA_OPT_WINDOW_GROUPING  1

// spectral signal analysis class
class SpecAnalyzer
{
private:

  // member variables
  uint16_t m_bandwidthOff[USAC_MAX_NUM_CHANNELS];
#if SA_IMPROVED_SFM_ESTIM
  uint8_t  m_magnCorrPrev[USAC_MAX_NUM_CHANNELS];
  uint32_t* m_magnSpectra[USAC_MAX_NUM_CHANNELS];
#endif
  uint32_t m_meanAbsValue[USAC_MAX_NUM_CHANNELS][1024 >> SA_BW_SHIFT];
  uint16_t m_numAnaBands [USAC_MAX_NUM_CHANNELS];
  short    m_parCorCoeffs[USAC_MAX_NUM_CHANNELS][MAX_PREDICTION_ORDER];
  uint32_t m_specAnaStats[USAC_MAX_NUM_CHANNELS];
  uint32_t m_tnsPredGains[USAC_MAX_NUM_CHANNELS];
  LinearPredictor* m_tnsPredictor;

public:

  // constructor
  SpecAnalyzer ();
  // destructor
#if SA_IMPROVED_SFM_ESTIM
  ~SpecAnalyzer () { for (unsigned ch = 0; ch < USAC_MAX_NUM_CHANNELS; ch++) MFREE (m_magnSpectra[ch]); }
#else
  ~SpecAnalyzer () { }
#endif
  // public functions
  unsigned getLinPredCoeffs (short parCorCoeffs[MAX_PREDICTION_ORDER], const unsigned channelIndex); // returns best filter order
  unsigned getMeanAbsValues (const int32_t* const mdctSignal, const int32_t* const mdstSignal, const unsigned nSamplesInFrame,
                             const unsigned channelIndex, const uint16_t* const bandStartOffsets, const unsigned nBands,
                             uint32_t* const meanBandValues);
  void getSpecAnalysisStats (uint32_t avgSpecAnaStats[USAC_MAX_NUM_CHANNELS], const unsigned nChannels);
  void getSpectralBandwidth (uint16_t bandwidthOffset[USAC_MAX_NUM_CHANNELS], const unsigned nChannels);
  unsigned initSigAnaMemory (LinearPredictor* const linPredictor, const unsigned nChannels, const unsigned maxTransfLength);
#if SA_OPT_WINDOW_GROUPING
  unsigned optimizeGrouping (const unsigned channelIndex, const unsigned preferredBandwidth, const unsigned preferredGrouping);
#endif
  unsigned spectralAnalysis (const int32_t* const mdctSignals[USAC_MAX_NUM_CHANNELS],
                             const int32_t* const mdstSignals[USAC_MAX_NUM_CHANNELS],
                             const unsigned nChannels, const unsigned nSamplesInFrame, const unsigned samplingRate,
                             const unsigned lfeChannelIndex = USAC_MAX_NUM_CHANNELS); // to skip an LFE channel
  int16_t stereoSigAnalysis (const int32_t* const mdctSignal1, const int32_t* const mdctSignal2,
                             const int32_t* const mdstSignal1, const int32_t* const mdstSignal2,
                             const unsigned nSamplesMax, const unsigned nSamplesInFrame, const bool shortTransforms,
                             uint8_t* const stereoCorrValue = nullptr); // per-band perceptual correlation data
}; // SpecAnalyzer

#endif // _SPEC_ANALYSIS_H_
