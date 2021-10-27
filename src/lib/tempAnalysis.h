/* tempAnalysis.h - header file for class providing temporal analysis of PCM signals
 * written by C. R. Helmrich, last modified in 2021 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under the exhale Copyright License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2021 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#ifndef _TEMP_ANALYSIS_H_
#define _TEMP_ANALYSIS_H_

#include "exhaleLibPch.h"

// constants, experimental macros
#define TA_EPS               4096

// temporal signal analysis class
class TempAnalyzer
{
private:

  // member variables
  unsigned m_avgAbsHpPrev[USAC_MAX_NUM_CHANNELS];
  unsigned m_maxAbsHpPrev[USAC_MAX_NUM_CHANNELS];
  int32_t  m_maxHfLevPrev[USAC_MAX_NUM_CHANNELS];
  unsigned m_maxIdxHpPrev[USAC_MAX_NUM_CHANNELS];
  unsigned m_pitchLagPrev[USAC_MAX_NUM_CHANNELS];
  int64_t  m_filtSampPrev[USAC_MAX_NUM_CHANNELS][6]; // for SBR subband calculation (NOTE: only approximate)
  uint32_t m_tempAnaStats[USAC_MAX_NUM_CHANNELS];
  int16_t  m_transientLoc[USAC_MAX_NUM_CHANNELS];

public:

  // constructor
  TempAnalyzer ();
  // destructor
  ~TempAnalyzer () { }
  // public functions
  void getTempAnalysisStats (uint32_t avgTempAnaStats[USAC_MAX_NUM_CHANNELS], const unsigned nChannels);
  void getTransientAndPitch (int16_t transIdxAndPitch[USAC_MAX_NUM_CHANNELS], const unsigned nChannels);
  uint8_t stereoPreAnalysis (const int32_t* const timeSignals[2], const uint8_t specFlatness[2], const unsigned nSamplesInSig);
  unsigned temporalAnalysis (const int32_t* const timeSignals[USAC_MAX_NUM_CHANNELS], const unsigned nChannels,
                             const int nSamplesInFrame, const unsigned lookaheadOffset, const uint8_t sbrShift,
                             int32_t* const lrCoreTimeSignals[USAC_MAX_NUM_CHANNELS] = nullptr, // if using SBR
                             const unsigned lfeChannelIndex = USAC_MAX_NUM_CHANNELS); // to skip an LFE channel
}; // TempAnalyzer

#endif // _TEMP_ANALYSIS_H_
