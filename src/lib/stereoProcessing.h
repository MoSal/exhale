/* stereoProcessing.h - header file for class providing M/S stereo coding functionality
 * written by C. R. Helmrich, last modified in 2022 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under the exhale Copyright License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2021 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#ifndef _STEREO_PROCESSING_H_
#define _STEREO_PROCESSING_H_

#include "exhaleLibPch.h"
#include "specAnalysis.h" // for SA_BW... constants
#include <random>

// constants, experimental macros
#define SP_0_DOT_1_16BIT     6554
#define SP_EPS                  1
#define SP_DIV (1.0 / 4294967296.0)

// joint-channel processing class
class StereoProcessor
{
private:

  // member variables
  int32_t m_originBandMdct1[320]; // i.e. 64 * 5 - NOTE: increase this when maximum grpLength > 5
  int32_t m_originBandMdct2[320];
  int32_t m_originBandMdst1[320];
  int32_t m_originBandMdst2[320];
  std::minstd_rand m_randomInt32;
  int32_t m_randomIntMemRe[1+MAX_NUM_SWB_LONG/2];
  int32_t m_randomIntMemIm[1+MAX_NUM_SWB_LONG/2];
  uint8_t m_stereoCorrValue[1024 >> SA_BW_SHIFT]; // one value for every 32 spectral coefficients

public:

  // constructor
  StereoProcessor ();
  // destructor
  ~StereoProcessor () { }
  // public functions
  unsigned applyPredJointStereo (int32_t* const mdctSpectrum1, int32_t* const mdctSpectrum2,
                                 int32_t* const mdstSpectrum1, int32_t* const mdstSpectrum2,
                                 SfbGroupData&  groupingData1, SfbGroupData&  groupingData2,
                                 const TnsData&   filterData1, const TnsData&   filterData2,
                                 const uint8_t    numSwbFrame, uint8_t* const sfbStereoData,
                                 const uint8_t    bitRateMode, const bool    useFullFrameMS,
                                 const bool    reversePredDir, const uint8_t realOnlyOffset,
                                 uint32_t* const sfbStepSize1, uint32_t* const sfbStepSize2);
}; // StereoProcessor

#endif // _STEREO_PROCESSING_H_
