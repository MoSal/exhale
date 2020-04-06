/* stereoProcessing.h - header file for class providing M/S stereo coding functionality
 * written by C. R. Helmrich, last modified in 2020 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under a Modified BSD-Style License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2020 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#ifndef _STEREO_PROCESSING_H_
#define _STEREO_PROCESSING_H_

#include "exhaleLibPch.h"
#include "specAnalysis.h" // for SA_BW... constants

// constants, experimental macros
#define SP_EPS                  1
#define SP_OPT_ALPHA_QUANT      1 // quantize alpha_q minimizing RMS distortion in louder channel

// joint-channel processing class
class StereoProcessor
{
private:

  // member variables
  uint8_t m_stereoCorrValue[1024 >> SA_BW_SHIFT]; // one value for every 32 spectral coefficients

public:

  // constructor
  StereoProcessor ();
  // destructor
  ~StereoProcessor () { }
  // public functions
  unsigned applyFullFrameMatrix (int32_t* const mdctSpectrum1, int32_t* const mdctSpectrum2,
                                 int32_t* const mdstSpectrum1, int32_t* const mdstSpectrum2,
                                 SfbGroupData&  groupingData1, SfbGroupData&  groupingData2,
                                 const TnsData&   filterData1, const TnsData&   filterData2,
                                 const uint8_t    numSwbFrame, uint8_t* const sfbStereoData,
                                 const uint8_t  useAltPredDir, const uint8_t useComplexCoef,
                                 uint32_t* const sfbStepSize1, uint32_t* const sfbStepSize2);
}; // StereoProcessor

#endif // _STEREO_PROCESSING_H_