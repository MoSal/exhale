/* stereoProcessing.cpp - source file for class providing M/S stereo coding functionality
 * written by C. R. Helmrich, last modified in 2020 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under a Modified BSD-Style License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2020 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#include "exhaleLibPch.h"
#include "stereoProcessing.h"

// static helper functions

// private helper function

// constructor
StereoProcessor::StereoProcessor ()
{
  for (unsigned ch = 0; ch < USAC_MAX_NUM_CHANNELS; ch++)
  {
    m_avgAbsHpPrev[ch] = 0;
    m_maxAbsHpPrev[ch] = 0;
    m_maxIdxHpPrev[ch] = 1;
    m_pitchLagPrev[ch] = 0;
    m_tempAnaStats[ch] = 0;
    m_transientLoc[ch] = -1;
  }
}

// public functions
unsigned StereoProcessor::applyFullFrameMatrix (int32_t* const mdctSpectrum1, int32_t* const mdctSpectrum2,
                                                int32_t* const mdstSpectrum1, int32_t* const mdstSpectrum2,
                                                SfbGroupData&  groupingData1, SfbGroupData&  groupingData2,
                                                const TnsData&   filterData1, const TnsData&   filterData2,
                                                const uint8_t    numSwbFrame, uint8_t* const sfbStereoData,
                                                uint32_t* const sfbStepSize1, uint32_t* const sfbStepSize2)
{
//const bool applyPredSte = (sfbStereoData != nullptr); // use real-valued predictive stereo
  const uint8_t maxSfbSte = __max (groupingData1.sfbsPerGroup, groupingData2.sfbsPerGroup);

  if ((mdctSpectrum1 == nullptr) || (mdctSpectrum2 == nullptr) || (groupingData1.numWindowGroups != groupingData2.numWindowGroups) ||
      (sfbStepSize1  == nullptr) || (sfbStepSize2  == nullptr) || (numSwbFrame < MIN_NUM_SWB_SHORT) || (numSwbFrame > MAX_NUM_SWB_LONG))
  {
    return 1; // invalid arguments error
  }

  for (uint16_t gr = 0; gr < groupingData1.numWindowGroups; gr++)
  {
    const bool realOnlyCalc = (filterData1.numFilters > 0 && gr == filterData1.filteredWindow) || (mdstSpectrum1 == nullptr) || 
                              (filterData2.numFilters > 0 && gr == filterData2.filteredWindow) || (mdstSpectrum2 == nullptr);
    const uint16_t*  grpOff = &groupingData1.sfbOffsets[numSwbFrame * gr];
    uint32_t* const grpRms1 = &groupingData1.sfbRmsValues[numSwbFrame * gr];
    uint32_t* const grpRms2 = &groupingData2.sfbRmsValues[numSwbFrame * gr];
    uint32_t* grpStepSizes1 = &sfbStepSize1[numSwbFrame * gr];
    uint32_t* grpStepSizes2 = &sfbStepSize2[numSwbFrame * gr];
    int32_t   prevReM = 0, prevReS = 0;

    if (realOnlyCalc) // preparation for first magnitude value
    {
      const uint16_t sPlus1 = grpOff[0] + 1;

      prevReM = int32_t (((int64_t) mdctSpectrum1[sPlus1] + (int64_t) mdctSpectrum2[sPlus1] + 1) >> 1);
      prevReS = int32_t (((int64_t) mdctSpectrum1[sPlus1] - (int64_t) mdctSpectrum2[sPlus1] + 1) >> 1);
    }

    for (uint16_t sfb = 0; sfb < maxSfbSte; sfb++)
    {
      const uint32_t sfbRmsL  = __max (SP_EPS, grpRms1[sfb]);
      const uint32_t sfbRmsR  = __max (SP_EPS, grpRms2[sfb]);
      const double   sfbFacLR = (sfbRmsL < (grpStepSizes1[sfb] >> 1) ? 1.0 : 2.0) * (sfbRmsR < (grpStepSizes2[sfb] >> 1) ? 1.0 : 2.0);
      const double   sfbRatLR = __min (1.0, grpStepSizes1[sfb] / (sfbRmsL * 2.0)) * __min (1.0, grpStepSizes2[sfb] / (sfbRmsR * 2.0)) * sfbFacLR;
      const uint16_t sfbStart = grpOff[sfb];
      const uint16_t sfbWidth = grpOff[sfb + 1] - sfbStart;
      int32_t* sfbMdct1 = &mdctSpectrum1[sfbStart];
      int32_t* sfbMdct2 = &mdctSpectrum2[sfbStart];
      double   sfbRmsMaxMS;
      uint64_t sumAbsValM = 0, sumAbsValS = 0;

      if (realOnlyCalc) // real data, only MDCTs are available
      {
        int32_t* sfbNext1 = &sfbMdct1[1];
        int32_t* sfbNext2 = &sfbMdct2[1];

        for (uint16_t s = sfbWidth - (sfb + 1 == numSwbFrame ? 1 : 0); s > 0; s--)
        {
          const int32_t dmixReM = int32_t (((int64_t) *sfbMdct1 + (int64_t) *sfbMdct2 + 1) >> 1);
          const int32_t dmixReS = int32_t (((int64_t) *sfbMdct1 - (int64_t) *sfbMdct2 + 1) >> 1);
          // TODO: improve the following lines since the calculation is partially redundant!
          const int32_t dmixImM = int32_t (((*sfbNext1 + (int64_t) *sfbNext2 + 1) >> 1) - (int64_t) prevReM) >> 1; // estimate, see also
          const int32_t dmixImS = int32_t (((*sfbNext1 - (int64_t) *sfbNext2 + 1) >> 1) - (int64_t) prevReS) >> 1; // getMeanAbsValues()

          const uint32_t absReM = abs (dmixReM);
          const uint32_t absReS = abs (dmixReS);   // Richard Lyons, 1997; en.wikipedia.org/
          const uint32_t absImM = abs (dmixImM);   // wiki/Alpha_max_plus_beta_min_algorithm
          const uint32_t absImS = abs (dmixImS);

          sumAbsValM += (absReM > absImM ? absReM + ((absImM * 3) >> 3) : absImM + ((absReM * 3) >> 3));
          sumAbsValS += (absReS > absImS ? absReS + ((absImS * 3) >> 3) : absImS + ((absReS * 3) >> 3));

          *(sfbMdct1++) = dmixReM;
          *(sfbMdct2++) = dmixReS;
          sfbNext1++; prevReM = dmixReM;
          sfbNext2++; prevReS = dmixReS;
        }
      }
      else // complex data, both MDCTs and MDSTs are available
      {
        int32_t* sfbMdst1 = &mdstSpectrum1[sfbStart];
        int32_t* sfbMdst2 = &mdstSpectrum2[sfbStart];

        for (uint16_t s = sfbWidth; s > 0; s--)
        {
          const int32_t dmixReM = int32_t (((int64_t) *sfbMdct1 + (int64_t) *sfbMdct2 + 1) >> 1);
          const int32_t dmixReS = int32_t (((int64_t) *sfbMdct1 - (int64_t) *sfbMdct2 + 1) >> 1);
          const int32_t dmixImM = int32_t (((int64_t) *sfbMdst1 + (int64_t) *sfbMdst2 + 1) >> 1);
          const int32_t dmixImS = int32_t (((int64_t) *sfbMdst1 - (int64_t) *sfbMdst2 + 1) >> 1);
#if SA_EXACT_COMPLEX_ABS
          const double cplxSqrM = (double) dmixReM * (double) dmixReM + (double) dmixImM * (double) dmixImM;
          const double cplxSqrS = (double) dmixReS * (double) dmixReS + (double) dmixImS * (double) dmixImS;

          sumAbsValM += uint64_t (sqrt (cplxSqrM) + 0.5);
          sumAbsValS += uint64_t (sqrt (cplxSqrS) + 0.5);
#else
          const uint32_t absReM = abs (dmixReM);
          const uint32_t absReS = abs (dmixReS);   // Richard Lyons, 1997; en.wikipedia.org/
          const uint32_t absImM = abs (dmixImM);   // wiki/Alpha_max_plus_beta_min_algorithm
          const uint32_t absImS = abs (dmixImS);

          sumAbsValM += (absReM > absImM ? absReM + ((absImM * 3) >> 3) : absImM + ((absReM * 3) >> 3));
          sumAbsValS += (absReS > absImS ? absReS + ((absImS * 3) >> 3) : absImS + ((absReS * 3) >> 3));
#endif
          *(sfbMdct1++) = dmixReM;
          *(sfbMdct2++) = dmixReS;
          *(sfbMdst1++) = dmixImM;
          *(sfbMdst2++) = dmixImS;
        }
      } // realOnlyCalc

      // average spectral sample magnitude across current band
      grpRms1[sfb] = uint32_t ((sumAbsValM + (sfbWidth >> 1)) / sfbWidth);
      grpRms2[sfb] = uint32_t ((sumAbsValS + (sfbWidth >> 1)) / sfbWidth);
      sfbRmsMaxMS  = __max (grpRms1[sfb], grpRms2[sfb]);

      if (sfbFacLR <= 1.0) // total simultaneous masking, no positive SNR in any channel SFB
      {
        double max = __max (sfbRmsL, sfbRmsR);
        grpStepSizes1[sfb] = grpStepSizes2[sfb] = uint32_t (__max (grpStepSizes1[sfb], grpStepSizes2[sfb]) * (sfbRmsMaxMS / max) + 0.5);
      }
      else // partial or no masking, redistribute positive SNR into at least one channel SFB
      {
        double min = __min (grpRms1[sfb], grpRms2[sfb]);
        grpStepSizes1[sfb] = grpStepSizes2[sfb] = uint32_t (__max (SP_EPS, (min > sfbRatLR * sfbRmsMaxMS ? sqrt (sfbRatLR * sfbRmsMaxMS *
                                                                            min) : __min (1.0/*TODO*/, sfbRatLR) * sfbRmsMaxMS)) + 0.5);
      }
    } // for sfb
  }

  return 0; // no error
}
