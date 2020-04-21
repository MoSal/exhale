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

// constructor
StereoProcessor::StereoProcessor ()
{
  memset (m_stereoCorrValue, 0, (1024 >> SA_BW_SHIFT) * sizeof (uint8_t));
}

// public functions
unsigned StereoProcessor::applyPredJointStereo (int32_t* const mdctSpectrum1, int32_t* const mdctSpectrum2,
                                                int32_t* const mdstSpectrum1, int32_t* const mdstSpectrum2,
                                                SfbGroupData&  groupingData1, SfbGroupData&  groupingData2,
                                                const TnsData&   filterData1, const TnsData&   filterData2,
                                                const uint8_t    numSwbFrame, uint8_t* const sfbStereoData,
                                                const bool    usePerCorrData, const bool    useFullFrameMS,
                                                const bool    reversePredDir, const bool    useComplexCoef,
                                                uint32_t* const sfbStepSize1, uint32_t* const sfbStepSize2)
{
  const bool applyPredSte = (sfbStereoData != nullptr); // use real-valued predictive stereo
  const bool alterPredDir = (applyPredSte && reversePredDir); // predict mid from side band?
  const SfbGroupData& grp = groupingData1;
  const bool  eightShorts = (grp.numWindowGroups > 1);
  const uint8_t maxSfbSte = (eightShorts ? __min (numSwbFrame, __max (grp.sfbsPerGroup, groupingData2.sfbsPerGroup) + 1) : numSwbFrame);
  const bool  perCorrData = (usePerCorrData && !eightShorts); // use perceptual correlation?
  uint32_t  numSfbPredSte = 0; // counter

  if ((mdctSpectrum1 == nullptr) || (mdctSpectrum2 == nullptr) || (numSwbFrame < maxSfbSte) || (grp.numWindowGroups != groupingData2.numWindowGroups) ||
      (sfbStepSize1  == nullptr) || (sfbStepSize2  == nullptr) || (numSwbFrame < MIN_NUM_SWB_SHORT) || (numSwbFrame > MAX_NUM_SWB_LONG))
  {
    return 1;  // invalid arguments error
  }
#if !SP_SFB_WISE_STEREO
  if (!useFullFrameMS)
  {
    if (applyPredSte) memset (sfbStereoData, 0, (MAX_NUM_SWB_SHORT * NUM_WINDOW_GROUPS) * sizeof (uint8_t));

    return 0; // zeroed ms_used, no pred.
  }
#endif

  if (applyPredSte && perCorrData) memcpy (m_stereoCorrValue, sfbStereoData, (grp.sfbOffsets[numSwbFrame] >> SA_BW_SHIFT) * sizeof (uint8_t));

  for (uint16_t gr = 0; gr < grp.numWindowGroups; gr++)
  {
    const bool realOnlyCalc = (filterData1.numFilters > 0 && gr == filterData1.filteredWindow) || (mdstSpectrum1 == nullptr) ||
                              (filterData2.numFilters > 0 && gr == filterData2.filteredWindow) || (mdstSpectrum2 == nullptr);
    const uint16_t*  grpOff = &grp.sfbOffsets[numSwbFrame * gr];
    uint32_t* const grpRms1 = &groupingData1.sfbRmsValues[numSwbFrame * gr];
    uint32_t* const grpRms2 = &groupingData2.sfbRmsValues[numSwbFrame * gr];
    uint32_t* grpStepSizes1 = &sfbStepSize1[numSwbFrame * gr];
    uint32_t* grpStepSizes2 = &sfbStepSize2[numSwbFrame * gr];
    int32_t  b = 0, prevReM = 0, prevReS = 0;
    uint32_t rmsSfbL[2] = {0, 0}, rmsSfbR[2] = {0, 0};

    if (realOnlyCalc) // preparation for first magnitude value
    {
      const uint16_t sPlus1 = grpOff[0] + 1;

      prevReM = int32_t (((int64_t) mdctSpectrum1[sPlus1] + (int64_t) mdctSpectrum2[sPlus1] + 1) >> 1);
      prevReS = int32_t (((int64_t) mdctSpectrum1[sPlus1] - (int64_t) mdctSpectrum2[sPlus1] + 1) >> 1);
    }

    for (uint16_t sfb = 0; sfb < maxSfbSte; sfb++)
    {
      const int32_t  sfbIsOdd = sfb & 1;
      const uint16_t sfbStart = grpOff[sfb];
      const uint16_t sfbWidth = grpOff[sfb + 1] - sfbStart;
      int32_t* sfbMdct1 = &mdctSpectrum1[sfbStart];
      int32_t* sfbMdct2 = &mdctSpectrum2[sfbStart];
#if SP_MDST_PRED
      int32_t* sfbMdst1 = &mdstSpectrum1[sfbStart];
      int32_t* sfbMdst2 = &mdstSpectrum2[sfbStart];
#endif
      uint64_t sumAbsValM = 0, sumAbsValS = 0;
      double   sfbTempVar;

#if SP_SFB_WISE_STEREO
      if ((sfbIsOdd == 0) && !useFullFrameMS) // save L/R data
      {
        const uint16_t cpyWidth = (grpOff[__min (maxSfbSte, sfb + 2)] - sfbStart) * sizeof (int32_t);

        memcpy (m_originBandMdct1, sfbMdct1, cpyWidth);
        memcpy (m_originBandMdct2, sfbMdct2, cpyWidth);
        memcpy (m_originBandMdst1, sfbMdst1, cpyWidth);
        memcpy (m_originBandMdst2, sfbMdst2, cpyWidth);
      }
#endif
      if (realOnlyCalc) // real data, only MDCTs are available
      {
        const int32_t* sfbNext1 = &sfbMdct1[1];
        const int32_t* sfbNext2 = &sfbMdct2[1];

        for (uint16_t s = sfbWidth - (sfb + 1 == numSwbFrame ? 1 : 0); s > 0; s--)
        {
          const int32_t dmixReM = int32_t (((int64_t) *sfbMdct1 + (int64_t) *sfbMdct2 + 1) >> 1);
          const int32_t dmixReS = int32_t (((int64_t) *sfbMdct1 - (int64_t) *sfbMdct2 + 1) >> 1);
          // TODO: improve the following lines since the calculation is partially redundant!
          const int32_t dmixImM = int32_t ((((*sfbNext1 + (int64_t) *sfbNext2 + 1) >> 1) - (int64_t) prevReM) >> 1); // estimate, see also
          const int32_t dmixImS = int32_t ((((*sfbNext1 - (int64_t) *sfbNext2 + 1) >> 1) - (int64_t) prevReS) >> 1); // getMeanAbsValues()

          const uint64_t absReM = abs (dmixReM);
          const uint64_t absReS = abs (dmixReS);   // Richard Lyons, 1997; en.wikipedia.org/
          const uint64_t absImM = abs (dmixImM);   // wiki/Alpha_max_plus_beta_min_algorithm
          const uint64_t absImS = abs (dmixImS);

          sumAbsValM += (absReM > absImM ? absReM + ((absImM * 3) >> 3) : absImM + ((absReM * 3) >> 3));
          sumAbsValS += (absReS > absImS ? absReS + ((absImS * 3) >> 3) : absImS + ((absReS * 3) >> 3));

          *(sfbMdct1++) = dmixReM;
          *(sfbMdct2++) = dmixReS;
#if SP_MDST_PRED
          *(sfbMdst1++) = dmixImM;
          *(sfbMdst2++) = dmixImS;
#endif
          sfbNext1++; prevReM = dmixReM;
          sfbNext2++; prevReS = dmixReS;
        }
        if (sfb + 1 == numSwbFrame) // handle remaining sample
        {
          const int32_t dmixReM = int32_t (((int64_t) *sfbMdct1 + (int64_t) *sfbMdct2 + 1) >> 1);
          const int32_t dmixReS = int32_t (((int64_t) *sfbMdct1 - (int64_t) *sfbMdct2 + 1) >> 1);

          sumAbsValM += abs (dmixReM);
          sumAbsValS += abs (dmixReS);

          *sfbMdct1 = dmixReM;
          *sfbMdct2 = dmixReS;
#if SP_MDST_PRED
          *sfbMdst1 = 0;
          *sfbMdst2 = 0;
#endif
        }
      }
      else // complex data, both MDCTs and MDSTs are available
      {
#if !SP_MDST_PRED
        int32_t* sfbMdst1 = &mdstSpectrum1[sfbStart];
        int32_t* sfbMdst2 = &mdstSpectrum2[sfbStart];
#endif
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
          const uint64_t absReM = abs (dmixReM);
          const uint64_t absReS = abs (dmixReS);   // Richard Lyons, 1997; en.wikipedia.org/
          const uint64_t absImM = abs (dmixImM);   // wiki/Alpha_max_plus_beta_min_algorithm
          const uint64_t absImS = abs (dmixImS);

          sumAbsValM += (absReM > absImM ? absReM + ((absImM * 3) >> 3) : absImM + ((absReM * 3) >> 3));
          sumAbsValS += (absReS > absImS ? absReS + ((absImS * 3) >> 3) : absImS + ((absReS * 3) >> 3));
#endif
          *(sfbMdct1++) = dmixReM;
          *(sfbMdct2++) = dmixReS;
          *(sfbMdst1++) = dmixImM;
          *(sfbMdst2++) = dmixImS;
        }
      } // realOnlyCalc

      rmsSfbL[sfbIsOdd] = grpRms1[sfb];
      rmsSfbR[sfbIsOdd] = grpRms2[sfb];
      // average spectral sample magnitude across current band
      grpRms1[sfb] = uint32_t ((sumAbsValM + (sfbWidth >> 1)) / sfbWidth);
      grpRms2[sfb] = uint32_t ((sumAbsValS + (sfbWidth >> 1)) / sfbWidth);

      if (applyPredSte) sfbStereoData[sfb + numSwbFrame * gr] = 16; // initialize to alpha=0

      if ((sfbIsOdd) || (sfb + 1 == maxSfbSte)) // finish pair
      {
        const uint16_t sfbEv = sfb & 0xFFFE; // even SFB index
        uint32_t  rmsSfbM[2] = {0, 0}, rmsSfbS[2] = {0, 0};
        bool nonZeroPredCoef = false;

        if (applyPredSte) // calc real-prediction coefficients
        {
          const uint16_t offEv = grpOff[sfbEv];
          const uint16_t width = grpOff[sfb + 1] - offEv;
          const int32_t* mdctA = (alterPredDir ? &mdctSpectrum2[offEv] : &mdctSpectrum1[offEv]);
          const int32_t* mdctB = (alterPredDir ? &mdctSpectrum1[offEv] : &mdctSpectrum2[offEv]);
#if SP_MDST_PRED
          const int32_t* mdstA = (alterPredDir ? &mdstSpectrum2[offEv] : &mdstSpectrum1[offEv]);
          const int32_t* mdstB = (alterPredDir ? &mdstSpectrum1[offEv] : &mdstSpectrum2[offEv]);
          int64_t sumPrdImAReB = 0, sumPrdImAImA = SP_EPS;
#endif
          int64_t sumPrdReAReB = 0, sumPrdReAReA = SP_EPS;  // stabilizes the division below
          double d, alphaLimit = 1.5; // max alpha_q magnitude

#if SP_MDST_PRED
          for (uint16_t s = width; s > 0; s--, mdctA++, mdctB++, mdstA++, mdstB++)
          {
            const int64_t prdReAReA = ((int64_t) *mdctA * (int64_t) *mdctA + SA_BW) >> (SA_BW_SHIFT + 1);
            const int64_t prdImAImA = ((int64_t) *mdstA * (int64_t) *mdstA + SA_BW) >> (SA_BW_SHIFT + 1);

            sumPrdReAReB += ((int64_t) *mdctA * (int64_t) *mdctB + SA_BW) >> (SA_BW_SHIFT + 1);
            sumPrdReAReA += prdReAReA;
            sumPrdImAReB += ((int64_t) *mdstA * (int64_t) *mdctB + SA_BW) >> (SA_BW_SHIFT + 1);
            sumPrdImAImA += prdImAImA;
            // add complex conjugate part, increases stability
            sumPrdReAReB += ((int64_t) *mdstA * (int64_t) *mdstB + SA_BW) >> (SA_BW_SHIFT + 1);
            sumPrdReAReA += prdImAImA;
            sumPrdImAReB -= ((int64_t) *mdctA * (int64_t) *mdstB + SA_BW) >> (SA_BW_SHIFT + 1);
            sumPrdImAImA += prdReAReA;
          }
#else
          for (uint16_t s = width; s > 0; s--, mdctA++, mdctB++)
          {
            sumPrdReAReB += ((int64_t) *mdctA * (int64_t) *mdctB + SA_BW) >> (SA_BW_SHIFT + 1);
            sumPrdReAReA += ((int64_t) *mdctA * (int64_t) *mdctA + SA_BW) >> (SA_BW_SHIFT + 1);
          }
          if (realOnlyCalc) // real data, only MDCTs available
          {
            const int32_t* nextA = (alterPredDir ? &mdctSpectrum2[offEv + 1] : &mdctSpectrum1[offEv + 1]);
            const int32_t* nextB = (alterPredDir ? &mdctSpectrum1[offEv + 1] : &mdctSpectrum2[offEv + 1]);

            if (sfbEv == 0) // exclude unavailable DC estimate
            {
              mdctA -= width; nextA++;
              mdctB -= width; nextB++;
            }
            else
            {
              mdctA -= (width + 1); // point to pre-
              mdctB -= (width + 1); // vious samples
            }
            for (uint16_t s = width - (sfbEv == 0 ? 2 : 1); s > 0; s--) // exclude unavailable final estimate
            {
              const int64_t mdstA = ((int64_t) *(nextA++) - (int64_t) *(mdctA++)) >> 1; // estimate, see also
              const int64_t mdstB = ((int64_t) *(nextB++) - (int64_t) *(mdctB++)) >> 1; // getMeanAbsValues()

              sumPrdReAReB += (mdstA * mdstB + SA_BW) >> (SA_BW_SHIFT + 1);
              sumPrdReAReA += (mdstA * mdstA + SA_BW) >> (SA_BW_SHIFT + 1);
            }
            if (sfb + 1 < numSwbFrame) // process final sample
            {
              const int64_t msSgn = (alterPredDir ? -1 : 1);
              const int64_t mdstA = ((((int64_t) *nextB + *nextA * msSgn + 1) >> 1) - (int64_t) *mdctA) >> 1;
              const int64_t mdstB = ((((int64_t) *nextA - *nextB * msSgn + 1) >> 1) - (int64_t) *mdctB) >> 1;

              sumPrdReAReB += (mdstA * mdstB + SA_BW) >> (SA_BW_SHIFT + 1);
              sumPrdReAReA += (mdstA * mdstA + SA_BW) >> (SA_BW_SHIFT + 1);
            }
          }
          else // complex data, both MDCTs and MDSTs available
          {
            const int32_t* mdstA = (alterPredDir ? &mdstSpectrum2[offEv] : &mdstSpectrum1[offEv]);
            const int32_t* mdstB = (alterPredDir ? &mdstSpectrum1[offEv] : &mdstSpectrum2[offEv]);

            for (uint16_t s = width; s > 0; s--, mdstA++, mdstB++)
            {
              sumPrdReAReB += ((int64_t) *mdstA * (int64_t) *mdstB + SA_BW) >> (SA_BW_SHIFT + 1);
              sumPrdReAReA += ((int64_t) *mdstA * (int64_t) *mdstA + SA_BW) >> (SA_BW_SHIFT + 1);
            }
          }
#endif
          for (b = sfbIsOdd; b >= 0; b--) // limit alpha_q to prevent residual RMS increases
          {
            const int idx = sfbEv + b;

            d = (alterPredDir ? (double) grpRms1[idx] / __max (SP_EPS, grpRms2[idx]) : (double) grpRms2[idx] / __max (SP_EPS, grpRms1[idx]));
            if (alphaLimit > d) alphaLimit = d;
          }
          sfbTempVar = CLIP_PM ((double) sumPrdReAReB / (double) sumPrdReAReA, alphaLimit);
#if SP_OPT_ALPHA_QUANT
          b = __max (512, 524 - int32_t (abs (10.0 * sfbTempVar))); // rounding optimization
          b = int32_t (10.0 * sfbTempVar + b * (sfbTempVar < 0 ? -0.0009765625 : 0.0009765625));
#else
          b = int32_t (10.0 * sfbTempVar + (sfbTempVar < 0 ? -0.5 : 0.5));// nearest integer
#endif
          sfbStereoData[sfbEv + numSwbFrame * gr] = uint8_t (b + 16); // finished alpha_q_re
#if SP_MDST_PRED
          alphaLimit = CLIP_PM ((double) sumPrdImAReB / (double) sumPrdImAImA, alphaLimit);
# if SP_OPT_ALPHA_QUANT
          b = __max (512, 524 - int32_t (abs (10.0 * alphaLimit))); // rounding optimization
          b = int32_t (10.0 * alphaLimit + b * (alphaLimit < 0 ? -0.0009765625 : 0.0009765625));
# else
          b = int32_t (10.0 * alphaLimit + (alphaLimit < 0 ? -0.5 : 0.5));// nearest integer
# endif
          if (sfbEv + 1 < numSwbFrame)
          sfbStereoData[sfbEv + 1 + numSwbFrame * gr] = uint8_t (b + 16); // init alpha_q_im
#endif

          if (perCorrData && ((offEv & (SA_BW - 1)) == 0) && ((width & (SA_BW - 1)) == 0))
          {
            const uint8_t* const perCorr = &m_stereoCorrValue[offEv >> SA_BW_SHIFT];

            // perceptual correlation data available from previous call to stereoSigAnalysis
            b = (width == SA_BW ? perCorr[0] : ((int32_t) perCorr[0] + (int32_t) perCorr[1] + 1) >> 1);
          }
          else b = UCHAR_MAX; // previous correlation data unavailable, assume maximum value

          if ((b > SCHAR_MAX && sfbStereoData[sfbEv + numSwbFrame * gr] != 16) ||  // signi-
              (2 <= abs ( (int) sfbStereoData[sfbEv + numSwbFrame * gr] - 16)))   // ficant?
          {
            nonZeroPredCoef = true;
          }
          sfbTempVar *= sfbTempVar;  // account for residual RMS reduction due to prediction
#if SP_MDST_PRED
          sfbTempVar += alphaLimit * alphaLimit; // including complex prediction by alpha_im
#endif
          for (b = sfbIsOdd; b >= 0; b--)
          {
            const int idx = sfbEv + b;

            if (alterPredDir)
            {
              d = (double) grpRms1[idx] * grpRms1[idx] - sfbTempVar * (double) grpRms2[idx] * grpRms2[idx];
              // consider discarding prediction if gain (residual RMS loss) is below -0.9 dB
              if ((double) grpRms1[idx] * grpRms1[idx] * 0.8125 < d) nonZeroPredCoef = false;
              rmsSfbM[b] = uint32_t (sqrt (__max (0.0, d)) + 0.5);
              rmsSfbS[b] = grpRms2[idx];
            }
            else // mid>side
            {
              d = (double) grpRms2[idx] * grpRms2[idx] - sfbTempVar * (double) grpRms1[idx] * grpRms1[idx];
              // consider discarding prediction if gain (residual RMS loss) is below -0.9 dB
              if ((double) grpRms2[idx] * grpRms2[idx] * 0.8125 < d) nonZeroPredCoef = false;
              rmsSfbS[b] = uint32_t (sqrt (__max (0.0, d)) + 0.5);
              rmsSfbM[b] = grpRms1[idx];
            }
          }
        } // if applyPredSte

#if SP_SFB_WISE_STEREO
        if (!useFullFrameMS) // test M/S compaction gain, revert to L/R if it's insufficient
        {
          const uint64_t bandSum1 = (sfbIsOdd > 0 ? (uint64_t) grpRms1[sfbEv] + (uint64_t) grpRms1[sfbEv + 1] : grpRms1[sfbEv]);
          const uint64_t bandSum2 = (sfbIsOdd > 0 ? (uint64_t) grpRms2[sfbEv] + (uint64_t) grpRms2[sfbEv + 1] : grpRms2[sfbEv]);
          const uint64_t bandSumL = (sfbIsOdd > 0 ? (uint64_t) rmsSfbL[0] + (uint64_t) rmsSfbL[1] : rmsSfbL[0]) >> 1;
          const uint64_t bandSumR = (sfbIsOdd > 0 ? (uint64_t) rmsSfbR[0] + (uint64_t) rmsSfbR[1] : rmsSfbR[0]) >> 1;
          const uint64_t bandSumM = (applyPredSte ? (uint64_t) rmsSfbM[0] + (uint64_t) rmsSfbM[1] : bandSum1) >> 1;
          const uint64_t bandSumS = (applyPredSte ? (uint64_t) rmsSfbS[0] + (uint64_t) rmsSfbS[1] : bandSum2) >> 1;

          if ((__min (bandSumM, bandSumS) * __max (bandSumL, bandSumR) >= __min (bandSumL, bandSumR) * __max (bandSumM, bandSumS)) ||
              (nonZeroPredCoef && (abs ( (int) sfbStereoData[sfbEv + numSwbFrame * gr] - 16) >= 10)))
          {
            const uint16_t sfbOffEv = grpOff[sfbEv];
            const uint16_t cpyWidth = (grpOff[sfb + 1] - sfbOffEv) * sizeof (int32_t);

            memcpy (&mdctSpectrum1[sfbOffEv], m_originBandMdct1, cpyWidth); // revert to L/R
            memcpy (&mdctSpectrum2[sfbOffEv], m_originBandMdct2, cpyWidth);
            memcpy (&mdstSpectrum1[sfbOffEv], m_originBandMdst1, cpyWidth);
            memcpy (&mdstSpectrum2[sfbOffEv], m_originBandMdst2, cpyWidth);

            for (b = sfbIsOdd; b >= 0; b--)
            {
              const int idx = sfbEv + b;

              grpRms1[idx] = rmsSfbL[b];
              grpRms2[idx] = rmsSfbR[b];
              if (applyPredSte) sfbStereoData[idx + numSwbFrame * gr] = 0; // zeroed ms_used
            }
            continue; // M/S is not used
          }
        }
#endif
        if (nonZeroPredCoef) numSfbPredSte++;  // a perceptually significant prediction band

        for (b = sfbIsOdd; b >= 0; b--)
        {
          const int idx = sfbEv + b;
          const uint32_t sfbRmsL = __max (SP_EPS, rmsSfbL[b]);
          const uint32_t sfbRmsR = __max (SP_EPS, rmsSfbR[b]);
          const double  sfbFacLR = (sfbRmsL < (grpStepSizes1[idx] >> 1) ? 1.0 : 2.0) * (sfbRmsR < (grpStepSizes2[idx] >> 1) ? 1.0 : 2.0);

          sfbTempVar = (applyPredSte ? __max (rmsSfbM[b], rmsSfbS[b]) : __max (grpRms1[idx], grpRms2[idx]));

          if ((grpStepSizes1[idx] == 0) || (grpStepSizes2[idx] == 0)) // HF noise filled SFB
          {
            grpStepSizes1[idx] = grpStepSizes2[idx] = 0;
          }
          else if (sfbFacLR <= 1.0)  // simultaneous masking - no positive SNR in either SFB
          {
            const double max = __max (sfbRmsL, sfbRmsR);

            grpStepSizes1[idx] = grpStepSizes2[idx] = uint32_t (__max (grpStepSizes1[idx], grpStepSizes2[idx]) * (sfbTempVar / max) + 0.5);
          }
          else // partial/no masking - redistribute positive SNR into at least 1 channel SFB
          {
            const double min = (applyPredSte ? __min (rmsSfbM[b], rmsSfbS[b]) : __min (grpRms1[idx], grpRms2[idx]));
            const double rat = __min (1.0, grpStepSizes1[idx] / (sfbRmsL * 2.0)) * __min (1.0, grpStepSizes2[idx] / (sfbRmsR * 2.0)) * sfbFacLR;

            grpStepSizes1[idx] = grpStepSizes2[idx] = uint32_t (__max (SP_EPS, (min > rat * sfbTempVar ? sqrt (rat * sfbTempVar * min) :
                                                                                __min (1.0, rat) * sfbTempVar)) + 0.5);
          }
        }
      } // if pair completed
    }
  } // for gr

  if (numSfbPredSte == 0) // discard prediction coefficients and stay with legacy M/S stereo
  {
    if (applyPredSte)
    for (uint16_t gr = 0; gr < grp.numWindowGroups; gr++)
    {
      uint8_t* const grpSData = &sfbStereoData[numSwbFrame * gr];

      for (uint16_t sfb = 0; sfb < maxSfbSte; sfb++)
      {
        if (grpSData[sfb] > 0) grpSData[sfb] = 16;
      }
      if (numSwbFrame > maxSfbSte) memset (&grpSData[maxSfbSte], (useFullFrameMS ? 16 : 0), (numSwbFrame - maxSfbSte) * sizeof (uint8_t));
    }
  }
  else // at least one "significant" prediction band, apply prediction and update RMS values
  {
    for (uint16_t gr = 0; gr < grp.numWindowGroups; gr++)
    {
      const bool realOnlyCalc = (filterData1.numFilters > 0 && gr == filterData1.filteredWindow) || (mdstSpectrum1 == nullptr) ||
                                (filterData2.numFilters > 0 && gr == filterData2.filteredWindow) || (mdstSpectrum2 == nullptr);
      const uint16_t*  grpOff = &grp.sfbOffsets[numSwbFrame * gr];
      uint32_t* const grpRms1 = &groupingData1.sfbRmsValues[numSwbFrame * gr];
      uint32_t* const grpRms2 = &groupingData2.sfbRmsValues[numSwbFrame * gr];
      uint8_t* const grpSData = &sfbStereoData[numSwbFrame * gr];
      int32_t prevResi = 0;

      if (realOnlyCalc) // preparation of res. magnitude value
      {
        const int64_t alphaRe = (grpSData[0] > 0 ? (int) grpSData[0] - 16 : 0) * SP_0_DOT_1_16BIT;
        const uint16_t sPlus1 = grpOff[0] + 1;

        prevResi = (alterPredDir ? mdctSpectrum1[sPlus1] - int32_t ((mdctSpectrum2[sPlus1] * alphaRe - SHRT_MIN) >> 16)
                                 : mdctSpectrum2[sPlus1] - int32_t ((mdctSpectrum1[sPlus1] * alphaRe - SHRT_MIN) >> 16));
      }

      for (uint16_t sfb = 0; sfb < maxSfbSte; sfb++)
      {
        const uint16_t sfbEv = sfb & 0xFFFE; // even SFB index
        const uint16_t sfbStart = grpOff[sfb];
        const uint16_t sfbWidth = grpOff[sfb + 1] - sfbStart;
        const int64_t   alphaRe = (grpSData[sfbEv] > 0 ? (int) grpSData[sfbEv] - 16 : 0) * SP_0_DOT_1_16BIT;
        int32_t* sfbMdctD = (alterPredDir ? &mdctSpectrum2[sfbStart] : &mdctSpectrum1[sfbStart]);
        int32_t* sfbMdctR = (alterPredDir ? &mdctSpectrum1[sfbStart] : &mdctSpectrum2[sfbStart]);
        uint64_t sumAbsValR = 0;

        if (alphaRe == 0)
        {
          if (realOnlyCalc) // update previous magnitude value
          {
            sfbMdctR += sfbWidth - 1;
            prevResi = (grpSData[sfbEv] > 0 ? *sfbMdctR : int32_t (((int64_t) sfbMdctD[sfbWidth - 1] +
                                          (alterPredDir ? 1 : -1) * (int64_t) *sfbMdctR + 1) >> 1));
          }
          continue; // nothing more to do, i.e., no prediction
        }

        if (realOnlyCalc) // real data, only MDCT is available
        {
          const int32_t* sfbNextD = &sfbMdctD[1];
          const int32_t* sfbNextR = &sfbMdctR[1];

          for (uint16_t s = sfbWidth - (sfb + 1 == numSwbFrame ? 1 : 0); s > 0; s--)
          {
            const int32_t  resiRe = *sfbMdctR - int32_t ((*sfbMdctD * alphaRe - SHRT_MIN) >> 16);
            // TODO: improve the following line since the calculation is partially redundant
            //       Also, in the final s index of this band, the wrong alphaRe may be used!
            const int32_t  resiIm = int32_t (((*sfbNextR - ((*sfbNextD * alphaRe - SHRT_MIN) >> 16)) - (int64_t) prevResi) >> 1);

            const uint64_t absReR = abs (resiRe);  // Richard Lyons, 1997; en.wikipedia.org/
            const uint64_t absImR = abs (resiIm);  // wiki/Alpha_max_plus_beta_min_algorithm

            sumAbsValR += (absReR > absImR ? absReR + ((absImR * 3) >> 3) : absImR + ((absReR * 3) >> 3));

            sfbMdctD++;
            *(sfbMdctR++) = resiRe;
            sfbNextD++;
            sfbNextR++; prevResi = resiRe;
          }
          if (sfb + 1 == numSwbFrame)  // process final sample
          {
            const int32_t resiRe = *sfbMdctR - int32_t ((*sfbMdctD * alphaRe - SHRT_MIN) >> 16);

            sumAbsValR += abs (resiRe);

            *sfbMdctR = resiRe;
          }
        }
        else  // complex data, both MDCT and MDST is available
        {
          int32_t* sfbMdstD = (alterPredDir ? &mdstSpectrum2[sfbStart] : &mdstSpectrum1[sfbStart]);
          int32_t* sfbMdstR = (alterPredDir ? &mdstSpectrum1[sfbStart] : &mdstSpectrum2[sfbStart]);

          for (uint16_t s = sfbWidth; s > 0; s--)
          {
            const int32_t  resiRe = *sfbMdctR - int32_t ((*sfbMdctD * alphaRe - SHRT_MIN) >> 16);
            const int32_t  resiIm = *sfbMdstR - int32_t ((*sfbMdstD * alphaRe - SHRT_MIN) >> 16);
#if SA_EXACT_COMPLEX_ABS
            const double cplxSqrR = (double) resiRe * (double) resiRe + (double) resiIm * (double) resiIm;

            sumAbsValR += uint64_t (sqrt (cplxSqrR) + 0.5);
#else
            const uint64_t absReR = abs (resiRe);  // Richard Lyons, 1997; en.wikipedia.org/
            const uint64_t absImR = abs (resiIm);  // wiki/Alpha_max_plus_beta_min_algorithm

            sumAbsValR += (absReR > absImR ? absReR + ((absImR * 3) >> 3) : absImR + ((absReR * 3) >> 3));
#endif
            sfbMdctD++;
            *(sfbMdctR++) = resiRe;
            sfbMdstD++;
            *(sfbMdstR++) = resiIm;
          }
        } // realOnlyCalc

        // average spectral res. magnitude across current band
        sumAbsValR = (sumAbsValR + (sfbWidth >> 1)) / sfbWidth;
        if (alterPredDir) grpRms1[sfb] = (uint32_t) sumAbsValR; else grpRms2[sfb] = (uint32_t) sumAbsValR;
      }
      if (numSwbFrame > maxSfbSte) memset (&grpSData[maxSfbSte], (useFullFrameMS ? 16 : 0), (numSwbFrame - maxSfbSte) * sizeof (uint8_t));

      if (alterPredDir) // swap channel data when pred_dir = 1
      {
        for (uint16_t sfb = 0; sfb < maxSfbSte; sfb++)
        {
          const uint16_t sfbStart = grpOff[sfb];
          int32_t* sfbMdct1 = &mdctSpectrum1[sfbStart];
          int32_t* sfbMdct2 = &mdctSpectrum2[sfbStart];

          if (grpSData[sfb] == 0) continue; // M/S is not used

          for (uint16_t s = grpOff[sfb + 1] - sfbStart; s > 0; s--)
          {
            const int32_t i = *sfbMdct1;
            *(sfbMdct1++)   = *sfbMdct2;
            *(sfbMdct2++)   = i;
          }
          numSfbPredSte = grpRms1[sfb];
          grpRms1[sfb]  = grpRms2[sfb];
          grpRms2[sfb]  = numSfbPredSte;
        }
      }
    } // for gr

    numSfbPredSte = 2;
  }

  return (numSfbPredSte); // no error
}
