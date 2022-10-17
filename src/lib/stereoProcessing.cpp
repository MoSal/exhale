/* stereoProcessing.cpp - source file for class providing M/S stereo coding functionality
 * written by C. R. Helmrich, last modified in 2022 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under the exhale Copyright License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2021 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#include "exhaleLibPch.h"
#include "stereoProcessing.h"
#include "bitAllocation.h" // define BA_MORE_CBR (more constant bit-rate, experimental!)

// static helper functions
static inline uint64_t complexAbsMS (const int32_t realPart, const int32_t imagPart)
{
#if SA_EXACT_COMPLEX_ABS
  const double complexSqr = (double) realPart * (double) realPart + (double) imagPart * (double) imagPart;

  return uint64_t (sqrt (complexSqr) + 0.5);
#else
  const uint64_t absReal  = abs (realPart); // Richard Lyons, 1997; en.wikipedia.org/
  const uint64_t absImag  = abs (imagPart); // wiki/Alpha_max_plus_beta_min_algorithm

  return (absReal > absImag ? absReal + ((absImag * 3) >> 3) : absImag + ((absReal * 3) >> 3));
#endif
}

static inline int32_t  getMidSample (const int64_t value1, const int64_t value2)
{
  return int32_t ((value1 + value2 + 1) >> 1);
}

static inline int32_t getSideSample (const int64_t value1, const int64_t value2)
{
  return int32_t ((value1 - value2 + 1) >> 1);
}

static inline int32_t getResiSample (const int64_t valueR, const int64_t valueD, const int64_t alphaQ)
{
  return int32_t (valueR - ((valueD * alphaQ - SHRT_MIN) >> 16));
}

static inline void   setStepSizesMS (const uint32_t* const rmsSfbL, const uint32_t* const rmsSfbR,
                                     const uint32_t* const rmsSfbM, const uint32_t* const rmsSfbS,
                                     const uint32_t* const grpRms1, const uint32_t* const grpRms2,
                                     uint32_t* const grpStepSizes1, uint32_t* const grpStepSizes2,
                                     const uint16_t sfb, const uint16_t b, const bool applyPredSte)
{
  const uint16_t     idx = sfb + b;
  const uint32_t sfbRmsL = __max (SP_EPS, rmsSfbL[b]);
  const uint32_t sfbRmsR = __max (SP_EPS, rmsSfbR[b]);
  const double  sfbFacLR = (sfbRmsL < (grpStepSizes1[idx] >> 1) ? 1.0 : 2.0) * (sfbRmsR < (grpStepSizes2[idx] >> 1) ? 1.0 : 2.0);
  const double  sfbMaxMS = (applyPredSte ? __max (rmsSfbM[b], rmsSfbS[b]) : __max (grpRms1[idx], grpRms2[idx]));

  if ((grpStepSizes1[idx] == 0) || (grpStepSizes2[idx] == 0))  // HF noise filled SFB
  {
    grpStepSizes1[idx] = grpStepSizes2[idx] = 0;
  }
  else if (sfbFacLR <= 1.0) // simultaneous masking, so no positive SNR in either SFB
  {
    const double max = __max (sfbRmsL, sfbRmsR);

    grpStepSizes1[idx] = grpStepSizes2[idx] = uint32_t (__max (grpStepSizes1[idx], grpStepSizes2[idx]) * (sfbMaxMS / max) + 0.5);
  }
  else // partial/no masking, redistribute positive SNR into at least one channel SFB
  {
    const double min = (applyPredSte ? __min (rmsSfbM[b], rmsSfbS[b]) : __min (grpRms1[idx], grpRms2[idx]));
    const double rat = __min (1.0, grpStepSizes1[idx] / (sfbRmsL * 2.0)) * __min (1.0, grpStepSizes2[idx] / (sfbRmsR * 2.0)) * sfbFacLR;

    grpStepSizes1[idx] = grpStepSizes2[idx] = uint32_t (__max (SP_EPS, (min > rat * sfbMaxMS ? sqrt (rat * sfbMaxMS * min) :
                                                                        __max (1.0/2048.0, __min (1.0, rat)) * sfbMaxMS)) + 0.5);
  }
}

// constructor
StereoProcessor::StereoProcessor ()
{
  memset (m_randomIntMemRe, 0, (1+MAX_NUM_SWB_LONG/2) * sizeof (int32_t));
  memset (m_randomIntMemIm, 0, (1+MAX_NUM_SWB_LONG/2) * sizeof (int32_t));
  memset (m_stereoCorrValue, 0, (1024 >> SA_BW_SHIFT) * sizeof (uint8_t));
}

// public functions
unsigned StereoProcessor::applyPredJointStereo (int32_t* const mdctSpectrum1, int32_t* const mdctSpectrum2,
                                                int32_t* const mdstSpectrum1, int32_t* const mdstSpectrum2,
                                                SfbGroupData&  groupingData1, SfbGroupData&  groupingData2,
                                                const TnsData&   filterData1, const TnsData&   filterData2,
                                                const uint8_t    numSwbFrame, uint8_t* const sfbStereoData,
                                                const uint8_t    bitRateMode, const bool    useFullFrameMS,
                                                const bool    reversePredDir, const uint8_t realOnlyOffset,
                                                uint32_t* const sfbStepSize1, uint32_t* const sfbStepSize2)
{
  const bool applyPredSte = (sfbStereoData != nullptr); // use real-valued predictive stereo
  const SfbGroupData& grp = groupingData1;
  const bool  eightShorts = (grp.numWindowGroups > 1);
  const uint8_t maxSfbSte = (eightShorts ? __min (numSwbFrame, __max (grp.sfbsPerGroup, groupingData2.sfbsPerGroup) + 1) : numSwbFrame);
  const bool  perCorrData = ((bitRateMode <= 5) && !eightShorts); // perceptual correlation?
  const bool  quantDither = ((bitRateMode >= 4) && !eightShorts); // quantization dithering?
  bool       alterPredDir = (applyPredSte && reversePredDir); // predict mid from side band?
  uint32_t rmsSfbL[2] = {0, 0}, rmsSfbR[2] = {0, 0};
  uint32_t  numSfbPredSte = 0; // counter
  uint16_t  numSfbNoMsSte = 0, idxSfbNoMsSte = 0, nNoMS = 0;
  uint32_t rms1NoMsSte[2] = {0, 0}, rms2NoMsSte[2] = {0, 0};
  uint32_t rmsMNoMsSte[2] = {0, 0}, rmsSNoMsSte[2] = {0, 0};
  uint8_t  dataNoMsSte[2] = {0, 0};
  bool nonZeroPredNoMsSte = false;

  if ((mdctSpectrum1 == nullptr) || (mdctSpectrum2 == nullptr) || (numSwbFrame < maxSfbSte) || (grp.numWindowGroups != groupingData2.numWindowGroups) ||
      (sfbStepSize1  == nullptr) || (sfbStepSize2  == nullptr) || (numSwbFrame < MIN_NUM_SWB_SHORT) || (numSwbFrame > MAX_NUM_SWB_LONG))
  {
    return 1;  // invalid arguments error
  }

  if (applyPredSte && (bitRateMode > 5) && !eightShorts && !reversePredDir) // pred_dir test
  {
    uint64_t sumRealM = 0, sumRealS = 0;

    for (uint16_t s = grp.sfbOffsets[numSwbFrame] - 1; s > 0; s--)
    {
      sumRealM += abs (mdctSpectrum1[s] + mdctSpectrum2[s]); // i.e., 2*mid,
      sumRealS += abs (mdctSpectrum1[s] - mdctSpectrum2[s]); // i.e., 2*side
    }
    alterPredDir = (sumRealS * 2 > sumRealM * 3);
  }

  if (applyPredSte && perCorrData) memcpy (m_stereoCorrValue, sfbStereoData, (grp.sfbOffsets[numSwbFrame] >> SA_BW_SHIFT) * sizeof (uint8_t));

  if ((bitRateMode >= 4) && eightShorts) // reset quantizer dither memory in short transform
  {
    for (uint16_t sfb = 0; sfb <= MAX_NUM_SWB_LONG / 2; sfb++)
    {
      m_randomIntMemRe[sfb] = (1 << 30);
      m_randomIntMemIm[sfb] = (1 << 30);
    }
  }

  for (uint16_t n = 0, gr = 0; gr < grp.numWindowGroups; gr++)
  {
    const uint16_t grOffset = numSwbFrame * gr;
    const uint8_t grpLength = grp.windowGroupLength[gr];
    const bool realOnlyCalc = ((grpLength == 1) && (filterData1.numFilters[n] > 0 || filterData2.numFilters[n] > 0)) || !mdstSpectrum1 || !mdstSpectrum2;
    const uint16_t*  grpOff = &grp.sfbOffsets[grOffset];
    uint32_t* const grpRms1 = &groupingData1.sfbRmsValues[grOffset];
    uint32_t* const grpRms2 = &groupingData2.sfbRmsValues[grOffset];
    uint32_t* grpStepSizes1 = &sfbStepSize1[grOffset];
    uint32_t* grpStepSizes2 = &sfbStepSize2[grOffset];
    int32_t  b = 0, prevReM = 0, prevReS = 0;

    if (realOnlyCalc) // preparation for first magnitude value
    {
      const uint16_t sIndex = grpOff[realOnlyOffset] + (grpOff[realOnlyOffset] > 0 ? -1 : 1);

      prevReM = getMidSample (mdctSpectrum1[sIndex], mdctSpectrum2[sIndex]);
      prevReS = getSideSample(mdctSpectrum1[sIndex], mdctSpectrum2[sIndex]);
    }

    for (uint16_t sfb = 0; sfb < maxSfbSte; sfb++)
    {
      const int32_t  sfbIsOdd = sfb & 1;
      const uint16_t sfbStart = grpOff[sfb];
      const uint16_t sfbWidth = grpOff[sfb + 1] - sfbStart;
      int32_t* sfbMdct1 = &mdctSpectrum1[sfbStart];
      int32_t* sfbMdct2 = &mdctSpectrum2[sfbStart];
      int32_t* sfbMdst1 = &mdstSpectrum1[sfbStart];
      int32_t* sfbMdst2 = &mdstSpectrum2[sfbStart];
      uint64_t sumAbsValM = 0, sumAbsValS = 0;
      double   sfbTempVar;

      if ((sfbIsOdd == 0) && !useFullFrameMS) // save L/R data
      {
        const uint16_t cpyWidth = (grpOff[__min (maxSfbSte, sfb + 2)] - sfbStart) * sizeof (int32_t);

        memcpy (m_originBandMdct1, sfbMdct1, cpyWidth);
        memcpy (m_originBandMdct2, sfbMdct2, cpyWidth);
        memcpy (m_originBandMdst1, sfbMdst1, cpyWidth);
        memcpy (m_originBandMdst2, sfbMdst2, cpyWidth);
      }

      if (realOnlyCalc && (sfb >= realOnlyOffset)) // real-valued data, only MDCTs available
      {
        const int32_t* sfbNext1 = &sfbMdct1[1];
        const int32_t* sfbNext2 = &sfbMdct2[1];

        for (uint16_t s = sfbWidth - (sfb + 1 == numSwbFrame ? 1 : 0); s > 0; s--)
        {
          const int32_t dmixReM = getMidSample (*sfbMdct1, *sfbMdct2);
          const int32_t dmixReS = getSideSample(*sfbMdct1, *sfbMdct2);
          // TODO: improve the following lines since the calculation is partially redundant!
          const int32_t dmixImM = int32_t ((getMidSample (*sfbNext1, *sfbNext2) - (int64_t) prevReM) >> 1); // estimate, see also
          const int32_t dmixImS = int32_t ((getSideSample(*sfbNext1, *sfbNext2) - (int64_t) prevReS) >> 1); // getMeanAbsValues()

          sumAbsValM += complexAbsMS (dmixReM, dmixImM);
          sumAbsValS += complexAbsMS (dmixReS, dmixImS);

          *(sfbMdct1++) = dmixReM;
          *(sfbMdct2++) = dmixReS;
          *(sfbMdst1++) = dmixImM;
          *(sfbMdst2++) = dmixImS;
          sfbNext1++; prevReM = dmixReM;
          sfbNext2++; prevReS = dmixReS;
        }
        if (sfb + 1 == numSwbFrame) // process the last sample
        {
          const int32_t dmixReM = getMidSample (*sfbMdct1, *sfbMdct2);
          const int32_t dmixReS = getSideSample(*sfbMdct1, *sfbMdct2);

          sumAbsValM += abs (dmixReM);
          sumAbsValS += abs (dmixReS);

          *sfbMdct1 = dmixReM;
          *sfbMdct2 = dmixReS;
          *sfbMdst1 = 0;
          *sfbMdst2 = 0;
        }
      }
      else // complex data, both MDCTs and MDSTs are available
      {
        for (uint16_t s = sfbWidth; s > 0; s--)
        {
          const int32_t dmixReM = getMidSample (*sfbMdct1, *sfbMdct2);
          const int32_t dmixReS = getSideSample(*sfbMdct1, *sfbMdct2);
          const int32_t dmixImM = getMidSample (*sfbMdst1, *sfbMdst2);
          const int32_t dmixImS = getSideSample(*sfbMdst1, *sfbMdst2);

          sumAbsValM += complexAbsMS (dmixReM, dmixImM);
          sumAbsValS += complexAbsMS (dmixReS, dmixImS);

          *(sfbMdct1++) = dmixReM;
          *(sfbMdct2++) = dmixReS;
          *(sfbMdst1++) = dmixImM;
          *(sfbMdst2++) = dmixImS;
        }
      } // realOnlyCalc && sfb >= realOnlyOffset

      rmsSfbL[sfbIsOdd] = grpRms1[sfb];
      rmsSfbR[sfbIsOdd] = grpRms2[sfb];
      // average spectral sample magnitude across current band
      grpRms1[sfb] = uint32_t ((sumAbsValM + (sfbWidth >> 1)) / sfbWidth);
      grpRms2[sfb] = uint32_t ((sumAbsValS + (sfbWidth >> 1)) / sfbWidth);

      if (applyPredSte) sfbStereoData[sfb + grOffset] = 16; // initialize alpha_q_.. to zero

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
          const int32_t* mdstA = (alterPredDir ? &mdstSpectrum2[offEv] : &mdstSpectrum1[offEv]);
          const int32_t* mdstB = (alterPredDir ? &mdstSpectrum1[offEv] : &mdstSpectrum2[offEv]);
          int64_t sumPrdReAReB = 0, sumPrdImAReB = 0, sumPrdReAReA = SP_EPS; // to stabilize
          double d, alphaLimit = 1.5; // max alpha_q magnitude

          for (uint16_t s = width; s > 0; s--, mdctA++, mdctB++, mdstA++, mdstB++)
          {
            const int64_t prdReAReA = ((int64_t) *mdctA * (int64_t) *mdctA + SA_BW) >> (SA_BW_SHIFT + 1);
            const int64_t prdImAImA = ((int64_t) *mdstA * (int64_t) *mdstA + SA_BW) >> (SA_BW_SHIFT + 1);

            sumPrdReAReB += ((int64_t) *mdctA * (int64_t) *mdctB + SA_BW) >> (SA_BW_SHIFT + 1);
            sumPrdReAReA += prdReAReA;
            sumPrdImAReB += ((int64_t) *mdstA * (int64_t) *mdctB + SA_BW) >> (SA_BW_SHIFT + 1);
            // add complex conjugate part, increases stability
            sumPrdReAReB += ((int64_t) *mdstA * (int64_t) *mdstB + SA_BW) >> (SA_BW_SHIFT + 1);
            sumPrdReAReA += prdImAImA;
            sumPrdImAReB -= ((int64_t) *mdctA * (int64_t) *mdstB + SA_BW) >> (SA_BW_SHIFT + 1);
          }
          for (b = sfbIsOdd; b >= 0; b--) // limit alpha_q to prevent residual RMS increases
          {
            const int idx = sfbEv + b;

            d = (alterPredDir ? (double) grpRms1[idx] / __max (SP_EPS, grpRms2[idx]) : (double) grpRms2[idx] / __max (SP_EPS, grpRms1[idx]));
            if (alphaLimit > d) alphaLimit = d;
          }
          sfbTempVar = CLIP_PM ((double) sumPrdReAReB / (double) sumPrdReAReA, alphaLimit);

          b = __max (512, 524 - int32_t (abs (10.0 * sfbTempVar))); // rounding optimization
          if (quantDither)
          {
            const int32_t r = (int32_t) m_randomInt32 ();
            const double dr = 10.0 * sfbTempVar + (r - m_randomIntMemRe[sfbEv >> 1]) * SP_DIV;

            b = int32_t (dr + b * (dr < 0.0 ? -0.0009765625 : 0.0009765625));
            m_randomIntMemRe[sfbEv >> 1] = r;
          }
          else b = int32_t (10.0 * sfbTempVar + b * (sfbTempVar < 0.0 ? -0.0009765625 : 0.0009765625));

          sfbStereoData[sfbEv + grOffset] = uint8_t (b + 16); // save SFB's final alpha_q_re
          alphaLimit = CLIP_PM ((double) sumPrdImAReB / (double) sumPrdReAReA, alphaLimit);

          b = __max (512, 524 - int32_t (abs (10.0 * alphaLimit))); // rounding optimization
          if (quantDither)
          {
            const int32_t r = (int32_t) m_randomInt32 ();
            const double dr = 10.0 * alphaLimit + (r - m_randomIntMemIm[sfbEv >> 1]) * SP_DIV;

            b = int32_t (dr + b * (dr < 0.0 ? -0.0009765625 : 0.0009765625));
            m_randomIntMemIm[sfbEv >> 1] = r;
          }
          else b = int32_t (10.0 * alphaLimit + b * (alphaLimit < 0.0 ? -0.0009765625 : 0.0009765625));

          if (sfbEv + 1 < numSwbFrame)
          sfbStereoData[sfbEv + 1 + grOffset] = uint8_t (b + 16); // save initial alpha_q_im

          if (perCorrData && ((offEv & (SA_BW - 1)) == 0) && ((width & (SA_BW - 1)) == 0))
          {
            const uint8_t* const perCorr = &m_stereoCorrValue[offEv >> SA_BW_SHIFT];

            // perceptual correlation data available from previous call to stereoSigAnalysis
            b = (width == SA_BW ? perCorr[0] : ((int32_t) perCorr[0] + (int32_t) perCorr[1] + 1) >> 1);
          }
          else b = UCHAR_MAX; // previous correlation data unavailable, assume maximum value

          if ((b > SCHAR_MAX && sfbStereoData[sfbEv + grOffset] != 16) || // if perceptually
              (2 <= abs ( (int) sfbStereoData[sfbEv + grOffset] - 16))) // significant pred.
          {
            nonZeroPredCoef = true;
          }
          sfbTempVar *= sfbTempVar;  // account for residual RMS reduction due to prediction
#if !BA_MORE_CBR
          if (bitRateMode > 0) sfbTempVar += alphaLimit * alphaLimit;  // including alpha_im
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

        if (!useFullFrameMS) // test M/S compaction gain, revert to L/R if it's insufficient
        {
          const uint64_t bandSum1 = (sfbIsOdd > 0 ? (uint64_t) grpRms1[sfbEv] + (uint64_t) grpRms1[sfbEv + 1] : grpRms1[sfbEv]);
          const uint64_t bandSum2 = (sfbIsOdd > 0 ? (uint64_t) grpRms2[sfbEv] + (uint64_t) grpRms2[sfbEv + 1] : grpRms2[sfbEv]);
          const uint64_t bandSumL = (sfbIsOdd > 0 ? (uint64_t) rmsSfbL[0] + (uint64_t) rmsSfbL[1] : rmsSfbL[0]) >> 1;
          const uint64_t bandSumR = (sfbIsOdd > 0 ? (uint64_t) rmsSfbR[0] + (uint64_t) rmsSfbR[1] : rmsSfbR[0]) >> 1;
          const uint64_t bandSumM = (applyPredSte ? (uint64_t) rmsSfbM[0] + (uint64_t) rmsSfbM[1] : bandSum1) >> 1;
          const uint64_t bandSumS = (applyPredSte ? (uint64_t) rmsSfbS[0] + (uint64_t) rmsSfbS[1] : bandSum2) >> 1;

          if ((__min (bandSumM, bandSumS) * __max (bandSumL, bandSumR) >= __min (bandSumL, bandSumR) * __max (bandSumM, bandSumS)) ||
              (nonZeroPredCoef && (abs ( (int) sfbStereoData[sfbEv + grOffset] - 16) >= 10)))
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

              if (numSfbNoMsSte == 0)
              {
                rms1NoMsSte[b] = grpRms1[idx];
                rms2NoMsSte[b] = grpRms2[idx];
              }
              grpRms1[idx] = rmsSfbL[b]; // recover left/right band energies
              grpRms2[idx] = rmsSfbR[b];
              if (applyPredSte)
              {
                if (numSfbNoMsSte == 0)
                {
                  rmsMNoMsSte[b] = rmsSfbM[b];
                  rmsSNoMsSte[b] = rmsSfbS[b];
                  dataNoMsSte[b] = sfbStereoData[idx + grOffset];
                }
                sfbStereoData[idx + grOffset] = 0;  // set ms_used flag to 0
              }
            }
            numSfbNoMsSte++;  nNoMS = n;
            idxSfbNoMsSte = sfbEv + grOffset;
            nonZeroPredNoMsSte = nonZeroPredCoef;

            continue; // M/S is not used
          }
        } // if !useFullFrameMS

        for (b = sfbIsOdd; b >= 0; b--) setStepSizesMS (rmsSfbL, rmsSfbR, rmsSfbM, rmsSfbS, grpRms1, grpRms2,
                                                        grpStepSizes1, grpStepSizes2, sfbEv, (uint16_t) b, applyPredSte);
        if (nonZeroPredCoef) numSfbPredSte++; // if perceptually significant prediction band
      } // if pair completed
    }
    if (grpLength == 1) n++;
  } // for gr

  if (numSfbNoMsSte == 1) // upgrade single L/R to M/S band to reduce M/S signaling overhead
  {
    const uint16_t   grNoMS = idxSfbNoMsSte / numSwbFrame;
    const uint16_t  offNoMS = numSwbFrame * grNoMS;
    const uint16_t  sfbNoMS = idxSfbNoMsSte - offNoMS;
    const uint8_t grpLength = grp.windowGroupLength[grNoMS];
    const bool realOnlyCalc = ((grpLength == 1) && (filterData1.numFilters[nNoMS] > 0 || filterData2.numFilters[nNoMS] > 0)) || !mdstSpectrum1 || !mdstSpectrum2;
    const uint16_t*  grpOff = &grp.sfbOffsets[offNoMS];
    uint32_t* const grpRms1 = &groupingData1.sfbRmsValues[offNoMS];
    uint32_t* const grpRms2 = &groupingData2.sfbRmsValues[offNoMS];

    for (int32_t b = (sfbNoMS + 1 < maxSfbSte ? 1 : 0); b >= 0; b--)
    {
      const int idx = sfbNoMS + b; // sfbNoMS = even SFB index
      const uint16_t sfbStart = grpOff[idx];
      const uint16_t sfbWidth = grpOff[idx + 1] - sfbStart;
      int32_t* sfbMdct1 = &mdctSpectrum1[sfbStart];
      int32_t* sfbMdct2 = &mdctSpectrum2[sfbStart];

      rmsSfbL[b] = grpRms1[idx];  // save left/right band RMSs
      rmsSfbR[b] = grpRms2[idx];
      grpRms1[idx] = rms1NoMsSte[b];  // recover M/S band RMSs
      grpRms2[idx] = rms2NoMsSte[b];
      if (applyPredSte) sfbStereoData[idx + offNoMS] = dataNoMsSte[b];

      if (realOnlyCalc && (idx >= realOnlyOffset)) // real-valued data, only MDCTs available
      {
        for (uint16_t s = sfbWidth; s > 0; s--)
        {
          const int32_t dmixReM = getMidSample (*sfbMdct1, *sfbMdct2);
          const int32_t dmixReS = getSideSample(*sfbMdct1, *sfbMdct2);

          *(sfbMdct1++) = dmixReM;
          *(sfbMdct2++) = dmixReS;
        }
      }
      else // complex data, both MDCTs and MDSTs are available
      {
        int32_t* sfbMdst1 = &mdstSpectrum1[sfbStart];
        int32_t* sfbMdst2 = &mdstSpectrum2[sfbStart];

        for (uint16_t s = sfbWidth; s > 0; s--)
        {
          const int32_t dmixReM = getMidSample (*sfbMdct1, *sfbMdct2);
          const int32_t dmixReS = getSideSample(*sfbMdct1, *sfbMdct2);
          const int32_t dmixImM = getMidSample (*sfbMdst1, *sfbMdst2);
          const int32_t dmixImS = getSideSample(*sfbMdst1, *sfbMdst2);

          *(sfbMdct1++) = dmixReM;
          *(sfbMdct2++) = dmixReS;
          *(sfbMdst1++) = dmixImM;
          *(sfbMdst2++) = dmixImS;
        }
      } // realOnlyCalc && idx >= realOnlyOffset

      setStepSizesMS (rmsSfbL, rmsSfbR, rmsMNoMsSte, rmsSNoMsSte, grpRms1, grpRms2,
                      &sfbStepSize1[offNoMS], &sfbStepSize2[offNoMS], sfbNoMS, (uint16_t) b, applyPredSte);
    }
    if (nonZeroPredNoMsSte) numSfbPredSte++; // was perceptually significant prediction band
  } // if numSfbNoMsSte == 1

  if (numSfbPredSte == 0) // discard prediction coefficients and stay with legacy M/S stereo
  {
    if (applyPredSte)
    {
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
  }
  else // at least one "significant" prediction band, apply prediction and update RMS values
  {
    for (uint16_t n = 0, gr = 0; gr < grp.numWindowGroups; gr++)
    {
      const uint16_t grOffset = numSwbFrame * gr;
      const uint8_t grpLength = grp.windowGroupLength[gr];
      const bool realOnlyCalc = ((grpLength == 1) && (filterData1.numFilters[n] > 0 || filterData2.numFilters[n] > 0)) || !mdstSpectrum1 || !mdstSpectrum2;
      const uint16_t*  grpOff = &grp.sfbOffsets[grOffset];
      uint32_t* const grpRms1 = &groupingData1.sfbRmsValues[grOffset];
      uint32_t* const grpRms2 = &groupingData2.sfbRmsValues[grOffset];
      uint8_t* const grpSData = &sfbStereoData[grOffset];
      int32_t prevResi = 0;

      if (realOnlyCalc) // preparation of res. magnitude value
      {
        const int64_t alphaRe = (grpSData[realOnlyOffset & 0xFE] > 0 ? (int) grpSData[realOnlyOffset & 0xFE] - 16 : 0) * SP_0_DOT_1_16BIT;
        const uint16_t sIndex = grpOff[realOnlyOffset] + (grpOff[realOnlyOffset] > 0 ? -1 : 1);

        prevResi = (alterPredDir ? getResiSample (mdctSpectrum1[sIndex], mdctSpectrum2[sIndex], alphaRe)
                                 : getResiSample (mdctSpectrum2[sIndex], mdctSpectrum1[sIndex], alphaRe));
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
          if (realOnlyCalc && (sfb >= realOnlyOffset)) // update previous residual MDCT data
          {
            sfbMdctR += sfbWidth - 1;
            prevResi = (grpSData[sfbEv] > 0 ? *sfbMdctR : int32_t (((int64_t) sfbMdctD[sfbWidth - 1] +
                                          (alterPredDir ? 1 : -1) * (int64_t) *sfbMdctR + 1) >> 1));
          }
          continue; // nothing more to do, i.e., no prediction
        }

        if (realOnlyCalc && (sfb >= realOnlyOffset))  // real-valued, only MDCT is available
        {
          const int32_t* sfbNextD = &sfbMdctD[1];
          const int32_t* sfbNextR = &sfbMdctR[1];

          for (uint16_t s = sfbWidth - (sfb + 1 == numSwbFrame ? 1 : 0); s > 0; s--)
          {
            const int32_t  resiRe = getResiSample (*sfbMdctR, *sfbMdctD, alphaRe);
            // TODO: improve the following line since the calculation is partially redundant
            const int32_t  resiIm = int32_t ((getResiSample (*sfbNextR, *sfbNextD, alphaRe) - (int64_t) prevResi) >> 1);

            sumAbsValR += complexAbsMS (resiRe, resiIm);

            sfbMdctD++;
            *(sfbMdctR++) = resiRe;
            sfbNextD++;
            sfbNextR++; prevResi = resiRe;
          }
          if (sfb + 1 == numSwbFrame)  // process final sample
          {
            const int32_t  resiRe = getResiSample (*sfbMdctR, *sfbMdctD, alphaRe);

            sumAbsValR += abs (resiRe);

            *sfbMdctR = resiRe;
          }
        }
        else // complex data, both MDCT and MDST are available
        {
          int32_t* sfbMdstD = (alterPredDir ? &mdstSpectrum2[sfbStart] : &mdstSpectrum1[sfbStart]);
          int32_t* sfbMdstR = (alterPredDir ? &mdstSpectrum1[sfbStart] : &mdstSpectrum2[sfbStart]);

          for (uint16_t s = sfbWidth; s > 0; s--)
          {
            const int32_t  resiRe = getResiSample (*sfbMdctR, *sfbMdctD, alphaRe);
            const int32_t  resiIm = getResiSample (*sfbMdstR, *sfbMdstD, alphaRe);

            sumAbsValR += complexAbsMS (resiRe, resiIm);

            sfbMdctD++;
            *(sfbMdctR++) = resiRe;
            sfbMdstD++;
            *(sfbMdstR++) = resiIm;
          }
        } // realOnlyCalc && sfb >= realOnlyOffset

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

          if (grpSData[sfb & 0xFFFE] == 0) continue; // no M/S

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
      if (grpLength == 1) n++;
    } // for gr

    numSfbPredSte = (applyPredSte && (alterPredDir != reversePredDir) ? 4 /*pred_dir=1*/ : 2);
  }

  return numSfbPredSte; // no error
}
