/* bitAllocation.cpp - source file for class needed for psychoacoustic bit-allocation
 * written by C. R. Helmrich, last modified in 2020 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under a Modified BSD-Style License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2020 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#include "exhaleLibPch.h"
#include "bitAllocation.h"

// static helper functions
static inline uint32_t jndModel (const uint32_t val, const uint32_t mean,
                                 const unsigned expTimes512, const unsigned mulTimes512)
{
  const double exp = (double) expTimes512 / 512.0; // exponent
  const double mul = (double) mulTimes512 / 512.0; // a factor
  const double res = pow (mul * (double) val, exp) * pow ((double) mean, 1.0 - exp);

  return uint32_t (__min ((double) UINT_MAX, res + 0.5));
}

static void jndPowerLawAndPeakSmoothing (uint32_t* const  stepSizes, const unsigned nStepSizes,
                                         const uint32_t avgStepSize, const uint8_t sfm, const uint8_t tfm)
{
  const unsigned  expTimes512 = 512u - sfm; // 1.0 - sfm / 2.0
  const unsigned  mulTimes512 = __min (expTimes512, 512u - tfm);
  uint32_t         stepSizeM3 = 0, stepSizeM2 = 0, stepSizeM1 = BA_EPS; // hearing threshold around zero Hz
  unsigned b;

  for (b = 0; b < __min (2, nStepSizes); b++)
  {
    stepSizeM3 = stepSizeM2;
    stepSizeM2 = stepSizeM1;
    stepSizeM1 = stepSizes[b] = jndModel (stepSizes[b], avgStepSize, expTimes512, mulTimes512);
  }
  stepSizes[0] = __min (stepSizeM1, stepSizes[0]); // `- becomes --
  for (/*b*/; b < nStepSizes; b++)
  {
    const uint32_t stepSizeB = jndModel (stepSizes[b], avgStepSize, expTimes512, mulTimes512);

    if ((stepSizeM3 <= stepSizeM2) && (stepSizeM3 <= stepSizeM1) && (stepSizeB <= stepSizeM2) && (stepSizeB <= stepSizeM1))
    {
      const uint32_t maxM3M0 = __max (stepSizeM3, stepSizeB); // smoothen local spectral peak of _´`- shape

      stepSizes[b - 2] = __min (maxM3M0, stepSizes[b - 2]); // _-`-
      stepSizes[b - 1] = __min (maxM3M0, stepSizes[b - 1]); // _---
    }
    stepSizeM3 = stepSizeM2;
    stepSizeM2 = stepSizeM1;
    stepSizeM1 = (stepSizes[b] = stepSizeB); // modified step-size may be smoothened in next loop iteration
  }
}

// constructor
BitAllocator::BitAllocator ()
{
  for (unsigned ch = 0; ch < USAC_MAX_NUM_CHANNELS; ch++)
  {
    m_avgStepSize[ch] = 0;
    m_avgSpecFlat[ch] = 0;
    m_avgTempFlat[ch] = 0;
  }
}

// public functions
void BitAllocator::getChAverageSpecFlat (uint8_t meanSpecFlatInCh[USAC_MAX_NUM_CHANNELS], const unsigned nChannels)
{
  if ((meanSpecFlatInCh == nullptr) || (nChannels > USAC_MAX_NUM_CHANNELS))
  {
    return;
  }
  memcpy (meanSpecFlatInCh, m_avgSpecFlat, nChannels * sizeof (uint8_t));
}

void BitAllocator::getChAverageTempFlat (uint8_t meanTempFlatInCh[USAC_MAX_NUM_CHANNELS], const unsigned nChannels)
{
  if ((meanTempFlatInCh == nullptr) || (nChannels > USAC_MAX_NUM_CHANNELS))
  {
    return;
  }
  memcpy (meanTempFlatInCh, m_avgTempFlat, nChannels * sizeof (uint8_t));
}

uint8_t BitAllocator::getScaleFac (const uint32_t sfbStepSize, const int32_t* const sfbSignal, const uint8_t sfbWidth,
                                   const uint32_t sfbRmsValue)
{
  uint8_t sf;
  uint32_t u;
#if !RESTRICT_TO_AAC
  uint64_t meanSpecLoudness = 0;
  double d;
#endif

  if ((sfbSignal == nullptr) || (sfbWidth == 0) || (sfbRmsValue < 46))
  {
    return 0; // use lowest scale factor
  }
#if RESTRICT_TO_AAC
  u = 0;
  for (sf = 0; sf < sfbWidth; sf++)
  {
    u += uint32_t (0.5 + sqrt (abs ((double) sfbSignal[sf])));
  }
  u = uint32_t ((u * 16384ui64 + (sfbWidth >> 1)) / sfbWidth);
  u = uint32_t (0.5 + sqrt ((double) u) * 128.0);

  if (u < 42567) return 0;

  u = uint32_t ((sfbStepSize * 42567ui64 + (u >> 1)) / u);
  sf = (u > 1 ? uint8_t (0.5 + 17.7169498394 * log10 ((double) u)) : 4);
#else
  for (sf = 0; sf < sfbWidth; sf++) // simple, low-complexity derivation method for USAC's arithmetic coder
  {
    const int64_t temp = ((int64_t) sfbSignal[sf] + 8) >> 4; // avoid overflow

    meanSpecLoudness += temp * temp;
  }
  meanSpecLoudness = uint64_t (0.5 + pow (256.0 * (double) meanSpecLoudness / sfbWidth, 0.25));

  u = uint32_t (0.5 + pow ((double) sfbRmsValue, 0.75) * 256.0);    // u = 2^8 * (sfbRmsValue^0.75)
  u = uint32_t ((meanSpecLoudness * sfbStepSize * 665 + (u >> 1)) / u); // u = sqrt(6.75) * m*thr/u
  d =  (u > 1 ? log10 ((double) u) : 0.25);

  u = uint32_t (0.5 + pow ((double) sfbRmsValue, 0.25) * 16384.0); // u = 2^14 * (sfbRmsValue^0.25)
  u = uint32_t (((uint64_t) sfbStepSize * 42567 + (u >> 1)) / u);         // u = sqrt(6.75) * thr/u
  d += (u > 1 ? log10 ((double) u) : 0.25);

  sf = uint8_t (0.5 + 8.8584749197 * d);  // sf = (8/3) * log2(u1*u2) = (8/3) * (log2(u1)+log2(u2))
#endif

  return __min (SCHAR_MAX, sf);
}

unsigned BitAllocator::initSfbStepSizes (const SfbGroupData* const groupData[USAC_MAX_NUM_CHANNELS], const uint8_t numSwbShort,
                                         const uint32_t specAnaStats[USAC_MAX_NUM_CHANNELS],
                                         const uint32_t tempAnaStats[USAC_MAX_NUM_CHANNELS],
                                         const unsigned nChannels, const unsigned samplingRate, uint32_t* const sfbStepSizes,
                                         const unsigned lfeChannelIndex /*= USAC_MAX_NUM_CHANNELS*/, const bool tnsDisabled /*= false*/)
{
  // equal-loudness weighting based on data from: K. Kurakata, T. Mizunami, and K. Matsushita, "Percentiles
  // of Normal Hearing-Threshold Distribution Under Free-Field Listening Conditions in Numerical Form," Ac.
  // Sci. Tech, vol. 26, no. 5, pp. 447-449, Jan. 2005, https://www.researchgate.net/publication/239433096.
  const unsigned HF/*idx*/= ((123456 - samplingRate) >> 11) + (samplingRate <= 34150 ? 2 : 0); // start SFB
  const unsigned LF/*idx*/= 9;
  const unsigned MF/*idx*/= (samplingRate < 28800 ? HF : __min (HF, 30u));
  const unsigned msShift  = (samplingRate + 36736) >> 15; // TODO: 768 smp
  const unsigned msOffset = 1 << (msShift - 1);
  uint32_t nMeans = 0, sumMeans = 0;

  if ((groupData == nullptr) || (specAnaStats == nullptr) || (tempAnaStats == nullptr) || (sfbStepSizes == nullptr) ||
      (numSwbShort < MIN_NUM_SWB_SHORT) || (numSwbShort > MAX_NUM_SWB_SHORT) || (nChannels > USAC_MAX_NUM_CHANNELS) ||
      (samplingRate < 7350) || (samplingRate > 96000) || (lfeChannelIndex > USAC_MAX_NUM_CHANNELS))
  {
    return 1; // invalid arguments error
  }

  for (unsigned ch = 0; ch < nChannels; ch++)
  {
    const SfbGroupData&   grpData = *groupData[ch];
    const uint32_t maxSfbInCh = grpData.sfbsPerGroup;
    const uint32_t nBandsInCh = grpData.numWindowGroups * maxSfbInCh;
    const uint32_t*   rms = grpData.sfbRmsValues;
    uint32_t*   stepSizes = &sfbStepSizes[ch * numSwbShort * NUM_WINDOW_GROUPS];
// --- apply INTRA-channel simultaneous masking, equal-loudness weighting, and thresholding to SFB RMS data
    uint32_t maskingSlope = 0, b, elw = 58254; // 8/9
    uint32_t rmsEqualLoud = 0;
    uint32_t sumStepSizes = 0;

    m_avgStepSize[ch] = 0;

    b = ((specAnaStats[ch] >> 16) & UCHAR_MAX); // start with squared spec. flatness from spectral analysis
    b = __max (b * b, (tempAnaStats[ch] >> 24) * (tempAnaStats[ch] >> 24)); // ..and from temporal analysis
    m_avgSpecFlat[ch] = uint8_t ((b + (1 << 7)) >> 8); // normalized maximum

    b = ((tempAnaStats[ch] >> 16) & UCHAR_MAX); // now derive squared temp. flatness from temporal analysis
    b = __max (b * b, (specAnaStats[ch] >> 24) * (specAnaStats[ch] >> 24)); // ..and from spectral analysis
    m_avgTempFlat[ch] = uint8_t ((b + (1 << 7)) >> 8); // normalized maximum

    if ((nBandsInCh == 0) || (grpData.numWindowGroups > NUM_WINDOW_GROUPS))
    {
      continue;
    }
    if ((ch == lfeChannelIndex) || (grpData.numWindowGroups != 1)) // LFE, SHORT windows: no masking or ELW
    {
      for (unsigned gr = 0; gr < grpData.numWindowGroups; gr++)
      {
        const uint32_t* gRms = &rms[numSwbShort * gr];
        uint32_t* gStepSizes = &stepSizes[numSwbShort * gr];

        for (b = numSwbShort - 1; b >= maxSfbInCh; b--)
        {
          gStepSizes[b] = 0;
        }
        for (/*b*/; b > 0; b--)
        {
          gStepSizes[b] = __max (gRms[b], BA_EPS);
          sumStepSizes += unsigned (0.5 + sqrt ((double) gStepSizes[b]));
        }
        gStepSizes[0]   = __max (gRms[0], BA_EPS);
        sumStepSizes   += unsigned (0.5 + sqrt ((double) gStepSizes[0]));
      } // for gr

      if (ch != lfeChannelIndex)
      {
// --- SHORT windows: apply perceptual just noticeable difference (JND) model and local band-peak smoothing
        nMeans++;
        m_avgStepSize[ch] = __min (USHRT_MAX, uint32_t ((sumStepSizes + (nBandsInCh >> 1)) / nBandsInCh));
        sumMeans += m_avgStepSize[ch];
        m_avgStepSize[ch] *= m_avgStepSize[ch];

        for (unsigned gr = 0; gr < grpData.numWindowGroups; gr++) // separate peak smoothing for each group
        {
          jndPowerLawAndPeakSmoothing (&stepSizes[numSwbShort * gr], maxSfbInCh, m_avgStepSize[ch], m_avgSpecFlat[ch], 0);
        }
      }
      continue;
    }

    stepSizes[0]   = __max (rms[0], BA_EPS);
    for (b = 1; b < __min (LF, maxSfbInCh); b++) // apply steeper low-frequency simultaneous masking slopes
    {
      maskingSlope = (stepSizes[b - 1] + (msOffset << (9u - b))) >> (msShift + 9u - b);
      stepSizes[b] = __max (rms[b], maskingSlope + BA_EPS);
    }
    for (/*b*/; b < __min (MF, maxSfbInCh); b++) // apply typical mid-frequency simultaneous masking slopes
    {
      maskingSlope = (stepSizes[b - 1] + msOffset) >> msShift;
      stepSizes[b] = __max (rms[b], maskingSlope + BA_EPS);
    }
    if ((samplingRate >= 28800) && (samplingRate <= 64000))
    {
      for (/*b*/; b < __min (HF, maxSfbInCh); b++) // compensate high-frequency slopes for linear SFB width
      {
        maskingSlope = ((uint64_t) stepSizes[b - 1] * (9u + b - MF) + (msOffset << 3u)) >> (msShift + 3u);
        stepSizes[b] = __max (rms[b], maskingSlope + BA_EPS);
      }
      for (/*b = HF region*/; b < maxSfbInCh; b++) // apply extra high-frequency equal-loudness attenuation
      {
        for (unsigned d = b - HF; d > 0; d--)
        {
          elw = (elw * 52430 - SHRT_MIN) >> 16; // elw *= 4/5
        }
        rmsEqualLoud = uint32_t (((uint64_t) rms[b] * elw - SHRT_MIN) >> 16);   // equal loudness weighting
        maskingSlope = ((uint64_t) stepSizes[b - 1] * (9u + b - MF) + (msOffset << 3u)) >> (msShift + 3u);
        stepSizes[b] = __max (rmsEqualLoud, maskingSlope + BA_EPS);
      }
    }
    else // no equal-loudness weighting for low or high rates
    {
      for (/*b = MF region*/; b < maxSfbInCh; b++) // compensate high-frequency slopes for linear SFB width
      {
        maskingSlope = ((uint64_t) stepSizes[b - 1] * (9u + b - MF) + (msOffset << 3u)) >> (msShift + 3u);
        stepSizes[b] = __max (rms[b], maskingSlope + BA_EPS);
      }
    }
    for (b -= 1; b > __min (MF, maxSfbInCh); b--) // complete simultaneous masking by reversing the pattern
    {
      sumStepSizes += unsigned (0.5 + sqrt ((double) stepSizes[b]));
      maskingSlope     = ((uint64_t) stepSizes[b] * (8u + b - MF) + (msOffset << 3u)) >> (msShift + 3u);
      stepSizes[b - 1] = __max (stepSizes[b - 1], maskingSlope);
    }
    for (/*b*/; b > __min (LF, maxSfbInCh); b--)  // typical reversed mid-freq. simultaneous masking slopes
    {
      sumStepSizes += unsigned (0.5 + sqrt ((double) stepSizes[b]));
      maskingSlope     = (stepSizes[b] + msOffset) >> msShift;
      stepSizes[b - 1] = __max (stepSizes[b - 1], maskingSlope);
    }
    for (/*b = min (9, maxSfbInCh)*/; b > 0; b--) // steeper reversed low-freq. simultaneous masking slopes
    {
      sumStepSizes += unsigned (0.5 + sqrt ((double) stepSizes[b]));
      maskingSlope     = (stepSizes[b] + (msOffset << (10u - b))) >> (msShift + 10u - b);
      stepSizes[b - 1] = __max (stepSizes[b - 1], maskingSlope);
    }
    sumStepSizes   += unsigned (0.5 + sqrt ((double) stepSizes[0]));

// --- LONG window: apply perceptual JND model and local band-peak smoothing, undo equal-loudness weighting
    nMeans++;
    m_avgStepSize[ch] = __min (USHRT_MAX, uint32_t ((sumStepSizes + (nBandsInCh >> 1)) / nBandsInCh));
    sumMeans += m_avgStepSize[ch];
    m_avgStepSize[ch] *= m_avgStepSize[ch];

    jndPowerLawAndPeakSmoothing (stepSizes, maxSfbInCh, m_avgStepSize[ch], m_avgSpecFlat[ch], tnsDisabled ? m_avgTempFlat[ch] : 0);

    if ((samplingRate >= 28800) && (samplingRate <= 64000))
    {
      elw = 36; // 36/32 = 9/8
      for (b = HF; b < grpData.sfbsPerGroup; b++)  // undo additional high-freq. equal-loudness attenuation
      {
        for (unsigned d = b - HF; d > 0; d--)
        {
          elw = (16u + elw * 40) >> 5; // inverse elw *= 5/4. NOTE: this may overflow for 64 kHz, that's OK
        }
        if (elw == 138 || elw >= 1024) elw--;
        rmsEqualLoud = uint32_t (__min (UINT_MAX, (16u + (uint64_t) stepSizes[b] * elw) >> 5));
        stepSizes[b] = rmsEqualLoud;
      }
    }
  } // for ch

  if ((nMeans < 2) || (sumMeans <= nMeans * BA_EPS)) // in case of one channel or low-RMS input, we're done
  {
    return 0; // no error
  }

  sumMeans = (sumMeans + (nMeans >> 1)) / nMeans;
  sumMeans *= sumMeans;  // since we've averaged square-roots
#if BA_INTER_CHAN_SIM_MASK
  if (nMeans > 3)
  {
    // TODO: cross-channel simultaneous masking for 4.0 - 7.1
  }
#endif

  for (unsigned ch = 0; ch < nChannels; ch++)
  {
    const SfbGroupData&   grpData = *groupData[ch];
    const uint32_t maxSfbInCh = grpData.sfbsPerGroup;
    const uint32_t nBandsInCh = grpData.numWindowGroups * maxSfbInCh;
    const uint32_t chStepSize = m_avgStepSize[ch];
    uint32_t*   stepSizes = &sfbStepSizes[ch * numSwbShort * NUM_WINDOW_GROUPS];
// --- apply INTER-channel simultaneous masking and JND modeling to calculated INTRA-channel step-size data
    uint64_t mAvgStepSize; // modified and averaged step-size

    if ((nBandsInCh == 0) || (grpData.numWindowGroups > NUM_WINDOW_GROUPS) || (ch == lfeChannelIndex))
    {
      continue;
    }
    mAvgStepSize = jndModel (chStepSize, sumMeans, 7 << 6 /*7/8*/, 512);

    for (unsigned gr = 0; gr < grpData.numWindowGroups; gr++)
    {
      uint32_t* gStepSizes = &stepSizes[numSwbShort * gr];

      for (unsigned b = 0; b < maxSfbInCh; b++)
      {
        gStepSizes[b] = uint32_t (__min (UINT_MAX, (mAvgStepSize * gStepSizes[b] + (chStepSize >> 1)) / chStepSize));
      }
    } // for gr

    m_avgStepSize[ch] = (uint32_t) mAvgStepSize;
  } // for ch

  return 0; // no error
}
