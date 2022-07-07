/* specGapFilling.cpp - source file for class with spectral gap filling coding methods
 * written by C. R. Helmrich, last modified in 2020 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under the exhale Copyright License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2021 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#include "exhaleLibPch.h"
#include "specGapFilling.h"

// ISO/IEC 23003-3, Table 109
static const uint16_t noiseFillingStartOffset[2 /*long/short*/][2 /*768/1024*/] = {{120, 160}, {15, 20}};

// constructor
SpecGapFiller::SpecGapFiller ()
{
  m_1stGapFillSfb = 0;
  memset (m_1stNonZeroSfb, 0, sizeof (m_1stNonZeroSfb));
}

#if SGF_SF_PEAK_SMOOTHING
static const unsigned smallDeltaHuffBitCount[15] = {8, 7, 6, 6, 5, 4, 3, 1, 4, 4, 5, 6, 6, 7, 7};

static inline unsigned huffBitCountEstimate (const int scaleFactorDelta)
{
  if (abs (scaleFactorDelta) < 8) return smallDeltaHuffBitCount[scaleFactorDelta + 7];
  return (abs (scaleFactorDelta) >> 1) + 4;
}
#endif

// public functions
uint8_t SpecGapFiller::getSpecGapFillParams (const SfbQuantizer& sfbQuantizer, const uint8_t* const quantMagn,
                                             const uint8_t numSwbShort, SfbGroupData& grpData /*modified*/,
                                             const unsigned nSamplesInFrame, const unsigned samplingRate,
                                             const unsigned sampRateBitSave, const uint8_t specFlat)
{
  const unsigned* const coeffMagn = sfbQuantizer.getCoeffMagnPtr ();
  const double* const  sfNormFacs = sfbQuantizer.getSfNormTabPtr ();
  const uint16_t       sfbsPerGrp = grpData.sfbsPerGroup;
  const uint16_t       windowNfso = noiseFillingStartOffset[grpData.numWindowGroups == 1 ? 0 : 1][nSamplesInFrame >> 10];
  const bool saveRate = (samplingRate >= sampRateBitSave);
  uint8_t scaleFacLim = 0;
  uint16_t u = 0;
  short diff = 0, s = 0;
  double    magnSum = 0.0;
#if SGF_OPT_SHORT_WIN_CALC
  double minGrpMean = (double) UINT_MAX;
  double sumGrpMean = 0.0; // for shorts
#endif

  if ((coeffMagn == nullptr) || (sfNormFacs == nullptr) || (quantMagn == nullptr) ||
      (numSwbShort < MIN_NUM_SWB_SHORT) || (numSwbShort > MAX_NUM_SWB_SHORT) || (nSamplesInFrame > 1024))
  {
    return 1; // invalid arguments error
  }

// --- determine noise_level as mean of all coeff magnitudes at zero-quantized coeff indices
  m_1stGapFillSfb = 0;
  memset (m_1stNonZeroSfb, -1, sizeof (m_1stNonZeroSfb));

  for (uint16_t gr = 0; gr < grpData.numWindowGroups; gr++)
  {
    const uint16_t*   grpOff = &grpData.sfbOffsets[numSwbShort * gr];
    const uint32_t*   grpRms = &grpData.sfbRmsValues[numSwbShort * gr]; // quant/coder stats
    const uint8_t* grpScFacs = &grpData.scaleFactors[numSwbShort * gr];
    const uint16_t grpLength = grpData.windowGroupLength[gr];
    const uint16_t   grpNfso = grpOff[0] + grpLength * windowNfso;
    const uint16_t  sfbLimit = (grpData.numWindowGroups == 1 ? sfbsPerGrp - (grpOff[sfbsPerGrp] >= nSamplesInFrame ? 1 : 0)
                                                             : __min (sfbsPerGrp, numSwbShort - 1)); // no high frequencies
#if SGF_OPT_SHORT_WIN_CALC
    uint16_t tempNum = u;
    double   tempSum = magnSum;
#endif
    for (uint16_t b = 0; b < sfbLimit; b++)  // determine first gap-fill SFB and noise_level
    {
      const uint16_t sfbStart = grpOff[b];
      const uint16_t sfbWidth = grpOff[b + 1] - sfbStart;
      const unsigned* const sfbMagn = &coeffMagn[sfbStart];
      const uint8_t* sfbQuant = &quantMagn[sfbStart];
      const uint8_t scaleFacB = grpScFacs[b];

      if (sfbStart < grpNfso) // SFBs below noiseFillingStartOffset
      {
        if ((grpRms[b] >> 16) > 0) // the SFB is non-zero quantized
        {
          if (m_1stNonZeroSfb[gr] < 0) m_1stNonZeroSfb[gr] = b;
          if (scaleFacLim < scaleFacB) scaleFacLim = scaleFacB;
        }
      }
      else // sfbStart >= grpNfso, so above noiseFillingStartOffset
      {
        if (m_1stNonZeroSfb[gr] < 0) m_1stNonZeroSfb[gr] = b;
        if (m_1stGapFillSfb == 0)    m_1stGapFillSfb = b;

        if ((grpRms[b] >> 16) > 0) // the SFB is non-zero quantized
        {
          uint64_t sfbMagnSum = 0;

          if (scaleFacLim < scaleFacB) scaleFacLim = scaleFacB;
#if SGF_OPT_SHORT_WIN_CALC
          if (grpLength > 1) // eight-short windows: SFB ungrouping
          {
            const uint32_t* sfbMagnPtr = sfbMagn;
            const uint8_t* sfbQuantPtr = sfbQuant;
            const int swbLength = (sfbWidth * oneTwentyEightOver[grpLength]) >> 7; // sfbWidth / grpLength
            unsigned sfbMagnMin = USHRT_MAX;
            uint16_t uMin = 0;

            for (uint16_t w = 0; w < grpLength; w++)
            {
              unsigned sfbMagnWin = 0;
              uint16_t uWin = 0;

              for (int i = swbLength - 1; i >= 0; i--, sfbMagnPtr++, sfbQuantPtr++)
              {
                if ((*sfbQuantPtr == 0) && (i == 0 || i == swbLength - 1 || *(sfbQuantPtr-1) + *(sfbQuantPtr+1) < 2))
                {
                  sfbMagnWin += *sfbMagnPtr;
                  uWin++;
                }
              }
              if (sfbMagnWin * (uint64_t) uMin < sfbMagnMin * (uint64_t) uWin) // new minimum
              {
                sfbMagnMin = sfbMagnWin;
                uMin = uWin;
              }
            } // for w

            sfbMagnSum += sfbMagnMin * grpLength; // scaled minimum
            u += uMin * grpLength;
          }
          else
#endif
          for (int i = sfbWidth - 1; i >= 0; i--)
          {
            if ((sfbQuant[i] == 0) && (sfbQuant[i - 1] + sfbQuant[i + 1] < 2))
            {
              sfbMagnSum += sfbMagn[i];
              u++;
            }
          }
          magnSum += sfbMagnSum * sfNormFacs[scaleFacB];
        }
      }
    } // for b

    // clip to non-negative value for get function and memset below
    if (m_1stNonZeroSfb[gr] < 0) m_1stNonZeroSfb[gr] = 0;
#if SGF_OPT_SHORT_WIN_CALC
    if ((grpData.numWindowGroups > 1) && (u > tempNum))
    {
      tempSum = (magnSum - tempSum) / double (u - tempNum);
      if (minGrpMean > tempSum) minGrpMean = tempSum;
      sumGrpMean += tempSum;  s++;
    }
#endif
  } // for gr

  // determine quantized noise_level from normalized mean magnitude
  if ((u < 4) || (magnSum * 359.0 < u * 16.0))
  {
    if (sfbsPerGrp <= m_1stGapFillSfb) return 0; // silent, level 0

    magnSum = 1.0;  u = 4; // max. level
  }
#if SGF_OPT_SHORT_WIN_CALC
  if ((s > 1) && (sumGrpMean > 0.0))
  {
    magnSum *= sqrt ((minGrpMean * s) / sumGrpMean);  // Robots fix
    if (magnSum * 64.0 < u * 3.0) // .05
    {
      magnSum = 3.0;  u = 64; // ensure noise_level remains nonzero
    }
  }
  s = 0;
#endif
  u = __min (7 + (specFlat >> 5), uint16_t (14.47118288 + 9.965784285 * log10 (magnSum / (double) u)));
  u = __max (1, u - int (specFlat >> 5)); // SFM-adaptive reduction

  magnSum = pow (2.0, (14 - u) / 3.0); // noiseVal^-1, 23003-3, 7.2
  magnSum *= 1.25 - specFlat * 0.0009765625;

// --- calculate gap-fill scale factors for zero quantized SFBs, then determine noise_offset
  u <<= 5;  // left-shift for bit-stream
  if (scaleFacLim < SGF_LIMIT) scaleFacLim = SGF_LIMIT;

  for (uint16_t gr = 0; gr < grpData.numWindowGroups; gr++)
  {
    const uint16_t*   grpOff = &grpData.sfbOffsets[numSwbShort * gr];
    const uint32_t*   grpRms = &grpData.sfbRmsValues[numSwbShort * gr]; // quant/coder stats
    uint8_t* const grpScFacs = &grpData.scaleFactors[numSwbShort * gr];
#if SGF_SF_PEAK_SMOOTHING
    uint16_t  lastNonZeroSfb = 0;
#endif
    for (uint16_t b = m_1stGapFillSfb; b < sfbsPerGrp; b++)  // get noise-fill scale factors
    {
      if ((grpRms[b] >> 16) == 0)  // the SFB is all-zero quantized
      {
        if (grpScFacs[b] > 0)
        {
          const uint16_t  sfbStart = grpOff[b];
          const int16_t sfbWidthM1 = grpOff[b + 1] - sfbStart - 1;
          const unsigned*  sfbMagn = &coeffMagn[sfbStart];
          unsigned sfbMagnMax = 0;
          uint64_t sfbMagnSum = 0;

          for (int i = sfbWidthM1; i >= 0; i--)
          {
            sfbMagnSum += sfbMagn[i];
            if (sfbMagnMax < sfbMagn[i]) sfbMagnMax = sfbMagn[i];  // sum up without maximum
          }
          grpScFacs[b] = sfbQuantizer.getScaleFacOffset (((sfbMagnSum - sfbMagnMax) * magnSum) / (double) sfbWidthM1);

          if ((samplingRate <= 32000) && (b < m_1stGapFillSfb + 4)) // lower mid-freq. noise
          {
            grpScFacs[b] = __max (0, grpScFacs[b] - int ((m_1stGapFillSfb + 4 - b) * (grpScFacs[b] >> 6)));
          }
          if (grpScFacs[b] > scaleFacLim) grpScFacs[b] = scaleFacLim;
        }
#if SGF_SF_PEAK_SMOOTHING
        // save delta-code bits by smoothing scale factor peaks in zero quantized SFB ranges
        if ((b > m_1stGapFillSfb) && ((grpRms[b - 1] >> 16) == 0) && ((grpRms[b - 2] >> 16) == 0))
        {
          const uint16_t next = grpScFacs[b];
          const uint16_t prev = grpScFacs[b - 2];
          uint8_t&       curr = grpScFacs[b - 1];

          if ((next | prev) && (curr > next) && (curr > prev)) curr = (curr + __max (next, prev)) >> 1;
          else if (saveRate && (curr < next) && (curr < prev)) curr = (curr + __min (next, prev) + 1) >> 1;
        }
#endif
      }
#if SGF_SF_PEAK_SMOOTHING
      else if (saveRate) lastNonZeroSfb = b;
#endif

      if ((b > m_1stGapFillSfb) && (((grpRms[b - 1] >> 16) > 0) ^ ((grpRms[b - 2] >> 16) > 0)))
      {
        diff += (int) grpScFacs[b - 1] - (int) grpScFacs[b - 2]; // sum up transition deltas
        s++;
      }
    } // for b
#if SGF_SF_PEAK_SMOOTHING
    if ((lastNonZeroSfb > 0) && (lastNonZeroSfb + 4 < sfbsPerGrp)) // HF factor line-fitting
    {
      const int32_t start = lastNonZeroSfb + 1;
      const int32_t size  = sfbsPerGrp - start - 1;
      const int32_t xSum  = (size * (size + 1)) >> 1;
      int32_t ySum = 0, a = 0, b = 0, y = 0;
      uint16_t x;

      for (x = start + 1; x < sfbsPerGrp; x++) ySum += grpScFacs[x]; // size * (mean factor)

      for (x = start + 1; x < sfbsPerGrp; x++)
      {
        const int32_t xZ  = size * (x - start) - xSum; // zero-mean
        const int32_t yZ  = size * grpScFacs[x] - ySum;

        a += xZ * xZ;
        b += xZ * yZ;
        y += yZ * yZ;
      }

      if ((a > 0) && (b * b > ((a * y) >> 3)))  // factor smoothing
      {
        unsigned countOld = 0, countNew = 0;

        b = CLIP_PM (((b * (1 << 8)) + (a >> 1)) / a, SHRT_MAX);
        a = ((ySum << 8) - b * xSum + (size >> 1)) / size;

        ySum = grpScFacs[start];
        for (x = start + 1; x < sfbsPerGrp; x++)
        {
          y = CLIP_UCHAR ((a + b * (x - start) - SCHAR_MIN) >> 8);

          countOld += huffBitCountEstimate ((int) grpScFacs[x] - grpScFacs[x - 1]);
          countNew += huffBitCountEstimate (y - ySum);
          ySum = y;
        }
        if (countNew < countOld)
        {
          for (x = start + 1; x < sfbsPerGrp; x++) grpScFacs[x] = CLIP_UCHAR ((a + b * (x - start) - SCHAR_MIN) >> 8);
        }
      }
    }
#endif
  } // for gr

  if (s > 0)
  {
    diff = (diff + (s >> 1)*(diff < 0 ? -1 : 1)) / s; // mean delta
    if (diff < -16) diff = -16;
    else
    if (diff >= 16) diff = 15;
  }
  s = __max (-diff, (short) scaleFacLim - SGF_LIMIT); // limit diff

  for (uint16_t gr = 0; gr < grpData.numWindowGroups; gr++)
  {
    const uint32_t*   grpRms = &grpData.sfbRmsValues[numSwbShort * gr]; // quant/coder stats
    uint8_t* const grpScFacs = &grpData.scaleFactors[numSwbShort * gr];

    for (uint16_t b = m_1stGapFillSfb; b < sfbsPerGrp; b++)  // account for the noise_offset
    {
      if ((grpRms[b] >> 16) == 0)  // the SFB is all-zero quantized
      {
        grpScFacs[b] = (uint8_t) __max (s, grpScFacs[b] - diff);

        if (grpScFacs[b] > scaleFacLim) grpScFacs[b] = scaleFacLim;
      }
    }

    // repeat first significant scale factor downwards to save bits
    memset (grpScFacs, grpScFacs[m_1stNonZeroSfb[gr]], m_1stNonZeroSfb[gr] * sizeof (uint8_t));
  } // for gr

  return CLIP_UCHAR (u | (diff + 16)); // combined level and offset
}
