/* specGapFilling.cpp - source file for class with spectral gap filling coding methods
 * written by C. R. Helmrich, last modified in 2019 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under a Modified BSD-Style License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2020 Christian R. Helmrich, project ecodis. All rights reserved.
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

// public functions
unsigned char SpecGapFiller::getSpecGapFillParams (const SfbQuantizer& sfbQuantizer, const unsigned char* const quantMagn,
                                                   const unsigned char numSwbShort,  SfbGroupData& grpData /*modified*/,
                                                   const unsigned nSamplesInFrame /*= 1024*/)
{
  const unsigned* const coeffMagn = sfbQuantizer.getCoeffMagnPtr ();
  const double* const  sfNormFacs = sfbQuantizer.getSfNormTabPtr ();
  const uint16_t       sfbsPerGrp = grpData.sfbsPerGroup;
  const uint16_t       windowNfso = noiseFillingStartOffset[grpData.numWindowGroups == 1 ? 0 : 1][nSamplesInFrame >> 10];
  unsigned char  scaleFactorLimit = 0;
  uint16_t u = 0;
  short diff = 0, s = 0;
  double  magnSum = 0.0;

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
    const uint16_t* grpOff = &grpData.sfbOffsets[numSwbShort * gr];
    const uint32_t* grpRms = &grpData.sfbRmsValues[numSwbShort * gr]; // quant./coding stats
    const unsigned char* grpScFacs = &grpData.scaleFactors[numSwbShort * gr];
    const uint16_t grpLength = grpData.windowGroupLength[gr];
    const uint16_t   grpNfso = grpOff[0] + grpLength * windowNfso;
    const uint16_t  sfbLimit = (grpData.numWindowGroups == 1 ? sfbsPerGrp - (grpOff[sfbsPerGrp] >= nSamplesInFrame ? 1 : 0)
                                                                   : __min (sfbsPerGrp, numSwbShort - 1)); // no high frequencies
    for (uint16_t b = 0; b < sfbLimit; b++) // find first gap-fill SFB and noise_level
    {
      const uint16_t sfbStart = grpOff[b];
      const uint16_t sfbWidth = grpOff[b + 1] - sfbStart;
      const unsigned* const sfbMagn = &coeffMagn[sfbStart];
      const unsigned char* sfbQuant = &quantMagn[sfbStart];
      const unsigned char      sFac = grpScFacs[b];

      if (sfbStart < grpNfso) // SFBs below noiseFillingStartOffset
      {
        if ((grpRms[b] >> 16) > 0) // the SFB is non-zero quantized
        {
          if (m_1stNonZeroSfb[gr] < 0) m_1stNonZeroSfb[gr] = b;
          if (scaleFactorLimit < sFac) scaleFactorLimit = sFac;
        }
      }
      else // sfbStart >= grpNfso, so above noiseFillingStartOffset
      {
        if (m_1stNonZeroSfb[gr] < 0) m_1stNonZeroSfb[gr] = b;
        if (m_1stGapFillSfb == 0)    m_1stGapFillSfb = b;

        if ((grpRms[b] >> 16) > 0) // the SFB is non-zero quantized
        {
          unsigned sfbMagnSum = 0; // NOTE: may overflow, but unlikely, and 32 bit is faster

          if (scaleFactorLimit < sFac) scaleFactorLimit = sFac;
#if SGF_OPT_SHORT_WIN_CALC
          if (grpLength > 1) // eight-short windows: SFB ungrouping
          {
            const uint32_t*   sfbMagnPtr = sfbMagn;
            const unsigned char* sfbQuantPtr = sfbQuant;
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
          magnSum += sfbMagnSum * sfNormFacs[grpScFacs[b]];
        }
      }
    } // for b

    // clip to non-negative value for get function and memset below
    if (m_1stNonZeroSfb[gr] < 0) m_1stNonZeroSfb[gr] = 0;
  } // for gr

  // determine quantized noise_level from normalized mean magnitude
  if ((u < 4) || (magnSum * 359.0 < u * 16.0))
  {
    if (sfbsPerGrp <= m_1stGapFillSfb) return 0; // silent, level 0

    magnSum = 1.0;  u = 4; // max. level
  }
  u = __min (7, uint16_t (14.47118288 + 9.965784285 * log10 (magnSum / (double) u)));

  magnSum = pow (2.0, (14 - u) / 3.0); // noiseVal^-1, 23003-3, 7.2

// --- calculate gap-fill scale factors for zero quantized SFBs, then determine noise_offset
  u <<= 5;  // left-shift for bit-stream
  if (scaleFactorLimit < SGF_LIMIT) scaleFactorLimit = SGF_LIMIT;

  for (uint16_t gr = 0; gr < grpData.numWindowGroups; gr++)
  {
    const uint16_t* grpOff = &grpData.sfbOffsets[numSwbShort * gr];
    const uint32_t* grpRms = &grpData.sfbRmsValues[numSwbShort * gr]; // quant./coding stats
    unsigned char* const grpScFacs = &grpData.scaleFactors[numSwbShort * gr];

    for (uint16_t b = m_1stGapFillSfb; b < sfbsPerGrp; b++) // calculate scale factors
    {
      if ((grpRms[b] >> 16) == 0)  // the SFB is all-zero quantized
      {
        if (grpScFacs[b] > 0)
        {
          const uint16_t  sfbStart = grpOff[b];
          const int16_t sfbWidthM1 = grpOff[b + 1] - sfbStart - 1;
          const unsigned*  sfbMagn = &coeffMagn[sfbStart];
          unsigned sfbMagnMax = 0;
          unsigned sfbMagnSum = 0; // NOTE: may overflow, but unlikely, and 32 bit is faster

          for (int i = sfbWidthM1; i >= 0; i--)
          {
            sfbMagnSum += sfbMagn[i];
            if (sfbMagnMax < sfbMagn[i]) sfbMagnMax = sfbMagn[i];  // sum up without maximum
          }
          grpScFacs[b] = sfbQuantizer.getScaleFacOffset (((sfbMagnSum - sfbMagnMax) * magnSum) / (double) sfbWidthM1);

          if (grpScFacs[b] > scaleFactorLimit) grpScFacs[b] = scaleFactorLimit;
        }
#if SGF_SF_PEAK_SMOOTHING
        // save delta-code bits by smoothing scale factor peaks in zero quantized SFB ranges
        if ((b >  m_1stGapFillSfb) && ((grpRms[b - 1] >> 16) == 0) && ((grpRms[b - 2] >> 16) == 0) &&
            (grpScFacs[b - 1] > grpScFacs[b]) && (grpScFacs[b - 1] > grpScFacs[b - 2]))
        {
          grpScFacs[b - 1] = (grpScFacs[b - 1] + __max (grpScFacs[b], grpScFacs[b - 2])) >> 1;
        }
#endif
      }

      if ((b > m_1stGapFillSfb) && (((grpRms[b - 1] >> 16) > 0) ^ ((grpRms[b - 2] >> 16) > 0)))
      {
        diff += (int) grpScFacs[b - 1] - (int) grpScFacs[b - 2]; // sum up transition deltas
        s++;
      }
    } // for b
  } // for gr

  if (s > 0)
  {
    diff = (diff + (s >> 1)*(diff < 0 ? -1 : 1)) / s; // mean delta
    if (diff < -16) diff = -16;
    else
    if (diff >= 16) diff = 15;
  }
  s = __max (-diff, (short) scaleFactorLimit - SGF_LIMIT); // limit

  for (uint16_t gr = 0; gr < grpData.numWindowGroups; gr++)
  {
    const uint32_t* grpRms = &grpData.sfbRmsValues[numSwbShort * gr]; // quant./coding stats
    unsigned char* const grpScFacs = &grpData.scaleFactors[numSwbShort * gr];

    for (uint16_t b = m_1stGapFillSfb; b < sfbsPerGrp; b++) // account f. noise_offset
    {
      if ((grpRms[b] >> 16) == 0)  // the SFB is all-zero quantized
      {
        grpScFacs[b] = (unsigned char) __max (s, grpScFacs[b] - diff);

        if (grpScFacs[b] > scaleFactorLimit) grpScFacs[b] = scaleFactorLimit;
      }
    } // for b

    // repeat first significant scale factor downwards to save bits
    memset (grpScFacs, grpScFacs[m_1stNonZeroSfb[gr]], m_1stNonZeroSfb[gr] * sizeof (unsigned char));
  } // for gr

  return CLIP_UCHAR (u | (diff + 16)); // combined level and offset
}
