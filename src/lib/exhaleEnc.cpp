/* exhaleEnc.cpp - source file for class providing Extended HE-AAC encoding capability
 * written by C. R. Helmrich, last modified in 2021 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under the exhale Copyright License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2021 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#include "exhaleLibPch.h"
#include "exhaleEnc.h"

// static helper functions
static double modifiedBesselFunctionOfFirstKind (const double x)
{
  const double xOver2 = x * 0.5;
  double d = 1.0, sum = 1.0;
  int    i = 0;

  do
  {
    const double x2di = xOver2 / double (++i);

    d *= (x2di * x2di);
    sum += d;
  }
  while (d > sum * 1.2e-38); // FLT_MIN

  return sum;
}

static int32_t* initWindowHalfCoeffs (const USAC_WSHP windowShape, const unsigned frameLength)
{
  int32_t* windowBuf = nullptr;
  unsigned u;

  if ((windowBuf = (int32_t*) malloc (frameLength * sizeof (int32_t))) == nullptr)
  {
    return nullptr; // allocation error
  }

  if (windowShape == WINDOW_SINE)
  {
    const double dNorm = 3.141592653589793 / (2.0 * frameLength);
    // MLT sine window half
    for (u = 0; u < frameLength; u++)
    {
      windowBuf[u] = int32_t (sin (dNorm * (u + 0.5)) * WIN_SCALE + 0.5);
    }
  }
  else  // if windowShape == WINDOW_KBD
  {
    const double alpha = 3.141592653589793 * (frameLength > 256 ? 4.0 : 6.0);
    const double dBeta = 1.0 / modifiedBesselFunctionOfFirstKind (alpha /*sqrt (1.0)*/);
    const double dNorm = 4.0 / (2.0 * frameLength);
    const double iScal = double (1u << 30);
    const double dScal = 1.0 / iScal;
    double d, sum = 0.0;
    // create Kaiser-Bessel window half
    for (u = 0; u < frameLength; u++)
    {
      const double du1 = dNorm * u - 1.0;

      d = dBeta * modifiedBesselFunctionOfFirstKind (alpha * sqrt (1.0 - du1 * du1));
      sum += d;
      windowBuf[u] = int32_t (d * iScal + 0.5);
    }
    d = 1.0 / sum; // normalized to sum
    sum = 0.0;
    // KBD window half
    for (u = 0; u < frameLength; u++)
    {
      sum += dScal * windowBuf[u];
      windowBuf[u] = int32_t (sqrt (d * sum /*cumulative sum*/) * WIN_SCALE + 0.5);
    }
  }
  return windowBuf;
}

static uint32_t quantizeSfbWithMinSnr (const unsigned* const coeffMagn, const uint16_t* const sfbOffset, const unsigned b,
                                       const uint8_t groupLength, uint8_t* const quantMagn, char* const arithTuples, const bool nonZeroSnr = false)
{
  const uint16_t sfbStart = sfbOffset[b];
  const uint16_t sfbWidth = sfbOffset[b + 1] - sfbStart;
  const unsigned* sfbMagn = &coeffMagn[sfbStart];
  uint32_t maxIndex = 0, maxLevel = sfbMagn[0];

  for (uint16_t s = sfbWidth - 1; s > 0; s--)
  {
    if (maxLevel < sfbMagn[s])  // find largest-level magn. in SFB
    {
      maxLevel = sfbMagn[s];
      maxIndex = s;
    }
  }
  if (quantMagn != nullptr)  // update quantized sample magnitudes
  {
    memset (&quantMagn[sfbStart], 0, sfbWidth * sizeof (uint8_t));

    if (nonZeroSnr) quantMagn[sfbStart + maxIndex] = 1; // magn. 1
  }

  if (arithTuples != nullptr)  // update entropy coding two-tuples
  {
    const uint16_t swbStart = ((sfbStart - sfbOffset[0]) * oneTwentyEightOver[groupLength]) >> 7;

    memset (&arithTuples[swbStart >> 1], 1, ((sfbWidth * oneTwentyEightOver[groupLength]) >> 8) * sizeof (char));

    if (nonZeroSnr && (groupLength == 1)) // max. two-tuple is 1+1
    {
      arithTuples[(swbStart + maxIndex) >> 1] = 2;
    }
  }

  return maxLevel;
}

// inline helper functions
static inline void applyStereoPreProcessingCplx (int32_t* mdctSample1, int32_t* mdctSample2,
                                                 int32_t* mdstSample1, int32_t* mdstSample2,
                                                 const int64_t factIn, const int64_t factDe, const int64_t sign)
{
  const int32_t  valI1 = *mdstSample1;
  const int32_t  valI2 = *mdstSample2;
  const int32_t  valR1 = *mdctSample1;
  const int32_t  valR2 = *mdctSample2;
  const int64_t  absR1 = abs (valR1);
  const int64_t  absR2 = abs (valR2);
  int64_t dmxI1, dmxR1 = valR1 * factDe + sign * valR2 * factIn; // cross
  int64_t dmxI2, dmxR2 = valR2 * factDe + sign * valR1 * factIn; // -talk
  double n, d;

  if (abs (dmxR1) < absR1 + absR2) // avoid destructive summations
  {
    if (absR1 * factDe < absR2 * factIn)
    {
      dmxR1 = valR2 * factIn - sign * valR1 * factDe;
      dmxI1 = valI2 * factIn - sign * valI1 * factDe;
    }
    else
    {
      dmxR1 = valR1 * factDe - sign * valR2 * factIn;
      dmxI1 = valI1 * factDe - sign * valI2 * factIn;
    }
  }
  else dmxI1 = valI1 * factDe + sign * valI2 * factIn;

  if (abs (dmxR2) < absR1 + absR2) // avoid destructive summations
  {
    if (absR1 * factIn < absR2 * factDe)
    {
      dmxR2 = valR2 * factDe - sign * valR1 * factIn;
      dmxI2 = valI2 * factDe - sign * valI1 * factIn;
    }
    else
    {
      dmxR2 = valR1 * factIn - sign * valR2 * factDe;
      dmxI2 = valI1 * factIn - sign * valI2 * factDe;
    }
  }
  else dmxI2 = valI2 * factDe + sign * valI1 * factIn;

  n = (double) valR1 * (double) valR1 + (double) valI1 * (double) valI1;
  d = (double) dmxR1 * (double) dmxR1 + (double) dmxI1 * (double) dmxI1;
  d = sqrt (n / __max (1.0, d));
  *mdctSample1 = int32_t (dmxR1 * d + (dmxR1 < 0 ? -0.5 : 0.5));
  *mdstSample1 = int32_t (dmxI1 * d + (dmxI1 < 0 ? -0.5 : 0.5));

  n = (double) valR2 * (double) valR2 + (double) valI2 * (double) valI2;
  d = (double) dmxR2 * (double) dmxR2 + (double) dmxI2 * (double) dmxI2;
  d = sqrt (n / __max (1.0, d));
  *mdctSample2 = int32_t (dmxR2 * d + (dmxR2 < 0 ? -0.5 : 0.5));
  *mdstSample2 = int32_t (dmxI2 * d + (dmxI2 < 0 ? -0.5 : 0.5));
}

static inline void applyStereoPreProcessingReal (int32_t* mdctSample1, int32_t* mdctSample2,
                                                 int32_t* prevSample1, int32_t* prevSample2,
                                                 const int64_t factIn, const int64_t factDe, const int64_t sign)
{
  const int64_t  valI1 = (*(mdctSample1 + 1) - (int64_t) *prevSample1) >> 1; // estimate, see also
  const int64_t  valI2 = (*(mdctSample2 + 1) - (int64_t) *prevSample2) >> 1; // getMeanAbsValues()
  const int32_t  valR1 = (*prevSample1 = *mdctSample1);
  const int32_t  valR2 = (*prevSample2 = *mdctSample2);
  const int64_t  absR1 = abs (valR1);
  const int64_t  absR2 = abs (valR2);
  int64_t dmxI1, dmxR1 = valR1 * factDe + sign * valR2 * factIn; // cross
  int64_t dmxI2, dmxR2 = valR2 * factDe + sign * valR1 * factIn; // -talk
  double n, d;

  if (abs (dmxR1) < absR1 + absR2) // avoid destructive summations
  {
    if (absR1 * factDe < absR2 * factIn)
    {
      dmxR1 = valR2 * factIn - sign * valR1 * factDe;
      dmxI1 = valI2 * factIn - sign * valI1 * factDe;
    }
    else
    {
      dmxR1 = valR1 * factDe - sign * valR2 * factIn;
      dmxI1 = valI1 * factDe - sign * valI2 * factIn;
    }
  }
  else dmxI1 = valI1 * factDe + sign * valI2 * factIn;

  if (abs (dmxR2) < absR1 + absR2) // avoid destructive summations
  {
    if (absR1 * factIn < absR2 * factDe)
    {
      dmxR2 = valR2 * factDe - sign * valR1 * factIn;
      dmxI2 = valI2 * factDe - sign * valI1 * factIn;
    }
    else
    {
      dmxR2 = valR1 * factIn - sign * valR2 * factDe;
      dmxI2 = valI1 * factIn - sign * valI2 * factDe;
    }
  }
  else dmxI2 = valI2 * factDe + sign * valI1 * factIn;

  n = (double) valR1 * (double) valR1 + (double) valI1 * (double) valI1;
  d = (double) dmxR1 * (double) dmxR1 + (double) dmxI1 * (double) dmxI1;
  *mdctSample1 = int32_t (dmxR1 * sqrt (n / __max (1.0, d)) + (dmxR1 < 0 ? -0.5 : 0.5));

  n = (double) valR2 * (double) valR2 + (double) valI2 * (double) valI2;
  d = (double) dmxR2 * (double) dmxR2 + (double) dmxI2 * (double) dmxI2;
  *mdctSample2 = int32_t (dmxR2 * sqrt (n / __max (1.0, d)) + (dmxR2 < 0 ? -0.5 : 0.5));
}

static inline void applyTnsCoeff2ChannelSynch (LinearPredictor& predictor, TnsData& tnsData1, TnsData& tnsData2,
                                               const uint16_t maxTnsOrder, const unsigned  n, bool* const commonFlag)
{
  int16_t* const parCor1 = tnsData1.coeffParCor[n];
  int16_t* const parCor2 = tnsData2.coeffParCor[n];

  for (uint16_t s = 0; s < maxTnsOrder; s++) // synchronize coeffs
  {
    parCor1[s] = (parCor1[s] + parCor2[s] + 1) >> 1;
  }
  tnsData1.coeffResLow[n] = false; // optimize synchronized coeffs
  tnsData1.filterOrder[n] = predictor.calcOptTnsCoeffs (parCor1, tnsData1.coeff[n], &tnsData1.coeffResLow[n],
                                                        maxTnsOrder, UCHAR_MAX /*max pred gain*/, 0, LP_DEPTH);
  tnsData1.numFilters[n] = (tnsData1.filterOrder[n] > 0 ? 1 : 0);
  memcpy (&tnsData2, &tnsData1, sizeof (TnsData));  // synchronize

  if (commonFlag != nullptr) *commonFlag &= (true);
}

static inline void applyTnsCoeffPreProcessing (LinearPredictor& predictor, TnsData& tnsData1, TnsData& tnsData2,
                                               const uint16_t maxTnsOrder, const unsigned  n, bool* const commonFlag, const int16_t fact)
{
  const int32_t  weightI = __min (64, fact); // crosstalk constant
  const int32_t  weightD = 128 - weightI; // (1 - crosstalk) * 128
  int16_t* const parCor1 = tnsData1.coeffParCor[n];
  int16_t* const parCor2 = tnsData2.coeffParCor[n];

  for (uint16_t s = 0; s < maxTnsOrder; s++) // apply crosstalking
  {
    const int16_t coeff1 = parCor1[s];

    parCor1[s] = int16_t ((coeff1 * weightD + parCor2[s] * weightI + 64) >> 7);
    parCor2[s] = int16_t ((coeff1 * weightI + parCor2[s] * weightD + 64) >> 7);
  }
  tnsData1.coeffResLow[n] = false; // optimize coeffs of channel 1
  tnsData1.filterOrder[n] = predictor.calcOptTnsCoeffs (parCor1, tnsData1.coeff[n], &tnsData1.coeffResLow[n],
                                                        maxTnsOrder, UCHAR_MAX /*max pred gain*/, 0, LP_DEPTH);
  tnsData1.numFilters[n] = (tnsData1.filterOrder[n] > 0 ? 1 : 0);

  tnsData2.coeffResLow[n] = false; // optimize coeffs of channel 2
  tnsData2.filterOrder[n] = predictor.calcOptTnsCoeffs (parCor2, tnsData2.coeff[n], &tnsData2.coeffResLow[n],
                                                        maxTnsOrder, UCHAR_MAX /*max pred gain*/, 0, LP_DEPTH);
  tnsData2.numFilters[n] = (tnsData2.filterOrder[n] > 0 ? 1 : 0);

  if (commonFlag != nullptr) *commonFlag &= (tnsData1.coeffResLow[n] == tnsData2.coeffResLow[n] && tnsData1.filterOrder[n] == tnsData2.filterOrder[n]);
  if (commonFlag != nullptr && *commonFlag)
  {
    const int32_t* coeff1 = (int32_t*) tnsData1.coeff[n];  // fast
    const int32_t* coeff2 = (int32_t*) tnsData2.coeff[n];  // comp

    *commonFlag &= (*coeff1 == *coeff2); // might not be portable!
  }
}

static inline uint8_t brModeAndFsToMaxSfbLong (const unsigned bitRateMode, const unsigned samplingRate)
{
  // max. for fs of 44 kHz: band 47 (19.3 kHz), 48 kHz: 45 (19.5 kHz), 64 kHz: 39 (22.0 kHz)
  return __max (39, (0x20A000 + (samplingRate >> 1)) / samplingRate) - 9 + bitRateMode - (samplingRate < 48000 ? bitRateMode >> 3 : 0);
}

static inline uint8_t brModeAndFsToMaxSfbShort(const unsigned bitRateMode, const unsigned samplingRate)
{
  // max. for fs of 44 kHz: band 13 (19.3 kHz), 48 kHz: 13 (21.0 kHz), 64 kHz: 11 (23.0 kHz)
  return (samplingRate > 51200 ? 11 : 13) - 2 + (bitRateMode >> 2);
}

static inline void findActualBandwidthShort (uint8_t* const maxSfbShort, const uint16_t* sfbOffsets,
                                             const int32_t* mdctSignals, const int32_t* mdstSignals, const unsigned nSamplesInShort)
{
  const uint16_t b = sfbOffsets[1];  // beginning of search region
  uint8_t   maxSfb = __max (1, *maxSfbShort);
  uint16_t sfbOffs = sfbOffsets[maxSfb - 1];

  for (uint16_t e = sfbOffsets[maxSfb] - 1; e >= b; e--) // search
  {
    int32_t maxAbs = abs (mdctSignals[e]);

    if (mdstSignals != nullptr)
    {
      maxAbs = __max (maxAbs, abs (mdstSignals[e]));
      for (uint16_t w = 7; w > 0; w--)
      {
        maxAbs = __max (maxAbs, abs (mdctSignals[e + w * nSamplesInShort]));
        maxAbs = __max (maxAbs, abs (mdstSignals[e + w * nSamplesInShort]));
      }
    }
    else
    {
      for (uint16_t w = 7; w > 0; w--)
      {
        maxAbs = __max (maxAbs, abs (mdctSignals[e + w * nSamplesInShort]));
      }
    }
    if (maxAbs > maxSfb * (SA_EPS >> 1)) break;

    if (e == sfbOffs) sfbOffs = sfbOffsets[(--maxSfb) - 1];
  }

  if (*maxSfbShort > maxSfb) *maxSfbShort = maxSfb;
}

static inline uint8_t stereoCorrGrouping (const SfbGroupData& grpData, const unsigned nSamplesInFrame, uint8_t* stereoCorrData)
{
  const uint16_t numWinGroup = grpData.numWindowGroups;
  const uint16_t numBandsWin = nSamplesInFrame >> (SA_BW_SHIFT + 3);
  uint32_t m = 0, w;

  for (uint16_t gr = 0; gr < numWinGroup; gr++)
  {
    const uint16_t grpLength = grpData.windowGroupLength[gr];
    const uint16_t grpLenFac = oneTwentyEightOver[grpLength]; // for grpStereoCorr/grpLength
    const uint16_t grpLenOff = ((grpLenFac & (grpLenFac - 1)) > 0 ? 0 : 64); // for rounding

    for (uint16_t b = 0; b < numBandsWin; b++)
    {
      uint32_t grpStereoCorr = 0;

      for (w = 0; w < grpLength; w++) grpStereoCorr += stereoCorrData[b + w * numBandsWin];

      if (b == 0) m += grpStereoCorr;
      grpStereoCorr = (grpStereoCorr * grpLenFac + grpLenOff) >> 7;

      for (w = 0; w < grpLength; w++) stereoCorrData[b + w * numBandsWin] = grpStereoCorr;
    }
    stereoCorrData += grpLength * numBandsWin;
  }

  return uint8_t ((m + 4) >> 3); // mean low-band correlation value
}

// ISO/IEC 23003-3, Table 75
static inline unsigned toFrameLength (const USAC_CCFL coreCoderFrameLength)
{
  return (unsigned) coreCoderFrameLength;
}

// ISO/IEC 23003-3, Table 73
static const uint8_t numberOfChannels[USAC_MAX_NUM_ELCONFIGS] = {0, 1, 2, 3, 4, 5, 6, 8, 2, 3, 4, 7, 8};

static inline unsigned toNumChannels (const USAC_CCI chConfigurationIndex)
{
  return numberOfChannels[__max (0, (signed char) chConfigurationIndex)];
}

// ISO/IEC 14496-3, Table 4.140
static const uint16_t sfbOffsetL0[42] = { // 88.2 and 96 kHz
    0,   4,   8,  12,  16,  20,  24,  28,  32,  36,  40,  44,  48,  52,  56,  64,  72,  80,  88,  96, 108,
  120, 132, 144, 156, 172, 188, 212, 240, 276, 320, 384, 448, 512, 576, 640, 704, 768, 832, 896, 960, 1024
};
// ISO/IEC 14496-3, Table 4.141
static const uint16_t sfbOffsetS0[13] = {
  0, 4, 8, 12, 16, 20, 24, 32, 40, 48, 64, 92, 128
};

// ISO/IEC 14496-3, Table 4.138
static const uint16_t sfbOffsetL1[48] = { // 64 kHz
    0,   4,   8,  12,  16,  20,  24,  28,  32,  36,  40,  44,  48,  52,  56,  64,  72,  80,  88, 100, 112, 124, 140, 156,
  172, 192, 216, 240, 268, 304, 344, 384, 424, 464, 504, 544, 584, 624, 664, 704, 744, 784, 824, 864, 904, 944, 984, 1024
};
// ISO/IEC 14496-3, Table 4.139
static const uint16_t sfbOffsetS1[13] = {
  0, 4, 8, 12, 16, 20, 24, 32, 40, 48, 64, 92, 128
};

// ISO/IEC 14496-3, Table 4.131
static const uint16_t sfbOffsetL2[52] = { // 32, 44.1, and 48 kHz
    0,   4,   8,  12,  16,  20,  24,  28,  32,  36,  40,  48,  56,  64,  72,  80,  88,  96, 108, 120, 132, 144, 160, 176, 196, 216, 240,
  264, 292, 320, 352, 384, 416, 448, 480, 512, 544, 576, 608, 640, 672, 704, 736, 768, 800, 832, 864, 896, 928, 960/*!*/, 992/*!*/, 1024
};
// ISO/IEC 14496-3, Table 4.130
static const uint16_t sfbOffsetS2[15] = {
  0, 4, 8, 12, 16, 20, 28, 36, 44, 56, 68, 80, 96, 112, 128
};

// ISO/IEC 14496-3, Table 4.136
static const uint16_t sfbOffsetL3[48] = { // 22.05 and 24 kHz
    0,   4,   8,  12,  16,  20,  24,  28,  32,  36,  40,  44,  52,  60,  68,  76,  84,  92, 100, 108, 116, 124, 136, 148,
  160, 172, 188, 204, 220, 240, 260, 284, 308, 336, 364, 396, 432, 468, 508, 552, 600, 652, 704, 768, 832, 896, 960, 1024
};
// ISO/IEC 14496-3, Table 4.137
static const uint16_t sfbOffsetS3[16] = {
  0, 4, 8, 12, 16, 20, 24, 28, 36, 44, 52, 64, 76, 92, 108, 128
};

// ISO/IEC 14496-3, Table 4.134
static const uint16_t sfbOffsetL4[44] = { // 11.025, 12, and 16 kHz
    0,   8,  16,  24,  32,  40,  48,  56,  64,  72,  80,  88, 100, 112, 124, 136, 148, 160, 172, 184, 196, 212,
  228, 244, 260, 280, 300, 320, 344, 368, 396, 424, 456, 492, 532, 572, 616, 664, 716, 772, 832, 896, 960, 1024
};
// ISO/IEC 14496-3, Table 4.135
static const uint16_t sfbOffsetS4[16] = {
  0, 4, 8, 12, 16, 20, 24, 28, 32, 40, 48, 60, 72, 88, 108, 128
};

// ISO/IEC 14496-3, Table 4.132
static const uint16_t sfbOffsetL5[41] = { // 8 kHz
    0,  12,  24,  36,  48,  60,  72,  84,  96, 108, 120, 132, 144, 156, 172, 188, 204, 220, 236, 252, 268,
  288, 308, 328, 348, 372, 396, 420, 448, 476, 508, 544, 580, 620, 664, 712, 764, 820, 880, 944, 1024
};
// ISO/IEC 14496-3, Table 4.133
static const uint16_t sfbOffsetS5[16] = {
  0, 4, 8, 12, 16, 20, 24, 28, 36, 44, 52, 60, 72, 88, 108, 128
};

// long-window SFB offset tables
static const uint16_t* swbOffsetsL[USAC_NUM_FREQ_TABLES] = {
  sfbOffsetL0, sfbOffsetL1, sfbOffsetL2, sfbOffsetL3, sfbOffsetL4, sfbOffsetL5
};
static const uint8_t numSwbOffsetL[USAC_NUM_FREQ_TABLES] = {42, 48, 52, 48, 44, 41};

// short-window SFB offset tables
static const uint16_t* swbOffsetsS[USAC_NUM_FREQ_TABLES] = {
  sfbOffsetS0, sfbOffsetS1, sfbOffsetS2, sfbOffsetS3, sfbOffsetS4, sfbOffsetS5
};
static const uint8_t numSwbOffsetS[USAC_NUM_FREQ_TABLES] = {13, 13, 15, 16, 16, 16};

// ISO/IEC 23003-3, Table 79
static const uint8_t freqIdxToSwbTableIdxAAC[USAC_NUM_SAMPLE_RATES + 2] = {
  /*96000*/ 0, 0, 1, 2, 2, 2,/*24000*/ 3, 3, 4, 4, 4, 5, 5, // AAC
  255, 255, 1, 2, 2, 2, 2, 2,/*25600*/ 3, 3, 3, 4, 4, 4, 4 // USAC
};
#if !RESTRICT_TO_AAC
static const uint8_t freqIdxToSwbTableIdx768[USAC_NUM_SAMPLE_RATES + 2] = {
  /*96000*/ 0, 0, 0, 1, 1, 2,/*24000*/ 2, 2, 3, 4, 4, 4, 4, // AAC
  255, 255, 0, 1, 2, 2, 2, 2,/*25600*/ 2, 3, 3, 3, 3, 4, 4 // USAC
};
#endif

// ISO/IEC 23003-3, Table 131
static const uint8_t tnsScaleFactorBandLimit[2 /*long/short*/][USAC_NUM_FREQ_TABLES] = { // TNS_MAX_BANDS
#if 0 // RESTRICT_TO_AAC
  {31, 34, 51 /*to be corrected to 42 (44.1) and 40 (48 kHz)!*/, 46, 42, 39}, {9, 10, 14, 14, 14, 14}
#else
  {31, 34, 51 /*to be corrected to 42 (44.1) and 40 (48 kHz)!*/, 47, 43, 40}, {9, 10, 14, 15, 15, 15}
#endif
};

static const uint8_t sbrRateOffset[10] = {7, 6, 6, 8, 7, 8, 8, 8, 8, 8}; // used for scaleSBR

// scale_factor_grouping map
// group lengths based on transient location:  1133, 1115, 2114, 3113, 4112, 5111, 3311, 1331
static const uint8_t scaleFactorGrouping[8] = {0x1B, 0x0F, 0x47, 0x63, 0x71, 0x78, 0x6C, 0x36};

static const uint8_t windowGroupingTable[8][NUM_WINDOW_GROUPS] = { // for window_group_length
  {1, 1, 3, 3}, {1, 1, 1, 5}, {2, 1, 1, 4}, {3, 1, 1, 3}, {4, 1, 1, 2}, {5, 1, 1, 1}, {3, 3, 1, 1}, {1, 3, 3, 1}
};

// window_sequence equalizer
static const USAC_WSEQ windowSequenceSynch[5][5] = {  // 1st: chan index 0, 2nd: chan index 1
  {ONLY_LONG,   LONG_START,  EIGHT_SHORT, LONG_STOP,   STOP_START }, // left: ONLY_LONG
#if RESTRICT_TO_AAC
  {LONG_START,  LONG_START,  EIGHT_SHORT, EIGHT_SHORT, STOP_START }, // Left: LONG_START
#else
  {LONG_START,  LONG_START,  EIGHT_SHORT, STOP_START,  STOP_START }, // Left: LONG_START
#endif
  {EIGHT_SHORT, EIGHT_SHORT, EIGHT_SHORT, EIGHT_SHORT, EIGHT_SHORT}, // Left: EIGHT_SHORT
#if RESTRICT_TO_AAC
  {LONG_STOP,   EIGHT_SHORT, EIGHT_SHORT, LONG_STOP,   STOP_START }, // Left: LONG_STOP
#else
  {LONG_STOP,   STOP_START,  EIGHT_SHORT, LONG_STOP,   STOP_START }, // Left: LONG_STOP
#endif
  {STOP_START,  STOP_START,  EIGHT_SHORT, STOP_START,  STOP_START }  // Left: STOP_START
};

// private helper functions
unsigned ExhaleEncoder::applyTnsToWinGroup (SfbGroupData& grpData, const uint8_t grpIndex, const uint8_t maxSfb, TnsData& tnsData,
                                            const unsigned channelIndex, const unsigned n, const bool realOnlyCalc)
{
  const uint16_t filtOrder = tnsData.filterOrder[n];
  const uint16_t*    grpSO = &grpData.sfbOffsets[m_numSwbShort * grpIndex];
  const bool   eightShorts = (grpData.numWindowGroups > 1);
  unsigned errorValue = 0; // no error

  if ((grpIndex >= NUM_WINDOW_GROUPS) || (maxSfb > (eightShorts ? MAX_NUM_SWB_SHORT : MAX_NUM_SWB_LONG)) || (channelIndex >= USAC_MAX_NUM_CHANNELS))
  {
    return 1; // invalid arguments error
  }

  if (filtOrder > 0) // determine TNS filter length in SFBs and apply TNS analysis filtering
  {
    const int numSwbWin = (eightShorts ? m_numSwbShort : m_numSwbLong);
    uint8_t tnsMaxBands = tnsScaleFactorBandLimit[eightShorts ? 1 : 0][m_swbTableIdx];
    int     tnsStartSfb = 3 + 32000 / toSamplingRate (m_frequencyIdx);  // 8-short TNS start

    if (!eightShorts)
    {
      const unsigned samplingRate = toSamplingRate (m_frequencyIdx); // refine TNS_MAX_BANDS
      const unsigned tnsStartOffs = (m_specAnaCurr[channelIndex] & 31) << SA_BW_SHIFT;

      if ((samplingRate >= 46009) && (samplingRate < 55426)) tnsMaxBands = 40; // for 48 kHz
      else
      if ((samplingRate >= 37566) && (samplingRate < 46009)) tnsMaxBands = 42; // & 44.1 kHz

      while (grpSO[tnsStartSfb] < tnsStartOffs) tnsStartSfb++;  // start band for TNS filter
    }
    if ((tnsMaxBands = __min (tnsMaxBands, maxSfb)) <= tnsStartSfb) tnsStartSfb = numSwbWin;

    if ((tnsData.filterLength[n] = __max (0, numSwbWin - tnsStartSfb)) > 0)
    {
      int32_t* const signal = m_mdctSignals[channelIndex];
      const short offs = grpSO[tnsStartSfb];
      uint16_t       s = grpSO[tnsMaxBands] - offs;
      short filterC[MAX_PREDICTION_ORDER] = {0, 0, 0, 0};
      int32_t* predSig = &signal[grpSO[tnsMaxBands]]; // end of signal region to be filtered

      errorValue |= m_linPredictor.quantTnsToLpCoeffs (tnsData.coeff[n], filtOrder, tnsData.coeffResLow[n], tnsData.coeffParCor[n], filterC);

      // back up the leading MDCT samples
      memcpy (m_tempIntBuf, &signal[offs - MAX_PREDICTION_ORDER], MAX_PREDICTION_ORDER * sizeof (int32_t));
      // TNS compliance: set them to zero
      memset (&signal[offs - MAX_PREDICTION_ORDER], 0, MAX_PREDICTION_ORDER * sizeof (int32_t));

      if (filtOrder >= 4) // max. order 4
      {
        for (predSig--; s > 0; s--)
        {
          const int64_t predSample = *(predSig - 1) * (int64_t) filterC[0] + *(predSig - 2) * (int64_t) filterC[1] +
                                     *(predSig - 3) * (int64_t) filterC[2] + *(predSig - 4) * (int64_t) filterC[3];
          *(predSig--) += int32_t ((predSample + (1 << (LP_DEPTH - 2))) >> (LP_DEPTH - 1));
        }
      }
      else if (filtOrder == 3) // order 3
      {
        for (predSig--; s > 0; s--)
        {
          const int64_t predSample = *(predSig - 1) * (int64_t) filterC[0] + *(predSig - 2) * (int64_t) filterC[1] +
                                     *(predSig - 3) * (int64_t) filterC[2];
          *(predSig--) += int32_t ((predSample + (1 << (LP_DEPTH - 2))) >> (LP_DEPTH - 1));
        }
      }
      else // save 1-2 MACs, order 2 or 1
      {
        for (predSig--; s > 0; s--)
        {
          const int64_t predSample = *(predSig - 1) * (int64_t) filterC[0] + *(predSig - 2) * (int64_t) filterC[1];

          *(predSig--) += int32_t ((predSample + (1 << (LP_DEPTH - 2))) >> (LP_DEPTH - 1));
        }
      }
      // restore the leading MDCT samples
      memcpy (&signal[offs - MAX_PREDICTION_ORDER], m_tempIntBuf, MAX_PREDICTION_ORDER * sizeof (int32_t));

      // compute RMS data after filtering
      errorValue |= m_specAnalyzer.getMeanAbsValues (signal, realOnlyCalc ? nullptr : m_mdstSignals[channelIndex],
                                                     grpSO[grpData.sfbsPerGroup], (eightShorts ? USAC_MAX_NUM_CHANNELS : channelIndex),
                                                     grpSO /*below TNS*/, __min ((int) grpData.sfbsPerGroup, tnsStartSfb),
                                                     &grpData.sfbRmsValues[m_numSwbShort * grpIndex]);
      errorValue |= m_specAnalyzer.getMeanAbsValues (signal, nullptr /*no TNS on MDST*/, grpSO[grpData.sfbsPerGroup], channelIndex,
                                                     &grpSO[tnsStartSfb], __max (0, grpData.sfbsPerGroup - tnsStartSfb),
                                                     &grpData.sfbRmsValues[m_numSwbShort * grpIndex + tnsStartSfb]);
    }
    else tnsData.filterOrder[n] = tnsData.numFilters[n] = 0; // disable length-0 TNS filters
  } // if order > 0

  return errorValue;
}

unsigned ExhaleEncoder::eightShortGrouping (SfbGroupData& grpData, uint16_t* const grpOffsets,
                                            int32_t* const mdctSignal, int32_t* const mdstSignal)
{
  const unsigned nSamplesInFrame = toFrameLength (m_frameLength);
  const unsigned nSamplesInShort = nSamplesInFrame >> 3;
  int32_t* const tempIntBuf/*2*/ = m_timeSignals[1]; // NOTE: requires at least stereo input
  unsigned grpStartLine = nSamplesInFrame;

  if ((grpOffsets == nullptr) || (mdctSignal == nullptr))
  {
    return 1; // invalid arguments error
  }

  for (short gr = grpData.numWindowGroups - 1; gr >= 0; gr--) // grouping, 14496-3 Fig. 4.24
  {
    const unsigned   grpLength = grpData.windowGroupLength[gr];
    uint16_t* const  grpOffset = &grpOffsets[m_numSwbShort * gr];
    int32_t* const  grpMdctSig = &mdctSignal[grpStartLine -= nSamplesInShort * grpLength];
    int32_t* const  grpMdstSig = (mdstSignal != nullptr ? &mdstSignal[grpStartLine] : nullptr);

    for (uint16_t b = 0; b < m_numSwbShort; b++)
    {
      const unsigned swbOffset = grpOffsets[b];
      const unsigned numCoeffs = __min (grpOffsets[b + 1], nSamplesInShort) - swbOffset;

      // adjust scale factor band offsets
      grpOffset[b] = uint16_t (grpStartLine + swbOffset * grpLength);
      // interleave spectral coefficients
      for (uint16_t w = 0; w < grpLength; w++)
      {
        memcpy (&m_tempIntBuf[grpOffset[b] + w * numCoeffs], &grpMdctSig[swbOffset + w * nSamplesInShort], numCoeffs * sizeof (int32_t));
        if (grpMdstSig != nullptr)
        {
          memcpy (&tempIntBuf[grpOffset[b] + w * numCoeffs], &grpMdstSig[swbOffset + w * nSamplesInShort], numCoeffs * sizeof (int32_t));
        }
      }
    }
    grpOffset[m_numSwbShort] = uint16_t (grpStartLine + nSamplesInShort * grpLength);
  } // for gr

  memcpy (mdctSignal, m_tempIntBuf, nSamplesInFrame * sizeof (int32_t));
  if (mdstSignal != nullptr)
  {
    memcpy (mdstSignal, tempIntBuf, nSamplesInFrame * sizeof (int32_t));
  }

  return 0; // no error
}

unsigned ExhaleEncoder::getOptParCorCoeffs (const SfbGroupData& grpData, const uint8_t maxSfb, TnsData& tnsData,
                                            const unsigned channelIndex, const uint8_t firstGroupIndexToTest /*= 0*/)
{
  const unsigned nSamplesInFrame = toFrameLength (m_frameLength);
  const unsigned tnsStartSfb = 3 + 32000 / toSamplingRate (m_frequencyIdx); // 8-short start
  uint32_t temp, predGainMax = 0;

  if ((maxSfb <= tnsStartSfb) || (channelIndex >= USAC_MAX_NUM_CHANNELS))
  {
    return 0; // invalid arguments error
  }

  if (grpData.numWindowGroups == 1) // LONG window: use ParCor coeffs from spectral analyzer
  {
    tnsData.coeffResLow[0] = false;
    tnsData.filterDownward[0] = false; // enforce direction = 0 for now, detection difficult
    tnsData.filterOrder[0] = (uint8_t) m_specAnalyzer.getLinPredCoeffs (tnsData.coeffParCor[0], channelIndex);
    tnsData.firstTnsWindow = 0;

    if (tnsData.filterOrder[0] > 0) // try to reduce TNS start band as long as SNR increases
    {
      const uint16_t filtOrder = tnsData.filterOrder[0];
      uint16_t s = 0,b = __min ((m_specAnaCurr[channelIndex] & 31) + 2, (nSamplesInFrame - filtOrder) >> SA_BW_SHIFT);
      short filterC[MAX_PREDICTION_ORDER] = {0, 0, 0, 0};
      int32_t* predSig = &m_mdctSignals[channelIndex][b << SA_BW_SHIFT]; // TNS start offset

      m_linPredictor.parCorToLpCoeffs (tnsData.coeffParCor[0], filtOrder, filterC);

      for (b--, predSig--; b > 0; b--) // start a bit higher; b is in spectr. analysis units
      {
        uint64_t sumAbsOrg = 0, sumAbsTns = 0;

        if (filtOrder >= 4) // max. order 4
        {
          for (s = SA_BW; s > 0; s--) // get 4th-order TNS residual
          {
            const int64_t predSample = *(predSig - 1) * (int64_t) filterC[0] + *(predSig - 2) * (int64_t) filterC[1] +
                                       *(predSig - 3) * (int64_t) filterC[2] + *(predSig - 4) * (int64_t) filterC[3];
            const int64_t mdctSample = *(predSig--);
            const int64_t resiSample = mdctSample + ((predSample + (1 << 8)) >> 9);

            sumAbsOrg += abs (mdctSample);  sumAbsTns += abs (resiSample);
          }
        }
        else if (filtOrder == 3) // order 3
        {
          for (s = SA_BW; s > 0; s--) // get 3rd-order TNS residual
          {
            const int64_t predSample = *(predSig - 1) * (int64_t) filterC[0] + *(predSig - 2) * (int64_t) filterC[1] +
                                       *(predSig - 3) * (int64_t) filterC[2];
            const int64_t mdctSample = *(predSig--);
            const int64_t resiSample = mdctSample + ((predSample + (1 << 8)) >> 9);

            sumAbsOrg += abs (mdctSample);  sumAbsTns += abs (resiSample);
          }
        }
        else // save 1-2 MACs, order 2 or 1
        {
          for (s = SA_BW; s > 0; s--) // get 2nd-order TNS residual
          {
            const int64_t predSample = *(predSig - 1) * (int64_t) filterC[0] + *(predSig - 2) * (int64_t) filterC[1];
            const int64_t mdctSample = *(predSig--);
            const int64_t resiSample = mdctSample + ((predSample + (1 << 8)) >> 9);

            sumAbsOrg += abs (mdctSample);  sumAbsTns += abs (resiSample);
          }
        }
        if (sumAbsOrg * 17 <= sumAbsTns * 16) break; // band SNR reduced by more than 0.5 dB
      }
      m_specAnaCurr[channelIndex] = (m_specAnaCurr[channelIndex] & (UINT_MAX - 31)) | (b + 1);
    } // if order > 0

    return (m_specAnaCurr[channelIndex] >> 24) & UCHAR_MAX; // spectral analyzer's pred gain
  }

  // SHORT window: for each length-1 group, get TNS filter, then determine best filter order
  tnsData.firstTnsWindow = UCHAR_MAX;
  for (uint8_t n = 0, gr = 0; gr < grpData.numWindowGroups; gr++)
  {
    if (grpData.windowGroupLength[gr] == 1)
    {
      tnsData.coeffResLow[n] = false;
      tnsData.filterDownward[n] = false; // force direction = 0 for now, detection difficult
      tnsData.filterOrder[n] = 0;
      if (tnsData.firstTnsWindow == UCHAR_MAX) tnsData.firstTnsWindow = gr;

      if (gr < firstGroupIndexToTest)
      {
        memset (tnsData.coeffParCor[n], 0, MAX_PREDICTION_ORDER * sizeof (int16_t));
      }
      else // first length-one group tested
      {
        const int32_t* signal = m_mdctSignals[channelIndex];
        const uint16_t* grpSO = &grpData.sfbOffsets[m_numSwbShort * gr];
        uint32_t predGainCurr, predGainPrev, bestOrder = MAX_PREDICTION_ORDER;

        temp = m_linPredictor.calcParCorCoeffs (&signal[grpSO[tnsStartSfb]], grpSO[maxSfb] - grpSO[tnsStartSfb], bestOrder, tnsData.coeffParCor[n]);
        if (predGainMax < temp) predGainMax = temp;  // maximum pred gain of filtered groups

        predGainCurr = (temp >> 24) & UCHAR_MAX;
        predGainPrev = (temp >> 16) & UCHAR_MAX;
        while ((bestOrder > 1) && (predGainPrev >= predGainCurr))  // lowest-order gain max.
        {
          bestOrder--;
          predGainCurr = predGainPrev;
          predGainPrev = (temp >> (8 * bestOrder - 16)) & UCHAR_MAX;
        }
        tnsData.filterOrder[n] = uint8_t ((bestOrder == 1) && (tnsData.coeffParCor[n] == 0) ? 0 : bestOrder);
      }
      n++;
    }
  } // for gr

  return (predGainMax >> 24) & UCHAR_MAX; // max pred gain of all orders and length-1 groups
}

#ifndef NO_DTX_MODE
uint32_t ExhaleEncoder::getThr (const unsigned channelIndex, const unsigned sfbIndex)
{
  const uint16_t* const sfbLoudMem = m_sfbLoudMem[channelIndex][sfbIndex];
  uint32_t sumSfbLoud = 0;

  for (int16_t s = 31; s >= 0; s--) sumSfbLoud += sfbLoudMem[s];
  sumSfbLoud = (sumSfbLoud + 32) >> 6; // -6 dB

  return sumSfbLoud * (sumSfbLoud >> (toSamplingRate (m_frequencyIdx) >> 13)); // scaled SMR
}
#endif

unsigned ExhaleEncoder::psychBitAllocation () // perceptual bit-allocation via scale factors
{
  const unsigned nChannels       = toNumChannels (m_channelConf);
  const unsigned nSamplesInFrame = toFrameLength (m_frameLength);
  const unsigned samplingRate    = toSamplingRate (m_frequencyIdx);
  const unsigned lfeChannelIndex = (m_channelConf >= CCI_6_CH ? __max (5, nChannels - 1) : USAC_MAX_NUM_CHANNELS);
  const bool     useMaxBandwidth = (samplingRate < 37566 || m_shiftValSBR > 0);
  const uint8_t  maxSfbLong      = (useMaxBandwidth ? m_numSwbLong : brModeAndFsToMaxSfbLong (m_bitRateMode, samplingRate));
  const uint16_t scaleSBR        = (m_shiftValSBR > 0 || m_nonMpegExt ? sbrRateOffset[m_bitRateMode] : 0); // -25% rate
  const uint64_t scaleSr         = (samplingRate < 27713 ? (samplingRate < 23004 ? 32 : 34) - __min (3 << m_shiftValSBR, m_bitRateMode)
                                                         : (samplingRate < 37566 && m_bitRateMode != 3u ? 36 : 37)) - (nChannels >> 1);
  const uint64_t scaleBr         = (m_bitRateMode == 0 || m_frameCount <= 1 ? __min (32, 17u + (((samplingRate + (1 << 11)) >> 12) << 1) - (nChannels >> 1))
                                   : scaleSr - eightTimesSqrt256Minus[256 - m_bitRateMode] - __min (3, (m_bitRateMode - 1) >> 1)) + scaleSBR;
  uint32_t* sfbStepSizes = (uint32_t*) m_tempIntBuf;
  uint8_t  meanSpecFlat[USAC_MAX_NUM_CHANNELS];
//uint8_t  meanTempFlat[USAC_MAX_NUM_CHANNELS];
  unsigned ci = 0, s; // running index
  unsigned errorValue = 0; // no error

  // psychoacoustic processing of SFB RMS values yielding masking thresholds in m_tempIntBuf
  errorValue |= m_bitAllocator.initSfbStepSizes (m_scaleFacData, m_numSwbShort, m_specAnaCurr, m_tempAnaCurr,
                                                 nChannels, samplingRate, sfbStepSizes, lfeChannelIndex);

  // get means of spectral and temporal flatness for every channel
  m_bitAllocator.getChAverageSpecFlat (meanSpecFlat, nChannels);
//m_bitAllocator.getChAverageTempFlat (meanTempFlat, nChannels);

  for (unsigned el = 0; el < m_numElements; el++)  // element loop
  {
    CoreCoderData& coreConfig = *m_elementData[el];
    const unsigned nrChannels = (coreConfig.elementType & 1) + 1; // for UsacCoreCoderData()

    if (coreConfig.elementType >= ID_USAC_LFE) // LFE/EXT elements
    {
      SfbGroupData& grpData = coreConfig.groupingData[0];
      uint32_t*   stepSizes = &sfbStepSizes[ci * m_numSwbShort * NUM_WINDOW_GROUPS];
      const uint16_t*   off = grpData.sfbOffsets;
      const uint32_t*   rms = grpData.sfbRmsValues;
      uint8_t* scaleFactors = grpData.scaleFactors;

      for (uint16_t b = 0; b < grpData.sfbsPerGroup; b++)
      {
        const unsigned lfConst = (samplingRate < 27713 ? 1 : 2);
        const unsigned lfAtten = 4 + b * lfConst; // LF SNR boost, cf my M.Sc. thesis, p. 54
        const uint8_t sfbWidth = off[b + 1] - off[b];
        const uint64_t   scale = scaleBr * __min (32, lfAtten); // rate control part 1 (SFB)

        // scale step-sizes according to VBR mode, then derive scale factors from step-sizes
        stepSizes[b] = uint32_t (__max (BA_EPS, ((1u << 9) + stepSizes[b] * scale) >> 10));

        scaleFactors[b] = m_bitAllocator.getScaleFac (stepSizes[b], &m_mdctSignals[ci][off[b]], sfbWidth, rms[b]);
      }
      ci++;
    }
    else // SCE or CPE: bandwidth-to-max_sfb mapping, short-window grouping for each channel
    {
      const bool  eightShorts0 = (coreConfig.icsInfoCurr[0].windowSequence == EIGHT_SHORT);
      const TnsData&  tnsData0 = coreConfig.tnsData[0];
      const TnsData&  tnsData1 = coreConfig.tnsData[1];
      uint8_t realOnlyStartSfb = (eightShorts0 ? m_numSwbShort : m_numSwbLong) - __max (tnsData0.filterLength[0], tnsData1.filterLength[0]);

      if (coreConfig.commonWindow && (coreConfig.stereoMode == 0) && (m_perCorrHCurr[el] > SCHAR_MAX || m_perCorrLCurr[el] > (UCHAR_MAX * 5) / 8))
      {
        coreConfig.stereoMode = 1;
      }
      if (m_perCorrHCurr[el] > 128) // execute stereo pre-processing to increase correlation
      {
        const int16_t chanCorrSign = (coreConfig.stereoConfig & 2 ? -1 : 1);
        const uint16_t nSamplesMax = (useMaxBandwidth ? nSamplesInFrame : swbOffsetsL[m_swbTableIdx][__min (m_numSwbLong, maxSfbLong + 1)]);
        const uint8_t steppFadeLen = (eightShorts0 ? 4 : (coreConfig.tnsActive ? 32 : 64));
        const uint8_t steppFadeOff = ((m_bitRateMode + 77000 / samplingRate) & 6) << (eightShorts0 ? 2 : 5);
        const int64_t steppWeightI = __min (64, m_perCorrHCurr[el] - 128) >> (eightShorts0 || coreConfig.tnsActive ? 1 : 0);
        const int64_t steppWeightD = 128 - steppWeightI; // decrement, (1 - crosstalk) * 128

        for (uint16_t n = 0, gr = 0; gr < coreConfig.groupingData[0].numWindowGroups; gr++)
        {
          const uint8_t grpLength = coreConfig.groupingData[0].windowGroupLength[gr];
          const uint16_t*  grpOff = &coreConfig.groupingData[0].sfbOffsets[m_numSwbShort * gr];
          const uint16_t grpStart = grpOff[0] + steppFadeOff * grpLength;
          int32_t* sigR0 = &m_mdctSignals[ci][grpStart];
          int32_t* sigR1 = &m_mdctSignals[ci + 1][grpStart];
          int64_t xTalkI = 0, xTalkD = 0; // weights for crosstalk

          if ((grpLength == 1) && (tnsData0.numFilters[n] > 0 || tnsData1.numFilters[n] > 0))
          {
            const uint16_t maxLen = (eightShorts0 ? grpOff[m_numSwbShort] - 1 : __min (nSamplesInFrame - 1u, nSamplesMax)) - grpStart;
            int32_t prevR0 = 0; // NOTE: functions also on grouped
            int32_t prevR1 = 0; // MDCT spectra, but not properly!

            for (uint16_t w = 0; w < grpLength; w++) // sub-window
            {
              prevR0 = *(sigR0++); prevR1 = *(sigR1++); // processing starts at offset of 1!
              xTalkI = steppWeightI;
              xTalkD = steppWeightD * (2 * steppFadeLen - 1);

              for (s = steppFadeLen - 1; s > 0; s--, sigR0++, sigR1++) // start with fade-in
              {
                applyStereoPreProcessingReal (sigR0, sigR1, &prevR0, &prevR1, xTalkI, xTalkD, chanCorrSign);
                xTalkI += steppWeightI;
                xTalkD -= steppWeightD;
              }
            }
            for (s = maxLen - steppFadeLen * grpLength; s > 0; s--, sigR0++, sigR1++) // end
            {
              applyStereoPreProcessingReal (sigR0, sigR1, &prevR0, &prevR1, xTalkI, xTalkD, chanCorrSign);
            }
            if (eightShorts0 || (nSamplesMax >= nSamplesInFrame)) *sigR0 = *sigR1 = 0;

            realOnlyStartSfb = __min (realOnlyStartSfb, __min ((eightShorts0 ? 5 : 24), steppFadeOff / (eightShorts0 ? 4 : 7)));
          }
          else // TNS inactive, both MDCTs and MDSTs are available
          {
            const uint16_t maxLen = (eightShorts0 ? grpOff[m_numSwbShort] : nSamplesMax) - grpStart;
            int32_t* sigI0 = &m_mdstSignals[ci][grpStart]; // imag
            int32_t* sigI1 = &m_mdstSignals[ci + 1][grpStart];

            for (uint16_t w = 0; w < grpLength; w++) // sub-window
            {
              sigR0++;  sigR1++;  sigI0++;  sigI1++; // processing starts at an offset of 1!
              xTalkI = steppWeightI;
              xTalkD = steppWeightD * (2 * steppFadeLen - 1);

              for (s = steppFadeLen - 1; s > 0; s--, sigR0++, sigR1++, sigI0++, sigI1++)
              {
                applyStereoPreProcessingCplx (sigR0, sigR1, sigI0, sigI1, xTalkI, xTalkD, chanCorrSign);
                xTalkI += steppWeightI;
                xTalkD -= steppWeightD;
              }
            }
            for (s = maxLen - steppFadeLen * grpLength; s > 0; s--, sigR0++, sigR1++, sigI0++, sigI1++)
            {
              applyStereoPreProcessingCplx (sigR0, sigR1, sigI0, sigI1, xTalkI, xTalkD, chanCorrSign);
            }
          }
          if (grpLength == 1) n++;
        }
      } // if m_perCorrHCurr[el] > 128

      if ((errorValue == 0) && (coreConfig.stereoMode > 0)) // perform M/S, synch statistics
      {
        const uint8_t   numSwbFrame = (eightShorts0 ? m_numSwbShort : __min (m_numSwbLong, maxSfbLong + 1));
        const uint32_t peakIndexSte = __max ((m_specAnaCurr[ci] >> 5) & 2047, (m_specAnaCurr[ci + 1] >> 5) & 2047) << 5;

        errorValue = m_stereoCoder.applyPredJointStereo (m_mdctSignals[ci], m_mdctSignals[ci + 1],
                                                         m_mdstSignals[ci], m_mdstSignals[ci + 1],
                                                         coreConfig.groupingData[0], coreConfig.groupingData[1],
                                                         coreConfig.tnsData[0], coreConfig.tnsData[1],
                                                         numSwbFrame, coreConfig.stereoDataCurr,
                                                         m_bitRateMode, coreConfig.stereoMode > 1,
                                                         (coreConfig.stereoConfig & 2) > 0, realOnlyStartSfb,
                                                         &sfbStepSizes[m_numSwbShort * NUM_WINDOW_GROUPS *  ci],
                                                         &sfbStepSizes[m_numSwbShort * NUM_WINDOW_GROUPS * (ci + 1)]);
        if (errorValue >= 2) // signal M/S with complex prediction
        {
#if SP_MDST_PRED
          coreConfig.stereoConfig |= errorValue - 2; // cplx coefs
#endif
          coreConfig.stereoMode += 2; errorValue = 0;
        }
        m_specAnaCurr[ci    ] = (m_specAnaCurr[ci    ] & (UINT_MAX - 65504)) | peakIndexSte;
        m_specAnaCurr[ci + 1] = (m_specAnaCurr[ci + 1] & (UINT_MAX - 65504)) | peakIndexSte;
        meanSpecFlat[ci] = meanSpecFlat[ci + 1] = ((uint16_t) meanSpecFlat[ci] + (uint16_t) meanSpecFlat[ci + 1]) >> 1;
     // meanTempFlat[ci] = meanTempFlat[ci + 1] = ((uint16_t) meanTempFlat[ci] + (uint16_t) meanTempFlat[ci + 1]) >> 1;
      }
      else memset (coreConfig.stereoDataCurr, 0, (MAX_NUM_SWB_SHORT * NUM_WINDOW_GROUPS) * sizeof (uint8_t));

      errorValue |= m_bitAllocator.imprSfbStepSizes (m_scaleFacData, m_numSwbShort, m_mdctSignals, nSamplesInFrame, nrChannels,
                                                     ((32 + 5 * m_shiftValSBR) * samplingRate) >> 5, sfbStepSizes, ci, meanSpecFlat,
                                                     coreConfig.commonWindow, coreConfig.stereoDataCurr, coreConfig.stereoConfig);

      for (unsigned ch = 0; ch < nrChannels; ch++) // channel loop
      {
        SfbGroupData&  grpData = coreConfig.groupingData[ch];
        const bool eightShorts = (coreConfig.icsInfoCurr[ch].windowSequence == EIGHT_SHORT);
        const uint8_t maxSfbCh = grpData.sfbsPerGroup;
#if !RESTRICT_TO_AAC
        const uint8_t numSwbCh = (eightShorts ? m_numSwbShort : m_numSwbLong);
#endif
        const uint16_t mSfmSqr = (m_bitRateMode < 2 && samplingRate >= 27713 ? ((uint16_t) meanSpecFlat[ci] * meanSpecFlat[ci]) >> m_bitRateMode : 0);
        const uint16_t mSfmFac = 256u - (((32u + m_bitRateMode) * ((uint32_t) meanSpecFlat[ci] << 4) - mSfmSqr + (1u << 9)) >> 10);
        uint32_t*    stepSizes = &sfbStepSizes[ci * m_numSwbShort * NUM_WINDOW_GROUPS];

        memset (grpData.scaleFactors, 0, (MAX_NUM_SWB_SHORT * NUM_WINDOW_GROUPS) * sizeof (uint8_t));

        for (uint16_t gr = 0; gr < grpData.numWindowGroups; gr++)
        {
          const uint16_t* grpOff = &grpData.sfbOffsets[m_numSwbShort * gr];
          const uint32_t* grpRms = &grpData.sfbRmsValues[m_numSwbShort * gr];
          uint8_t*  grpScaleFacs = &grpData.scaleFactors[m_numSwbShort * gr];
          uint32_t* grpStepSizes = &stepSizes[m_numSwbShort * gr];

          for (unsigned b = 0; b < maxSfbCh; b++)
          {
            const unsigned lfConst = (samplingRate < 27713 && !eightShorts ? 1 : 2); // lfAtten: LF SNR boost, as in my M.Sc. thesis
            const unsigned lfAtten = (b <= 5 ? (eightShorts ? 1 : 4) + b * lfConst : 5 * lfConst - 1 + b + ((b + 5) >> 4));
            const uint8_t sfbWidth = grpOff[b + 1] - grpOff[b];
            const uint64_t   scale = scaleBr * mSfmFac * __min (32, lfAtten * grpData.numWindowGroups); // rate control part 1 (SFB)

            // scale step-sizes according to VBR mode & derive scale factors from step-sizes
            grpStepSizes[b] = uint32_t (__max (BA_EPS, ((1u << 17) + grpStepSizes[b] * scale) >> 18));

            grpScaleFacs[b] = m_bitAllocator.getScaleFac (grpStepSizes[b], &m_mdctSignals[ci][grpOff[b]], sfbWidth, grpRms[b]);
          }
        } // for gr

#if !RESTRICT_TO_AAC
        if ((maxSfbCh > 0) && m_noiseFilling[el] && (m_bitRateMode <= 3 || !eightShorts))
        {
          const uint32_t maxSfbCurr = (eightShorts ? (useMaxBandwidth ? 17 - (samplingRate >> 13) // was 14, good for 22.05 - 32 kHz
                                                                      : brModeAndFsToMaxSfbShort (m_bitRateMode, samplingRate)) : maxSfbLong);
          const bool keepMaxSfbCurr = ((samplingRate < 37566) || (samplingRate >= 46009 && samplingRate < 55426 && eightShorts));
          const uint8_t numSwbFrame = __min ((numSwbCh * ((maxSfbCh == maxSfbCurr) || (m_bitRateMode <= 2) || (m_shiftValSBR > 0) ? 4u : 3u)) >> 2,
                                      (eightShorts ? maxSfbCh : maxSfbLong) + (m_bitRateMode < 2 || m_bitRateMode > 3 || keepMaxSfbCurr ? 0u : 1u));
#ifndef NO_DTX_MODE
          if ((m_bitRateMode < 1) && (m_numElements == 1) && (samplingRate < 27713) && eightShorts)
          {
            for (s = 0; s < 26; s++) m_sfbLoudMem[ch][s][m_frameCount & 31] = uint16_t (sqrt (double (getThr (ch, s) << (samplingRate >> 13))));
          }
          if ((maxSfbCh < numSwbFrame) || (m_bitRateMode <= 2)) // increase coding bandwidth
#else
          if (maxSfbCh < numSwbFrame) // increase coding bandwidth
#endif
          {
            for (uint16_t gr = 0; gr < grpData.numWindowGroups; gr++)
            {
#ifndef NO_DTX_MODE
              const uint32_t*  grpRms = &grpData.sfbRmsValues[m_numSwbShort * gr];

              if ((m_bitRateMode < 1) && (m_numElements == 1) && (samplingRate < 27713))
              {
                const uint32_t*  refRms = &coreConfig.groupingData[1 - ch].sfbRmsValues[m_numSwbShort * gr];
                uint8_t*  grpStereoData = &coreConfig.stereoDataCurr[m_numSwbShort * gr];
                const unsigned sfbStart = __max (samplingRate < 18783 ? 17 : 24, m_specGapFiller.getFirstGapFillSfb ());

                for (s = sfbStart; s < maxSfbCh; s++)
                {
                  const double rmsValue = double (grpStereoData[s] > 0 ? (grpRms[s] + (uint64_t) refRms[s] + 1) >> 1 : grpRms[s]);
                  const unsigned sfbIdx = s - sfbStart;

                  m_sfbLoudMem[ch][sfbIdx][m_frameCount & 31] = __max (BA_EPS, uint16_t (sqrt (rmsValue)));
                  if (grpRms[s] < getThr (ch, sfbIdx)) grpData.scaleFactors[s + m_numSwbShort * gr] = 0;
                }
              }
              else
              if ((m_bitRateMode <= 4) && (meanSpecFlat[ci] <= (SCHAR_MAX >> 1))) // low-RMS
              {
                for (s = __max (samplingRate < 27713 ? (samplingRate < 18783 ? 17 : 24) : 22, m_specGapFiller.getFirstGapFillSfb ()); s < maxSfbCh; s++)
                {
                  if (grpRms[s] < ((3 * TA_EPS) >> 1)) grpData.scaleFactors[s + m_numSwbShort * gr] = 0;
                }
              }
#endif
              memset (&grpData.scaleFactors[maxSfbCh + m_numSwbShort * gr], 0, (numSwbFrame - maxSfbCh) * sizeof (uint8_t));
            }
            grpData.sfbsPerGroup = coreConfig.icsInfoCurr[ch].maxSfb = numSwbFrame;
          }
          if (ch > 0 && coreConfig.commonWindow) // resynchronize the two max_sfb for stereo
          {
            uint8_t& maxSfb0 = coreConfig.icsInfoCurr[0].maxSfb;
            uint8_t& maxSfb1 = coreConfig.icsInfoCurr[1].maxSfb;

            if (coreConfig.stereoMode > 0)
            {
              maxSfb0 = maxSfb1 = coreConfig.groupingData[0].sfbsPerGroup = grpData.sfbsPerGroup = __max (maxSfb0, maxSfb1);
            }
            coreConfig.commonMaxSfb = (maxSfb0 == maxSfb1);
          }
        }
#ifndef NO_DTX_MODE
        else if (m_noiseFilling[el] && (m_bitRateMode < 1) && (m_numElements == 1) && (samplingRate < 27713))
        {
          for (s = 0; s < 26; s++) m_sfbLoudMem[ch][s][m_frameCount & 31] = BA_EPS;
        }
#endif
#endif // !RESTRICT_TO_AAC
        ci++;
      } // for ch

      for (unsigned ch = 0; ch < nrChannels; ch++) // channel loop
      {
        SfbGroupData& grpData = coreConfig.groupingData[ch];
        TnsData&      tnsData = coreConfig.tnsData[ch];

        if (tnsData.numFilters[0] + tnsData.numFilters[1] + tnsData.numFilters[2] > 0)
        {
          s = tnsData.firstTnsWindow = 0; // store length-1 group map for bit-stream writing
          for (uint16_t gr = 0; gr < grpData.numWindowGroups; s += grpData.windowGroupLength[gr++])
          {
            if (grpData.windowGroupLength[gr] == 1) tnsData.firstTnsWindow |= (1u << s);
          }
        }
      } // for ch
    }
  } // for el

  return errorValue;
}

unsigned ExhaleEncoder::quantizationCoding ()  // apply MDCT quantization and entropy coding
{
  const unsigned nChannels        = toNumChannels (m_channelConf);
  const unsigned nSamplesInFrame  = toFrameLength (m_frameLength);
  const unsigned samplingRate     = toSamplingRate (m_frequencyIdx);
  const unsigned nSamplesTempAna  = (nSamplesInFrame * 25) >> 4; // pre-delay for look-ahead
  const bool     useMaxBandwidth  = (samplingRate < 37566 || m_shiftValSBR > 0);
  const unsigned* const coeffMagn = m_sfbQuantizer.getCoeffMagnPtr ();
  uint8_t  meanSpecFlat[USAC_MAX_NUM_CHANNELS];
  uint8_t  meanTempFlat[USAC_MAX_NUM_CHANNELS] = {208, 208, 208, 208, 208, 208, 208, 208};
  unsigned ci = 0, s; // running index
  unsigned errorValue = (coeffMagn == nullptr ? 1 : 0);

  // get means of spectral and temporal flatness for every channel
  m_bitAllocator.getChAverageSpecFlat (meanSpecFlat, nChannels);
  if ((m_bitRateMode < (2u >> m_shiftValSBR)) && (samplingRate >= 23004) && (samplingRate < 37566))
  {
    m_bitAllocator.getChAverageTempFlat (meanTempFlat, nChannels);
  }

  for (unsigned el = 0; el < m_numElements; el++)  // element loop
  {
    CoreCoderData& coreConfig = *m_elementData[el];
    const unsigned nrChannels = (coreConfig.elementType & 1) + 1; // for UsacCoreCoderData()

    if ((coreConfig.elementType < ID_USAC_LFE) && (coreConfig.stereoMode > 0)) // synch SFMs
    {
      meanSpecFlat[ci] = meanSpecFlat[ci + 1] = ((uint16_t) meanSpecFlat[ci] + (uint16_t) meanSpecFlat[ci + 1]) >> 1;
      meanTempFlat[ci] = meanTempFlat[ci + 1] = ((uint16_t) meanTempFlat[ci] + (uint16_t) meanTempFlat[ci + 1]) >> 1;
    }

    for (unsigned ch = 0; ch < nrChannels; ch++)   // channel loop
    {
      EntropyCoder& entrCoder = m_entropyCoder[ci];
      SfbGroupData&   grpData = coreConfig.groupingData[ch];
      const bool shortWinCurr = (coreConfig.icsInfoCurr[ch].windowSequence == EIGHT_SHORT);
      const bool shortWinPrev = (coreConfig.icsInfoPrev[ch].windowSequence == EIGHT_SHORT);
      char* const arithTuples = entrCoder.arithGetTuplePtr ();
      uint8_t sfIdxPred = UCHAR_MAX;

      if ((errorValue > 0) || (arithTuples == nullptr))
      {
        return 0; // an internal error
      }

      // back up entropy coder memory for use by bit-stream writer
      memcpy (m_tempIntBuf, arithTuples, (nSamplesInFrame >> 1) * sizeof (char));
      errorValue |= (entrCoder.getIsShortWindow () != shortWinPrev ? 1 : 0); // sanity check

      memset (m_mdctQuantMag[ci], 0, nSamplesInFrame * sizeof (uint8_t));  // initialization

      for (uint16_t gr = 0; gr < grpData.numWindowGroups; gr++)
      {
        const uint8_t grpLength = grpData.windowGroupLength[gr];
        const uint16_t*  grpOff = &grpData.sfbOffsets[m_numSwbShort * gr];
        uint32_t* const  grpRms = &grpData.sfbRmsValues[m_numSwbShort * gr]; // coding stats
        uint8_t*   grpScaleFacs = &grpData.scaleFactors[m_numSwbShort * gr];
        uint32_t estimBitCount = 0;
        unsigned lastSfb = 0, lastSOff = 0;

        errorValue |= entrCoder.initWindowCoding (m_indepFlag && (gr == 0), shortWinCurr);
        s = 0;

        for (uint16_t b = 0; b < grpData.sfbsPerGroup; b++)
        {
          // partial SFB ungrouping for entropy coding setup below
          const uint16_t swbSize = ((grpOff[b + 1] - grpOff[b]) * oneTwentyEightOver[grpLength]) >> 7; // sfbWidth / grpLength
          uint8_t* const swbMagn = &m_mdctQuantMag[ci][grpOff[b + 1] - swbSize];

          grpScaleFacs[b] = m_sfbQuantizer.quantizeSpecSfb (entrCoder, m_mdctSignals[ci], grpLength, grpOff, grpRms,
                                                            b, grpScaleFacs[b], sfIdxPred, m_mdctQuantMag[ci]);
          if ((b > 0) && (grpScaleFacs[b] < UCHAR_MAX) && (sfIdxPred == UCHAR_MAX))
          {
            // back-propagate first nonzero-SFB scale factor index
            memset (grpScaleFacs, grpScaleFacs[b], b * sizeof (uint8_t));
          }
          sfIdxPred = grpScaleFacs[b];

          // correct previous scale factor if the delta exceeds 60
          if ((b > 0) && (grpScaleFacs[b] > grpScaleFacs[b - 1] + INDEX_OFFSET))
          {
            const uint16_t sfbM1Start = grpOff[b - 1];
            const uint16_t sfbM1Width = grpOff[b] - sfbM1Start;
            const uint16_t swbM1Size  = (sfbM1Width * oneTwentyEightOver[grpLength]) >> 7; // sfbM1Width / grpLength

            grpScaleFacs[b - 1] = grpScaleFacs[b] - (b > 1 ? INDEX_OFFSET : 0);  // zero-out
            memset (&m_mdctQuantMag[ci][sfbM1Start], 0, sfbM1Width * sizeof (uint8_t));

            // correct SFB statistics with some bit count estimate
            grpRms[b - 1] = 1 + (sfbM1Width >> 3) + entrCoder.indexGetBitCount (b > 1 ? (int) grpScaleFacs[b - 1] - grpScaleFacs[b - 2] : 0);
            // correct entropy coding 2-tuples for the next window
            memset (&arithTuples[lastSOff], 1, (swbM1Size >> 1) * sizeof (char));
          }

          if (b > 0)
          {
            if ((grpRms[b - 1] >> 16) > 0) lastSfb = b - 1;
            estimBitCount += grpRms[b - 1] & USHRT_MAX;
          }
          // set up entropy coding 2-tuples for next SFB or window
          lastSOff = s;
          for (uint16_t c = 0; c < swbSize; c += 2)
          {
            arithTuples[s++] = __min (0xF, swbMagn[c] + swbMagn[c + 1] + 1); // 23003-3, 7.4
          }
        } // for b

        if (grpData.sfbsPerGroup > 0) // rate control part 2 to reach constrained VBR (CVBR)
        {
          const uint8_t maxSfbLong  = (useMaxBandwidth ? 54 - (samplingRate >> 13) : brModeAndFsToMaxSfbLong (m_bitRateMode, samplingRate));
          const uint8_t maxSfbShort = (useMaxBandwidth ? 19 - (samplingRate >> 13) : brModeAndFsToMaxSfbShort(m_bitRateMode, samplingRate));
          const uint16_t peakIndex  = (shortWinCurr ? 0 : (m_specAnaCurr[ci] >> 5) & 2047);
          const unsigned sfmBasedSfbStart = (shortWinCurr ? maxSfbShort - 2 + (meanSpecFlat[ci] >> 6) : maxSfbLong  - 6 + (meanSpecFlat[ci] >> 5)) +
                                            (shortWinCurr ? -3 + (((1 << 5) + meanTempFlat[ci]) >> 6) : -7 + (((1 << 4) + meanTempFlat[ci]) >> 5));
          const unsigned targetBitCount25 = ((60000 + 20000 * ((m_bitRateMode + m_shiftValSBR) >> (m_frameCount <= 1 ? 2 : 0))) * nSamplesInFrame) /
                                            (samplingRate * ((grpData.numWindowGroups + 1) >> 1));
          unsigned b = grpData.sfbsPerGroup - 1;

          if ((grpRms[b] >> 16) > 0) lastSfb = b;
          estimBitCount += grpRms[b] & USHRT_MAX;

#if EC_TRELLIS_OPT_CODING
          if (grpLength == 1) // finalize bit count estimate, RDOC
          {
            estimBitCount = m_sfbQuantizer.quantizeSpecRDOC (entrCoder, grpScaleFacs, estimBitCount + 2u,
                                                             grpOff, grpRms, grpData.sfbsPerGroup, m_mdctQuantMag[ci]);
            for (b = 1; b < grpData.sfbsPerGroup; b++)
            {
              // correct previous scale factor if delta exceeds 60
              if (grpScaleFacs[b] > grpScaleFacs[b - 1] + INDEX_OFFSET)
              {
                const uint16_t sfbM1Start = grpOff[b - 1];
                const uint16_t sfbM1Width = grpOff[b] - sfbM1Start;

                grpScaleFacs[b - 1] = grpScaleFacs[b] - (b > 1 ? INDEX_OFFSET : 0); // 0-out
                memset (&m_mdctQuantMag[ci][sfbM1Start], 0, sfbM1Width * sizeof (uint8_t));

                // correct statistics with some bit count estimate
                grpRms[b - 1] = 1 + (sfbM1Width >> 3) + entrCoder.indexGetBitCount (b > 1 ? (int) grpScaleFacs[b - 1] - grpScaleFacs[b - 2] : 0);
                // correct entropy coding 2-tuples for next window
                memset (&arithTuples[(sfbM1Start - grpOff[0]) >> 1], 1, (sfbM1Width >> 1) * sizeof (char));
              }
            }
          }
#endif
          // coarse-quantize near-Nyquist SFB with SBR @ 48-64 kHz
          b = 40 + (samplingRate >> 12);
          if ((m_shiftValSBR == 0) || (samplingRate < 23004) || shortWinCurr || (b > lastSfb)) b = lastSfb;

          while ((b >= sfmBasedSfbStart + (m_bitRateMode >> 1)) && (grpOff[b] > peakIndex) && ((grpRms[b] >> 16) <= 1) /*coarse quantization*/ &&
                 ((estimBitCount * 5 > targetBitCount25 * 2) || (grpLength > 1 /*no accurate bit count estim. available for grouped spectrum*/)))
          {
            b--; // search first coarsely quantized high-freq. SFB
          }
          lastSOff = b;

          for (b++; b <= lastSfb; b++)
          {
            if ((grpRms[b] >> 16) > 0) // re-quantize nonzero band
            {
#if RESTRICT_TO_AAC
              uint32_t maxVal = 1;
#else
              uint32_t maxVal = (shortWinCurr || !m_noiseFilling[el] ? 1 : (m_specAnaCurr[ci] >> 23) & 1); // 1 or 0
#endif
              estimBitCount -= grpRms[b] & USHRT_MAX;
              grpRms[b] = (maxVal << 16) + maxVal; // bit estimate
              maxVal = quantizeSfbWithMinSnr (coeffMagn, grpOff, b, grpLength, m_mdctQuantMag[ci], arithTuples, maxVal > 0);

              grpScaleFacs[b] = __min (SCHAR_MAX, m_sfbQuantizer.getScaleFacOffset ((double) maxVal));

              // correct SFB statistics with estimate of bit count
              grpRms[b] += 3 + entrCoder.indexGetBitCount ((int) grpScaleFacs[b] - grpScaleFacs[b - 1]);
              estimBitCount += grpRms[b] & USHRT_MAX;
            }
            else // re-repeat scale factor for zero quantized band
            {
              grpScaleFacs[b] = grpScaleFacs[b - 1];
            }
          }

          if (estimBitCount > targetBitCount25) // too many bits!!
          {
            for (b = lastSOff; b > 0; b--)
            {
              if ((grpRms[b] >> 16) > 0) // emergency re-quantizer
              {
#if RESTRICT_TO_AAC
                uint32_t maxVal = 1;
#else
                uint32_t maxVal = (shortWinCurr || !m_noiseFilling[el] ? 1 : (m_specAnaCurr[ci] >> 23) & 1); // 1 or 0
#endif
                estimBitCount -= grpRms[b] & USHRT_MAX;
                grpRms[b] = (maxVal << 16) + maxVal; // bit estim.
                maxVal = quantizeSfbWithMinSnr (coeffMagn, grpOff, b, grpLength, m_mdctQuantMag[ci], arithTuples, maxVal > 0);

                grpScaleFacs[b] = __min (SCHAR_MAX, m_sfbQuantizer.getScaleFacOffset ((double) maxVal));

                // correct SFB statistics with estimated bit count
                grpRms[b] += 3 + entrCoder.indexGetBitCount ((int) grpScaleFacs[b] - grpScaleFacs[b - 1]);
                estimBitCount += grpRms[b] & USHRT_MAX;
              }
              if (estimBitCount <= targetBitCount25) break;
            }

            for (b++; b <= lastSfb; b++) // re-repeat scale factor
            {
              if ((grpRms[b] >> 16) == 0) // a zero quantized band
              {
                grpScaleFacs[b] = grpScaleFacs[b - 1];
              }
            }
          } // if estimBitCount > targetBitCount25

          for (b = lastSfb + 1; b < grpData.sfbsPerGroup; b++)
          {
            if ((grpRms[b] >> 16) == 0) // HF zero quantized bands
            {
              grpScaleFacs[b] = grpScaleFacs[b - 1];
            }
          }

          if ((grpScaleFacs[0] == UCHAR_MAX) &&
#if !RESTRICT_TO_AAC
              !m_noiseFilling[el] &&
#endif
              (lastSfb == 0))  // ensure all scale factors are set
          {
            memset (grpScaleFacs, (gr == 1 ? grpData.scaleFactors[grpData.sfbsPerGroup - 1] : 0), grpData.sfbsPerGroup * sizeof (uint8_t));
          }
        }
      } // for gr

      // restore entropy coder memory for use by bit-stream writer
      memcpy (arithTuples, m_tempIntBuf, (nSamplesInFrame >> 1) * sizeof (char));
      entrCoder.setIsShortWindow (shortWinPrev);
#if !RESTRICT_TO_AAC
      s = 22050 + 7350 * m_bitRateMode; // compute channel-wise noise_level and noise_offset
      coreConfig.specFillData[ch] = (!m_noiseFilling[el] ? 0 : m_specGapFiller.getSpecGapFillParams (m_sfbQuantizer, m_mdctQuantMag[ci], m_numSwbShort,
                                                                                                     grpData, nSamplesInFrame, samplingRate, s,
                                                                                                     shortWinCurr ? 0 : meanSpecFlat[ci]));
      // NOTE: gap-filling SFB bit count might be inaccurate now since scale factors changed
      if (coreConfig.specFillData[ch] == 1) errorValue |= 1;
#endif
      s = ci + nrChannels - 1 - 2 * ch; // other channel in stereo
      if ((coreConfig.elementType < ID_USAC_LFE) && (m_shiftValSBR > 0)) // collect SBR data
      {
        const uint8_t msfVal = (shortWinPrev ? 31 : __max (2, __max (m_meanSpecPrev[ci], meanSpecFlat[ci]) >> 3));
        const uint8_t msfSte = (coreConfig.stereoMode == 0 ? 0 : (coreConfig.icsInfoPrev[s + ch - ci].windowSequence ==
                                 EIGHT_SHORT ? 31 : __max (2, __max (m_meanSpecPrev[s ], meanSpecFlat[s ]) >> 3)));
        int32_t  tmpValSynch = 0;

        memset (m_coreSignals[ci], 0, 10 * sizeof (int32_t));
#if ENABLE_INTERTES
        m_coreSignals[ci][0] = 0x40000000; // bs_interTes = 1 for all frames
#endif
        m_coreSignals[ci][0] |= 4 - int32_t (sqrt (0.75 * msfVal)); // filter mode, 0 = none

        if (ch > 0 && coreConfig.stereoMode > 0) // synch. sbr_grid(), sbr_invf() for stereo
        {
          tmpValSynch = (m_coreSignals[s][0] >> 21) & 3; // nEnv, bits 23-22
          m_coreSignals[ci][0] |= m_coreSignals[s][0] & 0x10000F; // bits 21
          m_coreSignals[s][0] |= m_coreSignals[ci][0] & 0x10000F; // and 4-1
        }
        m_coreSignals[ci][0] |= getSbrEnvelopeAndNoise (&m_coreSignals[ci][nSamplesTempAna - 64 + nSamplesInFrame], msfVal,
                                                        __max (m_meanTempPrev[ci], meanTempFlat[ci]) >> 3, m_bitRateMode == 0,
                                                        m_indepFlag, msfSte, tmpValSynch, nSamplesInFrame, &m_coreSignals[ci][1]);
        if (ch + 1 == nrChannels) // update the flatness histories
        {
          m_meanSpecPrev[ci] = meanSpecFlat[ci];  m_meanSpecPrev[s] = meanSpecFlat[s];
          m_meanTempPrev[ci] = meanTempFlat[ci];  m_meanTempPrev[s] = meanTempFlat[s];
        }
      }
      ci++;
    }
  } // for el

  return (errorValue > 0 ? 0 : m_outStream.createAudioFrame (m_elementData, m_entropyCoder, m_mdctSignals, m_mdctQuantMag, m_indepFlag,
                                                             m_numElements, m_numSwbShort, (uint8_t* const) m_tempIntBuf,
#if !RESTRICT_TO_AAC
                                                             m_timeWarping, m_noiseFilling, m_frameCount - 1u, m_indepPeriod,
#endif
                                                             m_shiftValSBR, m_coreSignals, m_outAuData, nSamplesInFrame)); // returns AU size
}

unsigned ExhaleEncoder::spectralProcessing ()  // complete ics_info(), calc TNS and SFB data
{
  const unsigned nChannels       = toNumChannels (m_channelConf);
  const unsigned nSamplesInFrame = toFrameLength (m_frameLength);
  const unsigned nSamplesInShort = nSamplesInFrame >> 3;
  const unsigned samplingRate    = toSamplingRate (m_frequencyIdx);
  const unsigned lfeChannelIndex = (m_channelConf >= CCI_6_CH ? __max (5, nChannels - 1) : USAC_MAX_NUM_CHANNELS);
  const bool     useMaxBandwidth = (samplingRate < 37566 || m_shiftValSBR > 0);
  unsigned ci = 0, s; // running index
  unsigned errorValue = 0; // no error

  // get spectral channel statistics for last frame, used for input bandwidth (BW) detection
  m_specAnalyzer.getSpectralBandwidth (m_bandwidPrev, nChannels);

  // spectral analysis for current MCLT signal (windowed time-samples for the current frame)
  errorValue |= m_specAnalyzer.spectralAnalysis (m_mdctSignals, m_mdstSignals, nChannels, nSamplesInFrame, samplingRate, lfeChannelIndex);

  // get spectral channel statistics for this frame, used for perceptual model & BW detector
  m_specAnalyzer.getSpecAnalysisStats (m_specAnaCurr, nChannels);
  m_specAnalyzer.getSpectralBandwidth (m_bandwidCurr, nChannels);

  for (unsigned el = 0; el < m_numElements; el++)  // element loop
  {
    CoreCoderData& coreConfig = *m_elementData[el];
    const unsigned nrChannels = (coreConfig.elementType & 1) + 1; // for UsacCoreCoderData()

    coreConfig.commonMaxSfb   = false;
    coreConfig.commonTnsData  = false;
    coreConfig.tnsActive      = false;
    coreConfig.tnsOnLeftRight = true;  // enforce tns_on_lr = 1 for now, detection difficult
    memset (coreConfig.tnsData, 0, nrChannels * sizeof (TnsData));

    if (coreConfig.elementType >= ID_USAC_LFE) // LFE/EXT elements
    {
      SfbGroupData& grpData = coreConfig.groupingData[0];
      uint16_t*  grpSO = grpData.sfbOffsets;
      IcsInfo& icsCurr = coreConfig.icsInfoCurr[0];

      memcpy (grpSO, swbOffsetsL[m_swbTableIdx], numSwbOffsetL[m_swbTableIdx] * sizeof (uint16_t));

      icsCurr.maxSfb = MAX_NUM_SWB_LFE;
      while (grpSO[icsCurr.maxSfb] > LFE_MAX) icsCurr.maxSfb--; // limit coefficients in LFE
      grpData.sfbsPerGroup = icsCurr.maxSfb;
      ci++;
    }
    else // SCE or CPE: bandwidth-to-max_sfb mapping, short-window grouping for each channel
    {
      coreConfig.stereoConfig = coreConfig.stereoMode = 0;

      if (coreConfig.commonWindow && (m_bitRateMode <= 4)) // stereo pre-processing analysis
      {
        const bool     eightShorts = (coreConfig.icsInfoCurr[0].windowSequence == EIGHT_SHORT);
        const uint8_t meanSpecFlat = (((m_specAnaCurr[ci] >> 16) & UCHAR_MAX) + ((m_specAnaCurr[ci + 1] >> 16) & UCHAR_MAX) + 1) >> 1;
        const uint16_t* const swbo = swbOffsetsL[m_swbTableIdx];
        const uint16_t nSamplesMax = (useMaxBandwidth ? nSamplesInFrame : swbo[brModeAndFsToMaxSfbLong (m_bitRateMode, samplingRate)]);
        const int16_t  steAnaStats = m_specAnalyzer.stereoSigAnalysis (m_mdctSignals[ci], m_mdctSignals[ci + 1],
                                                                       m_mdstSignals[ci], m_mdstSignals[ci + 1], nSamplesMax,
                                                                       nSamplesInFrame, eightShorts, coreConfig.stereoDataCurr);
        if (steAnaStats == SHRT_MIN) errorValue = 1;

        if ((s = abs (steAnaStats)) * m_perCorrHCurr[el] == 0) // transition to/from silence
        {
          m_perCorrHCurr[el] = uint8_t ((32 + s * __min (64, eightTimesSqrt256Minus[meanSpecFlat])) >> 6);
        }
        else // gentle overlap length dependent temporal smoothing
        {
          const int16_t allowedDiff = (coreConfig.icsInfoCurr[0].windowSequence < EIGHT_SHORT ? 16 : 32);
          const int16_t prevPerCorr = __max (128, __min (192, m_perCorrHCurr[el]));
          const int16_t currPerCorr = (32 + s * __min (64, eightTimesSqrt256Minus[meanSpecFlat])) >> 6;

          m_perCorrHCurr[el] = (uint8_t) __max (prevPerCorr - allowedDiff, __min (prevPerCorr + allowedDiff, currPerCorr));
        }
        m_perCorrLCurr[el] = coreConfig.stereoDataCurr[0];

        if ((int) s == steAnaStats * -1) coreConfig.stereoConfig = 2;  // 2: S>M, pred_dir=1
        if (s > (UCHAR_MAX * 3) / 4) coreConfig.stereoMode = 2; // 2: all, ms_mask_present=2
        if (s >= UCHAR_MAX - 1) coreConfig.stereoConfig |= 8; // true: mono-in-stereo signal
      }
      else if (nrChannels > 1) m_perCorrHCurr[el] = m_perCorrLCurr[el] = 128; // "mid" value

      for (unsigned ch = 0; ch < nrChannels; ch++) // channel loop
      {
        SfbGroupData& grpData = coreConfig.groupingData[ch];
        uint16_t*  grpSO = grpData.sfbOffsets;
        IcsInfo& icsCurr = coreConfig.icsInfoCurr[ch];
        TnsData& tnsData = coreConfig.tnsData[ch];

        memset (grpSO, 0, (1 + MAX_NUM_SWB_SHORT * NUM_WINDOW_GROUPS) * sizeof (uint16_t));

        if (icsCurr.windowSequence != EIGHT_SHORT)
        {
          memcpy (grpSO, swbOffsetsL[m_swbTableIdx], numSwbOffsetL[m_swbTableIdx] * sizeof (uint16_t));

          icsCurr.maxSfb = 0;
          while (grpSO[icsCurr.maxSfb] < nSamplesInFrame) icsCurr.maxSfb++;  // num_swb_long
          grpSO[icsCurr.maxSfb] = (uint16_t) nSamplesInFrame;
          grpData.sfbsPerGroup = m_numSwbLong = icsCurr.maxSfb;  // changed to max_sfb later

          if (samplingRate > 32000) // set max_sfb based on VBR mode and bandwidth detection
          {
            if (icsCurr.maxSfb > 49) // may still be 51 for 32 kHz
            {
              grpData.sfbsPerGroup = m_numSwbLong = icsCurr.maxSfb = 49; // fix 44.1, 48 kHz
            }
            icsCurr.maxSfb = __min (icsCurr.maxSfb, brModeAndFsToMaxSfbLong (m_bitRateMode, samplingRate));
          }
          while (grpSO[icsCurr.maxSfb] > __max (m_bandwidCurr[ci], m_bandwidPrev[ci])) icsCurr.maxSfb--; // BW detector
        }
        else // icsCurr.windowSequence == EIGHT_SHORT
        {
          memcpy (grpSO, swbOffsetsS[m_swbTableIdx], numSwbOffsetS[m_swbTableIdx] * sizeof (uint16_t));

          icsCurr.maxSfb = 0;
          while (grpSO[icsCurr.maxSfb] < nSamplesInShort) icsCurr.maxSfb++; // num_swb_short
          grpSO[icsCurr.maxSfb] = (uint16_t) nSamplesInShort;
          grpData.sfbsPerGroup = m_numSwbShort = icsCurr.maxSfb; // changed to max_sfb later

          if (samplingRate > 32000) // set max_sfb based on VBR mode and zero-ness detection
          {
            icsCurr.maxSfb = __min (icsCurr.maxSfb, brModeAndFsToMaxSfbShort (m_bitRateMode, samplingRate));
          }
#if SA_OPT_WINDOW_GROUPING
          if (ch > 0 && coreConfig.commonWindow)  // resynchronize the scale_factor_grouping
          {
            if (icsCurr.windowGrouping != coreConfig.icsInfoCurr[0].windowGrouping)
            {
              icsCurr.windowGrouping = coreConfig.icsInfoCurr[0].windowGrouping;
            }
          }
          else // first element channel or not common_window, optimize scale_factor_grouping
          {
            if ((s = m_specAnalyzer.optimizeGrouping (ci, grpSO[icsCurr.maxSfb] << 3, icsCurr.windowGrouping)) < 8)
            {
              icsCurr.windowGrouping = (uint8_t) s;
            }
          }
          memcpy (grpData.windowGroupLength, windowGroupingTable[icsCurr.windowGrouping], NUM_WINDOW_GROUPS * sizeof (uint8_t));
#endif
          findActualBandwidthShort (&icsCurr.maxSfb, grpSO, m_mdctSignals[ci], nChannels < 2 ? nullptr : m_mdstSignals[ci], nSamplesInShort);

          errorValue |= eightShortGrouping (grpData, grpSO, m_mdctSignals[ci], nChannels < 2 ? nullptr : m_mdstSignals[ci]);
        } // if EIGHT_SHORT

        // compute and quantize optimal TNS coefficients, then find optimal TNS filter order
        s /*max. pred gain*/ = getOptParCorCoeffs (grpData, icsCurr.maxSfb, tnsData, ci,
                                                   ch > 0 && coreConfig.commonWindow ? coreConfig.tnsData[0].firstTnsWindow : 0);
        for (uint16_t n = 0, gr = 0; gr < grpData.numWindowGroups; gr++)
        {
          if (grpData.windowGroupLength[gr] == 1)
          {
            const uint8_t tonality = (m_specAnaCurr[ci] >> 16) & UCHAR_MAX;

            tnsData.filterOrder[n] = m_linPredictor.calcOptTnsCoeffs (tnsData.coeffParCor[n], tnsData.coeff[n], &tnsData.coeffResLow[n],
                                                                      tnsData.filterOrder[n], s, tonality >> (m_tempFlatPrev[ci] >> 5));
            tnsData.numFilters[n] = (tnsData.filterOrder[n] > 0 ? 1 : 0);
            if ((ch == 0) && (icsCurr.windowSequence == EIGHT_SHORT) && (tnsData.numFilters[n] == 0) && (tnsData.firstTnsWindow == gr))
            {
              tnsData.firstTnsWindow++; // simplify TNS stereo synching in eight-short frame
            }
            n++;
          }
        }
        m_tempFlatPrev[ci++] = (uint8_t) s;
      } // for ch

      if (coreConfig.commonWindow) // synchronization of all StereoCoreToolInfo() components
      {
        uint8_t& maxSfb0 = coreConfig.icsInfoCurr[0].maxSfb;
        uint8_t& maxSfb1 = coreConfig.icsInfoCurr[1].maxSfb;
        const uint8_t maxSfbSte = __max (maxSfb0, maxSfb1);   // max_sfb_ste, as in Table 24

        if ((maxSfb0 > 0) && (maxSfb1 > 0) && (maxSfbSte - __min (maxSfb0, maxSfb1) <= 1 || coreConfig.stereoMode > 0))
        {
          uint32_t& sac0 = m_specAnaCurr[ci-2];
          uint32_t& sac1 = m_specAnaCurr[ci-1];
          TnsData&  tnsData0 = coreConfig.tnsData[0];
          TnsData&  tnsData1 = coreConfig.tnsData[1];
          const int specFlat[2] = {int (sac0 >> 16) & UCHAR_MAX, int (sac1 >> 16) & UCHAR_MAX};
          const int tnsStart[2] = {int (sac0 & 31), int (sac1 & 31)}; // long TNS start band

          if ((abs (specFlat[0] - specFlat[1]) <= (UCHAR_MAX >> 3)) &&
              (abs (tnsStart[0] - tnsStart[1]) <= (UCHAR_MAX >> 4)))  // TNS synchronization
          {
            coreConfig.commonTnsData = true;
            for (uint16_t n = 0; n < 3; n++)
            {
              if ((s = __max (tnsData0.filterOrder[n], tnsData1.filterOrder[n])) == 0) continue;

              if ((coreConfig.stereoMode > 0) || m_linPredictor.similarParCorCoeffs (tnsData0.coeffParCor[n], tnsData1.coeffParCor[n], s, LP_DEPTH))
              {
                applyTnsCoeff2ChannelSynch (m_linPredictor, tnsData0, tnsData1, s, n, &coreConfig.commonTnsData);
              }
              else if ((m_bitRateMode <= 4) && (m_perCorrHCurr[el] > 128))
              {
                applyTnsCoeffPreProcessing (m_linPredictor, tnsData0, tnsData1, s, n, &coreConfig.commonTnsData, m_perCorrHCurr[el] - 128);
              }
              else coreConfig.commonTnsData = false;
            }

            if (coreConfig.commonTnsData || (abs (tnsStart[0] - tnsStart[1]) <= (UCHAR_MAX >> 5)))
            {
              const uint32_t avgTnsStart = (tnsStart[0] + tnsStart[1]) >> 1;  // synch start

              sac0 = (sac0 & (UINT_MAX - 31)) | avgTnsStart; // used by applyTnsToWinGroup()
              sac1 = (sac1 & (UINT_MAX - 31)) | avgTnsStart;
            }
          }
          maxSfb0 = maxSfb1 = maxSfbSte;

          if ((m_bitRateMode <= 4) && (coreConfig.icsInfoCurr[0].windowSequence == EIGHT_SHORT))
          {
            m_perCorrLCurr[el] = stereoCorrGrouping (coreConfig.groupingData[0], nSamplesInFrame, coreConfig.stereoDataCurr);
          }
        }
        else coreConfig.stereoMode = 0;  // since a max_sfb is 0

        coreConfig.commonMaxSfb = (maxSfb0 == maxSfb1); // synch
      } // if coreConfig.commonWindow
    }

    ci -= nrChannels; // zero frequency coefficients above num_swb for all channels, windows

    for (unsigned ch = 0; ch < nrChannels; ch++) // channel loop
    {
      SfbGroupData&  grpData = coreConfig.groupingData[ch];
      const uint16_t*  grpSO = grpData.sfbOffsets;
      const IcsInfo& icsCurr = coreConfig.icsInfoCurr[ch];
      const bool eightShorts = (icsCurr.windowSequence == EIGHT_SHORT);

      if (eightShorts) // map grouping table idx to scale_factor_grouping idx for bit-stream
      {
        coreConfig.icsInfoCurr[ch].windowGrouping = scaleFactorGrouping[icsCurr.windowGrouping];
      }
      s = 0;
      for (uint16_t gr = 0; gr < grpData.numWindowGroups; gr++)
      {
        const unsigned grMax = grpSO[grpData.sfbsPerGroup + m_numSwbShort * gr];

        s += (eightShorts ? nSamplesInShort : nSamplesInFrame) * grpData.windowGroupLength[gr];
        memset (&m_mdctSignals[ci][grMax], 0, (s - grMax) * sizeof (int32_t));
        memset (&m_mdstSignals[ci][grMax], 0, (s - grMax) * sizeof (int32_t));
      }
      memset (grpData.sfbRmsValues, 0, (MAX_NUM_SWB_SHORT * NUM_WINDOW_GROUPS) * sizeof (uint32_t));

      if (icsCurr.maxSfb > 0)
      {
        // use MCLTs for LONG but only MDCTs for SHORT windows when the MDSTs aren't grouped
        const uint8_t* nFilters = coreConfig.tnsData[ch].numFilters;
        const bool realOnlyCalc = (eightShorts && nChannels < 2);

        for (uint8_t n = 0, gr = 0; gr < grpData.numWindowGroups; gr++)
        {
          if (grpData.windowGroupLength[gr] == 1)
          {
            errorValue |= applyTnsToWinGroup (grpData, gr, grpData.sfbsPerGroup, coreConfig.tnsData[ch], ci, n, realOnlyCalc);
            coreConfig.tnsActive |= (nFilters[n++] > 0); // set tns_data_present, tns_active
          }
          if ((grpData.windowGroupLength[gr] > 1) || (nFilters[n - 1] == 0))
          {
            s = m_numSwbShort * gr;
            errorValue |= m_specAnalyzer.getMeanAbsValues (m_mdctSignals[ci], realOnlyCalc ? nullptr : m_mdstSignals[ci],
                                                           grpSO[grpData.sfbsPerGroup + s], (eightShorts ? USAC_MAX_NUM_CHANNELS : ci),
                                                           &grpSO[s], grpData.sfbsPerGroup, &grpData.sfbRmsValues[s]);
          }
        }
      }

      grpData.sfbsPerGroup = icsCurr.maxSfb; // change num_swb to max_sfb for coding process
      ci++;
    }
  } // for el

  return errorValue;
}

unsigned ExhaleEncoder::temporalProcessing () // determine time-domain aspects of ics_info()
{
  const unsigned nChannels       = toNumChannels (m_channelConf);
  const unsigned nSamplesInFrame = toFrameLength (m_frameLength) << m_shiftValSBR;
  const unsigned nSamplesTempAna = (nSamplesInFrame * 25) >> 4;  // pre-delay for look-ahead
  const unsigned lfeChannelIndex = (m_channelConf >= CCI_6_CH ? __max (5, nChannels - 1) : USAC_MAX_NUM_CHANNELS);
  unsigned ci = 0; // running ch index
  unsigned errorValue = 0; // no error

  // get temporal channel statistics for this frame, used for spectral grouping/quantization
  m_tempAnalyzer.getTempAnalysisStats (m_tempAnaCurr, nChannels);
  m_tempAnalyzer.getTransientAndPitch (m_tranLocCurr, nChannels);

  // temporal analysis for look-ahead signal (central nSamplesInFrame samples of next frame)
  errorValue |= m_tempAnalyzer.temporalAnalysis (m_timeSignals, nChannels, nSamplesInFrame, nSamplesTempAna,
                                                 m_shiftValSBR, m_coreSignals, lfeChannelIndex);
  // get temporal channel statistics for next frame, used for window length/overlap decision
  m_tempAnalyzer.getTempAnalysisStats (m_tempAnaNext, nChannels);
  m_tempAnalyzer.getTransientAndPitch (m_tranLocNext, nChannels);

#ifdef NO_PREROLL_DATA
  m_indepFlag = (((m_frameCount++) % m_indepPeriod) == 0); // configure usacIndependencyFlag
#else
  m_indepFlag = (((m_frameCount++) % m_indepPeriod) <= 1); // configure usacIndependencyFlag
#endif

  for (unsigned el = 0; el < m_numElements; el++)  // element loop
  {
    CoreCoderData& coreConfig = *m_elementData[el];
    const unsigned nrChannels = (coreConfig.elementType & 1) + 1; // for UsacCoreCoderData()

    coreConfig.commonWindow   = false;
    coreConfig.icsInfoPrev[0] = coreConfig.icsInfoCurr[0];
    coreConfig.icsInfoPrev[1] = coreConfig.icsInfoCurr[1];

    if (coreConfig.elementType >= ID_USAC_LFE) // LFE/EXT elements
    {
      IcsInfo& icsCurr = coreConfig.icsInfoCurr[0];

      icsCurr.windowGrouping  = 0;
      icsCurr.windowSequence  = ONLY_LONG;
#if RESTRICT_TO_AAC
      icsCurr.windowShape     = WINDOW_SINE;
#else
      icsCurr.windowShape     = WINDOW_KBD;
#endif
      ci++;
    }
    else // SCE or CPE: short-window, low-overlap, and sine-shape detection for each channel
    {
      unsigned tsCurr[2]; // save temporal stationarity values
      unsigned tsNext[2]; // for common_window decision in CPE

      for (unsigned ch = 0; ch < nrChannels; ch++) // channel loop
      {
        const IcsInfo& icsPrev = coreConfig.icsInfoPrev[ch];
              IcsInfo& icsCurr = coreConfig.icsInfoCurr[ch];
        const USAC_WSEQ wsPrev = icsPrev.windowSequence;
             USAC_WSEQ& wsCurr = icsCurr.windowSequence;
        // get temporal signal statistics, then determine overlap config. for the next frame
        const unsigned  plCurr = abs (m_tranLocCurr[ci]) & ((1024 << m_shiftValSBR) - 1);
        const unsigned  sfCurr = (m_tempAnaCurr[ci] >> 24) & UCHAR_MAX;
        const unsigned  tfCurr = (m_tempAnaCurr[ci] >> 16) & UCHAR_MAX;
        const unsigned  plNext = abs (m_tranLocNext[ci]) & ((1024 << m_shiftValSBR) - 1);
        const unsigned  sfNext = (m_tempAnaNext[ci] >> 24) & UCHAR_MAX;
        const unsigned  tfNext = (m_tempAnaNext[ci] >> 16) & UCHAR_MAX;
        const unsigned tThresh = UCHAR_MAX * (__max (plCurr, plNext) < 614 /*0.6 * 1024*/ ? 16 : 15 - (m_bitRateMode >> 3));

        tsCurr[ch] = (m_tempAnaCurr[ci] /*R*/) & UCHAR_MAX;
        tsNext[ch] = (m_tempAnaNext[ci] >>  8) & UCHAR_MAX;
        // save maximum spectral flatness of current and neighboring frames for quantization
        m_tempAnaCurr [ci] = (m_tempAnaCurr[ci] & 0xFFFFFF) | (__max (sfCurr, __max (m_specFlatPrev[ci], sfNext)) << 24);
        m_specFlatPrev[ci] = (uint8_t) sfCurr;

        const bool lowOlapNext = (m_tranLocNext[ci] >= 0) || (sfNext <= UCHAR_MAX / 4 && tfNext > (UCHAR_MAX * 13) / 16) ||
                                 (tsCurr[ch] > (UCHAR_MAX * 5) / 8) || (tsNext[ch] > (UCHAR_MAX * 5) / 8);
        const bool sineWinCurr = (sfCurr >= 170) && (sfNext >= 170) && (sfCurr < 221) && (sfNext < 221) && (tsCurr[ch] < 20) &&
                                 (tfCurr >= 153) && (tfNext >= 153) && (tfCurr < 184) && (tfNext < 184) && (tsNext[ch] < 20);
        // set window_sequence
        if ((wsPrev == ONLY_LONG) || (wsPrev == LONG_STOP)) // 1st window half - max overlap
        {
          wsCurr = (lowOlapNext ? LONG_START : ONLY_LONG);
        }
        else // LONG_START_SEQUENCE, STOP_START_SEQUENCE, EIGHT_SHORT_SEQUENCE - min overlap
        {
          wsCurr = (m_tranLocCurr[ci] >= 0) || (tsCurr[ch] > (UCHAR_MAX * 5) / 8) || (tfCurr > tThresh / 16) ? EIGHT_SHORT :
#if RESTRICT_TO_AAC
                   (lowOlapNext ? EIGHT_SHORT : LONG_STOP);
#else
                   (lowOlapNext ? STOP_START : LONG_STOP);
#endif
        }

        // set window_shape
        if ((wsCurr == ONLY_LONG) || (wsCurr == LONG_STOP)) // 2nd window half - max overlap
        {
          icsCurr.windowShape  = (sineWinCurr ? WINDOW_SINE : WINDOW_KBD);
        }
        else // LONG_START_SEQUENCE, STOP_START_SEQUENCE, EIGHT_SHORT_SEQUENCE - min overlap
        {
          icsCurr.windowShape  = (m_tranLocCurr[ci] >= 0) ? WINDOW_KBD :
                                 (sineWinCurr ? WINDOW_SINE : WINDOW_KBD);
        }

        // set scale_factor_grouping
        icsCurr.windowGrouping = (wsCurr == EIGHT_SHORT ? __max (0, m_tranLocCurr[ci]) / (2 * nSamplesInFrame) : 0);
        ci++;
      } // for ch

      if (nrChannels > 1) // common_window element detection for use in StereoCoreToolInfo()
      {
        IcsInfo&  icsInfo0 = coreConfig.icsInfoCurr[0];
        IcsInfo&  icsInfo1 = coreConfig.icsInfoCurr[1];
        USAC_WSEQ& winSeq0 = icsInfo0.windowSequence;
        USAC_WSEQ& winSeq1 = icsInfo1.windowSequence;

        if (winSeq0 != winSeq1) // try to synch window_sequences
        {
          const USAC_WSEQ initialWs0 = winSeq0;
          const USAC_WSEQ initialWs1 = winSeq1;

          winSeq0 = winSeq1 = windowSequenceSynch[initialWs0][initialWs1];   // equalization
          if ((winSeq0 != initialWs0) && (winSeq0 == EIGHT_SHORT))
          {
#if !RESTRICT_TO_AAC
            if ((tsCurr[0] * 7 < tsCurr[1] * 2) && (tsNext[0] * 7 < tsNext[1] * 2) &&
                (abs (m_specFlatPrev[ci - 1] - (int) m_specFlatPrev[ci - 2]) > UCHAR_MAX / 4))
            {
              winSeq0 = STOP_START; // don't synchronize to EIGHT_SHORT but keep low overlap
            }
            else
#endif
            icsInfo0.windowGrouping = icsInfo1.windowGrouping;
          }
          if ((winSeq1 != initialWs1) && (winSeq1 == EIGHT_SHORT))
          {
#if !RESTRICT_TO_AAC
            if ((tsCurr[1] * 7 < tsCurr[0] * 2) && (tsNext[1] * 7 < tsNext[0] * 2) &&
                (abs (m_specFlatPrev[ci - 1] - (int) m_specFlatPrev[ci - 2]) > UCHAR_MAX / 4))
            {
              winSeq1 = STOP_START; // don't synchronize to EIGHT_SHORT but keep low overlap
            }
            else
#endif
            icsInfo1.windowGrouping = icsInfo0.windowGrouping;
          }
        }
        else if (winSeq0 == EIGHT_SHORT) // resynchronize scale_factor_grouping if necessary
        {
          const int16_t tranLocSynch = __min (m_tranLocCurr[ci - 2], m_tranLocCurr[ci - 1]);

          icsInfo0.windowGrouping = icsInfo1.windowGrouping = __max (0, tranLocSynch) / (2 * nSamplesInFrame);
        }

        if ((icsInfo0.windowShape != WINDOW_SINE) || (icsInfo1.windowShape != WINDOW_SINE))
        {
          icsInfo0.windowShape = WINDOW_KBD; // always synchronize window_shapes in order to
          icsInfo1.windowShape = WINDOW_KBD; // encourage synch in next frame; KBD dominates
        }
        coreConfig.commonWindow = (winSeq0 == winSeq1); // synch

        memset (coreConfig.stereoDataPrev, 0, (MAX_NUM_SWB_LONG + 1) * sizeof (uint8_t));

        if (((winSeq0 == EIGHT_SHORT) == (coreConfig.icsInfoPrev[0].windowSequence == EIGHT_SHORT)) && !m_indepFlag &&
            ((winSeq1 == EIGHT_SHORT) == (coreConfig.icsInfoPrev[1].windowSequence == EIGHT_SHORT)) && (coreConfig.stereoMode > 0))
        {
          const unsigned lastGrpOffset = (coreConfig.icsInfoPrev[0].windowSequence == EIGHT_SHORT ? m_numSwbShort * (NUM_WINDOW_GROUPS - 1) : 0);
          const unsigned maxSfbStePrev = __max (coreConfig.icsInfoPrev[0].maxSfb, coreConfig.icsInfoPrev[1].maxSfb) + 1u; // for safety

          memcpy (coreConfig.stereoDataPrev, &coreConfig.stereoDataCurr[lastGrpOffset], __min (60 - lastGrpOffset, maxSfbStePrev) * sizeof (uint8_t));
        }
      } // if nrChannels > 1
    }

    ci -= nrChannels; // modulated complex lapped transform (MCLT) for all channels, windows

    for (unsigned ch = 0; ch < nrChannels; ch++) // channel loop
    {
      const IcsInfo& icsPrev = coreConfig.icsInfoPrev[ch];
      const IcsInfo& icsCurr = coreConfig.icsInfoCurr[ch];
      const int32_t* timeSig = (m_shiftValSBR > 0 ? m_coreSignals[ci] : m_timeSignals[ci]);
      const USAC_WSEQ wsCurr = icsCurr.windowSequence;
      const bool eightShorts = (wsCurr == EIGHT_SHORT);
      SfbGroupData&  grpData = coreConfig.groupingData[ch];

      grpData.numWindowGroups = (eightShorts ? NUM_WINDOW_GROUPS : 1);  // fill groupingData
      memcpy (grpData.windowGroupLength, windowGroupingTable[icsCurr.windowGrouping], NUM_WINDOW_GROUPS * sizeof (uint8_t));

      errorValue |= m_transform.applyMCLT (timeSig, eightShorts, icsPrev.windowShape != WINDOW_SINE, icsCurr.windowShape != WINDOW_SINE,
                                           wsCurr > LONG_START /*lOL*/, (wsCurr % 3) != ONLY_LONG /*lOR*/, m_mdctSignals[ci], m_mdstSignals[ci]);
      m_scaleFacData[ci++] = &grpData;
    }
  } // for el

  return errorValue;
}

// constructor
ExhaleEncoder::ExhaleEncoder (int32_t* const inputPcmData,           unsigned char* const outputAuData,
                              const unsigned sampleRate /*= 44100*/, const unsigned numChannels /*= 2*/,
                              const unsigned frameLength /*= 1024*/, const unsigned indepPeriod /*= 45*/,
                              const unsigned varBitRateMode /*= 3*/
#if !RESTRICT_TO_AAC
                            , const bool useNoiseFilling /*= true*/, const bool useEcodisExt /*= false*/
#endif
                              )
{
  // adopt basic coding parameters
  m_bitRateMode  = __min (9, varBitRateMode);
  m_channelConf  = (numChannels >= 7 ? CCI_UNDEF : (USAC_CCI) numChannels); // see 23003-3, Tables 73 & 161
  if (m_channelConf == CCI_CONF)
  {
    m_channelConf = CCI_2_CHM; // passing numChannels = 0 to ExhaleEncoder is interpreted as 2-ch dual-mono
  }
  m_numElements  = elementCountConfig[m_channelConf % USAC_MAX_NUM_ELCONFIGS]; // used in UsacDecoderConfig
  m_shiftValSBR  = (frameLength >= 1536 ? 1 : 0);
  m_frameCount   = 0;
  m_frameLength  = USAC_CCFL (frameLength >> m_shiftValSBR); // ccfl signaled using coreSbrFrameLengthIndex
  m_frequencyIdx = toSamplingFrequencyIndex (sampleRate >> m_shiftValSBR); // as usacSamplingFrequencyIndex
  m_indepFlag    = true; // usacIndependencyFlag in UsacFrame(), will be set per frame, true in first frame
  m_indepPeriod  = (indepPeriod == 0 ? USHRT_MAX : __min (USHRT_MAX, indepPeriod)); // random-access period
#if !RESTRICT_TO_AAC
  m_nonMpegExt   = useEcodisExt;
#endif
  m_numSwbLong   = MAX_NUM_SWB_LONG;
  m_numSwbShort  = MAX_NUM_SWB_SHORT;
  m_outAuData    = outputAuData;
  m_pcm24Data    = inputPcmData;
  m_tempIntBuf   = nullptr;

  // initialize all helper structs
  for (unsigned el = 0; el < USAC_MAX_NUM_ELEMENTS; el++)
  {
    const ELEM_TYPE et = elementTypeConfig[m_channelConf % USAC_MAX_NUM_ELCONFIGS][el];  // usacElementType

    m_elementData[el]  = nullptr;
    m_perCorrHCurr[el] = 0;
    m_perCorrLCurr[el] = 0;
#if !RESTRICT_TO_AAC
    m_noiseFilling[el] = (useNoiseFilling && (et < ID_USAC_LFE));
    m_timeWarping[el]  = (false /* N/A */ && (et < ID_USAC_LFE));
#endif
  }
  // initialize all signal buffers
  for (unsigned ch = 0; ch < USAC_MAX_NUM_CHANNELS; ch++)
  {
    m_bandwidCurr[ch]  = 0;
    m_bandwidPrev[ch]  = 0;
    m_coreSignals[ch]  = nullptr;
    m_mdctQuantMag[ch] = nullptr;
    m_mdctSignals[ch]  = nullptr;
    m_mdstSignals[ch]  = nullptr;
    m_meanSpecPrev[ch] = 0;
    m_meanTempPrev[ch] = 0;
    m_scaleFacData[ch] = nullptr;
    m_specAnaCurr[ch]  = 0;
    m_specFlatPrev[ch] = 0;
    m_tempAnaCurr[ch]  = 0;
    m_tempAnaNext[ch]  = 0;
    m_tempFlatPrev[ch] = 0;
    m_timeSignals[ch]  = nullptr;
    m_tranLocCurr[ch]  = -1;
    m_tranLocNext[ch]  = -1;
  }
  // initialize all window buffers
  for (unsigned ws = WINDOW_SINE; ws <= WINDOW_KBD; ws++)
  {
    m_timeWindowL[ws] = nullptr;
    m_timeWindowS[ws] = nullptr;
  }
}

// destructor
ExhaleEncoder::~ExhaleEncoder ()
{
  // free allocated helper structs
  for (unsigned el = 0; el < USAC_MAX_NUM_ELEMENTS; el++)
  {
    MFREE (m_elementData[el]);
  }
  // free allocated signal buffers
  for (unsigned ch = 0; ch < USAC_MAX_NUM_CHANNELS; ch++)
  {
    if (m_shiftValSBR > 0) MFREE (m_coreSignals[ch]);
    MFREE (m_mdctQuantMag[ch]);
    MFREE (m_mdctSignals[ch]);
    MFREE (m_mdstSignals[ch]);
    MFREE (m_timeSignals[ch]);
  }
  // free allocated window buffers
  for (unsigned ws = WINDOW_SINE; ws <= WINDOW_KBD; ws++)
  {
    MFREE (m_timeWindowL[ws]);
    MFREE (m_timeWindowS[ws]);
  }
  // execute sub-class destructors
}

// public functions
unsigned ExhaleEncoder::encodeLookahead ()
{
  const unsigned nChannels       = toNumChannels (m_channelConf);
  const unsigned nSamplesInFrame = toFrameLength (m_frameLength) << m_shiftValSBR;
  const unsigned nSamplesTempAna = (nSamplesInFrame * 25) >> 4;  // pre-delay for look-ahead
  const int32_t* chSig           = m_pcm24Data;
  unsigned ch, s;

  // copy nSamplesInFrame external channel-interleaved samples into internal channel buffers
  for (s = 0; s < nSamplesInFrame; s++) // sample loop
  {
    for (ch = 0; ch < nChannels; ch++) // channel loop
    {
      m_timeSignals[ch][nSamplesTempAna + s] = *(chSig++);
    }
  }

  // generate first nSamplesTempAna deinterleaved samples (previous frame data) by LP filter
  for (ch = 0; ch < nChannels; ch++)
  {
    short filterC[MAX_PREDICTION_ORDER] = {0, 0, 0, 0};
    short parCorC[MAX_PREDICTION_ORDER] = {0, 0, 0, 0};
    int32_t* predSig = &m_timeSignals[ch][nSamplesTempAna]; // end of signal to be predicted

    m_linPredictor.calcParCorCoeffs (predSig, uint16_t (nSamplesInFrame >> 1), MAX_PREDICTION_ORDER, parCorC);
    m_linPredictor.parCorToLpCoeffs (parCorC, MAX_PREDICTION_ORDER, filterC);

    for (s = nSamplesTempAna; s > 0; s--) // generate prediction signal without limit cycles
    {
      const int64_t predSample = *(predSig + 0) * (int64_t) filterC[0] + *(predSig + 1) * (int64_t) filterC[1] +
                                 *(predSig + 2) * (int64_t) filterC[2] + *(predSig + 3) * (int64_t) filterC[3];
      *(--predSig) = int32_t ((predSample > 0 ? -predSample + (1 << 9) - 1 : -predSample) >> 9);
    }
    if (m_shiftValSBR > 0) memset (m_coreSignals[ch], 0, ((nSamplesInFrame * 41) >> (4 + m_shiftValSBR)) * sizeof (int32_t));
  }

  // set initial temporal channel statistic to something meaningful before first coded frame
  m_tempAnalyzer.temporalAnalysis (m_timeSignals, nChannels, nSamplesInFrame, nSamplesTempAna - nSamplesInFrame,
                                   m_shiftValSBR, m_coreSignals); // default lfeChannelIndex
  if (temporalProcessing ()) // time domain: window length, overlap, grouping, and transform
  {
    return 2; // internal error in temporal processing
  }
  if (spectralProcessing ()) // MCLT domain: (common_)max_sfb, grouping 2, TNS, and SFB data
  {
    return 2; // internal error in spectral processing
  }
  if (psychBitAllocation ()) // SFB domain: psychoacoustic model and scale factor estimation
  {
    return 1; // internal error in bit-allocation code
  }

  return quantizationCoding (); // max(3, coded bytes)
}

unsigned ExhaleEncoder::encodeFrame ()
{
  const unsigned nChannels       = toNumChannels (m_channelConf);
  const unsigned nSamplesInFrame = toFrameLength (m_frameLength) << m_shiftValSBR;
  const unsigned nSamplesTempAna = (nSamplesInFrame * 25) >> 4;  // pre-delay for look-ahead
  const int32_t* chSig           = m_pcm24Data;
  unsigned ch, s;

  // move internal channel buffers nSamplesInFrame to the past to make room for next samples
  for (ch = 0; ch < nChannels; ch++)
  {
    memcpy (&m_timeSignals[ch][0], &m_timeSignals[ch][nSamplesInFrame], nSamplesInFrame * sizeof (int32_t));
    memcpy (&m_timeSignals[ch][nSamplesInFrame], &m_timeSignals[ch][2 * nSamplesInFrame], (nSamplesTempAna - nSamplesInFrame) * sizeof (int32_t));

    if (m_shiftValSBR > 0)
    {
      const unsigned nSmpInFrame = toFrameLength (m_frameLength); // core coder frame length

      memcpy (&m_coreSignals[ch][0], &m_coreSignals[ch][nSmpInFrame], nSmpInFrame * sizeof (int32_t));
      memcpy (&m_coreSignals[ch][nSmpInFrame], &m_coreSignals[ch][2 * nSmpInFrame], (nSamplesInFrame >> 2) * sizeof (int32_t));
    }
  }

  // copy nSamplesInFrame external channel-interleaved samples into internal channel buffers
  for (s = 0; s < nSamplesInFrame; s++) // sample loop
  {
    for (ch = 0; ch < nChannels; ch++) // channel loop
    {
      m_timeSignals[ch][nSamplesTempAna + s] = *(chSig++);
    }
  }

  if (temporalProcessing ()) // time domain: window length, overlap, grouping, and transform
  {
    return 2; // internal error in temporal processing
  }
  if (spectralProcessing ()) // MCLT domain: (common_)max_sfb, grouping 2, TNS, and SFB data
  {
    return 2; // internal error in spectral processing
  }
  if (psychBitAllocation ()) // SFB domain: psychoacoustic model and scale factor estimation
  {
    return 1; // internal error in bit-allocation code
  }

  return quantizationCoding (); // max(3, coded bytes)
}

unsigned ExhaleEncoder::initEncoder (unsigned char* const audioConfigBuffer, uint32_t* const audioConfigBytes /*= nullptr*/)
{
  const unsigned nChannels       = toNumChannels (m_channelConf);
  const unsigned nSamplesInFrame = toFrameLength (m_frameLength);
  const unsigned specSigBufSize  = nSamplesInFrame * sizeof (int32_t);
  const unsigned timeSigBufSize  = (((nSamplesInFrame << m_shiftValSBR) * 41) >> 4) * sizeof (int32_t); // core-codec delay*4
  const unsigned char chConf     = m_channelConf;
  unsigned errorValue = 0; // no error

  // check user's input parameters
#if RESTRICT_TO_AAC
  if ((m_channelConf <= CCI_CONF) || (m_channelConf > CCI_8_CH))
#else
  if ((m_channelConf <= CCI_CONF) || (m_channelConf > CCI_8_CHS))
#endif
  {
    errorValue |= 128;
  }
#if RESTRICT_TO_AAC
  if (m_frameLength != CCFL_1024)
#else
  if ((m_frameLength != CCFL_768) && (m_frameLength != CCFL_1024))
#endif
  {
    errorValue |=  64;
  }
  if ((m_frequencyIdx < 0) || (m_bitRateMode > (toSamplingRate (m_frequencyIdx) >> (m_shiftValSBR > 0 ? 11 : 12)) + 2))
  {
    errorValue |=  32;
  }
  if ((m_outAuData == nullptr) || (m_pcm24Data == nullptr))
  {
    errorValue |=  16;
  }
  if (errorValue > 0) return errorValue;

  // get window band table index
  errorValue = (unsigned) m_frequencyIdx; // for temporary storage
#if RESTRICT_TO_AAC
  m_swbTableIdx = freqIdxToSwbTableIdxAAC[errorValue];
#else
  m_swbTableIdx = (m_frameLength == CCFL_768 ? freqIdxToSwbTableIdx768[errorValue] : freqIdxToSwbTableIdxAAC[errorValue]);
#endif
  errorValue = 0;

  if (m_elementData[0] != nullptr) // initEncoder was called before, don't reallocate memory
  {
    if (audioConfigBuffer != nullptr) // recreate the UsacConfig()
    {
      errorValue = m_outStream.createAudioConfig (m_frequencyIdx, m_frameLength != CCFL_1024, chConf, m_numElements,
                                                  elementTypeConfig[chConf], audioConfigBytes ? *audioConfigBytes : 0,
#if !RESTRICT_TO_AAC
                                                  m_timeWarping, m_noiseFilling,
#endif
                                                  m_shiftValSBR, audioConfigBuffer);
      if (audioConfigBytes) *audioConfigBytes = errorValue; // size of UsacConfig() in bytes
      errorValue = (errorValue == 0 ? 1 : 0);
    }

    return errorValue;
  }

  // allocate all helper structs
  for (unsigned el = 0; el < m_numElements; el++)  // element loop
  {
    if ((m_elementData[el] = (CoreCoderData*) malloc (sizeof (CoreCoderData))) == nullptr)
    {
      errorValue |= 8;
    }
    else
    {
      memset (m_elementData[el], 0, sizeof (CoreCoderData));
      m_elementData[el]->elementType = elementTypeConfig[chConf][el]; // usacElementType[el]
    }
  }
#ifndef NO_DTX_MODE
  memset (m_sfbLoudMem, 1, 2 * 26 * 32 * sizeof (uint16_t));
#endif
  // allocate all signal buffers
  if (m_shiftValSBR > 0)
  {
    if (m_shiftValSBR > 1)
    {
      return (errorValue | 4); // >2:1 not supported at the moment
    }
    else for (unsigned ch = 0; ch < nChannels; ch++)
    {
      if ((m_coreSignals[ch] = (int32_t*) malloc (timeSigBufSize >> m_shiftValSBR)) == nullptr)
      {
        errorValue |= 4;
      }
    }
  }
  for (unsigned ch = 0; ch < nChannels; ch++)
  {
    if ((m_entropyCoder[ch].initCodingMemory (nSamplesInFrame) > 0) ||
        (m_mdctQuantMag[ch]= (uint8_t*) malloc (nSamplesInFrame * sizeof (uint8_t))) == nullptr ||
        (m_mdctSignals[ch] = (int32_t*) malloc (specSigBufSize)) == nullptr ||
        (m_mdstSignals[ch] = (int32_t*) malloc (specSigBufSize)) == nullptr ||
        (m_timeSignals[ch] = (int32_t*) malloc (timeSigBufSize)) == nullptr)
    {
      errorValue |= 4;
    }
  }
  // allocate all window buffers
  for (unsigned ws = WINDOW_SINE; ws <= WINDOW_KBD; ws++)
  {
    if ((m_timeWindowL[ws] = initWindowHalfCoeffs ((USAC_WSHP) ws, nSamplesInFrame)) == nullptr ||
        (m_timeWindowS[ws] = initWindowHalfCoeffs ((USAC_WSHP) ws, nSamplesInFrame >> 3)) == nullptr)
    {
      errorValue |= 2;
    }
  }
  if (errorValue > 0) return errorValue;

  // initialize coder class memory
  m_tempIntBuf = m_timeSignals[0];
  if (m_bitAllocator.initAllocMemory (&m_linPredictor, numSwbOffsetL[m_swbTableIdx] - 1, m_bitRateMode >> ((nChannels - 1) >> 2)) > 0 ||
#if EC_TRELLIS_OPT_CODING
      m_sfbQuantizer.initQuantMemory (nSamplesInFrame, numSwbOffsetL[m_swbTableIdx] - 1, m_bitRateMode, toSamplingRate (m_frequencyIdx)) > 0 ||
#else
      m_sfbQuantizer.initQuantMemory (nSamplesInFrame) > 0 ||
#endif
      m_specAnalyzer.initSigAnaMemory (&m_linPredictor, m_bitRateMode <= 4 ? nChannels : 0, nSamplesInFrame) > 0 ||
      m_transform.initConstants (m_tempIntBuf, m_timeWindowL, m_timeWindowS, nSamplesInFrame) > 0)
  {
    errorValue |= 1;
  }

  if ((errorValue == 0) && (audioConfigBuffer != nullptr)) // save UsacConfig() for writeout
  {
    const uint32_t loudnessInfo = (audioConfigBytes ? *audioConfigBytes : 0);
#if 1 //FULL_FRM_LOOKAHEAD
    if (*audioConfigBuffer > 0) m_frameCount--; // to skip 1 frame
#endif
    errorValue = m_outStream.createAudioConfig (m_frequencyIdx, m_frameLength != CCFL_1024, chConf, m_numElements,
                                                elementTypeConfig[chConf], loudnessInfo,
#if !RESTRICT_TO_AAC
                                                m_timeWarping, m_noiseFilling,
#endif
                                                m_shiftValSBR, audioConfigBuffer);
    if (audioConfigBytes) *audioConfigBytes = errorValue; // length of UsacConfig() in bytes
    errorValue = (errorValue == 0 ? 1 : 0);

    // NOTE: In the following, an error value of 256 is actually a warning, not an error. If
    // the exhale library is used for realtime encoding and a nonzero program loudness level
    // is provided before any frames have been encoded, this warning reminds the implementer
    // to apply short-term loudness normalization of the incoming live audio before encoding
    // each frame, preferably to a program level of -23 LUFS and as recommended in EBU R128.
    if ((m_frameCount == 0) && ((loudnessInfo & 16383) > 0)) errorValue |= 256;
  }

  return errorValue;
}
