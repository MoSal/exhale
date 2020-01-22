/* exhaleEnc.cpp - source file for class providing Extended HE-AAC encoding capability
 * written by C. R. Helmrich, last modified in 2020 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under a Modified BSD-Style License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2020 Christian R. Helmrich, project ecodis. All rights reserved.
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

static inline uint32_t getComplexRmsValue (const uint32_t rmsValue, const unsigned sfbGroup, const unsigned sfbIndex,
                                           const uint8_t numSwb, const TnsData& tnsData)
{
  // compensate for missing MDST coefficients in RMS calculation of SFBs where TNS is active
  return ((tnsData.numFilters > 0) && (sfbGroup == tnsData.filteredWindow) && (rmsValue <= UINT_MAX / 3) &&
          (tnsData.filterLength[0] + sfbIndex >= numSwb) ? (rmsValue * 3u) >> 1 : rmsValue);
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
  return numberOfChannels[__max (0, (char) chConfigurationIndex)];
}

// ISO/IEC 23003-3, Table 68
static const uint8_t  elementCountConfig[USAC_MAX_NUM_ELCONFIGS] = {0, 1, 1, 2, 3, 3, 4, 5, 2, 2, 2, 5, 5};

static const ELEM_TYPE elementTypeConfig[USAC_MAX_NUM_ELCONFIGS][USAC_MAX_NUM_ELEMENTS] = {
  {ID_EL_UNDEF, ID_EL_UNDEF, ID_EL_UNDEF, ID_EL_UNDEF, ID_EL_UNDEF}, // CCI_UNDEF
  {ID_USAC_SCE, ID_EL_UNDEF, ID_EL_UNDEF, ID_EL_UNDEF, ID_EL_UNDEF}, // CCI_1_CH
  {ID_USAC_CPE, ID_EL_UNDEF, ID_EL_UNDEF, ID_EL_UNDEF, ID_EL_UNDEF}, // CCI_2_CH
  {ID_USAC_SCE, ID_USAC_CPE, ID_EL_UNDEF, ID_EL_UNDEF, ID_EL_UNDEF}, // CCI_3_CH
  {ID_USAC_SCE, ID_USAC_CPE, ID_USAC_SCE, ID_EL_UNDEF, ID_EL_UNDEF}, // CCI_4_CH
  {ID_USAC_SCE, ID_USAC_CPE, ID_USAC_CPE, ID_EL_UNDEF, ID_EL_UNDEF}, // CCI_5_CH
  {ID_USAC_SCE, ID_USAC_CPE, ID_USAC_CPE, ID_USAC_LFE, ID_EL_UNDEF}, // CCI_6_CH
  {ID_USAC_SCE, ID_USAC_CPE, ID_USAC_CPE, ID_USAC_CPE, ID_USAC_LFE}, // CCI_8_CH
  {ID_USAC_SCE, ID_USAC_SCE, ID_EL_UNDEF, ID_EL_UNDEF, ID_EL_UNDEF}, // CCI_2_CHM
  {ID_USAC_CPE, ID_USAC_SCE, ID_EL_UNDEF, ID_EL_UNDEF, ID_EL_UNDEF}, // CCI_3_CHR
  {ID_USAC_CPE, ID_USAC_CPE, ID_EL_UNDEF, ID_EL_UNDEF, ID_EL_UNDEF}, // CCI_4_CHR
  {ID_USAC_SCE, ID_USAC_CPE, ID_USAC_CPE, ID_USAC_SCE, ID_USAC_LFE}, // CCI_7_CH
  {ID_USAC_SCE, ID_USAC_CPE, ID_USAC_CPE, ID_USAC_CPE, ID_USAC_LFE}  // CCI_8_CHM
};

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
#endif // !RESTRICT_TO_AAC

// ISO/IEC 23003-3, Table 131
static const uint8_t tnsScaleFactorBandLimit[2 /*long/short*/][USAC_NUM_FREQ_TABLES] = { // TNS_MAX_BANDS
#if RESTRICT_TO_AAC
  {31, 34, 51 /*to be corrected to 42 (44.1) and 40 (48 kHz)!*/, 46, 42, 39}, {9, 10, 14, 14, 14, 14}
#else
  {31, 34, 51 /*to be corrected to 42 (44.1) and 40 (48 kHz)!*/, 47, 43, 40}, {9, 10, 14, 15, 15, 15}
#endif
};

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
unsigned ExhaleEncoder::applyTnsToWinGroup (TnsData& tnsData, SfbGroupData& grpData, const bool eightShorts, const uint8_t maxSfb,
                                            const unsigned channelIndex)
{
  const uint16_t filtOrder = tnsData.filterOrder[0];
  const uint16_t*    grpSO = &grpData.sfbOffsets[m_numSwbShort * tnsData.filteredWindow];
  const unsigned nSamplesInFrame = toFrameLength (m_frameLength);
  unsigned errorValue = 0; // no error

  if ((maxSfb > (eightShorts ? 15 : 51)) || (channelIndex >= USAC_MAX_NUM_CHANNELS))
  {
    return 1; // invalid arguments error
  }

  if (filtOrder > 0) // determine TNS filter length in SFBs and apply TNS analysis filtering
  {
    uint8_t numSwbFrame = (eightShorts ? numSwbOffsetS[m_swbTableIdx] : numSwbOffsetL[m_swbTableIdx]) - 1;
    uint8_t tnsMaxBands = tnsScaleFactorBandLimit[eightShorts ? 1 : 0][m_swbTableIdx];
    uint8_t tnsStartSfb = 3 + 32000 / toSamplingRate (m_frequencyIdx);  // 8-short TNS start

    if (!eightShorts)
    {
      const unsigned samplingRate = toSamplingRate (m_frequencyIdx); // refine TNS_MAX_BANDS
      const unsigned tnsStartOffs = (m_specAnaCurr[channelIndex] & 31) << SA_BW_SHIFT;

      if ((samplingRate >= 46009) && (samplingRate < 55426)) // ~48kHz
      {
        numSwbFrame = 49;
        tnsMaxBands = 40;
      }
      else
      if ((samplingRate >= 37566) && (samplingRate < 46009)) // ~44kHz
      {
        numSwbFrame = 49;
        tnsMaxBands = 42;
      }
      while (grpSO[tnsStartSfb] < tnsStartOffs) tnsStartSfb++;  // start band for TNS filter
    }
    tnsMaxBands = __min (tnsMaxBands, maxSfb);

    if ((tnsData.filterLength[0] = __max (0, numSwbFrame - (int) tnsStartSfb)) > 0)
    {
      int32_t* const mdctSignal = m_mdctSignals[channelIndex];
      const short offs = grpSO[tnsStartSfb];
      uint16_t       s = grpSO[tnsMaxBands] - offs;
      short filterC[MAX_PREDICTION_ORDER] = {0, 0, 0, 0};
      int32_t* predSig = &mdctSignal[grpSO[tnsMaxBands]]; // end of spectrum to be predicted

      errorValue |= m_linPredictor.quantTnsToLpCoeffs (tnsData.coeff[0], filtOrder, tnsData.coeffResLow, tnsData.coeffParCor, filterC);

      // back up the leading MDCT samples
      memcpy (m_tempIntBuf, &mdctSignal[offs - MAX_PREDICTION_ORDER], MAX_PREDICTION_ORDER * sizeof (int32_t));
      // TNS compliance: set them to zero
      memset (&mdctSignal[offs - MAX_PREDICTION_ORDER], 0, MAX_PREDICTION_ORDER * sizeof (int32_t));

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
      memcpy (&mdctSignal[offs - MAX_PREDICTION_ORDER], m_tempIntBuf, MAX_PREDICTION_ORDER * sizeof (int32_t));

      // recalculate SFB RMS in TNS range
      errorValue |= m_specAnalyzer.getMeanAbsValues (mdctSignal, nullptr, nSamplesInFrame, 0, &grpSO[tnsStartSfb], __max (0, tnsMaxBands - (int) tnsStartSfb),
                                                     &grpData.sfbRmsValues[tnsStartSfb + m_numSwbShort * tnsData.filteredWindow]);
    }
    else tnsData.filterOrder[0] = tnsData.numFilters = 0; // disable zero-length TNS filters
  } // if order > 0

  return errorValue;
}

unsigned ExhaleEncoder::eightShortGrouping (SfbGroupData& grpData, uint16_t* const grpOffsets, int32_t* const mdctSignal)
{
  const unsigned nSamplesInFrame = toFrameLength (m_frameLength);
  const unsigned nSamplesInShort = nSamplesInFrame >> 3;
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
      }
    }
    grpOffset[m_numSwbShort] = uint16_t (grpStartLine + nSamplesInShort * grpLength);
  } // for gr

  memcpy (mdctSignal, m_tempIntBuf, nSamplesInFrame * sizeof (int32_t));

  return 0; // no error
}

unsigned ExhaleEncoder::getOptParCorCoeffs (const int32_t* const mdctSignal, const SfbGroupData& grpData, const uint8_t maxSfb,
                                            const unsigned channelIndex, TnsData& tnsData, const uint8_t firstGroupIndexToTest /*= 0*/)
{
  const unsigned nSamplesInFrame = toFrameLength (m_frameLength);
  const unsigned tnsStartSfb = 3 + 32000 / toSamplingRate (m_frequencyIdx); // 8-short start
  unsigned bestOrder = MAX_PREDICTION_ORDER, predGainCurr, predGainPrev, temp = 0;
  int16_t parCorBuffer[MAX_PREDICTION_ORDER];

  tnsData.filterOrder[0] = tnsData.filteredWindow = tnsData.numFilters = 0; // zero TNS data
  tnsData.filterDownward[0] = false;   // enforce direction = 0 for now, detection difficult

  if ((mdctSignal == nullptr) || (tnsData.coeffParCor == nullptr) || (maxSfb <= tnsStartSfb) || (channelIndex >= USAC_MAX_NUM_CHANNELS))
  {
    return 0; // invalid arguments error
  }

  if (grpData.numWindowGroups == 1) // LONG window: use ParCor coeffs from spectral analyzer
  {
    tnsData.filterOrder[0] = (uint8_t) m_specAnalyzer.getLinPredCoeffs (tnsData.coeffParCor, channelIndex);

#if EE_OPT_TNS_SPEC_RANGE
    if (tnsData.filterOrder[0] > 0) // try to reduce TNS start band as long as SNR increases
    {
      const uint16_t filtOrder = tnsData.filterOrder[0];
      uint16_t b = __min (m_specAnaCurr[channelIndex] & 31, (nSamplesInFrame - filtOrder) >> SA_BW_SHIFT);
      short filterC[MAX_PREDICTION_ORDER] = {0, 0, 0, 0};
      int32_t* predSig = &m_mdctSignals[channelIndex][b << SA_BW_SHIFT]; // TNS start offset

      m_linPredictor.parCorToLpCoeffs (tnsData.coeffParCor, filtOrder, filterC);

      for (b = (b > 0 ? b - 1 : 0), predSig--; b > 0; b--) // b is in spectr. analysis units
      {
        uint64_t sumAbsOrg = 0, sumAbsTns = 0;

        if (filtOrder >= 4) // max. order 4
        {
          for (uint16_t s = 1 << SA_BW_SHIFT; s > 0; s--) // produce the TNS filter residual
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
          for (uint16_t s = 1 << SA_BW_SHIFT; s > 0; s--) // produce the TNS filter residual
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
          for (uint16_t s = 1 << SA_BW_SHIFT; s > 0; s--) // produce the TNS filter residual
          {
            const int64_t predSample = *(predSig - 1) * (int64_t) filterC[0] + *(predSig - 2) * (int64_t) filterC[1];
            const int64_t mdctSample = *(predSig--);
            const int64_t resiSample = mdctSample + ((predSample + (1 << 8)) >> 9);

            sumAbsOrg += abs (mdctSample);  sumAbsTns += abs (resiSample);
          }
        }
        if (sumAbsOrg * 9 <= sumAbsTns * 8) break; // band SNR was reduced by more than 1 dB
      }
      m_specAnaCurr[channelIndex] = (m_specAnaCurr[channelIndex] & (UINT_MAX - 31)) | (b + 1);
    } // if order > 0
#endif // EE_OPT_TNS_SPEC_RANGE
    return (m_specAnaCurr[channelIndex] >> 24) & UCHAR_MAX; // spectral analyzer's pred gain
  }
  // SHORT window: find short group with maximum pred gain, then determine best filter order
  for (uint8_t gr = firstGroupIndexToTest; gr < grpData.numWindowGroups; gr++)
  {
    if (grpData.windowGroupLength[gr] == 1)
    {
      const uint16_t* grpSO = &grpData.sfbOffsets[m_numSwbShort * gr];

      predGainCurr = m_linPredictor.calcParCorCoeffs (&mdctSignal[grpSO[tnsStartSfb]], grpSO[maxSfb] - grpSO[tnsStartSfb], bestOrder, parCorBuffer);

      if (temp < predGainCurr)  // current pred gain set is "better" than best pred gain set
      {
        temp = predGainCurr;
        tnsData.filteredWindow = gr; // changed later
        memcpy (tnsData.coeffParCor, parCorBuffer, bestOrder * sizeof (int16_t));
      }
    }
  } // for gr

  predGainCurr = (temp >> 24) & UCHAR_MAX;
  predGainPrev = (temp >> 16) & UCHAR_MAX;
  while ((bestOrder > 1) && (predGainPrev >= predGainCurr)) // get lowest-order gain maximum
  {
    bestOrder--;
    predGainCurr = predGainPrev;
    predGainPrev = (temp >> (8 * bestOrder - 16)) & UCHAR_MAX;
  }
  tnsData.filterOrder[0] = ((bestOrder == 1) && (tnsData.coeffParCor[0] == 0) ? 0 : bestOrder);

  return predGainCurr;  // maximum pred gain of all filter orders and length-1 window groups
}

unsigned ExhaleEncoder::psychBitAllocation () // perceptual bit-allocation via scale factors
{
  const unsigned nChannels       = toNumChannels (m_channelConf);
  const unsigned nSamplesInFrame = toFrameLength (m_frameLength);
  const unsigned samplingRate    = toSamplingRate (m_frequencyIdx);
  const unsigned lfeChannelIndex = (m_channelConf >= CCI_6_CH ? __max (5, nChannels - 1) : USAC_MAX_NUM_CHANNELS);
  const uint32_t maxSfbLong      = (samplingRate < 37566 ? 51 /*32 kHz*/ : brModeAndFsToMaxSfbLong (m_bitRateMode, samplingRate));
#if EC_TRELLIS_OPT_CODING
  const uint64_t scaleSr         = (samplingRate < 27713 ? 37 - m_bitRateMode : 39 - (m_bitRateMode > 2 ? 1 : 0));
#else
  const uint64_t scaleSr         = (samplingRate < 27713 ? 36 - m_bitRateMode : 37);
#endif
  const uint64_t scaleBr         = (m_bitRateMode == 0 ? 32 : scaleSr - eightTimesSqrt256Minus[256 - m_bitRateMode] - ((m_bitRateMode - 1) >> 1));
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
        const unsigned lfAtten = 4 + b * 2; // LF SNR boost, cf my M.Sc. thesis, p. 54
        const uint8_t sfbWidth = off[b + 1] - off[b];
        const uint64_t   scale = scaleBr * __min (32, lfAtten); // rate control part 1

        // scale step-sizes according to VBR mode, then derive scale factors from step-sizes
        stepSizes[b] = uint32_t (__max (BA_EPS, ((1u << 9) + stepSizes[b] * scale) >> 10));

        scaleFactors[b] = m_bitAllocator.getScaleFac (stepSizes[b], &m_mdctSignals[ci][off[b]], sfbWidth, rms[b]);
      }
      ci++;
    }
    else // SCE or CPE: bandwidth-to-max_sfb mapping, short-window grouping for each channel
    {
      for (unsigned ch = 0; ch < nrChannels; ch++) // channel loop
      {
        SfbGroupData&  grpData = coreConfig.groupingData[ch];
        const bool eightShorts = (coreConfig.icsInfoCurr[ch].windowSequence == EIGHT_SHORT);
        const uint8_t  mSfmFac = eightTimesSqrt256Minus[meanSpecFlat[ci]];
        uint32_t*    stepSizes = &sfbStepSizes[ci * m_numSwbShort * NUM_WINDOW_GROUPS];
        uint8_t    numSwbFrame = (eightShorts ? numSwbOffsetS[m_swbTableIdx] : numSwbOffsetL[m_swbTableIdx]) - 1;

        if (!eightShorts && (samplingRate >= 37566) && (samplingRate < 55426)) // fix numSwb
        {
          numSwbFrame = 49;
        }
        memset (grpData.scaleFactors, 0, (MAX_NUM_SWB_SHORT * NUM_WINDOW_GROUPS) * sizeof (uint8_t));

        for (uint16_t gr = 0; gr < grpData.numWindowGroups; gr++)
        {
          const uint16_t* grpOff = &grpData.sfbOffsets[m_numSwbShort * gr];
          const uint32_t* grpRms = &grpData.sfbRmsValues[m_numSwbShort * gr];
          const uint32_t* refRms = &coreConfig.groupingData[1 - ch].sfbRmsValues[m_numSwbShort * gr];
          uint8_t*  grpScaleFacs = &grpData.scaleFactors[m_numSwbShort * gr];
          uint32_t* grpStepSizes = &stepSizes[m_numSwbShort * gr];
          uint32_t  b, grpRmsMin = INT_MAX; // min. RMS value, used for overcoding reduction

          // undercoding reduction for case where large number of coefs is quantized to zero
          s = (eightShorts ? (nSamplesInFrame * grpData.windowGroupLength[gr]) >> 1 : nSamplesInFrame << 2);
          for (b = 0; b < grpData.sfbsPerGroup; b++)
          {
            const uint32_t rmsComp = getComplexRmsValue (grpRms[b], gr, b, numSwbFrame, coreConfig.tnsData[ch]);
            const uint32_t rmsRef9 = (!coreConfig.commonWindow ? rmsComp :
                                     getComplexRmsValue (refRms[b], gr, b, numSwbFrame, coreConfig.tnsData[1 - ch]) >> 9);

            if (rmsComp < grpRmsMin) grpRmsMin = rmsComp;
            if (rmsComp >= rmsRef9 && (rmsComp < (grpStepSizes[b] >> 1)))  // zero-quantized
            {
              s -= ((grpOff[b + 1] - grpOff[b]) * 3 * __min (2 * SA_EPS, rmsComp) + SA_EPS) >> 11; // / (2 * SA_EPS)
            }
          }
          if ((samplingRate >= 27713) && (b < maxSfbLong) && !eightShorts)  // uncoded coefs
          {
            const uint32_t rmsComp = getComplexRmsValue (grpRms[b], gr, b, numSwbFrame, coreConfig.tnsData[ch]);
            const uint32_t rmsRef9 = (!coreConfig.commonWindow ? rmsComp :
                                     getComplexRmsValue (refRms[b], gr, b, numSwbFrame, coreConfig.tnsData[1 - ch]) >> 9);

            if (rmsComp >= rmsRef9) // check only first SFB above max_sfb for simplification
            {
              s -= ((grpOff[maxSfbLong] - grpOff[b]) * 3 * __min (2 * SA_EPS, rmsComp) + SA_EPS) >> 11; // / (2 * SA_EPS)
            }
          }
          s = (eightShorts ? s / ((nSamplesInFrame * grpData.windowGroupLength[gr]) >> 8) : s / (nSamplesInFrame >> 5));

          for (b = 0; b < grpData.sfbsPerGroup; b++)
          {
            const unsigned lfAtten = (b <= 5 ? (eightShorts ? 1 : 4) + b * 2 : 9 + b + ((b + 5) >> 4));  // LF SNR boost
            const uint8_t sfbWidth = grpOff[b + 1] - grpOff[b];
            const uint64_t rateFac = mSfmFac * s * __min (32, lfAtten * grpData.numWindowGroups); // rate control part 1
            const uint64_t sScaled = ((1u << 23) + __max (grpRmsMin, grpStepSizes[b]) * scaleBr * rateFac) >> 24;

            // scale step-sizes according to VBR mode & derive scale factors from step-sizes
            grpStepSizes[b] = uint32_t (__max (BA_EPS, __min (UINT_MAX, sScaled)));

            grpScaleFacs[b] = m_bitAllocator.getScaleFac (grpStepSizes[b], &m_mdctSignals[ci][grpOff[b]], sfbWidth, grpRms[b]);
          }
        } // for gr

#if !RESTRICT_TO_AAC
        if (grpData.sfbsPerGroup > 0 && m_noiseFilling[el] && !eightShorts) // HF noise-fill
        {
          numSwbFrame = __min (numSwbFrame, maxSfbLong); // bit-rate dependent max bandwidth

          if (grpData.sfbsPerGroup < numSwbFrame)
          {
            memset (&grpData.scaleFactors[grpData.sfbsPerGroup], 0, (numSwbFrame - grpData.sfbsPerGroup) * sizeof (uint8_t));
            grpData.sfbsPerGroup = coreConfig.icsInfoCurr[ch].maxSfb = numSwbFrame;
          }
          if (ch > 0) coreConfig.commonMaxSfb = (coreConfig.icsInfoCurr[0].maxSfb == coreConfig.icsInfoCurr[1].maxSfb);
        }
#endif
        ci++;
      } // for ch

      for (unsigned ch = 0; ch < nrChannels; ch++) // channel loop
      {
        SfbGroupData& grpData = coreConfig.groupingData[ch];
        TnsData&      tnsData = coreConfig.tnsData[ch];

        if (tnsData.numFilters > 0) // convert TNS group index to window index for write-out
        {
          s = 0;
          for (uint16_t gr = 0; gr < grpData.numWindowGroups; gr++)
          {
            if (gr == tnsData.filteredWindow)
            {
              tnsData.filteredWindow = (uint8_t) s;
              break;
            }
            s += grpData.windowGroupLength[gr];
          }
        }
      } // for ch
    }
  } // for el

  return errorValue;
}

unsigned ExhaleEncoder::quantizationCoding ()  // apply MDCT quantization and entropy coding
{
  const unsigned nSamplesInFrame  = toFrameLength (m_frameLength);
  const unsigned samplingRate     = toSamplingRate (m_frequencyIdx);
  const unsigned* const coeffMagn = m_sfbQuantizer.getCoeffMagnPtr ();
  unsigned ci = 0, s; // running index
  unsigned errorValue = (coeffMagn == nullptr ? 1 : 0);

  for (unsigned el = 0; el < m_numElements; el++)  // element loop
  {
    CoreCoderData& coreConfig = *m_elementData[el];
    const unsigned nrChannels = (coreConfig.elementType & 1) + 1; // for UsacCoreCoderData()

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

            grpScaleFacs[b - 1] = grpScaleFacs[b] - INDEX_OFFSET; // reset prev. SFB to zero
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
          const uint8_t maxSfbLong  = (samplingRate < 37566 ? 51 /*32 kHz*/ : brModeAndFsToMaxSfbLong (m_bitRateMode, samplingRate));
          const uint8_t maxSfbShort = (samplingRate < 37566 ? 14 /*32 kHz*/ : brModeAndFsToMaxSfbShort(m_bitRateMode, samplingRate));
          const uint16_t peakIndex  = (shortWinCurr ? 0 : (m_specAnaCurr[ci] >> 5) & 2047);
#if RESTRICT_TO_AAC
          const unsigned sfmBasedSfbStart = (shortWinCurr ? maxSfbShort : maxSfbLong) - 6 + (m_bitRateMode >> 1) + ((m_specAnaCurr[ci] >> 21) & 7);
#else
          const unsigned highFreqMinStart = (m_noiseFilling[el] ? 6 : 6 - (m_bitRateMode >> 1));
          const unsigned sfmBasedSfbStart = (shortWinCurr ? maxSfbShort : maxSfbLong) - highFreqMinStart + ((m_specAnaCurr[ci] >> 21) & 7);
#endif
          const unsigned targetBitCountX2 = ((48000 + 16000 * m_bitRateMode) * nSamplesInFrame) / (samplingRate * grpData.numWindowGroups);
          unsigned b = grpData.sfbsPerGroup - 1;

          if ((grpRms[b] >> 16) > 0) lastSfb = b;
          estimBitCount += grpRms[b] & USHRT_MAX;

          if (grpLength == 1) // finalize bit count estimate, RDOC
          {
            estimBitCount += ((entrCoder.arithGetCtxState () >> 17) & 31) + 2; // m_acBits+2
#if EC_TRELLIS_OPT_CODING
            estimBitCount = m_sfbQuantizer.quantizeSpecRDOC (entrCoder, grpScaleFacs, m_bitRateMode, // __min (estimBitCount, targetBitCountX2),
                                                             grpOff, grpRms, grpData.sfbsPerGroup, m_mdctQuantMag[ci]);
#endif
          }
          b = lastSfb;
          while ((b >= sfmBasedSfbStart) && (grpOff[b] > peakIndex) && ((grpRms[b] >> 16) <= 1) /*coarse quantization*/ &&
                 ((estimBitCount * 2 > targetBitCountX2) || (grpLength > 1 /*no accurate bit count estimate available*/)))
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

          if (estimBitCount > targetBitCountX2) // too many bits!!
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
              if (estimBitCount <= targetBitCountX2) break;
            }

            for (b++; b <= lastSfb; b++) // re-repeat scale factor
            {
              if ((grpRms[b] >> 16) == 0) // a zero quantized band
              {
                grpScaleFacs[b] = grpScaleFacs[b - 1];
              }
            }
          } // if (estimBitCount > targetBitCountX2)

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
      // obtain channel-wise noise_level and noise_offset for USAC
      coreConfig.specFillData[ch] = (!m_noiseFilling[el] ? 0 : m_specGapFiller.getSpecGapFillParams (m_sfbQuantizer, m_mdctQuantMag[ci],
                                                                                                     m_numSwbShort, grpData, nSamplesInFrame));
      // NOTE: gap-filling SFB bit count might be inaccurate now since scale factors changed
      if (coreConfig.specFillData[ch] == 1) errorValue |= 1;
#endif
      ci++;
    }
  } // for el

  return (errorValue > 0 ? 0 : m_outStream.createAudioFrame (m_elementData, m_entropyCoder, m_mdctSignals, m_mdctQuantMag, m_indepFlag,
                                                             m_numElements, m_numSwbShort, (uint8_t* const) m_tempIntBuf,
#if !RESTRICT_TO_AAC
                                                             m_timeWarping, m_noiseFilling,
#endif
                                                             m_outAuData, nSamplesInFrame)); // returns AU size
}

unsigned ExhaleEncoder::spectralProcessing ()  // complete ics_info(), calc TNS and SFB data
{
  const unsigned nChannels       = toNumChannels (m_channelConf);
  const unsigned nSamplesInFrame = toFrameLength (m_frameLength);
  const unsigned nSamplesInShort = nSamplesInFrame >> 3;
  const unsigned samplingRate    = toSamplingRate (m_frequencyIdx);
  const unsigned lfeChannelIndex = (m_channelConf >= CCI_6_CH ? __max (5, nChannels - 1) : USAC_MAX_NUM_CHANNELS);
  unsigned ci = 0, s; // running index
  unsigned errorValue = 0; // no error

  // get spectral channel statistics for last frame, used for input bandwidth (BW) detection
  m_specAnalyzer.getSpecAnalysisStats (m_specAnaPrev, nChannels);
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
      ci++;
    }
    else // SCE or CPE: bandwidth-to-max_sfb mapping, short-window grouping for each channel
    {
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
          grpData.sfbsPerGroup = icsCurr.maxSfb; // initialization, changed to max_sfb later

          if (samplingRate > 32000) // set max_sfb based on VBR mode and bandwidth detection
          {
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
          while (grpSO[icsCurr.maxSfb] > __max (m_bandwidCurr[ci], m_bandwidPrev[ci])) icsCurr.maxSfb--; // not a bug!!

          errorValue |= eightShortGrouping (grpData, grpSO, m_mdctSignals[ci]);
        } // if EIGHT_SHORT

        // compute and quantize optimal TNS coefficients, then find optimal TNS filter order
        s /*linear pred gain*/ = getOptParCorCoeffs (m_mdctSignals[ci], grpData, icsCurr.maxSfb, ci, tnsData,
                                                     ch > 0 && coreConfig.commonWindow ? coreConfig.tnsData[0].filteredWindow : 0);
        tnsData.filterOrder[0] = m_linPredictor.calcOptTnsCoeffs (tnsData.coeffParCor, tnsData.coeff[0], &tnsData.coeffResLow,
                                                                  tnsData.filterOrder[0], s, (m_specAnaCurr[ci] >> 16) & UCHAR_MAX);
        tnsData.numFilters = (tnsData.filterOrder[0] > 0 ? 1 : 0);
        ci++;
      } // for ch

      if (coreConfig.commonWindow) // synchronization of all StereoCoreToolInfo() components
      {
        uint8_t& maxSfb0 = coreConfig.icsInfoCurr[0].maxSfb;
        uint8_t& maxSfb1 = coreConfig.icsInfoCurr[1].maxSfb;
        const uint8_t maxSfbSte = __max (maxSfb0, maxSfb1);   // max_sfb_ste, as in Table 24

        if ((maxSfb0 > 0) && (maxSfb1 > 0) && (maxSfbSte - __min (maxSfb0, maxSfb1) <= 1))
        {
          uint32_t& sa0 = m_specAnaCurr[ci-2];
          uint32_t& sa1 = m_specAnaCurr[ci-1];
          const int specFlat[2] = {int (sa0 >> 16) & UCHAR_MAX, int (sa1 >> 16) & UCHAR_MAX};
          const int tnsStart[2] = {int (sa0 & 31), int (sa1 & 31)}; // long TNS start offset

          if ((coreConfig.tnsData[0].filteredWindow == coreConfig.tnsData[1].filteredWindow) &&
              (abs (specFlat[0] - specFlat[1]) <= (UCHAR_MAX >> 3)) &&
              (abs (tnsStart[0] - tnsStart[1]) <= (UCHAR_MAX >> 5)))  // TNS synchronization
          {
            const uint16_t maxTnsOrder = __max (coreConfig.tnsData[0].filterOrder[0], coreConfig.tnsData[1].filterOrder[0]);
            TnsData& tnsData0 = coreConfig.tnsData[0];
            TnsData& tnsData1 = coreConfig.tnsData[1];

            if (m_linPredictor.similarParCorCoeffs (tnsData0.coeffParCor, tnsData1.coeffParCor, maxTnsOrder, LP_DEPTH))
            {
              coreConfig.commonTnsData = true; // synch tns_data
              for (s = 0; s < maxTnsOrder; s++)
              {
                tnsData0.coeffParCor[s] = (tnsData0.coeffParCor[s] + tnsData1.coeffParCor[s] + 1) >> 1;
              }
              tnsData0.coeffResLow = false; // reoptimize coeffs
              tnsData0.filterOrder[0] = m_linPredictor.calcOptTnsCoeffs (tnsData0.coeffParCor, tnsData0.coeff[0], &tnsData0.coeffResLow,
                                                                         maxTnsOrder, UCHAR_MAX /*maximum pred gain*/, 0, LP_DEPTH);
              tnsData0.numFilters = (tnsData0.filterOrder[0] > 0 ? 1 : 0);
              memcpy (&tnsData1, &tnsData0, sizeof (TnsData));
            }
            else if ((maxTnsOrder > 0) && (tnsData0.coeffResLow == tnsData1.coeffResLow) && (tnsData0.filterOrder[0] == tnsData1.filterOrder[0]))
            {
              const int32_t* coeff0 = (int32_t*) tnsData0.coeff[0]; // fast comparison code,
              const int32_t* coeff1 = (int32_t*) tnsData1.coeff[0]; // might not be portable

              coreConfig.commonTnsData = (*coeff0 == *coeff1); // first four coeffs the same
            }
            if (coreConfig.commonTnsData) // synch TNS start SFB
            {
              const uint32_t avgTnsStart = (tnsStart[0] + tnsStart[1]) >> 1;  // mean offset

              sa0 = (sa0 & (UINT_MAX - 31)) | avgTnsStart;  // is used by applyTnsToWinGroup
              sa1 = (sa1 & (UINT_MAX - 31)) | avgTnsStart;
            }
          }
          maxSfb0 = maxSfb1 = maxSfbSte;
        }
        coreConfig.commonMaxSfb = (maxSfb0 == maxSfb1); // synch
        coreConfig.stereoConfig = coreConfig.stereoMode = 0;
      } // if coreConfig.commonWindow
    }

    ci -= nrChannels; // zero frequency coefficients above num_swb for all channels, windows

    for (unsigned ch = 0; ch < nrChannels; ch++) // channel loop
    {
      SfbGroupData&  grpData = coreConfig.groupingData[ch];
      const uint16_t*  grpSO = grpData.sfbOffsets;
      const IcsInfo& icsCurr = coreConfig.icsInfoCurr[ch];
      const bool eightShorts = (icsCurr.windowSequence == EIGHT_SHORT);
      unsigned grpEndLine = 0;

      if (eightShorts) // map grouping table idx to scale_factor_grouping idx for bit-stream
      {
        coreConfig.icsInfoCurr[ch].windowGrouping = scaleFactorGrouping[icsCurr.windowGrouping];
      }

      for (uint16_t gr = 0; gr < grpData.numWindowGroups; gr++)
      {
        const unsigned grpSOStart = grpSO[grpData.sfbsPerGroup + m_numSwbShort * gr];

        grpEndLine += (eightShorts ? nSamplesInShort : nSamplesInFrame) * grpData.windowGroupLength[gr];
        memset (&m_mdctSignals[ci][grpSOStart], 0, (grpEndLine - grpSOStart) * sizeof (int32_t));
        memset (&m_mdstSignals[ci][grpSOStart], 0, (grpEndLine - grpSOStart) * sizeof (int32_t));
      }
      memset (grpData.sfbRmsValues, 0, (MAX_NUM_SWB_SHORT * NUM_WINDOW_GROUPS) * sizeof (uint32_t));

      if (icsCurr.maxSfb > 0)
      {
        // use MCLTs for LONG but only MDCTs for SHORT windows (since MDSTs are not grouped)
        errorValue |= m_specAnalyzer.getMeanAbsValues (m_mdctSignals[ci], eightShorts ? nullptr : m_mdstSignals[ci], nSamplesInFrame,
                                                       ci, grpSO, grpData.sfbsPerGroup * grpData.numWindowGroups, grpData.sfbRmsValues);
        errorValue |= applyTnsToWinGroup (coreConfig.tnsData[ch], grpData, eightShorts, icsCurr.maxSfb, ci);
        coreConfig.tnsActive |= (coreConfig.tnsData[ch].numFilters > 0); // tns_data_present
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
  const unsigned nSamplesInFrame = toFrameLength (m_frameLength);
  const unsigned nSamplesTempAna = (nSamplesInFrame * 25) >> 4;  // pre-delay for look-ahead
  const unsigned lfeChannelIndex = (m_channelConf >= CCI_6_CH ? __max (5, nChannels - 1) : USAC_MAX_NUM_CHANNELS);
  unsigned ci = 0; // running ch index
  unsigned errorValue = 0; // no error

  // get temporal channel statistics for this frame, used for spectral grouping/quantization
  m_tempAnalyzer.getTempAnalysisStats (m_tempAnaCurr, nChannels);
  m_tempAnalyzer.getTransientLocation (m_tranLocCurr, nChannels);

  // temporal analysis for look-ahead signal (central nSamplesInFrame samples of next frame)
  errorValue |= m_tempAnalyzer.temporalAnalysis (m_timeSignals, nChannels, nSamplesInFrame, nSamplesTempAna, lfeChannelIndex);

  // get temporal channel statistics for next frame, used for window length/overlap decision
  m_tempAnalyzer.getTempAnalysisStats (m_tempAnaNext, nChannels);
  m_tempAnalyzer.getTransientLocation (m_tranLocNext, nChannels);

  m_indepFlag = (((m_frameCount++) % m_indepPeriod) == 0); // configure usacIndependencyFlag

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
        const unsigned  sfCurr = (m_tempAnaCurr[ci] >> 24) & UCHAR_MAX;
        const unsigned  tfCurr = (m_tempAnaCurr[ci] >> 16) & UCHAR_MAX;
        const unsigned  sfNext = (m_tempAnaNext[ci] >> 24) & UCHAR_MAX;
        const unsigned  tfNext = (m_tempAnaNext[ci] >> 16) & UCHAR_MAX;

        tsCurr[ch] = (m_tempAnaCurr[ci] /*R*/) & UCHAR_MAX;
        tsNext[ch] = (m_tempAnaNext[ci] >>  8) & UCHAR_MAX;

        const bool lowOlapNext = (m_tranLocNext[ci] >= 0) || (sfNext < 68 && tfNext >= 204) || (tsCurr[ch] >= 153) || (tsNext[ch] >= 153);
        const bool sineWinCurr = (sfCurr >= 170) && (sfNext >= 170) && (sfCurr < 221) && (sfNext < 221) && (tsCurr[ch] < 20) &&
                                 (tfCurr >= 153) && (tfNext >= 153) && (tfCurr < 184) && (tfNext < 184) && (tsNext[ch] < 20);
        // set window_sequence
        if ((wsPrev == ONLY_LONG) || (wsPrev == LONG_STOP)) // 1st window half - max overlap
        {
          wsCurr = (lowOlapNext ? LONG_START : ONLY_LONG);
        }
        else // LONG_START_SEQUENCE, STOP_START_SEQUENCE, EIGHT_SHORT_SEQUENCE - min overlap
        {
          wsCurr = (m_tranLocCurr[ci] >= 0) ? EIGHT_SHORT :
#if RESTRICT_TO_AAC
                   (lowOlapNext && (m_tranLocNext[ci] >= 0 || wsPrev != EIGHT_SHORT) ? EIGHT_SHORT : LONG_STOP);
#else
                   (lowOlapNext && (m_tranLocNext[ci] >= 0 || wsPrev != STOP_START) ? STOP_START : LONG_STOP);
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
        icsCurr.windowGrouping = (wsCurr == EIGHT_SHORT ? (m_tranLocCurr[ci] * 8) / (int16_t) nSamplesInFrame : 0);
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
          const USAC_WSEQ initialWs0 = (USAC_WSEQ) winSeq0;
          const USAC_WSEQ initialWs1 = (USAC_WSEQ) winSeq1;

          winSeq0 = winSeq1 = windowSequenceSynch[initialWs0][initialWs1];   // equalization
          if ((winSeq0 != initialWs0) && (winSeq0 == EIGHT_SHORT))
          {
#if !RESTRICT_TO_AAC
            if ((tsCurr[0] * 7 < tsCurr[1] * 2) && (tsNext[0] * 7 < tsNext[1] * 2))
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
            if ((tsCurr[1] * 7 < tsCurr[0] * 2) && (tsNext[1] * 7 < tsNext[0] * 2))
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

          icsInfo0.windowGrouping = icsInfo1.windowGrouping = (tranLocSynch * 8) / (int16_t) nSamplesInFrame;
        }

        if ((icsInfo0.windowShape != WINDOW_SINE) || (icsInfo1.windowShape != WINDOW_SINE))
        {
          icsInfo0.windowShape = WINDOW_KBD; // always synchronize window_shapes in order to
          icsInfo1.windowShape = WINDOW_KBD; // encourage synch in next frame; KBD dominates
        }
        coreConfig.commonWindow = (winSeq0 == winSeq1); // synch
      } // if nrChannels > 1
    }

    ci -= nrChannels; // modulated complex lapped transform (MCLT) for all channels, windows

    for (unsigned ch = 0; ch < nrChannels; ch++) // channel loop
    {
      const IcsInfo& icsPrev = coreConfig.icsInfoPrev[ch];
      const IcsInfo& icsCurr = coreConfig.icsInfoCurr[ch];
      const USAC_WSEQ wsCurr = icsCurr.windowSequence;
      const bool eightShorts = (wsCurr == EIGHT_SHORT);
      SfbGroupData&  grpData = coreConfig.groupingData[ch];

      grpData.numWindowGroups = (eightShorts ? NUM_WINDOW_GROUPS : 1);  // fill groupingData
      memcpy (grpData.windowGroupLength, windowGroupingTable[icsCurr.windowGrouping], NUM_WINDOW_GROUPS * sizeof (uint8_t));

      errorValue |= m_transform.applyMCLT (m_timeSignals[ci], eightShorts, icsPrev.windowShape != WINDOW_SINE, icsCurr.windowShape != WINDOW_SINE,
                                           wsCurr > LONG_START /*lOL*/, (wsCurr % 3) != ONLY_LONG /*lOR*/, m_mdctSignals[ci], m_mdstSignals[ci]);
      m_scaleFacData[ci] = &grpData;
      ci++;
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
  m_frameCount   = 0;
  m_frameLength  = (USAC_CCFL) frameLength; // coreCoderFrameLength, signaled using coreSbrFrameLengthIndex
  m_frequencyIdx = toSamplingFrequencyIndex (sampleRate);  // I/O sample rate as usacSamplingFrequencyIndex
  m_indepFlag    = true; // usacIndependencyFlag in UsacFrame(), will be set per frame, true in first frame
  m_indepPeriod  = (indepPeriod == 0 ? UINT_MAX : indepPeriod); // RAP, signaled using usacIndependencyFlag
#if !RESTRICT_TO_AAC
  m_nonMpegExt   = useEcodisExt;
#endif
  m_numSwbShort  = MAX_NUM_SWB_SHORT;
  m_outAuData    = outputAuData;
  m_pcm24Data    = inputPcmData;
  m_tempIntBuf   = nullptr;

  // initialize all helper structs
  for (unsigned el = 0; el < USAC_MAX_NUM_ELEMENTS; el++)
  {
    const ELEM_TYPE et = elementTypeConfig[m_channelConf % USAC_MAX_NUM_ELCONFIGS][el];  // usacElementType

    m_elementData[el]  = nullptr;
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
    m_mdctQuantMag[ch] = nullptr;
    m_mdctSignals[ch]  = nullptr;
    m_mdstSignals[ch]  = nullptr;
    m_scaleFacData[ch] = nullptr;
    m_specAnaCurr[ch]  = 0;
    m_specAnaPrev[ch]  = 0;
    m_tempAnaCurr[ch]  = 0;
    m_tempAnaNext[ch]  = 0;
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
  const unsigned nSamplesInFrame = toFrameLength (m_frameLength);
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
  }

  // set initial temporal channel statistic to something meaningful before first coded frame
  m_tempAnalyzer.temporalAnalysis (m_timeSignals, nChannels, nSamplesInFrame, nSamplesTempAna - nSamplesInFrame);

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
  const unsigned nSamplesInFrame = toFrameLength (m_frameLength);
  const unsigned nSamplesTempAna = (nSamplesInFrame * 25) >> 4;  // pre-delay for look-ahead
  const int32_t* chSig           = m_pcm24Data;
  unsigned ch, s;

  // move internal channel buffers nSamplesInFrame to the past to make room for next samples
  for (ch = 0; ch < nChannels; ch++)
  {
    memcpy (&m_timeSignals[ch][0], &m_timeSignals[ch][nSamplesInFrame], nSamplesInFrame * sizeof (int32_t));
    memcpy (&m_timeSignals[ch][nSamplesInFrame], &m_timeSignals[ch][2 * nSamplesInFrame], (nSamplesTempAna - nSamplesInFrame) * sizeof (int32_t));
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
  const unsigned timeSigBufSize  = ((nSamplesInFrame * 41) >> 4) * sizeof (int32_t); // core-codec delay*4
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
  if (m_frequencyIdx < 0)
  {
    errorValue |=  32;
  }
  if ((m_outAuData == nullptr) || (m_pcm24Data == nullptr))
  {
    errorValue |=  16;
  }
  if (errorValue > 0) return errorValue;

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
  // allocate all signal buffers
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
  errorValue = (unsigned) m_frequencyIdx; // for temporary storage
#if RESTRICT_TO_AAC
  m_swbTableIdx = freqIdxToSwbTableIdxAAC[errorValue];
#else
  m_swbTableIdx = (m_frameLength == CCFL_768 ? freqIdxToSwbTableIdx768[errorValue] : freqIdxToSwbTableIdxAAC[errorValue]);
#endif
  m_tempIntBuf  = m_timeSignals[0];
  errorValue = 0;
#if EC_TRELLIS_OPT_CODING
  if (m_sfbQuantizer.initQuantMemory (nSamplesInFrame, numSwbOffsetL[m_swbTableIdx] - 1, m_bitRateMode) > 0 ||
#else
  if (m_sfbQuantizer.initQuantMemory (nSamplesInFrame) > 0 ||
#endif
      m_specAnalyzer.initLinPredictor (&m_linPredictor) > 0 ||
      m_transform.initConstants (m_tempIntBuf, m_timeWindowL, m_timeWindowS, nSamplesInFrame) > 0)
  {
    errorValue |= 1;
  }

  if ((errorValue == 0) && (audioConfigBuffer != nullptr)) // save UsacConfig() for writeout
  {
    errorValue = m_outStream.createAudioConfig (m_frequencyIdx, m_frameLength != CCFL_1024, chConf, m_numElements,
                                                elementTypeConfig[chConf], false /*usacConfigExtensionPresent=0*/,
#if !RESTRICT_TO_AAC
                                                m_timeWarping, m_noiseFilling,
#endif
                                                audioConfigBuffer);
    if (audioConfigBytes) *audioConfigBytes = errorValue; // length of UsacConfig() in bytes
    errorValue = (errorValue == 0 ? 1 : 0);
  }

  return errorValue;
}
