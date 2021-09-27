/* bitStreamWriter.cpp - source file for class with basic bit-stream writing capability
 * written by C. R. Helmrich, last modified in 2021 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under the exhale Copyright License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2021 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#include "exhaleLibPch.h"
#include "bitStreamWriter.h"
#include "bitAllocation.h" // define BA_MORE_CBR (more constant bit-rate, experimental!)

#ifndef NO_PREROLL_DATA
static const uint8_t zeroAu[2][14] = { // single-element AUs incl. SBR for digital silence
  {132,   0,  2, 0, 8, 0, 0,  0, 0, 0, 0, 0, 0, 0}, // SCE, 8 bytes
  {132, 129, 16, 0, 8, 0, 0, 32, 0, 0, 0, 0, 0, 0} // CPE, 14 bytes
};
#endif

// static helper functions
static inline int getPredCoefPrevGrp (const uint8_t aqIdxPrevGrp)
{
  return int (aqIdxPrevGrp > 0 ? aqIdxPrevGrp & 31 : 16);
}

static uint32_t getDeltaCodeTimeFlag (const uint8_t* const alphaQCurr, const unsigned numWinGroups, const unsigned numSwbShort,
                                      const uint8_t* const alphaQPrev, const unsigned maxSfbSte, const EntropyCoder& entrCoder,
                                      const bool complexCoef)
{
  unsigned b, g, bitCountFreq = 0, bitCountTime = 0;

  if ((alphaQCurr == nullptr) || (alphaQPrev == nullptr)) return 0;

  for (g = 0; g < numWinGroups; g++)
  {
    const uint8_t* const aqReIdxPrvGrp = (g == 0 ? alphaQPrev : &alphaQCurr[numSwbShort * (g - 1)]);
    const uint8_t* const aqImIdxPrvGrp = &aqReIdxPrvGrp[1];
    const uint8_t* const gCplxPredUsed = &alphaQCurr[numSwbShort * g];
    int  aqReIdxPred = 16, aqImIdxPred = 16; // init alpha_q_.. = 0

    for (b = 0; b < maxSfbSte; b += SFB_PER_PRED_BAND)
    {
      if (gCplxPredUsed[b] > 0) // count dpcm_alpha_q_re/_q_im bits
      {
        int aqIdx = gCplxPredUsed[b] & 31; // range -15,...0,...,15

        bitCountFreq += entrCoder.indexGetBitCount (aqIdx - aqReIdxPred);
        bitCountTime += entrCoder.indexGetBitCount (aqIdx - getPredCoefPrevGrp (aqReIdxPrvGrp[b]));

        aqReIdxPred = aqIdx;

        if (complexCoef)
        {
          aqIdx = gCplxPredUsed[b + 1] & 31; // TODO: <32 kHz short

          bitCountFreq += entrCoder.indexGetBitCount (aqIdx - aqImIdxPred);
          bitCountTime += entrCoder.indexGetBitCount (aqIdx - getPredCoefPrevGrp (aqImIdxPrvGrp[b]));

          aqImIdxPred = aqIdx;
        }
      }
      else aqReIdxPred = aqImIdxPred = 16;
    }
  } // for g

  return (bitCountFreq > bitCountTime ? 1 : 0);
}

static uint8_t getNumEnvBits (const int32_t int32Value, const uint8_t maxBitCount, const uint8_t minBitCount)
{
  uint8_t bits = __max (1, maxBitCount);
  int32_t mask = 1 << (bits - 1);

  while (bits > minBitCount && (int32Value & mask) == 0) // get MSB
  {
    bits--;
    mask >>= 1;
  }
  return bits;
}

#ifndef NO_PREROLL_DATA
static unsigned getLowRatePreRollAU (uint8_t* const byteBuffer, const CoreCoderData& elData, EntropyCoder& entrCoder,
                                     const uint8_t* ipfAuState, const uint8_t sbrRatioShiftValue)
{
  const uint16_t et = elData.elementType & 1;
  const bool notLFE = (elData.elementType < ID_USAC_LFE);
  const bool useSbr = (sbrRatioShiftValue > 0 && notLFE);
  unsigned byteCount;

  if (ipfAuState[0] == 0) // create zero-spectrum AU
  {
    byteCount = (8 + et * 6) >> (useSbr ? 0 : 1);
    memcpy (byteBuffer, zeroAu[et], byteCount);

    if (notLFE) // write appropriate window_sequence
    {
      const USAC_WSEQ wsPrev0 = elData.icsInfoPrev[0].windowSequence;
      const uint8_t wsPreRoll = uint8_t (wsPrev0 == EIGHT_SHORT || wsPrev0 == STOP_START ? LONG_START : wsPrev0);
      // SCE/CPE: window_sequence @ offset 2/0, left-shifted by 2/0
      byteBuffer[2 - 2 * et] |= wsPreRoll << (2 - 2 * et);
    }
  }
  else // complete AU with only 1 non-zero MDCT line
  {
    OutputStream au;
    unsigned ci = 0;

    byteCount = ((unsigned) ipfAuState[0] << 1) | (ipfAuState[1] >> 7);
    while (ci < byteCount) au.write (byteBuffer[ci++], 8);
    au.heldBitCount = ipfAuState[1] & SCHAR_MAX;
    au.heldBitChunk = ipfAuState[2];

    if (ipfAuState[3] > 0)
    {
      const uint16_t li = ipfAuState[4]; // line idx
      const uint16_t lg = __min (li + 56, (uint16_t) ipfAuState[3] << 2);

      memset (byteBuffer, 0, lg); // MDCT line coder
      byteBuffer[li] = ipfAuState[5] >> 1;
      entrCoder.initWindowCoding (true);
      entrCoder.arithCodeSigMagn (byteBuffer, 0, lg, true, &au);
      if (byteBuffer[li]) au.write (ipfAuState[5] & 1, 1); // sign
    }
    au.write (0, 1);// fac_data_present, no fac_data

    if (useSbr) // UsacSbrData()
    {
      au.write (1, 7);// SbrInfo(), sbrUseDfltHeader
      if (et) au.write (0, 1);  // fix _coupling = 0
      au.write (0, 8 << et);
      au.write (0, 16 << et);
# if ENABLE_INTERTES
      au.write (0, et + 1);
# endif
    }

    if (au.heldBitCount > 0) au.stream.push_back (au.heldBitChunk);
    byteCount = (unsigned) au.stream.size ();
    memcpy (byteBuffer, &au.stream.front (), byteCount);
  }

  return byteCount;
}
#endif // !NO_PREROLL_DATA

static uint8_t getOptMsMaskModeValue (const uint8_t* const msUsed, const unsigned numWinGroups, const unsigned numSwbShort,
                                      const uint8_t    msMaskMode, const unsigned maxSfbSte)
{
  const unsigned sfbStep = (msMaskMode < 3 ? 1 : SFB_PER_PRED_BAND);
  unsigned b, g;

  if ((msUsed == nullptr) || ((msMaskMode & 1) == 0)) return msMaskMode;

  for (g = 0; g < numWinGroups; g++)
  {
    const uint8_t* const gMsUsed = &msUsed[numSwbShort * g];

    for (b = 0; b < maxSfbSte; b += sfbStep)
    {
      if (gMsUsed[b] == 0) return msMaskMode;  // M/S in some bands
    }
  } // for g

  return (msMaskMode + 1); // upgrade mask mode to M/S in all bands
}

// private helper functions
void BitStreamWriter::writeByteAlignment () // write '0' bits until stream is byte-aligned
{
  if (m_auBitStream.heldBitCount > 0)
  {
    m_auBitStream.stream.push_back (m_auBitStream.heldBitChunk);
    m_auBitStream.heldBitChunk = 0;
    m_auBitStream.heldBitCount = 0;
  }
}

unsigned BitStreamWriter::writeChannelWiseIcsInfo (const IcsInfo& icsInfo)   // ics_info()
{
#if RESTRICT_TO_AAC
  m_auBitStream.write ((unsigned) icsInfo.windowSequence, 2);
#else
  m_auBitStream.write (unsigned (icsInfo.windowSequence == STOP_START ? LONG_START : icsInfo.windowSequence), 2);
#endif
  m_auBitStream.write ((unsigned) icsInfo.windowShape, 1);
  if (icsInfo.windowSequence == EIGHT_SHORT)
  {
    m_auBitStream.write (icsInfo.maxSfb, 4);
    m_auBitStream.write (icsInfo.windowGrouping, 7); // scale_factor_grouping

    return 14;
  }
  m_auBitStream.write (icsInfo.maxSfb, 6);

  return 9;
}

unsigned BitStreamWriter::writeChannelWiseSbrData (const int32_t* const sbrDataCh0, const int32_t* const sbrDataCh1,
                                                   const bool indepFlag /*= false*/)
{
  const unsigned nb = (sbrDataCh0 != nullptr ? 2 * ((sbrDataCh0[0] >> 23) & 1) + 2 : 0); // noise bits/ch = 2 or 4
  const unsigned ob = (indepFlag ? 1 : 0);  // indepFlag dependent bit offset
  const bool stereo = (sbrDataCh1 != nullptr);
  const int32_t ch0 = (nb > 0 ? sbrDataCh0[0] : 0);
  const int32_t ch1 = (stereo ? sbrDataCh1[0] : 0);
#if ENABLE_INTERTES
  const bool issTes = (((ch0 >> 30) & 1) != 0);
  const int8_t  res = (ch0 >> 29) & 1; // bs_amp_res
#else
  const int16_t res = ch0 >> 29; // short bs_amp_res
#endif
  const bool couple = (((ch1 >> 23) & 1) != 0);
  unsigned bitCount = (stereo ? (couple ? 2 : 7 + nb) : 0) + 6 + nb;
  unsigned i, envCh0, envCh1; // bs_num_env[], 1 - 8

  if (nb == 0) return 0;

  envCh0 = 1 << ((ch0 >> 21) & 3);
  envCh1 = 1 << (((stereo && !couple ? ch1 : ch0) >> 21) & 3);

  if (stereo) m_auBitStream.write (couple ? 1 : 0, 1); // _coupling

  // sbr_grid(), assumes bs_frame_class[ch] == 0, i.e. class FIXFIX
  m_auBitStream.write ((ch0 >> 20) & 7, 5);
  if (stereo && !couple) m_auBitStream.write ((ch1 >> 20) & 7, 5);

  // sbr_dtdf(), assumes bs_pvc == 0, i.e. no PVC like rest of code
  for (i = ob; i < envCh0; i++) m_auBitStream.write ((ch0 >> (12 + i)) & 1, 1);
  bitCount += i - ob;
  for (i = ob; i < __min (2, envCh0); i++) m_auBitStream.write ((ch0 >> (4 + i)) & 1, 1);
  bitCount += i - ob;

  if (stereo)
  {
    for (i = ob; i < envCh1; i++) m_auBitStream.write ((ch1 >> (12 + i)) & 1, 1);
    bitCount += i - ob;
    for (i = ob; i < __min (2, envCh1); i++) m_auBitStream.write ((ch1 >> (4 + i)) & 1, 1);
    bitCount += i - ob;
  }

  // sbr_invf(), assumes dflt_noise_bands < 3, i.e. 1-2 noise bands
  i = (1 << nb) - 1;
  m_auBitStream.write (ch0 & i, nb); // 2- or 4-bit bs_invf_mode[0]
  if (stereo && !couple) m_auBitStream.write (ch1 & i, nb);

  // sbr_envelope() for mono/left channel, assumes bs_df_env[] == 0
  for (i = 1; i <= envCh0; i++) // dt loop
  {
    const uint8_t bits = ((ch0 & (1 << (11 + i))) != 0 ? 2 : (res > 0 && envCh0 > 1 ? 6 : 7));
    uint8_t nCodedBits = getNumEnvBits (sbrDataCh0[i], 8, bits);

    m_auBitStream.write (sbrDataCh0[i] & ((1 << nCodedBits) - 1), nCodedBits); // 1st env.
    bitCount += nCodedBits;
    nCodedBits = getNumEnvBits (sbrDataCh0[i], 32, 9) - 9;  // avoid writing MSB delimiter

    m_auBitStream.write ((sbrDataCh0[i] >> 8) & ((1 << nCodedBits) - 1), nCodedBits);
    bitCount += nCodedBits;
#if ENABLE_INTERTES
    if (issTes)
    {
      m_auBitStream.write ((sbrDataCh0[9] >> (i - 1)) & 1, 1); // bs_temp_shape[ch][env=i]
      bitCount++;
      if ((sbrDataCh0[9] >> (i - 1)) & 1)
      {
        m_auBitStream.write (GAMMA, 2); // bs_inter_temp_shape_mode
        bitCount += 2;
      }
    }
#endif
  }

  if (stereo && !couple)
  {
    for (i = 1; i <= envCh1; i++) // decoup. sbr_envelope() dt loop
    {
      const uint8_t bits = ((ch1 & (1 << (11 + i))) != 0 ? 2 : (res > 0 && envCh1 > 1 ? 6 : 7));
      uint8_t nCodedBits = getNumEnvBits (sbrDataCh1[i], 8, bits);

      m_auBitStream.write (sbrDataCh1[i] & ((1 << nCodedBits) - 1), nCodedBits);
      bitCount += nCodedBits;
      nCodedBits = getNumEnvBits (sbrDataCh1[i], 32, 9) - 9;

      m_auBitStream.write ((sbrDataCh1[i] >> 8) & ((1 << nCodedBits) - 1), nCodedBits);
      bitCount += nCodedBits;
#if ENABLE_INTERTES
      if (issTes)
      {
        m_auBitStream.write ((sbrDataCh1[9] >> (i - 1)) & 1, 1); // bs_temp_shape[ch][env]
        bitCount++;
        if ((sbrDataCh1[9] >> (i - 1)) & 1)
        {
          m_auBitStream.write (GAMMA, 2);
          bitCount += 2;
        }
      }
#endif
    }
  }

  // sbr_noise() for mono/left channel, assumes bs_df_noise[i] == 0
  for (i = 1; i <= __min (2, envCh0); i++) // dt loop
  {
    const uint8_t bits = ((ch0 & (1 << (3 + i))) != 0 ? 1 : 5);

    m_auBitStream.write ((sbrDataCh0[9] >> (13 * i)) & 31, bits); // _noise[]
    bitCount += bits;
    if (nb == 4)
    {
      m_auBitStream.write ((sbrDataCh0[9] >> (13 * i - 5)) & 31, 1);
      bitCount++;
    }
  }

  if (stereo)
  {
    if (couple)
    {
      for (i = 1; i <= envCh1; i++) // coup. sbr_envelope() dt loop
      {
        const uint8_t bits = ((ch1 & (1 << (11 + i))) != 0 ? 1 : (res > 0 && envCh1 > 1 ? 5 : 6));
        uint8_t nCodedBits = getNumEnvBits (sbrDataCh1[i], 8, bits);

        m_auBitStream.write (sbrDataCh1[i] & ((1 << nCodedBits) - 1), nCodedBits);
        bitCount += nCodedBits;
        nCodedBits = getNumEnvBits (sbrDataCh1[i], 32, 9) - 9;

        m_auBitStream.write ((sbrDataCh1[i] >> 8) & ((1 << nCodedBits) - 1), nCodedBits);
        bitCount += nCodedBits;
#if ENABLE_INTERTES
        if (issTes)
        {
          m_auBitStream.write ((sbrDataCh1[9] >> (i - 1)) & 1, 1); // bs_temp_shape[ch][i]
          bitCount++;
          if ((sbrDataCh1[9] >> (i - 1)) & 1)
          {
            m_auBitStream.write (GAMMA, 2);
            bitCount += 2;
          }
        }
#endif
      }
    }

    for (i = 1; i <= __min (2, envCh1); i++) // sbr_noise() dt loop
    {
      const uint8_t bits = ((ch1 & (1 << (3 + i))) != 0 ? 1 : 5);

      m_auBitStream.write ((sbrDataCh1[9] >> (13 * i)) & 31, bits);
      bitCount += bits;
      if (nb == 4)
      {
        m_auBitStream.write ((sbrDataCh1[9] >> (13 * i - 5)) & 31, 1);
        bitCount++;
      }
    }
  }

  m_auBitStream.write (0, 1);  // fixed bs_add_harmonic_flag[0] = 0
  if (stereo) m_auBitStream.write (0, 1);

  return bitCount;
}

unsigned BitStreamWriter::writeChannelWiseTnsData (const TnsData& tnsData, const bool eightShorts)
{
  const unsigned numWindows = (eightShorts ? 8 : 1);
  const unsigned offsetBits = (eightShorts ? 1 : 2);
  unsigned bitCount = 0, i;

  for (unsigned n = 0, w = 0; w < numWindows; w++)
  {
    bitCount += offsetBits;
    if ((n >= 3) || ((tnsData.firstTnsWindow & (1u << w)) == 0))
    {
      m_auBitStream.write (0/*n_filt[w] = 0*/, offsetBits);
    }
    else // first, second or third length-1 window group in frame and channel
    {
      const unsigned numFiltersInWindow = tnsData.numFilters[n];

      m_auBitStream.write (numFiltersInWindow, offsetBits);
      if (numFiltersInWindow > 0)
      {
        m_auBitStream.write (tnsData.coeffResLow[n] ? 0 : 1, 1);  // coef_res
        bitCount++;
        for (unsigned f = 0; f < numFiltersInWindow; f++)
        {
          const unsigned order = tnsData.filterOrder[n + f];

          m_auBitStream.write (tnsData.filterLength[n + f], 2 + offsetBits * 2);
          m_auBitStream.write (order, 2 + offsetBits);
          bitCount += 4 + offsetBits * 3;
          if (order > 0)
          {
            const int8_t* coeff = tnsData.coeff[n + f];
            unsigned   coefBits = (tnsData.coeffResLow[n] ? 3 : 4);
            int8_t coefMaxValue = (tnsData.coeffResLow[n] ? 2 : 4);
            bool   dontCompress = false;

            m_auBitStream.write (tnsData.filterDownward[n + f] ? 1 : 0, 1);
            for (i = 0; i < order; i++) // get coef_compress, then write coef
            {
              dontCompress |= ((coeff[i] < -coefMaxValue) || (coeff[i] >= coefMaxValue));
            }
            m_auBitStream.write (dontCompress ? 0 : 1, 1);
            coefMaxValue <<= 1;
            if (dontCompress) coefMaxValue <<= 1; else coefBits--;
            for (i = 0; i < order; i++)
            {
              m_auBitStream.write (unsigned (coeff[i] < 0 ? coefMaxValue + coeff[i] : coeff[i]), coefBits);
            }
            bitCount += 2 + order * coefBits;
          }
        }
      } // if n_filt[w] > 0
      n++;
    }
  } // for w

  return bitCount;
}

unsigned BitStreamWriter::writeFDChannelStream (const CoreCoderData& elData, EntropyCoder& entrCoder, const unsigned ch,
                                                const int32_t* const mdctSignal, const uint8_t* const mdctQuantMag,
#if !RESTRICT_TO_AAC
                                                const bool timeWarping, const bool noiseFilling, uint8_t* ipfAuState,
#endif
                                                const bool indepFlag /*= false*/)
{
  const IcsInfo&  icsInfo = elData.icsInfoCurr[ch];
  const TnsData&  tnsData = elData.tnsData[ch];
  const SfbGroupData& grp = elData.groupingData[ch];
  const unsigned   maxSfb = grp.sfbsPerGroup;
  const bool  eightShorts = icsInfo.windowSequence == EIGHT_SHORT;
  uint8_t* const sf = (uint8_t* const) grp.scaleFactors;
  uint8_t sfIdxPred = CLIP_UCHAR (sf[0] > SCHAR_MAX ? 0 : sf[0] + (eightShorts ? 68 : 80));
  unsigned bitCount = 8, g, b, i;

  m_auBitStream.write (sfIdxPred, 8);  // adjusted global_gain
#if !RESTRICT_TO_AAC
  if (noiseFilling)
  {
    m_auBitStream.write (elData.specFillData[ch], 8); // noise level | offset
    bitCount += 8;
  }
#endif
  if (!elData.commonWindow)
  {
    bitCount += writeChannelWiseIcsInfo (icsInfo); // ics_info
  }
#if !RESTRICT_TO_AAC
  if (timeWarping) // && (!common_tw)
  {
    m_auBitStream.write (0, 1); // enforce tw_data_present = 0
    bitCount++;
  }
#endif
  sfIdxPred = sf[0]; // scale factors
  for (g = 0; g < grp.numWindowGroups; g++)
  {
    uint8_t* const gSf = &sf[m_numSwbShort * g];

    for (b = 0; b < maxSfb; b++)
    {
      uint8_t sfIdx = gSf[b];

      if ((g + 1 < grp.numWindowGroups) && (b + 1 == maxSfb) && ((unsigned) sfIdx + INDEX_OFFSET < sf[m_numSwbShort * (g + 1)]))
      {
        // ugly, avoidable if each gr. had its own global_gain
        gSf[b] = sfIdx = sf[m_numSwbShort * (g + 1)] - INDEX_OFFSET;
      }
      if ((g > 0) || (b > 0))
      {
        int sfIdxDpcm = (int) sfIdx - sfIdxPred;
        unsigned sfBits;

        if (sfIdxDpcm > INDEX_OFFSET) // just as sanity checks
        {
          sfIdxDpcm =  INDEX_OFFSET;
          sfIdxPred += INDEX_OFFSET;
        }
        else if (sfIdxDpcm < -INDEX_OFFSET) // highly unlikely
        {
          sfIdxDpcm = -INDEX_OFFSET;
          sfIdxPred -= INDEX_OFFSET;
        }
        else // scale factor range OK
        {
          sfIdxPred = sfIdx;
        }
        sfBits = entrCoder.indexGetBitCount (sfIdxDpcm);
        m_auBitStream.write (entrCoder.indexGetHuffCode (sfIdxDpcm), sfBits);
        bitCount += sfBits;
      }
    }
  } // for g

  if (!elData.commonTnsData && (tnsData.numFilters[0] + tnsData.numFilters[1] + tnsData.numFilters[2] > 0))
  {
    bitCount += writeChannelWiseTnsData (tnsData, eightShorts);
  }

  bitCount += (indepFlag ? 1 : 2); // arith_reset_flag, fac_data_present bits

  if (maxSfb == 0) // zeroed spectrum
  {
    entrCoder.initWindowCoding (!eightShorts /*reset*/, eightShorts);

    if (!indepFlag) m_auBitStream.write (1, 1); // force reset
#ifndef NO_PREROLL_DATA
    if (ipfAuState) memset (ipfAuState, 0, 4);  // no spectrum
#endif
  }
  else // not zeroed, nasty since SFB ungrouping may be needed
  {
    const uint16_t* grpOff = grp.sfbOffsets;
    uint8_t grpLen = grp.windowGroupLength[0];
    uint8_t grpWin = 0;
    uint8_t swbSize[MAX_NUM_SWB_SHORT];
    const uint8_t* winMag = (grpLen > 1 ? m_uCharBuffer : mdctQuantMag);
    const uint16_t lg     = (grpLen > 1 ? grpOff[maxSfb] / grpLen : grpOff[maxSfb]);

    if (eightShorts || (grpLen > 1)) // ungroup the SFB widths
    {
      for (b = 0, i = oneTwentyEightOver[grpLen]; b < maxSfb; b++)
      {
        swbSize[b] = ((grpOff[b+1] - grpOff[b]) * i) >> 7; // sfbWidth/grpLen
      }
    }
    g = 0;
    for (int w = 0; w < (eightShorts ? 8 : 1); w++, grpWin++)  // window loop
    {
      if (grpWin >= grpLen) // next g
      {
        grpOff += m_numSwbShort;
        grpLen = grp.windowGroupLength[++g];
        grpWin = 0;
        winMag = (grpLen > 1 ? m_uCharBuffer : &mdctQuantMag[grpOff[0]]);
      }
      if (eightShorts && (grpLen > 1))
      {
        for (b = i = 0; b < maxSfb; b++) // ungroup magnitudes
        {
          memcpy (&m_uCharBuffer[i], &mdctQuantMag[grpOff[b] + grpWin * swbSize[b]], swbSize[b] * sizeof (uint8_t));
          i += swbSize[b];
        }
      }
      entrCoder.initWindowCoding (indepFlag && (w == 0), eightShorts);

      if (!indepFlag && (w == 0)) // optimize arith_reset_flag
      {
        if ((b = entrCoder.arithGetResetBit (winMag, 0, lg)) != 0)
        {
          entrCoder.arithResetMemory ();
          entrCoder.arithSetCodState (USHRT_MAX << 16);
          entrCoder.arithSetCtxState (0);
        }
        m_auBitStream.write (b, 1); // write adapted reset bit
      }
#ifndef NO_PREROLL_DATA
      if (ipfAuState && (w == 0))
      {
        b = (unsigned) m_auBitStream.stream.size ();

        if (eightShorts || (b > 511) || !indepFlag)
        {
          memset (ipfAuState, 0, 4); // grouped or no residual
        }
        else
        {
          const int32_t* const winSig = &mdctSignal[grpOff[0]];
          int32_t sigPk = 0;

          ipfAuState[0] = uint8_t (b >> 1);
          ipfAuState[1] = uint8_t ((b & 1) << 7) | m_auBitStream.heldBitCount;
          ipfAuState[2] = m_auBitStream.heldBitChunk;
          ipfAuState[3] = CLIP_UCHAR (lg >> 2);

          for (b = i = 0; i < __min (256u, lg); i++)
          {
            if ((winMag[i] != 0) && (abs (winSig[i]) > sigPk))
            {
              sigPk = abs (winSig[i]);
              b = i;
            }
          }
          ipfAuState[4] = (uint8_t) b;
          ipfAuState[5] = winMag[b] << 1;
          if (winSig[b] > 0) ipfAuState[5] |= 1; // store sign of single peak
        }
      }
#endif
      bitCount += entrCoder.arithCodeSigMagn (winMag, 0, lg, true, &m_auBitStream);

      if (eightShorts && (grpLen > 1))
      {
        for (b = i = 0; b < maxSfb; b++) // ungroup coef signs
        {
          const int32_t* const swbSig = &mdctSignal[grpOff[b] + grpWin * swbSize[b]];

          for (unsigned j = 0; j < swbSize[b]; j++, i++)
          {
            if (winMag[i] != 0)
            {
              m_auBitStream.write (swbSig[j] < 0 ? 0 : 1, 1); // - = 0, + = 1
              bitCount++;
            }
          }
        }
      }
      else // not grouped long window
      {
        const int32_t* const winSig = &mdctSignal[grpOff[0]];

        for (i = 0; i < lg; i++)
        {
          if (winMag[i] != 0)
          {
            m_auBitStream.write (winSig[i] < 0 ? 0 : 1, 1); // -1 = 0, +1 = 1
            bitCount++;
          }
        }
      }
    } // for w
  } // if maxSfb == 0

  m_auBitStream.write (0, 1); // fac_data_present, no fac_data

  return bitCount;
}

unsigned BitStreamWriter::writeStereoCoreToolInfo (const CoreCoderData& elData, EntropyCoder& entrCoder,
#if !RESTRICT_TO_AAC
                                                   const bool timeWarping, bool* const commonTnsFlag,
#endif
                                                   const bool indepFlag /*= false*/)
{
  const IcsInfo& icsInfo0 = elData.icsInfoCurr[0];
  const IcsInfo& icsInfo1 = elData.icsInfoCurr[1];
  const TnsData& tnsData0 = elData.tnsData[0];
  const TnsData& tnsData1 = elData.tnsData[1];
  const uint16_t nWinGrps = elData.groupingData[0].numWindowGroups;
  const bool eightShorts0 = (icsInfo0.windowSequence == EIGHT_SHORT);
  unsigned bitCount = 2, g, b;

  m_auBitStream.write (elData.tnsActive ? 1 : 0, 1); // tns_active
  m_auBitStream.write (elData.commonWindow ? 1 : 0, 1);
  if (elData.commonWindow)
  {
    const unsigned maxSfbSte = __max (icsInfo0.maxSfb, icsInfo1.maxSfb);
    const unsigned  sfb1Bits = (eightShorts0 ? 4 : 6);
    const uint8_t msMaskMode = getOptMsMaskModeValue (elData.stereoDataCurr, nWinGrps, m_numSwbShort, elData.stereoMode, maxSfbSte);

    bitCount += writeChannelWiseIcsInfo (icsInfo0);  // ics_info()
    m_auBitStream.write (elData.commonMaxSfb ? 1 : 0, 1);
    if (!elData.commonMaxSfb)
    {
      m_auBitStream.write (icsInfo1.maxSfb, sfb1Bits); // max_sfb1
      bitCount += sfb1Bits;
    }
    m_auBitStream.write (__min (3, msMaskMode), 2); // ms_mask_pr.
    bitCount += 3;
    if (msMaskMode == 1)  // some M/S, write SFB-wise ms_used flag
    {
      for (g = 0; g < nWinGrps; g++)
      {
        const uint8_t* const gMsUsed = &elData.stereoDataCurr[m_numSwbShort * g];

        for (b = 0; b < maxSfbSte; b++)
        {
          m_auBitStream.write (gMsUsed[b] > 0 ? 1 : 0, 1);
        }
      }
      bitCount += maxSfbSte * g;
    }
#if !RESTRICT_TO_AAC
    else if (msMaskMode >= 3) // pred. M/S, write cplx_pred_data()
    {
      const bool complexCoef = (elData.stereoConfig & 1);
      uint32_t deltaCodeTime = 0;

      m_auBitStream.write (msMaskMode - 3, 1);    // cplx_pred_all
      if (msMaskMode == 3)
      {
        for (g = 0; g < nWinGrps; g++)
        {
          const uint8_t* const gCplxPredUsed = &elData.stereoDataCurr[m_numSwbShort * g];

          for (b = 0; b < maxSfbSte; b += SFB_PER_PRED_BAND)
          {
            m_auBitStream.write (gCplxPredUsed[b] > 0 ? 1 : 0, 1);
          }
        }
        bitCount += ((maxSfbSte + 1) / SFB_PER_PRED_BAND) * g;
      }
      m_auBitStream.write (elData.stereoConfig & 3, 2);// pred_dir
      bitCount += 3;
      if (!indepFlag)  // write use_prev_frame and delta_code_time
      {
        if (complexCoef)
        {
          m_auBitStream.write (elData.stereoConfig & 4 ? 1 : 0, 1);
          bitCount++;
        }
#ifndef NO_WORKAROUND_FOR_APPLE_ISSUE_FB8928108
        if ((eightShorts0 && elData.icsInfoPrev[0].windowSequence != EIGHT_SHORT) ||  // first ch. in CPE
            (elData.icsInfoPrev[0].windowSequence == EIGHT_SHORT && !eightShorts0) ||
            (eightShorts0 && elData.icsInfoPrev[1].windowSequence != EIGHT_SHORT) || // second ch. in CPE
            (elData.icsInfoPrev[1].windowSequence == EIGHT_SHORT && !eightShorts0))
        {
          deltaCodeTime = 0;
        }
        else
#endif
        deltaCodeTime = getDeltaCodeTimeFlag (elData.stereoDataCurr, nWinGrps, m_numSwbShort, elData.stereoDataPrev, maxSfbSte, entrCoder, complexCoef);
        m_auBitStream.write (deltaCodeTime, 1);
        bitCount++;
      }

      for (g = 0; g < nWinGrps; g++)
      {
        const uint8_t* const aqReIdxPrvGrp = (g == 0 ? elData.stereoDataPrev : &elData.stereoDataCurr[m_numSwbShort * (g - 1)]);
        const uint8_t* const aqImIdxPrvGrp = &aqReIdxPrvGrp[1];
        const uint8_t* const gCplxPredUsed = &elData.stereoDataCurr[m_numSwbShort * g];
        int  aqReIdxPred = 16, aqImIdxPred = 16; // alpha_q_.. = 0

        for (b = 0; b < maxSfbSte; b += SFB_PER_PRED_BAND)
        {
          if (gCplxPredUsed[b] > 0) // write dpcm_alpha_q_re/_q_im
          {
            int aqIdx = gCplxPredUsed[b] & 31; // range -15,...,15
            int aqIdxDpcm = aqIdx - (deltaCodeTime > 0 ? getPredCoefPrevGrp (aqReIdxPrvGrp[b]) : aqReIdxPred);
            unsigned bits = entrCoder.indexGetBitCount (aqIdxDpcm);

            if (deltaCodeTime == 0) aqReIdxPred = aqIdx;
            m_auBitStream.write (entrCoder.indexGetHuffCode (aqIdxDpcm), bits);
            bitCount += bits;

            if (complexCoef)
            {
              aqIdx = gCplxPredUsed[b + 1] & 31; // <32 kHz short!
              aqIdxDpcm = aqIdx - (deltaCodeTime > 0 ? getPredCoefPrevGrp (aqImIdxPrvGrp[b]) : aqImIdxPred);
              bits = entrCoder.indexGetBitCount (aqIdxDpcm);

              if (deltaCodeTime == 0) aqImIdxPred = aqIdx;
              m_auBitStream.write (entrCoder.indexGetHuffCode (aqIdxDpcm), bits);
              bitCount += bits;
            }
          }
          else if (deltaCodeTime == 0) aqReIdxPred = aqImIdxPred = 16;
        }
      } // for g
    }
#endif
  } // common_window
#if !RESTRICT_TO_AAC
  if (timeWarping)
  {
    m_auBitStream.write (0, 1); // common_tw not needed in BL USAC
    bitCount++;
  } // tw_mdct
#endif
  if (elData.tnsActive)
  {
    bool commonTns = elData.commonTnsData;

    if (elData.commonWindow)
    {
#if !RESTRICT_TO_AAC
      if ((commonTnsFlag != nullptr) && !commonTns)  // common_tns
      {
        const uint8_t* data1 = (uint8_t*) &tnsData0; // fast comp.
        const uint8_t* data2 = (uint8_t*) &tnsData1; // portable??

        commonTns = true;
        for (b = 0; b < sizeof (TnsData); b++) commonTns &= (data1[b] == data2[b]);
        *commonTnsFlag = commonTns;
      }
#endif
      m_auBitStream.write (/*optim.*/commonTns ? 1 : 0, 1);
      bitCount++;
    }
    m_auBitStream.write (elData.tnsOnLeftRight ? 1 : 0, 1);
    bitCount++;
    if (commonTns)
    {
      bitCount += writeChannelWiseTnsData (tnsData0, eightShorts0);
    }
    else  // tns_present_both and tns_data_present[1]
    {
      const bool tnsPresentBoth = (tnsData0.numFilters[0] + tnsData0.numFilters[1] + tnsData0.numFilters[2] > 0) &&
                                  (tnsData1.numFilters[0] + tnsData1.numFilters[1] + tnsData1.numFilters[2] > 0);
      m_auBitStream.write (tnsPresentBoth ? 1 : 0, 1);
      bitCount++;
      if (!tnsPresentBoth)
      {
        m_auBitStream.write (tnsData1.numFilters[0] + tnsData1.numFilters[1] + tnsData1.numFilters[2] > 0 ? 1 : 0, 1);
        bitCount++;
      }
    }
  } // tns_active

  return bitCount;
}

// public functions
unsigned BitStreamWriter::createAudioConfig (const char samplingFrequencyIndex,  const bool shortFrameLength,
                                             const uint8_t chConfigurationIndex, const uint8_t numElements,
                                             const ELEM_TYPE* const elementType, const uint32_t loudnessInfo,
#if !RESTRICT_TO_AAC
                                             const bool* const tw_mdct /*N/A*/,  const bool* const noiseFilling,
#endif
                                             const uint8_t sbrRatioShiftValue,   unsigned char* const audioConfig)
{
  const uint8_t fli = (sbrRatioShiftValue == 0 ? 1 /*no SBR*/ : __min (2, sbrRatioShiftValue & 3) + 2);
  const int8_t usfi = __max (0, samplingFrequencyIndex - 3 * (sbrRatioShiftValue & 3)); // TODO: nonstandard rates
  unsigned bitCount = 37, auLen;
#ifndef NO_PREROLL_DATA
  unsigned ucOffset = (samplingFrequencyIndex < AAC_NUM_SAMPLE_RATES ? 2 : 5);
#endif

  if ((elementType == nullptr) || (audioConfig == nullptr) || (chConfigurationIndex >= USAC_MAX_NUM_ELCONFIGS) ||
#if !RESTRICT_TO_AAC
      (noiseFilling == nullptr) || (tw_mdct == nullptr) ||
#endif
      (numElements == 0) || (numElements > USAC_MAX_NUM_ELEMENTS) || (samplingFrequencyIndex < 0) || (samplingFrequencyIndex >= 0x1F))
  {
    return 0; // invalid arguments error
  }

  m_auBitStream.reset ();
// --- AudioSpecificConfig(): https://wiki.multimedia.cx/index.php/MPEG-4_Audio/
  m_auBitStream.write (0x7CA, 11); // audio object type (AOT) 32 (esc) + 10 = 42
  if (samplingFrequencyIndex < AAC_NUM_SAMPLE_RATES)
  {
    m_auBitStream.write (usfi, 4);
  }
  else
  {
    m_auBitStream.write (0xF, 4); // esc
    m_auBitStream.write (toSamplingRate (usfi), 24);
    bitCount += 24;
  }
  // for multichannel audio, refer to channel mapping of AotSpecificConfig below
  m_auBitStream.write (chConfigurationIndex > 2 ? 0 : chConfigurationIndex, 4);

// --- AotSpecificConfig(): UsacConfig()
  m_auBitStream.write (usfi, 5); // usacSamplingFrequencyIndex (after SBR dec.!)
  m_auBitStream.write (shortFrameLength ? 0 : fli, 3);// coreSbrFrameLengthIndex
  m_auBitStream.write (chConfigurationIndex, 5);    // channelConfigurationIndex
#ifdef NO_PREROLL_DATA
  m_auBitStream.write (numElements - 1, 4);  // numElements in UsacDecoderConfig
#else
  m_auBitStream.write (numElements, 4); // 4bit numElements in UsacDecoderConfig

  m_auBitStream.write (ID_USAC_EXT, 2); // usacElementType[0] = 3, for IPF stuff
  m_auBitStream.write (3, 4); // UsacExtElementConfig(), ID_EXT_ELE_AUDIOPREROLL
  m_auBitStream.write (0, 6); // usacExtElementConfigLength = 0, rest of config.
  bitCount += 12;
#endif
  for (unsigned el = 0; el < numElements; el++) // el element loop
  {
    m_auBitStream.write ((unsigned) elementType[el], 2);  // usacElementType[el]
    bitCount += 2;
    if (elementType[el] < ID_USAC_LFE) // SCE, CPE: UsacCoreConfig
    {
#if RESTRICT_TO_AAC
      m_auBitStream.write (0, 2);  // time warping and noise filling not allowed
#else
      m_auBitStream.write ((tw_mdct[el] ? 2 : 0) | (noiseFilling[el] ? 1 : 0), 2);
#endif
      bitCount += 2;
      if (sbrRatioShiftValue > 0)  // sbrRatioIndex > 0: SbrConfig
      {
        const uint32_t sf = (samplingFrequencyIndex == 6 || samplingFrequencyIndex < 5 ? 10 : (samplingFrequencyIndex < 8 ? 9 : 8)); // bs_stop_freq
#if ENABLE_INTERTES
        m_auBitStream.write (2, 3);  // bs_interTes = 1, harmonicSBR, bs_pvc = 0
#else
        m_auBitStream.write (0, 3);  // fix harmonicSBR, bs_interTes, bs_pvc = 0
#endif
        bitCount += 13; // incl. SbrDfltHeader following hereafter
        m_auBitStream.write (15 - (sbrRatioShiftValue / 4), 4); // bs_start_freq
        m_auBitStream.write (sf, 4); // 16193 @ 44.1, 18375 @ 48, 22500 @ 64 kHz
        if (loudnessInfo >> 30)
        {
          m_auBitStream.write (2, 2);// set dflt_header_extra1 = 1
          m_auBitStream.write (2 + (loudnessInfo >> 31), 2);
          m_auBitStream.write (4 | ((loudnessInfo >> 29) & 2), 3);
          bitCount += 5;
        }
        else m_auBitStream.write (0, 2); // dflt_header_extra* = 0

        if (elementType[el] == ID_USAC_CPE)
        {
          m_auBitStream.write (0, 2); // fix stereoConfigIndex = 0
          bitCount += 2;
        }
      }
    }
  } // for el

  m_auBitStream.write (loudnessInfo > 0 ? 1 : 0, 1); // ..ConfigExtensionPresent
  if (loudnessInfo > 0) // ISO 23003-4: loudnessInfo()
  {
    const unsigned methodDefinition = (loudnessInfo >> 14) & 0xF;
    const unsigned methodValueBits  = (methodDefinition == 7 ? 5 : (methodDefinition == 8 ? 2 : 8));

    m_auBitStream.write (0, 2); // numConfigExtensions
    m_auBitStream.write (2, 4); // ..EXT_LOUDNESS_INFO
    m_auBitStream.write (methodValueBits < 3 ? 7 : 8, 4); // usacConfigExtLength

    m_auBitStream.write (1, 12);// loudnessInfoCount=1
    m_auBitStream.write (1, 14);// samplePeakLevel..=1
    m_auBitStream.write ((loudnessInfo >> 18) & 0xFFF, 12); // bsSamplePeakLevel
    m_auBitStream.write (1, 5);  // measurementCount=1
    m_auBitStream.write (methodDefinition, 4);
    m_auBitStream.write ((loudnessInfo >> 6) & ((1 << methodValueBits) - 1), methodValueBits);
    m_auBitStream.write ((loudnessInfo >> 2) & 0xF, 4);     // measurementSystem
    m_auBitStream.write ((loudnessInfo & 0x3), 2);  // reliability, 3 = accurate

    m_auBitStream.write (0, 1);  // loudnessInfoSetExtPresent=0, payload padding
    bitCount += (methodValueBits < 3 ? 66 : 74);
    if (methodValueBits >= 3) m_auBitStream.write (0, 10 - methodValueBits);
  }

  bitCount += (8 - m_auBitStream.heldBitCount) & 7;
  writeByteAlignment ();  // flush bytes
  auLen = __min (18u + fli, bitCount >> 3);
#ifndef NO_PREROLL_DATA
  m_usacConfigLen = uint16_t (__max (15, auLen - ucOffset)); // excl ASC payload
  memcpy (m_usacConfig, &m_auBitStream.stream.at (ucOffset), auLen - ucOffset);
#endif
  memcpy (audioConfig,  &m_auBitStream.stream.front (), auLen);

  return (bitCount >> 3);  // byte count
}

unsigned BitStreamWriter::createAudioFrame (CoreCoderData** const elementData,  EntropyCoder* const entropyCoder,
                                            int32_t** const mdctSignals,        uint8_t** const mdctQuantMag,
                                            const bool usacIndependencyFlag,    const uint8_t numElements,
                                            const uint8_t numSwbShort,          uint8_t* const tempBuffer,
#if !RESTRICT_TO_AAC
                                            const bool* const tw_mdct /*N/A*/,  const bool* const noiseFilling,
                                            const uint32_t frameCount,          const uint32_t indepPeriod,  uint32_t* rate,
#endif
                                            const uint8_t sbrRatioShiftValue,   int32_t** const sbrInfoAndData,
                                            unsigned char* const accessUnit,    const unsigned nSamplesInFrame)
{
#ifndef NO_PREROLL_DATA
  const uint8_t ipf = (frameCount == 1 ? 2 : ((frameCount % (indepPeriod << 1)) == 1 ? 1 : 0));
#endif
#if !RESTRICT_TO_AAC
  uint8_t* ipfState = (frameCount > 0 && (frameCount % (indepPeriod << 1)) == 0 && numElements == 1 ? m_usacIpfState : nullptr);
#endif
  unsigned bitCount = 1, ci = 0;

  if ((elementData == nullptr) || (entropyCoder == nullptr) || (tempBuffer == nullptr) || (sbrInfoAndData == nullptr) ||
      (mdctSignals == nullptr) || (mdctQuantMag == nullptr) || (accessUnit == nullptr) || (nSamplesInFrame > 2048) ||
#if !RESTRICT_TO_AAC
      (noiseFilling == nullptr) || (tw_mdct == nullptr) ||
# ifndef NO_PREROLL_DATA
      (ipf && !usacIndependencyFlag) ||
# endif
#endif
      (numElements == 0) || (numElements > USAC_MAX_NUM_ELEMENTS) || (numSwbShort < MIN_NUM_SWB_SHORT) || (numSwbShort > MAX_NUM_SWB_SHORT))
  {
    return 0; // invalid arguments error
  }
#ifndef NO_PREROLL_DATA
  if (ipf)
  {
    bitCount = ((ipf == 2) || (ipf == 1 && (numElements > 1 || !noiseFilling[0]))
                ? __min (nSamplesInFrame << 2, (unsigned) m_auBitStream.stream.size ())
                : ((unsigned) m_usacIpfState[0] << 1) | (m_usacIpfState[1] >> 7));
    memcpy (tempBuffer, &m_auBitStream.stream.front (), bitCount); // prev fr AU
  }
#endif
  m_auBitStream.reset ();
  m_numSwbShort = numSwbShort;
  m_uCharBuffer = tempBuffer;
  m_auBitStream.write (usacIndependencyFlag ? 1 : 0, 1);

#ifndef NO_PREROLL_DATA
  m_auBitStream.write (ipf ? 1 : 0, 1); // UsacExtElement, usacExtElementPresent
  if (ipf)
  {
    const bool lowRatePreRollExt = (ipf == 1 && numElements == 1 && noiseFilling[0]);
    const unsigned   extraLength = (m_usacConfigLen > 14 ? 4 : 3) + m_usacConfigLen;
    const unsigned payloadLength = (lowRatePreRollExt ? getLowRatePreRollAU (tempBuffer, *elementData[0], entropyCoder[0],
                                    m_usacIpfState, sbrRatioShiftValue) : bitCount) + extraLength; // in bytes
    m_auBitStream.write (0, 1); // usacExtElementUseDefaultLength = 0 (variable)
    m_auBitStream.write (CLIP_UCHAR (payloadLength), 8);
    if (payloadLength > 254) m_auBitStream.write (payloadLength - 253, 16);

    m_auBitStream.write (__min (15, m_usacConfigLen), 4); // configLen (part #1)
    if (m_usacConfigLen > 14) m_auBitStream.write (m_usacConfigLen - 15, 4);

    m_auBitStream.write (m_usacConfig[ci++] & 31, 5); // 1st 3 bits are from ASC
    while (ci < m_usacConfigLen) m_auBitStream.write (m_usacConfig[ci++], 8);
    ci = 0;
    m_auBitStream.write (0, 8 - 5); // pad end of UsacConfig() data

    m_auBitStream.write (0, 2); // applyCrossfade = 0 and reserved = 0 (part #2)
    m_auBitStream.write (1, 2); // numPreRollFrames, only one supported for now!
    m_auBitStream.write (payloadLength - extraLength, 16); // auLen

    if (lowRatePreRollExt) bitCount = payloadLength - extraLength;
    while (ci < bitCount) m_auBitStream.write (tempBuffer[ci++], 8); // write AU
    ci = 0;
    if (m_usacConfigLen > 14) m_auBitStream.write (0, 4); // pad end of ext data

    bitCount = (payloadLength > 254 ? 26 : 10) + (payloadLength << 3); // for PR
  }
  bitCount++; // for ElementPresent flag
#endif // !NO_PREROLL_DATA
  for (unsigned el = 0; el < numElements; el++) // el element loop
  {
    const CoreCoderData* const elData = elementData[el];

    if (elData == nullptr)
    {
      return 0; // internal memory error
    }
    switch (elData->elementType)  // write out UsacCoreCoderData()
    {
      case ID_USAC_SCE: // UsacSingleChannelElement()
      {
        m_auBitStream.write (CORE_MODE_FD, 1);
        m_auBitStream.write (elData->tnsActive ? 1 : 0, 1);  // tns_data_present
        bitCount += 2;
        bitCount += writeFDChannelStream (*elData, entropyCoder[ci], 0,
                                          mdctSignals[ci], mdctQuantMag[ci],
#if !RESTRICT_TO_AAC
                                          tw_mdct[el], noiseFilling[el], ipfState,
#endif
                                          usacIndependencyFlag);
        if (sbrRatioShiftValue > 0) // UsacSbrData()
        {
          if (usacIndependencyFlag)
          {
            m_auBitStream.write ((sbrInfoAndData[ci][0] >> 24) & 63, 6);  // SbrInfo(), bs_pvc = 0
            m_auBitStream.write (1, 1);// fix sbrUseDfltHeader = 1
            bitCount += 7;
          }
          else
          {
            m_auBitStream.write (0, 1);  // fix sbrInfoPresent = 0
            bitCount++;
          }
          bitCount += writeChannelWiseSbrData (sbrInfoAndData[ci], nullptr, // L (mono) only, no R
                                               usacIndependencyFlag);
        }
        ci++;
        break;
      }
      case ID_USAC_CPE: // UsacChannelPairElement()
      {
        m_auBitStream.write (CORE_MODE_FD, 1); // L
        m_auBitStream.write (CORE_MODE_FD, 1); // R
        bitCount += 2;
        bitCount += writeStereoCoreToolInfo (*elData, entropyCoder[ci], // L
#if !RESTRICT_TO_AAC
                                             tw_mdct[el], &elementData[el]->commonTnsData,
#endif
                                             usacIndependencyFlag);
        bitCount += writeFDChannelStream (*elData, entropyCoder[ci], 0, // L
                                          mdctSignals[ci], mdctQuantMag[ci],
#if !RESTRICT_TO_AAC
                                          tw_mdct[el], noiseFilling[el], nullptr,
#endif
                                          usacIndependencyFlag);
        ci++;
        bitCount += writeFDChannelStream (*elData, entropyCoder[ci], 1, // R
                                          mdctSignals[ci], mdctQuantMag[ci],
#if !RESTRICT_TO_AAC
                                          tw_mdct[el], noiseFilling[el], ipfState,
#endif
                                          usacIndependencyFlag);
        if (sbrRatioShiftValue > 0) // UsacSbrData()
        {
          if (usacIndependencyFlag)
          {
            m_auBitStream.write ((sbrInfoAndData[ci][0] >> 24) & 63, 6);  // SbrInfo(), bs_pvc = 0
            m_auBitStream.write (1, 1);// fix sbrUseDfltHeader = 1
            bitCount += 7;
          }
          else
          {
            m_auBitStream.write (0, 1);  // fix sbrInfoPresent = 0
            bitCount++;
          }
          bitCount += writeChannelWiseSbrData (sbrInfoAndData[ci - 1], sbrInfoAndData[ci], // L, R
                                               usacIndependencyFlag);
        }
        ci++;
        break;
      }
      case ID_USAC_LFE: // UsacLfeElement()
      {
        bitCount += writeFDChannelStream (*elData, entropyCoder[ci], 0,
                                          mdctSignals[ci], mdctQuantMag[ci],
#if !RESTRICT_TO_AAC
                                          false, false, ipfState,
#endif
                                          usacIndependencyFlag);
        ci++;
        break;
      }
      default: break;
    }
  } // for el

  bitCount += (8 - m_auBitStream.heldBitCount) & 7;
  writeByteAlignment ();  // flush bytes

#if RESTRICT_TO_AAC || defined (NO_PREROLL_DATA)
  memcpy (accessUnit, &m_auBitStream.stream.front (), __min (768 * ci, bitCount >> 3));
#else
  m_auByteCount += bitCount >> 3;
  if (rate != nullptr)  // sampling rate
  {
    const double framesPerSec = (double) *rate / nSamplesInFrame;
    const unsigned targetRate = (4 - (sbrRatioShiftValue & 1)) * ci; // frame average for preset 1

    if (framesPerSec > 0.0 && targetRate > 0 && frameCount < UINT_MAX) // running overcoding ratio
    {
#if BA_MORE_CBR
      *rate = uint32_t (0.5 + (m_auByteCount * framesPerSec) / (__max (framesPerSec, (double) frameCount) * targetRate));
#else
      *rate = uint32_t (0.5 + (m_auByteCount * framesPerSec) / (__max (20.0 * framesPerSec, (double) frameCount) * targetRate));
#endif
    }
    else *rate = 0; // insufficient data
  }
  memcpy (accessUnit, &m_auBitStream.stream.front (), __min (ci * (ipf ? 1152 : 768), bitCount >> 3));
#endif
  return (bitCount >> 3);  // byte count
}
