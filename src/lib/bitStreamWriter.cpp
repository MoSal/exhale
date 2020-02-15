/* bitStreamWriter.cpp - source file for class with basic bit-stream writing capability
 * written by C. R. Helmrich, last modified in 2020 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under a Modified BSD-Style License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2020 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#include "exhaleLibPch.h"
#include "bitStreamWriter.h"

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

unsigned BitStreamWriter::writeChannelWiseTnsData (const TnsData& tnsData, const bool eightShorts)
{
  const unsigned numWindows = (eightShorts ? 8 : 1);
  const unsigned offsetBits = (eightShorts ? 1 : 2);
  unsigned bitCount = 0, i;

  for (unsigned w = 0; w < numWindows; w++)
  {
    bitCount += offsetBits;
    if (w != tnsData.filteredWindow)
    {
      m_auBitStream.write (0/*n_filt[w] = 0*/, offsetBits);
    }
    else
    {
      m_auBitStream.write (tnsData.numFilters, offsetBits);
      if (tnsData.numFilters > 0)
      {
        m_auBitStream.write (tnsData.coeffResLow ? 0 : 1, 1);  // coef_res[w]
        bitCount++;
        for (unsigned f = 0; f < tnsData.numFilters; f++)
        {
          const unsigned order = tnsData.filterOrder[f];

          m_auBitStream.write (tnsData.filterLength[f], 2 + offsetBits * 2);
          m_auBitStream.write (order, 2 + offsetBits);
          bitCount += 4 + offsetBits * 3;
          if (order > 0)
          {
            const int8_t* coeff = tnsData.coeff[f];
            unsigned   coefBits = (tnsData.coeffResLow ? 3 : 4);
            char   coefMaxValue = (tnsData.coeffResLow ? 2 : 4);
            bool   dontCompress = false;

            m_auBitStream.write (tnsData.filterDownward[f] ? 1 : 0, 1);
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
      } // if (n_filt[w])
    }
  } // for w

  return bitCount;
}

unsigned BitStreamWriter::writeFDChannelStream (const CoreCoderData& elData, EntropyCoder& entrCoder, const unsigned ch,
                                                const int32_t* const mdctSignal, const uint8_t* const mdctQuantMag,
#if !RESTRICT_TO_AAC
                                                const bool timeWarping, const bool noiseFilling,
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

  if (!elData.commonTnsData && (tnsData.numFilters > 0))
  {
    bitCount += writeChannelWiseTnsData (tnsData, eightShorts);
  }

  bitCount += (indepFlag ? 1 : 2); // arith_reset_flag, fac_data_present bits

  if (maxSfb == 0) // zeroed spectrum
  {
    entrCoder.initWindowCoding (!eightShorts /*reset*/, eightShorts);

    if (!indepFlag) m_auBitStream.write (1, 1); // force reset
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
  } // if (maxSfb == 0)

  m_auBitStream.write (0, 1); // fac_data_present, no fac_data

  return bitCount;
}

unsigned BitStreamWriter::writeStereoCoreToolInfo (const CoreCoderData& elData,
#if !RESTRICT_TO_AAC
                                                   const bool timeWarping,
#endif
                                                   const bool indepFlag /*= false*/)
{
  const IcsInfo& icsInfo0 = elData.icsInfoCurr[0];
  const IcsInfo& icsInfo1 = elData.icsInfoCurr[1];
  const TnsData& tnsData0 = elData.tnsData[0];
  const TnsData& tnsData1 = elData.tnsData[1];
  unsigned bitCount = 2, g, b;

  m_auBitStream.write (elData.tnsActive ? 1 : 0, 1); // tns_active
  m_auBitStream.write (elData.commonWindow ? 1 : 0, 1);
  if (elData.commonWindow)
  {
    const unsigned maxSfbSte = __max (icsInfo0.maxSfb, icsInfo1.maxSfb);
    const unsigned sfb1Bits  = icsInfo1.windowSequence == EIGHT_SHORT ? 4 : 6;

    bitCount += writeChannelWiseIcsInfo (icsInfo0);  // ics_info()
    m_auBitStream.write (elData.commonMaxSfb ? 1 : 0, 1);
    if (!elData.commonMaxSfb)
    {
      m_auBitStream.write (icsInfo1.maxSfb, sfb1Bits); // max_sfb1
      bitCount += sfb1Bits;
    }
    m_auBitStream.write (__min (3, elData.stereoMode), 2); // ms_mask_present
    bitCount += 3;
    if (elData.stereoMode == 1) // write SFB-wise ms_used[][] flag
    {
      for (g = 0; g < elData.groupingData[0].numWindowGroups; g++)
      {
        const char* const gMsUsed = &elData.stereoData[m_numSwbShort * g];

        for (b = 0; b < maxSfbSte; b++)
        {
          m_auBitStream.write (gMsUsed[b] > 0 ? 1 : 0, 1);
        }
      }
      bitCount += maxSfbSte * g;
    }
#if !RESTRICT_TO_AAC
    else if (elData.stereoMode >= 3)  // SFB-wise cplx_pred_data()
    {
      m_auBitStream.write (elData.stereoMode - 3, 1); // _pred_all
      if (elData.stereoMode == 3)
      {
        for (g = 0; g < elData.groupingData[0].numWindowGroups; g++)
        {
          const char* const gCplxPredUsed = &elData.stereoData[m_numSwbShort * g];

          for (b = 0; b < maxSfbSte; b += SFB_PER_PRED_BAND)
          {
            m_auBitStream.write (gCplxPredUsed[b] > 0 ? 1 : 0, 1);
          }
        }
        bitCount += ((maxSfbSte + 1) / SFB_PER_PRED_BAND) * g;
      }
      // pred_dir and complex_coef. TODO: rest of cplx_pred_data()
      m_auBitStream.write (elData.stereoConfig & 3, 2);
      bitCount += 3;
    }
#endif
  } // common_window
#if !RESTRICT_TO_AAC
  if (timeWarping)
  {
    m_auBitStream.write (0, 1); // common_tw not needed in xHE-AAC
    bitCount++;
  } // tw_mdct
#endif
  if (elData.tnsActive)
  {
    if (elData.commonWindow)
    {
      m_auBitStream.write (elData.commonTnsData ? 1 : 0, 1);
      bitCount++;
    }
    m_auBitStream.write (elData.tnsOnLeftRight ? 1 : 0, 1);
    bitCount++;
    if (elData.commonTnsData)
    {
      bitCount += writeChannelWiseTnsData (tnsData0, icsInfo0.windowSequence == EIGHT_SHORT);
    }
    else  // tns_present_both and tns_data_present[1]
    {
      const bool tnsPresentBoth = (tnsData0.numFilters > 0) && (tnsData1.numFilters > 0);

      m_auBitStream.write (tnsPresentBoth ? 1 : 0, 1);
      bitCount++;
      if (!tnsPresentBoth)
      {
        m_auBitStream.write (tnsData1.numFilters > 0 ? 1 : 0, 1);
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
                                             unsigned char* const audioConfig)
{
  unsigned bitCount = 37;

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
    m_auBitStream.write (samplingFrequencyIndex, 4);
  }
  else
  {
    m_auBitStream.write (0xF, 4); // esc
    m_auBitStream.write (toSamplingRate (samplingFrequencyIndex), 24);
    bitCount += 24;
  }
  // for multichannel audio, refer to channel mapping of AotSpecificConfig below
  m_auBitStream.write (chConfigurationIndex > 2 ? 0 : chConfigurationIndex, 4);

// --- AotSpecificConfig(): UsacConfig()
  m_auBitStream.write (samplingFrequencyIndex, 5); // usacSamplingFrequencyIndex
  m_auBitStream.write (shortFrameLength ? 0 : 1, 3);  // coreSbrFrameLengthIndex
  m_auBitStream.write (chConfigurationIndex, 5);    // channelConfigurationIndex
  m_auBitStream.write (numElements - 1, 4);  // numElements in UsacDecoderConfig

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
    }
  } // for el

  m_auBitStream.write (loudnessInfo > 0 ? 1 : 0, 1); // ..ConfigExtensionPresent
  if (loudnessInfo > 0) // ISO 23003-4: loudnessInfo()
  {
    const unsigned methodDefinition = (loudnessInfo >> 14) & 0xF;
    const unsigned methodValueBits  = (methodDefinition == 7 ? 5 : (methodDefinition == 8 ? 2 : 8)); 

    m_auBitStream.write (0, 2); // numConfigExtensions
    m_auBitStream.write (ID_EXT_LOUDNESS_INFO, 4);
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

  memcpy (audioConfig, &m_auBitStream.stream.front (), __min (16, bitCount >> 3));

  return (bitCount >> 3);  // byte count
}

unsigned BitStreamWriter::createAudioFrame (CoreCoderData** const elementData,  EntropyCoder* const entropyCoder,
                                            int32_t** const mdctSignals,        uint8_t** const mdctQuantMag,
                                            const bool usacIndependencyFlag,    const uint8_t numElements,
                                            const uint8_t numSwbShort,          uint8_t* const tempBuffer,
#if !RESTRICT_TO_AAC
                                            const bool* const tw_mdct /*N/A*/,  const bool* const noiseFilling,
#endif
                                            unsigned char* const accessUnit,    const unsigned nSamplesInFrame /*= 1024*/)
{
  unsigned bitCount = 1, ci = 0;

  if ((elementData == nullptr) || (entropyCoder == nullptr) || (tempBuffer == nullptr) ||
      (mdctSignals == nullptr) || (mdctQuantMag == nullptr) || (accessUnit == nullptr) || (nSamplesInFrame > 2048) ||
#if !RESTRICT_TO_AAC
      (noiseFilling == nullptr) || (tw_mdct == nullptr) ||
#endif
      (numElements == 0) || (numElements > USAC_MAX_NUM_ELEMENTS) || (numSwbShort < MIN_NUM_SWB_SHORT) || (numSwbShort > MAX_NUM_SWB_SHORT))
  {
    return 0; // invalid arguments error
  }

  m_auBitStream.reset ();
  m_frameLength = nSamplesInFrame;
  m_numSwbShort = numSwbShort;
  m_uCharBuffer = tempBuffer;
  m_auBitStream.write (usacIndependencyFlag ? 1 : 0, 1);

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
                                          tw_mdct[el], noiseFilling[el],
#endif
                                          usacIndependencyFlag);
        ci++;
        break;
      }
      case ID_USAC_CPE: // UsacChannelPairElement()
      {
        m_auBitStream.write (CORE_MODE_FD, 1); // L
        m_auBitStream.write (CORE_MODE_FD, 1); // R
        bitCount += 2;
        bitCount += writeStereoCoreToolInfo (*elData,
#if !RESTRICT_TO_AAC
                                             tw_mdct[el],
#endif
                                             usacIndependencyFlag);
        bitCount += writeFDChannelStream (*elData, entropyCoder[ci], 0, // L
                                          mdctSignals[ci], mdctQuantMag[ci],
#if !RESTRICT_TO_AAC
                                          tw_mdct[el], noiseFilling[el],
#endif
                                          usacIndependencyFlag);
        ci++;
        bitCount += writeFDChannelStream (*elData, entropyCoder[ci], 1, // R
                                          mdctSignals[ci], mdctQuantMag[ci],
#if !RESTRICT_TO_AAC
                                          tw_mdct[el], noiseFilling[el],
#endif
                                          usacIndependencyFlag);
        ci++;
        break;
      }
      case ID_USAC_LFE: // UsacLfeElement()
      {
        bitCount += writeFDChannelStream (*elData, entropyCoder[ci], 0,
                                          mdctSignals[ci], mdctQuantMag[ci],
#if !RESTRICT_TO_AAC
                                          false, false,
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

  memcpy (accessUnit, &m_auBitStream.stream.front (), __min (768 * ci, bitCount >> 3));

  return (bitCount >> 3);  // byte count
}
