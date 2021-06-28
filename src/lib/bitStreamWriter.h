/* bitStreamWriter.h - header file for class with basic bit-stream writing capability
 * written by C. R. Helmrich, last modified in 2021 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under the exhale Copyright License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2021 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#ifndef _BIT_STREAM_WRITER_H_
#define _BIT_STREAM_WRITER_H_

#include "exhaleLibPch.h"
#include "entropyCoding.h"

// constants, experimental macros
#define CORE_MODE_FD            0
#define SFB_PER_PRED_BAND       2

// output bit-stream writer class
class BitStreamWriter
{
private:

  // member variables
  OutputStream m_auBitStream; // access unit bit-stream to write
  uint64_t     m_auByteCount;
  uint8_t      m_numSwbShort; // max. SFB count in short windows
  uint8_t*     m_uCharBuffer; // temporary buffer for ungrouping
#ifndef NO_PREROLL_DATA
  uint8_t      m_usacConfig[20]; // buffer for UsacConfig in IPF
  uint16_t     m_usacConfigLen;
#endif
#if !RESTRICT_TO_AAC
  uint8_t      m_usacIpfState[6];
#endif
  // helper functions
  void     writeByteAlignment (); // write 0s for byte alignment
  unsigned writeChannelWiseIcsInfo (const IcsInfo& icsInfo);
  unsigned writeChannelWiseSbrData (const int32_t* const sbrDataCh0, const int32_t* const sbrDataCh1,
                                    const bool indepFlag = false);
  unsigned writeChannelWiseTnsData (const TnsData& tnsData, const bool eightShorts);
  unsigned writeFDChannelStream    (const CoreCoderData& elData, EntropyCoder& entrCoder, const unsigned ch,
                                    const int32_t* const mdctSignal, const uint8_t* const mdctQuantMag,
#if !RESTRICT_TO_AAC
                                    const bool timeWarping, const bool noiseFilling, uint8_t* ipfAuState,
#endif
                                    const bool indepFlag = false);
  unsigned writeStereoCoreToolInfo (const CoreCoderData& elData, EntropyCoder& entrCoder,
#if !RESTRICT_TO_AAC
                                    const bool timeWarping, bool* const commonTnsFlag,
#endif
                                    const bool indepFlag = false);

public:

  // constructor
  BitStreamWriter () { m_auBitStream.reset (); m_auByteCount = m_numSwbShort = 0; m_uCharBuffer = nullptr;
#ifndef NO_PREROLL_DATA
                       memset (m_usacConfig, 0, 20); m_usacConfigLen = 0; memset (m_usacIpfState, 0, 4);
#endif
    }
  // destructor
  ~BitStreamWriter() { m_auBitStream.reset (); }
  // public functions
  unsigned createAudioConfig (const char samplingFrequencyIndex,  const bool shortFrameLength,
                              const uint8_t chConfigurationIndex, const uint8_t numElements,
                              const ELEM_TYPE* const elementType, const uint32_t loudnessInfo,
#if !RESTRICT_TO_AAC
                              const bool* const tw_mdct /*N/A*/,  const bool* const noiseFilling,
#endif
                              const uint8_t sbrRatioShiftValue,   unsigned char* const audioConfig);
  unsigned createAudioFrame  (CoreCoderData** const elementData,  EntropyCoder* const entropyCoder,
                              int32_t** const mdctSignals,        uint8_t** const mdctQuantMag,
                              const bool usacIndependencyFlag,    const uint8_t numElements,
                              const uint8_t numSwbShort,          uint8_t* const tempBuffer,
#if !RESTRICT_TO_AAC
                              const bool* const tw_mdct /*N/A*/,  const bool* const noiseFilling,
                              const uint32_t frameCount,          const uint32_t indepPeriod,  uint32_t* rate,
#endif
                              const uint8_t sbrRatioShiftValue,   int32_t** const sbrInfoAndData,
                              unsigned char* const accessUnit,    const unsigned nSamplesInFrame);
}; // BitStreamWriter

#endif // _BIT_STREAM_WRITER_H_
