/* bitStreamWriter.h - header file for class with basic bit-stream writing capability
 * written by C. R. Helmrich, last modified in 2020 - see License.htm for legal notices
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
#define ID_EXT_LOUDNESS_INFO    2
#define ID_EXT_ELE_FILL         0
#define SFB_PER_PRED_BAND       2

// output bit-stream writer class
class BitStreamWriter
{
private:

  // member variables
  OutputStream m_auBitStream; // access unit bit-stream to write
  uint32_t     m_frameLength; // number of samples in full frame
  uint8_t      m_numSwbShort; // max. SFB count in short windows
  uint8_t*     m_uCharBuffer; // temporary buffer for ungrouping

  // helper functions
  void     writeByteAlignment (); // write 0s for byte alignment
  unsigned writeChannelWiseIcsInfo (const IcsInfo& icsInfo); // ics_info()
  unsigned writeChannelWiseSbrData (const int32_t* const sbrDataCh0, const int32_t* const sbrDataCh1,
                                    const bool indepFlag = false);
  unsigned writeChannelWiseTnsData (const TnsData& tnsData, const bool eightShorts);
  unsigned writeFDChannelStream    (const CoreCoderData& elData, EntropyCoder& entrCoder, const unsigned ch,
                                    const int32_t* const mdctSignal, const uint8_t* const mdctQuantMag,
#if !RESTRICT_TO_AAC
                                    const bool timeWarping, const bool noiseFilling,
#endif
                                    const bool indepFlag = false);
  unsigned writeStereoCoreToolInfo (const CoreCoderData& elData, EntropyCoder& entrCoder,
#if !RESTRICT_TO_AAC
                                    const bool timeWarping, bool* const commonTnsFlag,
#endif
                                    const bool indepFlag = false);

public:

  // constructor
  BitStreamWriter () { m_auBitStream.reset (); m_frameLength = 0; m_numSwbShort = 0; m_uCharBuffer = nullptr; }
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
                              const bool* const tw_mdct /*N/A*/,  const bool* const noiseFilling, const bool ipf,
#endif
                              const uint8_t sbrRatioShiftValue,   int32_t** const sbrInfoAndData,
                              unsigned char* const accessUnit,    const unsigned nSamplesInFrame = 1024);
}; // BitStreamWriter

#endif // _BIT_STREAM_WRITER_H_
