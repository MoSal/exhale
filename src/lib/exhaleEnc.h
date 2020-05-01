/* exhaleEnc.h - header file for class providing Extended HE-AAC encoding capability
 * written by C. R. Helmrich, last modified in 2020 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under a Modified BSD-Style License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2020 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#ifndef _EXHALE_ENC_H_
#define _EXHALE_ENC_H_

#include "exhaleDecl.h"
#include "exhaleLibPch.h"
#include "bitAllocation.h"
#include "bitStreamWriter.h"
#include "entropyCoding.h"
#include "lappedTransform.h"
#include "linearPrediction.h"
#include "quantization.h"
#include "specAnalysis.h"
#include "specGapFilling.h"
#include "stereoProcessing.h"
#include "tempAnalysis.h"

// constant and experimental macro
#define WIN_SCALE double (1 << 23)

// channelConfigurationIndex setup
typedef enum USAC_CCI : char
{
  CCI_UNDEF = -1,
  CCI_CONF  = 0,  // channel-to-speaker mapping defined in UsacChannelConfig() (not to be used here!)
  CCI_1_CH  = 1,  // 1.0: front-center
  CCI_2_CH  = 2,  // 2.0: front-left, front-right
  CCI_3_CH  = 3,  // 3.0: front-center, front-left, front-right
  CCI_4_CH  = 4,  // 4.0: front-center, front-left, front-right, back-center
  CCI_5_CH  = 5,  // 5.0: front-center, front-left, front-right, back-left, back-right
  CCI_6_CH  = 6,  // 5.1: front-center, front-left, front-right, back-left, back-right, LFE
  CCI_8_CH  = 7,  // 7.1: front-center, front-left, front-right, side-left, side-right, back-left, back-right, LFE
  CCI_2_CHM = 8,  // 2.0, dual-mono: channel1, channel2
  CCI_3_CHR = 9,  // 3.0, R-rotated: front-left, front-right, back-center
  CCI_4_CHR = 10, // 4.0, R-rotated: front-left, front-right, back-left, back-right
  CCI_7_CH  = 11, // 6.1: front-center, front-left, front-right, back-left, back-right, back-center, LFE
  CCI_8_CHS = 12  // 7.1, surround: front-center, front-L, front-R, surround-L, surround-R, back-L, back-R, LFE
} USAC_CCI;

// coreCoderFrameLength definition
typedef enum USAC_CCFL : short
{
  CCFL_UNDEF = -1,
#if !RESTRICT_TO_AAC
  CCFL_768   = 768, // LD
#endif
  CCFL_1024  = 1024 // LC
} USAC_CCFL;

// overall xHE-AAC encoding class
class ExhaleEncoder : public ExhaleEncAPI
{
private:

  // member variables
  uint16_t        m_bandwidCurr[USAC_MAX_NUM_CHANNELS];
  uint16_t        m_bandwidPrev[USAC_MAX_NUM_CHANNELS];
  BitAllocator    m_bitAllocator; // for scale factor init
  uint8_t         m_bitRateMode;
  USAC_CCI        m_channelConf;
  CoreCoderData*  m_elementData[USAC_MAX_NUM_ELEMENTS];
  EntropyCoder    m_entropyCoder[USAC_MAX_NUM_CHANNELS];
  uint32_t        m_frameCount;
  USAC_CCFL       m_frameLength;
  char            m_frequencyIdx;
  bool            m_indepFlag; // usacIndependencyFlag bit
  uint32_t        m_indepPeriod;
  LinearPredictor m_linPredictor; // for pre-roll est, TNS
  uint8_t*        m_mdctQuantMag[USAC_MAX_NUM_CHANNELS];
  int32_t*        m_mdctSignals[USAC_MAX_NUM_CHANNELS];
  int32_t*        m_mdstSignals[USAC_MAX_NUM_CHANNELS];
#if !RESTRICT_TO_AAC
  bool            m_noiseFilling[USAC_MAX_NUM_ELEMENTS];
  bool            m_nonMpegExt;
#endif
  uint8_t         m_numElements;
  uint8_t         m_numSwbLong;
  uint8_t         m_numSwbShort;
  unsigned char*  m_outAuData;
  BitStreamWriter m_outStream; // for access unit creation
  int32_t*        m_pcm24Data;
  uint8_t         m_perCorrHCurr[USAC_MAX_NUM_ELEMENTS];
  uint8_t         m_perCorrLCurr[USAC_MAX_NUM_ELEMENTS];
  SfbGroupData*   m_scaleFacData[USAC_MAX_NUM_CHANNELS];
  SfbQuantizer    m_sfbQuantizer; // powerlaw quantization
  SpecAnalyzer    m_specAnalyzer; // for spectral analysis
  uint32_t        m_specAnaCurr[USAC_MAX_NUM_CHANNELS];
  uint8_t         m_specFlatPrev[USAC_MAX_NUM_CHANNELS];
#if !RESTRICT_TO_AAC
  SpecGapFiller   m_specGapFiller;// for noise/gap filling
#endif
  StereoProcessor m_stereoCoder;  // for M/S stereo coding
  uint8_t         m_swbTableIdx;
  TempAnalyzer    m_tempAnalyzer; // for temporal analysis
  uint32_t        m_tempAnaCurr[USAC_MAX_NUM_CHANNELS];
  uint32_t        m_tempAnaNext[USAC_MAX_NUM_CHANNELS];
  int32_t*        m_tempIntBuf;  // temporary int32 buffer
  int32_t*        m_timeSignals[USAC_MAX_NUM_CHANNELS];
#if !RESTRICT_TO_AAC
  bool            m_timeWarping[USAC_MAX_NUM_ELEMENTS];
#endif
  int32_t*        m_timeWindowL[2];  // long window halves
  int32_t*        m_timeWindowS[2]; // short window halves
  int16_t         m_tranLocCurr[USAC_MAX_NUM_CHANNELS];
  int16_t         m_tranLocNext[USAC_MAX_NUM_CHANNELS];
  LappedTransform m_transform; // time-frequency transform

  // helper functions
  unsigned applyTnsToWinGroup (TnsData& tnsData, SfbGroupData& grpData, const bool eightShorts, const uint8_t maxSfb,
                               const unsigned channelIndex);
  unsigned eightShortGrouping (SfbGroupData& grpData, uint16_t* const grpOffsets, int32_t* const mdctSignal,
                               int32_t* const mdstSignal = nullptr);
  unsigned getOptParCorCoeffs (const int32_t* const mdctSignal, const SfbGroupData& grpData, const uint8_t maxSfb,
                               const unsigned channelIndex, TnsData& tnsData, const uint8_t firstGroupIndexToTest = 0);
  unsigned psychBitAllocation ();
  unsigned quantizationCoding ();
  unsigned spectralProcessing ();
  unsigned temporalProcessing ();

public:

  // constructor
  ExhaleEncoder (int32_t* const inputPcmData,       unsigned char* const outputAuData,
                 const unsigned sampleRate = 44100, const unsigned numChannels = 2,
                 const unsigned frameLength = 1024, const unsigned indepPeriod = 45,
                 const unsigned varBitRateMode = 3
#if !RESTRICT_TO_AAC
               , const bool useNoiseFilling = true, const bool useEcodisExt = false
#endif
                );
  // destructor
  virtual ~ExhaleEncoder ();
  // public functions
  unsigned encodeLookahead ();
  unsigned encodeFrame ();
  unsigned initEncoder (unsigned char* const audioConfigBuffer, uint32_t* const audioConfigBytes = nullptr);

}; // ExhaleEncoder

#ifdef EXHALE_DYN_LINK
// C constructor
extern "C" EXHALE_DECL ExhaleEncAPI* exhaleCreate (int32_t* const inputPcmData,       unsigned char* const outputAuData,
                                                   const unsigned sampleRate = 44100, const unsigned numChannels = 2,
                                                   const unsigned frameLength = 1024, const unsigned indepPeriod = 45,
                                                   const unsigned varBitRateMode = 3, const bool useNoiseFilling = true,
                                                   const bool useEcodisExt = false)
{
  return new ExhaleEncoder (inputPcmData, outputAuData, sampleRate, numChannels, frameLength, indepPeriod, varBitRateMode
#if !RESTRICT_TO_AAC
                          , useNoiseFilling, useEcodisExt
#endif
                            );
}

// C destructor
extern "C" EXHALE_DECL unsigned exhaleDelete (ExhaleEncAPI* exhaleEnc)
{
  if (exhaleEnc != NULL) { delete exhaleEnc;  return 0; }

  return USHRT_MAX; // error
}

// C initializer
extern "C" EXHALE_DECL unsigned exhaleInitEncoder (ExhaleEncAPI* exhaleEnc, unsigned char* const audioConfigBuffer,
                                                   uint32_t* const audioConfigBytes = nullptr)
{
  if (exhaleEnc != NULL) return exhaleEnc->initEncoder (audioConfigBuffer, audioConfigBytes);

  return USHRT_MAX; // error
}

// C lookahead encoder
extern "C" EXHALE_DECL unsigned exhaleEncodeLookahead (ExhaleEncAPI* exhaleEnc)
{
  if (exhaleEnc != NULL) return exhaleEnc->encodeLookahead ();

  return USHRT_MAX; // error
}

// C frame encoder
extern "C" EXHALE_DECL unsigned exhaleEncodeFrame (ExhaleEncAPI* exhaleEnc)
{
  if (exhaleEnc != NULL) return exhaleEnc->encodeFrame ();

  return USHRT_MAX; // error
}
#endif // EXHALE_DYN_LINK

#endif // _EXHALE_ENC_H_
