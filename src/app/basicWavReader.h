/* basicWavReader.h - header file for class with basic WAVE file reading capability
 * written by C. R. Helmrich, last modified in 2020 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under a Modified BSD-Style License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2020 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#ifndef _BASIC_WAV_READER_H_
#define _BASIC_WAV_READER_H_

#include "exhaleAppPch.h"

// constant data sizes & limits
#define BWR_BUFFERED_READ        1 // faster reader
#define BWR_READ_FRACT           5 // 2^-READ_FRACT
#define CHUNK_FORMAT_MAX        20
#define CHUNK_FORMAT_SIZE       16
#define CHUNK_HEADER_SIZE        8
#define FILE_HEADER_SIZE        12
#define MAX_VALUE_AUDIO24  8388607 // (1 << 23) - 1
#define MIN_VALUE_AUDIO24 -8388608 // (1 << 23) *-1

// WAVE data format definitions
typedef enum WAV_TYPE
{
  WAV_PCM = 0, // linear PCM
  WAV_ADPCM,   // ADPCM
  WAV_FLOAT    // IEEE float
} WAV_TYPE;

// data reader function pointer
typedef unsigned (*ReadFunc) (const int, int32_t*, const unsigned, const unsigned, void*);

// basic WAV audio reader class
class BasicWavReader
{
private:

  // member variables
  char*    m_byteBuffer;
  unsigned m_bytesRead;
  int64_t  m_bytesRemaining;
  int64_t  m_chunkLength;
  int      m_fileHandle;
  unsigned m_frameLimit;
  ReadFunc m_readDataFunc;
  int64_t  m_readOffset;
  unsigned m_waveBitDepth;
  unsigned m_waveBitRate;
  unsigned m_waveChannels;
  WAV_TYPE m_waveDataType;
  unsigned m_waveFrameRate;
  unsigned m_waveFrameSize;
  // private reader functions
  bool     readRiffHeader ();
  bool     readFormatChunk();
  bool     readDataHeader ();
  // private helper function
  bool     seekToChunkTag (uint8_t* const buf, const uint32_t tagName);
  // static reading functions
  static unsigned readDataFloat16 (const int fileHandle, int32_t* frameBuf, const unsigned frameCount,
                                   const unsigned chanCount, void* tempBuf);
  static unsigned readDataFloat32 (const int fileHandle, int32_t* frameBuf, const unsigned frameCount,
                                   const unsigned chanCount, void* tempBuf);
  static unsigned readDataLnPcm08 (const int fileHandle, int32_t* frameBuf, const unsigned frameCount,
                                   const unsigned chanCount, void* tempBuf);
  static unsigned readDataLnPcm16 (const int fileHandle, int32_t* frameBuf, const unsigned frameCount,
                                   const unsigned chanCount, void* tempBuf);
  static unsigned readDataLnPcm24 (const int fileHandle, int32_t* frameBuf, const unsigned frameCount,
                                   const unsigned chanCount, void* tempBuf);
  static unsigned readDataLnPcm32 (const int fileHandle, int32_t* frameBuf, const unsigned frameCount,
                                   const unsigned chanCount, void* tempBuf);
public:

  // constructor
  BasicWavReader () { m_fileHandle = -1;  reset (); }
  // destructor
  ~BasicWavReader() { if (m_byteBuffer != nullptr) free ((void*) m_byteBuffer); }
  // public functions
  int64_t  getDataBytesLeft () const { return m_bytesRemaining; }
  int64_t  getDataBytesRead () const { return m_chunkLength; }
  unsigned getBitDepth      () const { return m_waveBitDepth; }
  unsigned getNumChannels   () const { return m_waveChannels; }
  unsigned getSampleRate    () const { return m_waveFrameRate; }
  unsigned open  (const int wavFileHandle, const uint16_t maxFrameRead, const int64_t fileLength = LLONG_MAX /*for stdin*/);
  unsigned read  (int32_t* const frameBuf, const uint16_t frameCount);
  void     reset ();
}; // BasicWavReader

#endif // _BASIC_WAV_READER_H_
