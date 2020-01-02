/* basicWavReader.cpp - source file for class with basic WAVE file reading capability
 * written by C. R. Helmrich, last modified in 2019 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under a Modified BSD-Style License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2020 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#include "exhaleAppPch.h"
#include "basicWavReader.h"

// static helper functions
static unsigned reverseFourBytes (const uint8_t* b)
{
  return ((unsigned) b[3] << 24) | ((unsigned) b[2] << 16) | ((unsigned) b[1] << 8) | (unsigned) b[0];
}

static int64_t fourBytesToLength (const uint8_t* b, const int64_t lengthLimit)
{
  int64_t chunkLength = (int64_t) reverseFourBytes (b);

  chunkLength += chunkLength & 1;  // make sure it is even

  return __min (lengthLimit, chunkLength); // for security
}

// private reader functions
bool BasicWavReader::readRiffHeader ()
{
  uint8_t b[FILE_HEADER_SIZE] = {0};  // temp. byte buffer

  if ((m_bytesRead = _READ (m_fileHandle, b, FILE_HEADER_SIZE)) != FILE_HEADER_SIZE) return false; // error
  m_bytesRemaining -= m_bytesRead;
  m_chunkLength = fourBytesToLength (&b[4], m_bytesRemaining) - 4; // minus 4 bytes for WAVE tag

  return (b[0] == 'R' && b[1] == 'I' && b[2] == 'F' && b[3] == 'F' &&
          b[8] == 'W' && b[9] == 'A' && b[10]== 'V' && b[11]== 'E' &&
          m_bytesRemaining > 32);  // true: RIFF supported
}

bool BasicWavReader::readFormatChunk ()
{
  uint8_t b[CHUNK_FORMAT_MAX] = {0};  // temp. byte buffer

  if (!seekToChunkTag (b, 0x20746D66 /*fmt */) || (m_chunkLength < CHUNK_FORMAT_SIZE) || (m_chunkLength > CHUNK_FORMAT_MAX))
  {
    return false; // fmt_ chunk invalid or read incomplete
  }
  if ((m_bytesRead = _READ (m_fileHandle, b, (unsigned) m_chunkLength)) != m_chunkLength) return false; // error
  m_bytesRemaining -= m_bytesRead;

  m_waveDataType  = WAV_TYPE (b[0]-1); // 1: PCM, 3: float
  m_waveChannels  = b[2];  // only 1, 2, ..., 63 supported
  m_waveFrameRate = reverseFourBytes (&b[4]);  // frames/s
  m_waveBitRate   = reverseFourBytes (&b[8]) * 8; // bit/s
  m_waveFrameSize = b[12];  // bytes/s divided by frames/s
  m_waveBitDepth  = b[14]; // only 8, 16, 24, 32 supported

  return ((m_waveDataType == WAV_PCM || (m_waveDataType == WAV_FLOAT && (m_waveBitDepth & 15) == 0)) &&
          (m_waveChannels > 0 && m_waveChannels <= 63) && isSamplingRateSupported (m_waveFrameRate) &&
          (m_waveBitRate == 8 * m_waveFrameRate * m_waveFrameSize) && (b[ 1] == 0) && (b[ 3] == 0) &&
          (m_waveFrameSize * 8 == m_waveBitDepth * m_waveChannels) && (b[13] == 0) && (b[15] == 0) &&
          (m_waveBitDepth >= 8 && m_waveBitDepth <= 32 && (m_waveBitDepth & 7) == 0) &&
          m_bytesRemaining > 8); // true: format supported
}

bool BasicWavReader::readDataHeader ()
{
  uint8_t b[CHUNK_HEADER_SIZE] = {0}; // temp. byte buffer

  if (!seekToChunkTag (b, 0x61746164 /*data*/))
  {
    return false; // data chunk invalid or read incomplete
  }
  return (m_chunkLength > 0); // true: WAVE data available
}

// private helper function
bool BasicWavReader::seekToChunkTag (uint8_t* const buf, const uint32_t tagID)
{
  if ((m_bytesRead = _READ (m_fileHandle, buf, CHUNK_HEADER_SIZE)) != CHUNK_HEADER_SIZE) return false; // error
  m_bytesRemaining -= m_bytesRead;
  m_chunkLength = fourBytesToLength (&buf[4], m_bytesRemaining);

  while ((*((uint32_t* const) buf) != tagID) &&
         (m_bytesRemaining > 0)) // seek until tagID found
  {
    if ((m_readOffset = _SEEK (m_fileHandle, m_chunkLength, 1 /*SEEK_CUR*/)) == -1)
    {
      // for stdin compatibility, don't abort, try reading
      for (int64_t i = m_chunkLength >> 1; i > 0; i--)
      {
        _READ (m_fileHandle, buf, 2); // as length is even
      }
    }
    m_bytesRemaining -= m_chunkLength;
    if (m_bytesRemaining <= 0)
    {
      return false; // an error which should never happen!
    }
    if ((m_bytesRead = _READ (m_fileHandle, buf, CHUNK_HEADER_SIZE)) != CHUNK_HEADER_SIZE) return false; // error
    m_bytesRemaining -= m_bytesRead;
    m_chunkLength = fourBytesToLength (&buf[4], m_bytesRemaining);
  }
  return (m_bytesRemaining > 0);
}

// static reading functions
unsigned BasicWavReader::readDataFloat16 (const int fileHandle, int32_t* frameBuf, const unsigned frameCount,
                                          const unsigned chanCount, void* tempBuf)
{
#if BWR_BUFFERED_READ
  const int16_t* fBuf = (const int16_t*) tempBuf; // words
  const int bytesRead = _READ (fileHandle, tempBuf, frameCount * chanCount * 2);
  unsigned framesRead = __max (0, bytesRead / (chanCount * 2));

  for (unsigned i = framesRead * chanCount; i > 0; i--)
  {
    const int16_t i16 = *(fBuf++);
    const int32_t e = ((i16 & 0x7C00) >> 10) - 18; // exp.
    // an exponent e <= -12 will lead to zero-quantization
    *frameBuf = int32_t (e < 0 ? (1024 + (i16 & 0x03FF) + (1 << (-1 - e)) /*rounding offset*/) >> -e
                               : (e > 12 ? MAX_VALUE_AUDIO24 /*inf*/ : (1024 + (i16 & 0x03FF)) << e));
    if ((i16 & 0x8000) != 0) *frameBuf *= -1; // neg. sign
    frameBuf++;
  }
  if (framesRead < frameCount) // zero out missing samples
  {
    memset (frameBuf, 0, (frameCount - framesRead) * chanCount * sizeof (int32_t));
  }
  return framesRead;
#else
  unsigned bytesRead = 0;

  for (unsigned i = frameCount * chanCount; i > 0; i--)
  {
    int16_t i16 = 0;

    bytesRead += _READ (fileHandle, &i16, 2); // two bytes

    const int32_t e = ((i16 & 0x7C00) >> 10) - 18; // exp.
    // an exponent e <= -12 will lead to zero-quantization
    *frameBuf = int32_t (e < 0 ? (1024 + (i16 & 0x03FF) + (1 << (-1 - e)) /*rounding offset*/) >> -e
                               : (e > 12 ? MAX_VALUE_AUDIO24 /*inf*/ : (1024 + (i16 & 0x03FF)) << e));
    if ((i16 & 0x8000) != 0) *frameBuf *= -1; // neg. sign
    frameBuf++;
  }
  return bytesRead / (chanCount * 2);
#endif
}

unsigned BasicWavReader::readDataFloat32 (const int fileHandle, int32_t* frameBuf, const unsigned frameCount,
                                          const unsigned chanCount, void* tempBuf)
{
#if BWR_BUFFERED_READ
  const float*   fBuf = (const float*) tempBuf; // 4 bytes
  const int bytesRead = _READ (fileHandle, tempBuf, frameCount * chanCount * 4);
  unsigned framesRead = __max (0, bytesRead / (chanCount * 4));

  for (unsigned i = framesRead * chanCount; i > 0; i--)
  {
    const float   f32 = *fBuf * float (1 << 23); // * 2^23
    fBuf++;
    *frameBuf = int32_t (f32 + (f32 < 0.0 ? -0.5 : 0.5)); // rounding
    if (*frameBuf < MIN_VALUE_AUDIO24) *frameBuf = MIN_VALUE_AUDIO24;
    else
    if (*frameBuf > MAX_VALUE_AUDIO24) *frameBuf = MAX_VALUE_AUDIO24;
    frameBuf++;
  }
  if (framesRead < frameCount) // zero out missing samples
  {
    memset (frameBuf, 0, (frameCount - framesRead) * chanCount * sizeof (int32_t));
  }
  return framesRead;
#else
  unsigned bytesRead = 0;

  for (unsigned i = frameCount * chanCount; i > 0; i--)
  {
    float f32 = 0.0; // IEEE-754 normalized floating point

    bytesRead += _READ (fileHandle, &f32, 4);   // 4 bytes
    *frameBuf = int32_t (f32 * (1 << 23) + (f32 < 0.0 ? -0.5 : 0.5)); // * 2^23 with rounding
    if (*frameBuf < MIN_VALUE_AUDIO24) *frameBuf = MIN_VALUE_AUDIO24;
    else
    if (*frameBuf > MAX_VALUE_AUDIO24) *frameBuf = MAX_VALUE_AUDIO24;
    frameBuf++;
  }
  return bytesRead / (chanCount * 4);
#endif
}

unsigned BasicWavReader::readDataLnPcm08 (const int fileHandle, int32_t* frameBuf, const unsigned frameCount,
                                          const unsigned chanCount, void* tempBuf)
{
#if BWR_BUFFERED_READ
  const uint8_t* iBuf = (uint8_t*) tempBuf;
  const int bytesRead = _READ (fileHandle, tempBuf, frameCount * chanCount);
  unsigned framesRead = __max (0, bytesRead / chanCount);

  for (unsigned i = framesRead * chanCount; i > 0; i--)
  {
    *(frameBuf++) = ((int32_t) *(iBuf++) - 128) << 16; // * 2^16
  }
  if (framesRead < frameCount) // zero out missing samples
  {
    memset (frameBuf, 0, (frameCount - framesRead) * chanCount * sizeof (int32_t));
  }
  return framesRead;
#else
  unsigned bytesRead = 0;

  for (unsigned i = frameCount * chanCount; i > 0; i--)
  {
    uint8_t ui8 = 128;

    bytesRead += _READ (fileHandle, &ui8, 1);  // one byte
    *(frameBuf++) = ((int32_t) ui8 - 128) << 16; // * 2^16
  }
  return bytesRead / chanCount;
#endif
}

unsigned BasicWavReader::readDataLnPcm16 (const int fileHandle, int32_t* frameBuf, const unsigned frameCount,
                                          const unsigned chanCount, void* tempBuf)
{
#if BWR_BUFFERED_READ
  const int16_t* iBuf = (const int16_t*) tempBuf; // words
  const int bytesRead = _READ (fileHandle, tempBuf, frameCount * chanCount * 2);
  unsigned framesRead = __max (0, bytesRead / (chanCount * 2));

  for (unsigned i = framesRead * chanCount; i > 0; i--)
  {
    *(frameBuf++) = (int32_t) *(iBuf++) << 8; // * 2^8
  }
  if (framesRead < frameCount) // zero out missing samples
  {
    memset (frameBuf, 0, (frameCount - framesRead) * chanCount * sizeof (int32_t));
  }
  return framesRead;
#else
  unsigned bytesRead = 0;

  for (unsigned i = frameCount * chanCount; i > 0; i--)
  {
    int16_t i16 = 0;

    bytesRead += _READ (fileHandle, &i16, 2); // two bytes
    *(frameBuf++) = (int32_t) i16 << 8; // * 2^8
  }
  return bytesRead / (chanCount * 2);
#endif
}

unsigned BasicWavReader::readDataLnPcm24 (const int fileHandle, int32_t* frameBuf, const unsigned frameCount,
                                          const unsigned chanCount, void* tempBuf)
{
#if BWR_BUFFERED_READ
  const uint8_t* iBuf = (uint8_t*) tempBuf;
  const int bytesRead = _READ (fileHandle, tempBuf, frameCount * chanCount * 3);
  unsigned framesRead = __max (0, bytesRead / (chanCount * 3));

  for (unsigned i = framesRead * chanCount; i > 0; i--)
  {
    const int32_t i24 = (int32_t) iBuf[0] | ((int32_t) iBuf[1] << 8) | ((int32_t) iBuf[2] << 16);
    iBuf += 3;
    *(frameBuf++) = (i24 > MAX_VALUE_AUDIO24 ? i24 + 2 * MIN_VALUE_AUDIO24 : i24);
  }
  if (framesRead < frameCount) // zero out missing samples
  {
    memset (frameBuf, 0, (frameCount - framesRead) * chanCount * sizeof (int32_t));
  }
  return framesRead;
#else
  unsigned bytesRead = 0;

  for (unsigned i = frameCount * chanCount; i > 0; i--)
  {
    int32_t i24 = 0;

    bytesRead += _READ (fileHandle, &i24, 3);   // 3 bytes
    *(frameBuf++) = (i24 > MAX_VALUE_AUDIO24 ? i24 + 2 * MIN_VALUE_AUDIO24 : i24);
  }
  return bytesRead / (chanCount * 3);
#endif
}

unsigned BasicWavReader::readDataLnPcm32 (const int fileHandle, int32_t* frameBuf, const unsigned frameCount,
                                          const unsigned chanCount, void* tempBuf)
{
#if BWR_BUFFERED_READ
  const int32_t* iBuf = (const int32_t*) tempBuf; // dword
  const int bytesRead = _READ (fileHandle, tempBuf, frameCount * chanCount * 4);
  unsigned framesRead = __max (0, bytesRead / (chanCount * 4));

  for (unsigned i = framesRead * chanCount; i > 0; i--)
  {
    const int32_t i24 = ((*iBuf >> 1) + (1 << 6)) >> 7; // * 2^-8 with rounding, overflow-safe
    iBuf++;
    *(frameBuf++) = __min (MAX_VALUE_AUDIO24, i24);
  }
  if (framesRead < frameCount) // zero out missing samples
  {
    memset (frameBuf, 0, (frameCount - framesRead) * chanCount * sizeof (int32_t));
  }
  return framesRead;
#else
  unsigned bytesRead = 0;

  for (unsigned i = frameCount * chanCount; i > 0; i--)
  {
    int32_t i24 = 0;
    bytesRead += _READ (fileHandle, &i24, 4);   // 4 bytes
    i24 = ((i24 >> 1) + (1 << 6)) >> 7; // * 2^-8 with rounding, overflow-safe
    *(frameBuf++) = __min (MAX_VALUE_AUDIO24, i24);
  }
  return bytesRead / (chanCount * 4);
#endif
}

// public functions
unsigned BasicWavReader::open (const int wavFileHandle, const uint16_t maxFrameRead, const int64_t fileLength /*= LLONG_MAX*/)
{
  m_bytesRemaining = fileLength;
  m_fileHandle     = wavFileHandle;

  if ((m_fileHandle == -1) || (fileLength <= 44))
  {
    return 1; // file handle invalid or file too small
  }
#if defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64)
  if ((fileLength < LLONG_MAX) && (m_readOffset = _telli64 (m_fileHandle)) != 0)
#else // Linux, MacOS, Unix
  if ((fileLength < LLONG_MAX) && (m_readOffset = lseek (m_fileHandle, 0, 1 /*SEEK_CUR*/)) != 0)
#endif
  {
    m_readOffset = _SEEK (m_fileHandle, 0, 0 /*SEEK_SET*/);
  }
  if ((m_readOffset != 0) || !readRiffHeader ())
  {
    return 2; // file type invalid or file seek failed
  }
  if (!readFormatChunk ())
  {
    return 3; // audio format invalid or not supported
  }
  if (!readDataHeader ())
  {
    return 4; // WAVE data part invalid or unsupported
  }
  if ((m_byteBuffer = (char*) malloc (m_waveFrameSize * maxFrameRead)) == nullptr)
  {
    return 5; // read-in byte buffer allocation failed
  }
  m_frameLimit = maxFrameRead;

  // ready to read audio data: initialize byte counter
  if (m_bytesRemaining > m_chunkLength)
  {
    m_bytesRemaining = m_chunkLength;
  }
  m_chunkLength = 0;

  if (m_waveDataType == WAV_PCM) // & function pointer
  {
    switch (m_waveBitDepth)
    {
      case 8:
        m_readDataFunc = readDataLnPcm08; break;
      case 16:
        m_readDataFunc = readDataLnPcm16; break;
      case 24:
        m_readDataFunc = readDataLnPcm24; break;
      default:
        m_readDataFunc = readDataLnPcm32; break;
    }
  }
  else  m_readDataFunc = (m_waveBitDepth == 16 ? readDataFloat16 : readDataFloat32);

  return (m_readDataFunc == nullptr ? 6 : 0); // 0: OK
}

unsigned BasicWavReader::read (int32_t* const frameBuf, const uint16_t frameCount)
{
  unsigned framesRead;

  if ((frameBuf == nullptr) || (m_fileHandle == -1) || (__min (m_frameLimit, frameCount) == 0) || (m_byteBuffer == nullptr))
  {
    return 0; // invalid args or class not initialized
  }
  framesRead  = m_readDataFunc (m_fileHandle, frameBuf, __min (m_frameLimit, frameCount), m_waveChannels, m_byteBuffer);
  m_bytesRead = m_waveFrameSize * framesRead;
  m_bytesRemaining -= m_bytesRead;
  m_chunkLength    += m_bytesRead;

  return framesRead;
}

void BasicWavReader::reset ()
{
  m_byteBuffer     = nullptr;
  m_bytesRead      = 0;
  m_bytesRemaining = 0;
  m_chunkLength    = 0;
  m_frameLimit     = 0;
  m_readDataFunc   = nullptr;
  m_readOffset     = 0;
  m_waveBitDepth   = 0;
  m_waveChannels   = 0;
  m_waveFrameRate  = 0;

  if (m_fileHandle != -1) _SEEK (m_fileHandle, 0, 0 /*SEEK_SET*/);
}
