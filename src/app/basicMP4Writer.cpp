/* basicMP4Writer.cpp - source file for class with basic MPEG-4 file writing capability
 * written by C. R. Helmrich, last modified in 2019 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under a Modified BSD-Style License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2020 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#include "exhaleAppPch.h"
#include "basicMP4Writer.h"

#if 0 // DEBUG
static const uint8_t muLawHeader[44] = {
  0x52, 0x49, 0x46, 0x46, 0xF0, 0xFF, 0xFF, 0xFF, 0x57, 0x41, 0x56, 0x45, 0x66, 0x6D, 0x74, 0x20,
  0x10, 0x00, 0x00, 0x00, 0x07, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x08, 0x00, 0x64, 0x61, 0x74, 0x61, 0xF0, 0xFF, 0xFF, 0xFF
};
#endif

static const uint8_t staticHeaderTemplate[STAT_HEADER_SIZE] = {
  0x00, 0x00, 0x00, 0x18, 0x66, 0x74, 0x79, 0x70, 0x6D, 0x70, 0x34, 0x32, 0x00, 0x00, 0x00, 0x00, // ftyp
  0x6D, 0x70, 0x34, 0x32, 0x69, 0x73, 0x6F, 0x6D, 0x00, 0x00, MOOV_BSIZE, 0x6D, 0x6F, 0x6F, 0x76, // moov
  0x00, 0x00, 0x00, 0x6C, 0x6D, 0x76, 0x68, 0x64, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // mvhd
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00,
  0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x40, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x02, // end atom 2.1 (mvhd)
  0x00, 0x00, 0x00, 0x18, 0x69, 0x6F, 0x64, 0x73, 0x00, 0x00, 0x00, 0x00, 0x10, 0x80, 0x80, 0x80, // iods
  0x07, 0x00, 0x4F, 0xFF, 0xFF, 0x49, 0xFF, 0xFF, 0x00, 0x00, TRAK_BSIZE, // end atom 2.2 (iods)
  0x74, 0x72, 0x61, 0x6B, 0x00, 0x00, 0x00, 0x5C, 0x74, 0x6B, 0x68, 0x64, 0x00, 0x00, 0x00, 0x07, // tkhd
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x01, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x40, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x24, 0x65, 0x64, 0x74, 0x73, 0x00, 0x00, 0x00, 0x1C, 0x65, 0x6C, 0x73, 0x74, // elst
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x01, 0x00, 0x00, 0x00, 0x00, MDIA_BSIZE, 0x6D, 0x64, 0x69, 0x61, 0x00, 0x00, 0x00, 0x20, // mdhd
  0x6D, 0x64, 0x68, 0x64, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x55, 0xC4, 0x00, 0x00, 0x00, 0x00, 0x00, 0x24,
  0x68, 0x64, 0x6C, 0x72, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x73, 0x6F, 0x75, 0x6E, // hdlr
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x63, 0x72, 0x68, 0x00,
  0x00, 0x00, MINF_BSIZE, 0x6D, 0x69, 0x6E, 0x66, 0x00, 0x00, 0x00, 0x10, 0x73, 0x6D, 0x68, 0x64,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x24, 0x64, 0x69, 0x6E, 0x66, // dinf
  0x00, 0x00, 0x00, 0x1C, 0x64, 0x72, 0x65, 0x66, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01,
  0x00, 0x00, 0x00, 0x0C, 0x75, 0x72, 0x6C, 0x20, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, STBL_BSIZE,
  0x73, 0x74, 0x62, 0x6C, 0x00, 0x00, 0x00, 0x20, 0x73, 0x74, 0x74, 0x73, 0x00, 0x00, 0x00, 0x00, // stts
  0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, STSD_BSIZE, 0x73, 0x74, 0x73, 0x64, 0x00, 0x00, 0x00, 0x00, // stsd
  0x00, 0x00, 0x00, 0x01, 0x00, 0x00, MP4A_BSIZE, 0x6D, 0x70, 0x34, 0x61, 0x00, 0x00, 0x00, 0x00, // mp4a
  0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, ESDS_BSIZE, 0x65, 0x73, 0x64, 0x73, // esds
  0x00, 0x00, 0x00, 0x00, 0x03, 0x80, 0x80, 0x80, 0x25, 0x00, 0x00, 0x00, 0x04, 0x80, 0x80, 0x80, // tag4
  0x17, 0x40, 0x15, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x05, 0x80, // tag5
  0x80, 0x80, 0x05, 0x00, 0x00, 0x00, 0x00, 0x00 // ASC continued in m_dynamicHeader if >5 bytes
};

// static helper functions
static uint32_t toBigEndian (const unsigned ui) // Motorola endianness
{
  return ((ui & UCHAR_MAX) << 24) | (((ui >> 8) & UCHAR_MAX) << 16) | (((ui >> 16) & UCHAR_MAX) << 8) | ((ui >> 24) & UCHAR_MAX);
}

static uint16_t toUShortValue (const uint8_t hiByte, const uint8_t loByte)
{
  return ((uint16_t) hiByte << 8) | (uint16_t) loByte;
}

// public functions
int BasicMP4Writer::addFrameAU (const uint8_t* byteBuf, const uint32_t byteOffset, const uint32_t byteCount)
{
  if ((m_fileHandle == -1) || (m_m4aMdatSize > 0xFFFFFFF0u - byteCount))
  {
    return 1; // invalid file handle or file getting too big
  }

  // add frame byte-size, in Big Endian format, to frame size list (stsz)
  m_dynamicHeader.push_back ((byteCount >> 24) & UCHAR_MAX);
  m_dynamicHeader.push_back ((byteCount >> 16) & UCHAR_MAX);
  m_dynamicHeader.push_back ((byteCount >>  8) & UCHAR_MAX);
  m_dynamicHeader.push_back ( byteCount        & UCHAR_MAX);

  if (((m_frameCount++) % m_rndAccPeriod) == 0) // add RAP to list (stco)
  {
    m_rndAccOffsets.push_back (m_m4aMdatSize);
  }
  m_m4aMdatSize += byteCount;

  return _WRITE (m_fileHandle, byteBuf, byteCount);  // write access unit
}

int BasicMP4Writer::finishFile (const unsigned avgBitrate, const unsigned maxBitrate, const uint32_t audioLength,
                                const uint32_t modifTime /*= 0*/)
{
  const unsigned numFramesFirstPeriod = __min (m_frameCount, m_rndAccPeriod);
  const unsigned numFramesFinalPeriod = (m_frameCount <= m_rndAccPeriod ? 0 : m_frameCount % m_rndAccPeriod);
  const unsigned numSamplesFinalFrame = (audioLength + m_pregapLength) % m_frameLength;
  const uint32_t stszAtomSize = STSX_BSIZE + 4 /*bytes for sampleSize*/ + m_frameCount * 4;
  const uint32_t stscAtomSize = STSX_BSIZE + (numFramesFinalPeriod == 0 ? 12 : 24);
  const uint32_t stcoAtomSize = STSX_BSIZE + (uint32_t) m_rndAccOffsets.size () * 4;
  const uint32_t stblIncrSize = m_ascSizeM5 + stszAtomSize + stscAtomSize + stcoAtomSize;
  const uint32_t moovAtomSize = toBigEndian (toUShortValue (MOOV_BSIZE) + stblIncrSize);
  const uint32_t trakAtomSize = toBigEndian (toUShortValue (TRAK_BSIZE) + stblIncrSize);
  const uint32_t mdiaAtomSize = toBigEndian (toUShortValue (MDIA_BSIZE) + stblIncrSize);
  const uint32_t minfAtomSize = toBigEndian (toUShortValue (MINF_BSIZE) + stblIncrSize);
  const uint32_t stblAtomSize = toBigEndian (toUShortValue (STBL_BSIZE) + stblIncrSize);
  const uint32_t numSamplesBE = toBigEndian (audioLength);
  const uint32_t  timeStampBE = toBigEndian (modifTime);
  const uint32_t  headerBytes = STAT_HEADER_SIZE + (uint32_t) m_dynamicHeader.size () + stscAtomSize + stcoAtomSize;
  uint32_t* const header4Byte = (uint32_t* const) m_staticHeader;
  int bytesWritten = 0;

  if ((m_fileHandle == -1) || (m_m4aMdatSize > 0xFFFFFFF0u - headerBytes))
  {
    return 1; // invalid file handle or file getting too big
  }

  // finish setup of fixed-length part of MPEG-4 file header
  if (modifTime > 0)
  {
    header4Byte[ 48>>2] = timeStampBE; // mvhd
    header4Byte[188>>2] = timeStampBE; // tkhd
    header4Byte[324>>2] = timeStampBE; // mdhd
  }
  header4Byte[ 24>>2] = moovAtomSize;
  header4Byte[ 56>>2] = numSamplesBE;
  header4Byte[164>>2] = trakAtomSize;
  header4Byte[200>>2] = numSamplesBE;
  header4Byte[300>>2] = mdiaAtomSize;
  header4Byte[332>>2] = toBigEndian (audioLength + m_pregapLength);
  header4Byte[376>>2] = minfAtomSize;
  header4Byte[288>>2] = numSamplesBE;  // elst
  header4Byte[436>>2] = stblAtomSize;
  header4Byte[460>>2] = toBigEndian (m_frameCount - 1); // 2 entries used
  header4Byte[472>>2] = toBigEndian (numSamplesFinalFrame == 0 ? m_frameLength : numSamplesFinalFrame);

  m_staticHeader[558] = ((maxBitrate >> 24) & UCHAR_MAX);
  m_staticHeader[559] = ((maxBitrate >> 16) & UCHAR_MAX);
  m_staticHeader[560] = ((maxBitrate >>  8) & UCHAR_MAX);
  m_staticHeader[561] = ( maxBitrate        & UCHAR_MAX);
  m_staticHeader[562] = ((avgBitrate >> 24) & UCHAR_MAX);
  m_staticHeader[563] = ((avgBitrate >> 16) & UCHAR_MAX);
  m_staticHeader[564] = ((avgBitrate >>  8) & UCHAR_MAX);
  m_staticHeader[565] = ( avgBitrate        & UCHAR_MAX);

  // finish dynamically-sized 2nd part of MPEG-4 file header
  m_dynamicHeader.at (m_ascSizeM5 +  6) = ((stszAtomSize >> 24) & UCHAR_MAX);
  m_dynamicHeader.at (m_ascSizeM5 +  7) = ((stszAtomSize >> 16) & UCHAR_MAX);
  m_dynamicHeader.at (m_ascSizeM5 +  8) = ((stszAtomSize >>  8) & UCHAR_MAX);
  m_dynamicHeader.at (m_ascSizeM5 +  9) = ( stszAtomSize        & UCHAR_MAX);
  m_dynamicHeader.at (m_ascSizeM5 + 22) = ((m_frameCount >> 24) & UCHAR_MAX);
  m_dynamicHeader.at (m_ascSizeM5 + 23) = ((m_frameCount >> 16) & UCHAR_MAX);
  m_dynamicHeader.at (m_ascSizeM5 + 24) = ((m_frameCount >>  8) & UCHAR_MAX);
  m_dynamicHeader.at (m_ascSizeM5 + 25) = ( m_frameCount        & UCHAR_MAX);

  m_dynamicHeader.push_back (0x00); m_dynamicHeader.push_back (0x00);
  m_dynamicHeader.push_back (0x00); m_dynamicHeader.push_back (stscAtomSize);
  m_dynamicHeader.push_back (0x73); m_dynamicHeader.push_back (0x74);
  m_dynamicHeader.push_back (0x73); m_dynamicHeader.push_back (0x63); // stsc
  m_dynamicHeader.push_back (0x00); m_dynamicHeader.push_back (0x00);
  m_dynamicHeader.push_back (0x00); m_dynamicHeader.push_back (0x00);
  m_dynamicHeader.push_back (0x00); m_dynamicHeader.push_back (0x00);
  m_dynamicHeader.push_back (0x00); m_dynamicHeader.push_back (numFramesFinalPeriod == 0 ? 1 : 2);
  m_dynamicHeader.push_back (0x00); m_dynamicHeader.push_back (0x00);
  m_dynamicHeader.push_back (0x00); m_dynamicHeader.push_back (0x01);  // 1st
  m_dynamicHeader.push_back ((numFramesFirstPeriod >> 24) & UCHAR_MAX);
  m_dynamicHeader.push_back ((numFramesFirstPeriod >> 16) & UCHAR_MAX);
  m_dynamicHeader.push_back ((numFramesFirstPeriod >>  8) & UCHAR_MAX);
  m_dynamicHeader.push_back ( numFramesFirstPeriod        & UCHAR_MAX);
  m_dynamicHeader.push_back (0x00); m_dynamicHeader.push_back (0x00);
  m_dynamicHeader.push_back (0x00); m_dynamicHeader.push_back (0x01);  // idx

  if (numFramesFinalPeriod > 0)
  {
    m_dynamicHeader.push_back ((m_rndAccOffsets.size () >> 24) & UCHAR_MAX);
    m_dynamicHeader.push_back ((m_rndAccOffsets.size () >> 16) & UCHAR_MAX);
    m_dynamicHeader.push_back ((m_rndAccOffsets.size () >>  8) & UCHAR_MAX);
    m_dynamicHeader.push_back ( m_rndAccOffsets.size ()        & UCHAR_MAX);
    m_dynamicHeader.push_back ((numFramesFinalPeriod >> 24) & UCHAR_MAX);
    m_dynamicHeader.push_back ((numFramesFinalPeriod >> 16) & UCHAR_MAX);
    m_dynamicHeader.push_back ((numFramesFinalPeriod >>  8) & UCHAR_MAX);
    m_dynamicHeader.push_back ( numFramesFinalPeriod        & UCHAR_MAX);
    m_dynamicHeader.push_back (0x00); m_dynamicHeader.push_back (0x00);
    m_dynamicHeader.push_back (0x00); m_dynamicHeader.push_back (0x01);// idx
  }

  m_dynamicHeader.push_back ((stcoAtomSize >> 24) & UCHAR_MAX);
  m_dynamicHeader.push_back ((stcoAtomSize >> 16) & UCHAR_MAX);
  m_dynamicHeader.push_back ((stcoAtomSize >>  8) & UCHAR_MAX);
  m_dynamicHeader.push_back ( stcoAtomSize        & UCHAR_MAX);
  m_dynamicHeader.push_back (0x73); m_dynamicHeader.push_back (0x74);
  m_dynamicHeader.push_back (0x63); m_dynamicHeader.push_back (0x6F); // stco
  m_dynamicHeader.push_back (0x00); m_dynamicHeader.push_back (0x00);
  m_dynamicHeader.push_back (0x00); m_dynamicHeader.push_back (0x00);
  m_dynamicHeader.push_back ((m_rndAccOffsets.size () >> 24) & UCHAR_MAX);
  m_dynamicHeader.push_back ((m_rndAccOffsets.size () >> 16) & UCHAR_MAX);
  m_dynamicHeader.push_back ((m_rndAccOffsets.size () >>  8) & UCHAR_MAX);
  m_dynamicHeader.push_back ( m_rndAccOffsets.size ()        & UCHAR_MAX);

  // add header size corrected random-access offsets to file
  for (unsigned i = 0; i < m_rndAccOffsets.size (); i++)
  {
    const uint32_t rndAccOffset = m_rndAccOffsets.at (i) + headerBytes;

    m_dynamicHeader.push_back ((rndAccOffset >> 24) & UCHAR_MAX);
    m_dynamicHeader.push_back ((rndAccOffset >> 16) & UCHAR_MAX);
    m_dynamicHeader.push_back ((rndAccOffset >>  8) & UCHAR_MAX);
    m_dynamicHeader.push_back ( rndAccOffset        & UCHAR_MAX);
  }
  m_dynamicHeader.push_back ((m_m4aMdatSize >> 24) & UCHAR_MAX);
  m_dynamicHeader.push_back ((m_m4aMdatSize >> 16) & UCHAR_MAX);
  m_dynamicHeader.push_back ((m_m4aMdatSize >>  8) & UCHAR_MAX);
  m_dynamicHeader.push_back ( m_m4aMdatSize        & UCHAR_MAX);
  m_dynamicHeader.push_back (0x6D); m_dynamicHeader.push_back (0x64);
  m_dynamicHeader.push_back (0x61); m_dynamicHeader.push_back (0x74); // mdat

  _SEEK (m_fileHandle, 0, 0 /*SEEK_SET*/);  // back to start

  bytesWritten += _WRITE (m_fileHandle, m_staticHeader, STAT_HEADER_SIZE);
  bytesWritten += _WRITE (m_fileHandle, &m_dynamicHeader.front (), (unsigned) m_dynamicHeader.size ());

  return bytesWritten;
}

int BasicMP4Writer::initHeader (const uint32_t audioLength) // reserve bytes for header in file
{
#if 0 // DEBUG
  const uint8_t numChannels = m_staticHeader[517];

  // write basic �-Law WAVE header for testing
  memcpy (m_staticHeader, muLawHeader, 44 * sizeof (uint8_t));
  m_staticHeader[22] = numChannels;
  m_staticHeader[24] = uint8_t (m_sampleRate & 0xFF);
  m_staticHeader[25] = uint8_t (m_sampleRate >>  8u);
  m_staticHeader[26] = uint8_t (m_sampleRate >> 16u);
  m_staticHeader[28] = uint8_t ((m_sampleRate * numChannels) & 0xFF);
  m_staticHeader[29] = uint8_t ((m_sampleRate * numChannels) >>  8u);
  m_staticHeader[30] = uint8_t ((m_sampleRate * numChannels) >> 16u);
  m_staticHeader[32] = numChannels;  // byte count per frame

  return _WRITE (m_fileHandle, m_staticHeader, 44);
#else
  const bool flushFrameUsed = ((audioLength + m_pregapLength) % m_frameLength) > 0;
  const unsigned frameCount = ((audioLength + m_frameLength - 1) / m_frameLength) + (flushFrameUsed ? 2 : 1);
  const unsigned chunkCount = ((frameCount + m_rndAccPeriod - 1) / m_rndAccPeriod);
  const unsigned finalChunk = (frameCount <= m_rndAccPeriod ? 0 : frameCount % m_rndAccPeriod);
  const int estimHeaderSize = STAT_HEADER_SIZE + m_ascSizeM5 + 6+4 + frameCount * 4 /*stsz*/ + STSX_BSIZE * 3 +
                              (finalChunk == 0 ? 12 : 24) /*stsc*/ + chunkCount * 4 /*stco*/ + 8 /*mdat*/;
  int bytesWritten = 0;

  for (int i = estimHeaderSize; i > 0; i -= STAT_HEADER_SIZE)
  {
    bytesWritten += _WRITE (m_fileHandle, m_staticHeader, __min (i, STAT_HEADER_SIZE));
  }
  return bytesWritten;
#endif
}

unsigned BasicMP4Writer::open (const int mp4FileHandle, const unsigned sampleRate,  const unsigned numChannels,
                               const unsigned bitDepth, const unsigned frameLength, const unsigned pregapLength,
                               const unsigned raPeriod, const uint8_t* ascBuf,      const unsigned ascSize,
                               const uint32_t creatTime /*= 0*/, const char vbrQuality /*= 0*/)
{
  const uint32_t  frameSizeBE = toBigEndian (frameLength);
  const uint32_t pregapSizeBE = toBigEndian (pregapLength);
  const uint32_t sampleRateBE = toBigEndian (sampleRate);
  const uint32_t  timeStampBE = toBigEndian (creatTime);
  uint32_t* const header4Byte = (uint32_t* const) m_staticHeader;

  if ((mp4FileHandle == -1) || (frameLength == 0) || (sampleRate == 0) || (numChannels == 0) || (numChannels * 3 > UCHAR_MAX) ||
      (raPeriod == 0) || (ascBuf == nullptr) || (ascSize < 5) || (ascSize > 108) || (bitDepth == 0) || (bitDepth > UCHAR_MAX))
  {
    return 1; // invalid file handle or other input variable
  }

  m_fileHandle = mp4FileHandle;
  reset (frameLength, pregapLength, raPeriod);
#if 0 // DEBUG
  m_sampleRate = sampleRate;
#endif
  // create fixed-length 576-byte part of MPEG-4 file header
  memcpy (m_staticHeader, staticHeaderTemplate, STAT_HEADER_SIZE * sizeof (uint8_t));

  header4Byte[ 44>>2] = timeStampBE;
  header4Byte[ 48>>2] = timeStampBE;
  header4Byte[ 52>>2] = sampleRateBE;
  header4Byte[184>>2] = timeStampBE;
  header4Byte[188>>2] = timeStampBE;
  header4Byte[292>>2] = pregapSizeBE; // pregap size in elst
  header4Byte[320>>2] = timeStampBE;
  header4Byte[324>>2] = timeStampBE;
  header4Byte[328>>2] = sampleRateBE;
  header4Byte[332>>2] = pregapSizeBE; // +audio length later
  header4Byte[464>>2] = frameSizeBE;

  m_staticHeader[339] = vbrQuality;
  m_staticHeader[517] = (uint8_t) numChannels;
  m_staticHeader[519] = (uint8_t) bitDepth;
  m_staticHeader[523] = (sampleRate >> 16) & UCHAR_MAX; // ?
  m_staticHeader[524] = (sampleRate >>  8) & UCHAR_MAX;
  m_staticHeader[525] = sampleRate & UCHAR_MAX;
  m_staticHeader[556] = (uint8_t) numChannels * 3; // 6144 bit/chan

  memcpy (&m_staticHeader[571], ascBuf, 5 * sizeof (uint8_t));

  if (ascSize > 5) // increase atom byte-sizes
  {
    const uint8_t inc = m_ascSizeM5 = ascSize - 5;

    m_staticHeader[ 27] += inc;  // MOOV_BSIZE
    m_staticHeader[167] += inc;  // TRAK_BSIZE
    m_staticHeader[303] += inc;  // MDIA_BSIZE
    if (m_staticHeader[379] + m_ascSizeM5 > UCHAR_MAX) m_staticHeader[378]++;
    m_staticHeader[379] += inc;  // MINF_BSIZE
    m_staticHeader[439] += inc;  // STBL_BSIZE
    m_staticHeader[479] += inc;  // STSD_BSIZE
    m_staticHeader[495] += inc;  // MP4A_BSIZE
    m_staticHeader[531] += inc;  // ESDS_BSIZE
    m_staticHeader[544] += inc;  // esds tag 3
    m_staticHeader[552] += inc;  // esds tag 4
    m_staticHeader[570] += inc;  // esds tag 5

    for (unsigned i = 0; i < m_ascSizeM5; i++) m_dynamicHeader.push_back (ascBuf[5 + i]);
  }

  // prepare variable-length remainder of MPEG-4 file header
  m_dynamicHeader.push_back (0x06); m_dynamicHeader.push_back (0x80); // esds
  m_dynamicHeader.push_back (0x80); m_dynamicHeader.push_back (0x80);
  m_dynamicHeader.push_back (0x01); m_dynamicHeader.push_back (0x02);
  m_dynamicHeader.push_back (0x00); m_dynamicHeader.push_back (0x00); // + 4N
  m_dynamicHeader.push_back (0x00); m_dynamicHeader.push_back (STSX_BSIZE + 4);
  m_dynamicHeader.push_back (0x73); m_dynamicHeader.push_back (0x74);
  m_dynamicHeader.push_back (0x73); m_dynamicHeader.push_back (0x7A); // stsz
  m_dynamicHeader.push_back (0x00); m_dynamicHeader.push_back (0x00);
  m_dynamicHeader.push_back (0x00); m_dynamicHeader.push_back (0x00);
  m_dynamicHeader.push_back (0x00); m_dynamicHeader.push_back (0x00);
  m_dynamicHeader.push_back (0x00); m_dynamicHeader.push_back (0x00);
  m_dynamicHeader.push_back (0x00); m_dynamicHeader.push_back (0x00); // := N
  m_dynamicHeader.push_back (0x00); m_dynamicHeader.push_back (0x00);

  return 0; // correct operation
}

void BasicMP4Writer::reset (const unsigned frameLength /*= 0*/, const unsigned pregapLength /*= 0*/, const unsigned raPeriod /*= 0*/)
{
  m_ascSizeM5    = 0;
  m_frameCount   = 0;
  m_frameLength  = frameLength;
  m_m4aMdatSize  = 8; // bytes for mdat header
  m_pregapLength = pregapLength;
  m_rndAccPeriod = raPeriod;
  m_sampleRate   = 0;
  m_dynamicHeader.clear ();
  m_rndAccOffsets.clear ();

  if (m_fileHandle != -1) _SEEK (m_fileHandle, 0, 0 /*SEEK_SET*/);
}
