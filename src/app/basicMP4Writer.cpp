/* basicMP4Writer.cpp - source file for class with basic MPEG-4 file writing capability
 * written by C. R. Helmrich, last modified in 2021 - see License.htm for legal notices
 * pre-roll serializer and related code added by J. Calhoun in 2020, see merge request 4
 *
 * The copyright in this software is being made available under the exhale Copyright License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2021 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#include "exhaleAppPch.h"
#include "basicMP4Writer.h"
#include "version.h"

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
static uint32_t toBigEndian (const unsigned ui) // to Motorola endianness
{
  return ((ui & UCHAR_MAX) << 24) | (((ui >> 8) & UCHAR_MAX) << 16) | (((ui >> 16) & UCHAR_MAX) << 8) | ((ui >> 24) & UCHAR_MAX);
}

static uint16_t toUShortValue (const uint8_t hiByte, const uint8_t loByte)
{
  return ((uint16_t) hiByte << 8) | (uint16_t) loByte;
}

// private helper function
void BasicMP4Writer::push32BitValue (const uint32_t value) // push to dynamic header
{
  m_dynamicHeader.push_back ((value >> 24) & UCHAR_MAX);
  m_dynamicHeader.push_back ((value >> 16) & UCHAR_MAX);
  m_dynamicHeader.push_back ((value >>  8) & UCHAR_MAX);
  m_dynamicHeader.push_back ((value      ) & UCHAR_MAX);
}

// public functions
int BasicMP4Writer::addFrameAU (const uint8_t* byteBuf, const uint32_t byteCount)
{
  if ((m_fileHandle == -1) || (m_mediaSize > 0xFFFFFFF0u - m_mediaOffset - byteCount))
  {
    return 1; // invalid file handle or file getting too big
  }

  // add frame byte-size, in Big Endian format, to frame size list (stsz)
  push32BitValue (byteCount);

  if (((m_frameCount++) % m_rndAccPeriod) == 0) // add RAP to list (stco)
  {
    m_rndAccOffsets.push_back (m_mediaSize);
#ifndef NO_PREROLL_DATA
    if (((m_frameCount - 1u) % (m_rndAccPeriod << 1)) == 0)  // every 2nd
    {
      if ((byteBuf[0] & 0xE0) == 0xC0) // IPF?
      {
        // byte-offset of UsacConfig() in AU (excl. first 5 config bits!)
        m_ipfCfgOffsets.push_back (byteBuf[0] == 0xDF && (byteBuf[1] & 0xE0) == 0xE0 ? 5 : 3);
      }
      else // 3-bit ID is missing, not an IPF!
      {
        m_ipfCfgOffsets.push_back (0);
      }
    }
#endif
  }
  m_mediaSize += byteCount;

  return _WRITE (m_fileHandle, byteBuf, byteCount);  // write access unit
}

int BasicMP4Writer::finishFile (const unsigned avgBitrate, const unsigned maxBitrate, const uint32_t audioLength,
                                const uint32_t modifTime /*= 0*/, const uint8_t* ascBuf /*= nullptr*/)
{
  const unsigned numFramesFirstPeriod = __min (m_frameCount, m_rndAccPeriod);
  const unsigned numFramesFinalPeriod = (m_frameCount <= m_rndAccPeriod ? 0 : m_frameCount % m_rndAccPeriod);
  const unsigned numSamplesFinalFrame = (audioLength + m_preLength + m_postLength) % m_frameLength;
  const uint32_t raOffsetSize = (uint32_t) m_rndAccOffsets.size ();
  const uint32_t stszAtomSize = STSX_BSIZE + 4 /*bytes for sampleSize*/ + m_frameCount * 4;
  const uint32_t stscAtomSize = STSX_BSIZE + (numFramesFinalPeriod == 0 ? 12 : 24);
  const uint32_t stcoAtomSize = STSX_BSIZE + raOffsetSize * 4;
#ifdef NO_PREROLL_DATA
  const uint32_t stssAtomSize = STSX_BSIZE + 4;
#else
  const uint32_t preRollCount = (m_frameCount + (m_rndAccPeriod << 1) - 1) / (m_rndAccPeriod << 1);
  const uint32_t stssAtomSize = STSX_BSIZE + preRollCount * 4;
#endif
  /* The following code creates a 'prol' sample group with a repeating pattern of membership,
  indicating that the first sample in each increment of m_rndAccPeriod samples is a member and
  thus an independent frame (but not an immediate playout frame), and the remainder are not. */
  const uint32_t sgpdAtomSize = STSX_BSIZE + 4 /*defaultLength == 2*/ + 4 /*entryCount == 1*/ + 2 /*rollDistance*/;
  const unsigned patternEntryCount = __min (m_frameCount, m_rndAccPeriod << 1);
  const uint32_t patternLengthSize = (patternEntryCount > UINT8_MAX ? 2 : 1);
  const uint32_t compPatternLength = (patternEntryCount + 1) >> 1; // four bit per index with byte alignment
  const uint32_t sampleCountSize   = (m_frameCount > UINT16_MAX ? 4 : (m_frameCount > UINT8_MAX ? 2 : 1));
  const uint32_t csgpAtomSize = STSX_BSIZE + 4 /*patternCount == 1*/ + patternLengthSize + sampleCountSize + compPatternLength;
  const uint32_t stblIncrSize = m_ascSizeM5 + stszAtomSize + stscAtomSize + stcoAtomSize + stssAtomSize + sgpdAtomSize + csgpAtomSize;
  const uint32_t moovAtomSize = toBigEndian (toUShortValue (MOOV_BSIZE) + stblIncrSize + UDTA_BSIZE);
  const uint32_t trakAtomSize = toBigEndian (toUShortValue (TRAK_BSIZE) + stblIncrSize);
  const uint32_t mdiaAtomSize = toBigEndian (toUShortValue (MDIA_BSIZE) + stblIncrSize);
  const uint32_t minfAtomSize = toBigEndian (toUShortValue (MINF_BSIZE) + stblIncrSize);
  const uint32_t stblAtomSize = toBigEndian (toUShortValue (STBL_BSIZE) + stblIncrSize);
  const uint32_t numSamplesBE = toBigEndian (audioLength);
  const uint32_t  timeStampBE = toBigEndian (modifTime);
  uint32_t* const header4Byte = (uint32_t* const) m_staticHeader;
  int bytesWritten = 0;
  uint32_t i;

  if ((m_fileHandle == -1) || (m_mediaSize > 0xFFFFFFF0u - m_mediaOffset))
  {
    return 1; // invalid file handle or file getting too big
  }

  if (ascBuf != nullptr) // update ASC + UC data if required
  {
    memcpy (&m_staticHeader[571], ascBuf, 5 * sizeof (uint8_t));

    for (i = 0; i < m_ascSizeM5; i++) m_dynamicHeader.at (i) = ascBuf[5 + i];
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
  header4Byte[332>>2] = toBigEndian (audioLength + m_preLength);
  header4Byte[376>>2] = minfAtomSize;
  header4Byte[288>>2] = numSamplesBE;  // elst
  header4Byte[436>>2] = stblAtomSize;
  header4Byte[460>>2] = toBigEndian (m_frameCount - 1); // 2 entries used
  header4Byte[472>>2] = toBigEndian ((numSamplesFinalFrame == 0 ? m_frameLength : __max (m_postLength + 1u, numSamplesFinalFrame)) - m_postLength);

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

  push32BitValue (stscAtomSize);
  m_dynamicHeader.push_back (0x73); m_dynamicHeader.push_back (0x74);
  m_dynamicHeader.push_back (0x73); m_dynamicHeader.push_back (0x63); // stsc
  push32BitValue (0);
  push32BitValue (numFramesFinalPeriod == 0 ? 1 : 2);
  push32BitValue (1); // 1st
  push32BitValue (numFramesFirstPeriod);
  push32BitValue (1); // idx

  if (numFramesFinalPeriod > 0)
  {
    push32BitValue (raOffsetSize);
    push32BitValue (numFramesFinalPeriod);
    push32BitValue (1);
  }

  push32BitValue (stcoAtomSize);
  m_dynamicHeader.push_back (0x73); m_dynamicHeader.push_back (0x74);
  m_dynamicHeader.push_back (0x63); m_dynamicHeader.push_back (0x6F); // stco
  push32BitValue (0);
  push32BitValue (raOffsetSize);

  // add header size corrected random-access offsets to file
  for (i = 0; i < raOffsetSize; i++) push32BitValue (m_rndAccOffsets.at (i) + m_mediaOffset);

  push32BitValue (stssAtomSize);
  m_dynamicHeader.push_back (0x73); m_dynamicHeader.push_back (0x74);
  m_dynamicHeader.push_back (0x73); m_dynamicHeader.push_back (0x73); // stss
  push32BitValue (0);
#ifdef NO_PREROLL_DATA
  push32BitValue (1); // 1st
  push32BitValue (1); // AU!
#else
  push32BitValue (preRollCount);
  for (i = 1; i <= m_frameCount; i += (m_rndAccPeriod << 1)) push32BitValue (i);
#endif
  push32BitValue (sgpdAtomSize);
  m_dynamicHeader.push_back (0x73); m_dynamicHeader.push_back (0x67);
  m_dynamicHeader.push_back (0x70); m_dynamicHeader.push_back (0x64); // sgpd
  push32BitValue (1 << 24);
  m_dynamicHeader.push_back (0x70); m_dynamicHeader.push_back (0x72);
  m_dynamicHeader.push_back (0x6F); m_dynamicHeader.push_back (0x6C); // prol
  push32BitValue (2);
  push32BitValue (1);
  m_dynamicHeader.push_back (0);  m_dynamicHeader.push_back (m_frameLength > 1024 ? 2 : 1); // roll_distance
  push32BitValue (csgpAtomSize);
  m_dynamicHeader.push_back (0x63); m_dynamicHeader.push_back (0x73);
  m_dynamicHeader.push_back (0x67); m_dynamicHeader.push_back (0x70); // csgp
  push32BitValue ((__min (3, patternLengthSize) << 4) | (__min (3, sampleCountSize) << 2));
  m_dynamicHeader.push_back (0x70); m_dynamicHeader.push_back (0x72);
  m_dynamicHeader.push_back (0x6F); m_dynamicHeader.push_back (0x6C); // prol
  push32BitValue (1);
  if (patternLengthSize > 1)
  {
    m_dynamicHeader.push_back ((patternEntryCount >> 8) & UCHAR_MAX);
  }
  m_dynamicHeader.push_back   ( patternEntryCount       & UCHAR_MAX);
  if (sampleCountSize > 2)
  {
    m_dynamicHeader.push_back ((m_frameCount >> 24) & UCHAR_MAX);
    m_dynamicHeader.push_back ((m_frameCount >> 16) & UCHAR_MAX);
  }
  if (sampleCountSize > 1)
  {
    m_dynamicHeader.push_back ((m_frameCount >>  8) & UCHAR_MAX);
  }
  m_dynamicHeader.push_back   ( m_frameCount        & UCHAR_MAX);

  for (i = 0; i < (numFramesFirstPeriod >> 1); i++) m_dynamicHeader.push_back (0); // nonmembers, first part
  if (i < compPatternLength)
  {
    m_dynamicHeader.push_back (1 << ((numFramesFirstPeriod & 1) ? 0 : 4));
    for (i++; i < compPatternLength; i++) m_dynamicHeader.push_back (0);  // rest of nonmembers, second part
  }
#if UDTA_BSIZE
  const char ver[] = EXHALELIB_VERSION_MAJOR "." EXHALELIB_VERSION_MINOR EXHALELIB_VERSION_BUGFIX;

  push32BitValue (UDTA_BSIZE);
  m_dynamicHeader.push_back (0x75); m_dynamicHeader.push_back (0x64);
  m_dynamicHeader.push_back (0x74); m_dynamicHeader.push_back (0x61); // udta
  push32BitValue (UDTA_BSIZE - 8);
  m_dynamicHeader.push_back (0x6D); m_dynamicHeader.push_back (0x65);
  m_dynamicHeader.push_back (0x74); m_dynamicHeader.push_back (0x61); // meta
  push32BitValue (0);

  push32BitValue (33);
  m_dynamicHeader.push_back (0x68); m_dynamicHeader.push_back (0x64);
  m_dynamicHeader.push_back (0x6C); m_dynamicHeader.push_back (0x72); // hdlr
  push32BitValue (0);
  push32BitValue (0);
  m_dynamicHeader.push_back (0x6D); m_dynamicHeader.push_back (0x64);
  m_dynamicHeader.push_back (0x69); m_dynamicHeader.push_back (0x72); // mdir
  m_dynamicHeader.push_back (0x61); m_dynamicHeader.push_back (0x70);
  m_dynamicHeader.push_back (0x70); m_dynamicHeader.push_back (0x6C); // appl
  push32BitValue (0);
  push32BitValue (0);
  m_dynamicHeader.push_back (0);

  push32BitValue (UDTA_BSIZE - 53);
  m_dynamicHeader.push_back (0x69); m_dynamicHeader.push_back (0x6C);
  m_dynamicHeader.push_back (0x73); m_dynamicHeader.push_back (0x74); // ilst
  push32BitValue (UDTA_BSIZE - 53 - 8);
  m_dynamicHeader.push_back (0xA9); m_dynamicHeader.push_back (0x74);
  m_dynamicHeader.push_back (0x6F); m_dynamicHeader.push_back (0x6F); // ©too
  push32BitValue (UDTA_BSIZE - 53 - 16);
  m_dynamicHeader.push_back (0x64); m_dynamicHeader.push_back (0x61);
  m_dynamicHeader.push_back (0x74); m_dynamicHeader.push_back (0x61); // data
  push32BitValue (1);
  push32BitValue (0);
  m_dynamicHeader.push_back (0x65); m_dynamicHeader.push_back (0x78);
  m_dynamicHeader.push_back (0x68); m_dynamicHeader.push_back (0x61);
  m_dynamicHeader.push_back (0x6C); m_dynamicHeader.push_back (0x65); // exhale
  m_dynamicHeader.push_back (0x20);
  for (i = 0; i < 5; i++) m_dynamicHeader.push_back ((uint8_t) ver[i]);
#endif
  const uint32_t moovAndMdatOverhead = STAT_HEADER_SIZE + (uint32_t) m_dynamicHeader.size () + 8;
  const uint32_t headerPaddingLength = uint32_t (m_mediaOffset - moovAndMdatOverhead);

  if (moovAndMdatOverhead > m_mediaOffset) // header has grown to encroach upon the media data - fatal error
  {
    return 1;
  }
  else // start the 'mdat' atom early with a small amount of unused but informative data before the first AU
  {
    m_mediaSize += headerPaddingLength;
  }

  push32BitValue (m_mediaSize + 8);
  m_dynamicHeader.push_back (0x6D); m_dynamicHeader.push_back (0x64);
  m_dynamicHeader.push_back (0x61); m_dynamicHeader.push_back (0x74); // mdat
  for (uint32_t pNdx = 0; pNdx < headerPaddingLength; pNdx++)
  {
    if (pNdx == 0)  // add padding byte with library version
    {
#if !UDTA_BSIZE
      const char ver[] = EXHALELIB_VERSION_MAJOR "." EXHALELIB_VERSION_MINOR EXHALELIB_VERSION_BUGFIX;
#endif
      const int verInt = (ver[0] - 0x30) * 100 + (ver[2] - 0x30) * 10 + (ver[4] - 0x30);

      m_dynamicHeader.push_back (__max (0, __min (UCHAR_MAX, verInt)));
    }
    else if (pNdx == 1) // add 8-bit cyclic redundancy check
    {
      uint8_t crc8 = m_dynamicHeader.back(); // Baicheva '98

      for (i = 8; i > 0; i--) if (crc8 & 0x80) crc8 = (crc8 << 1) ^ 0x2F; else crc8 <<= 1;
      m_dynamicHeader.push_back (crc8); // add padding CRC-8
    }
    else
    {
      m_dynamicHeader.push_back (0); // add one padding byte
    }
  }

  _SEEK (m_fileHandle, 0, 0 /*SEEK_SET*/);  // back to start

  bytesWritten += _WRITE (m_fileHandle, m_staticHeader, STAT_HEADER_SIZE);
  bytesWritten += _WRITE (m_fileHandle, &m_dynamicHeader.front (), (uint32_t) m_dynamicHeader.size ());

  return bytesWritten;
}

int BasicMP4Writer::initHeader (const uint32_t audioLength, const unsigned extraDelay)
{
  const unsigned frameCount = (audioLength + m_preLength + m_postLength - extraDelay + m_frameLength - 1u) / m_frameLength;
  const unsigned chunkCount = ((frameCount + m_rndAccPeriod - 1) / m_rndAccPeriod);
  const unsigned numFramesFinalPeriod = (frameCount <= m_rndAccPeriod ? 0 : frameCount % m_rndAccPeriod);
  const unsigned patternEntryCount = __min (frameCount, m_rndAccPeriod << 1);
  const unsigned smpGrpSize = 10 /*sgpd*/ + (patternEntryCount > UINT8_MAX ? 10 : 9) + ((patternEntryCount + 1) >> 1) /*csgp*/;
  const int estimHeaderSize = STAT_HEADER_SIZE + m_ascSizeM5 + 6 + 4 + frameCount * 4 /*stsz*/ + STSX_BSIZE * 6 + smpGrpSize + chunkCount * 4 /*stco*/ +
#ifdef NO_PREROLL_DATA
                              4 /*minimum stss*/ + UDTA_BSIZE +
#else
                              ((chunkCount + 1) >> 1) * 4 /*stss*/ + UDTA_BSIZE +
#endif
                              (numFramesFinalPeriod == 0 ? (frameCount > m_rndAccPeriod && m_frameLength == 2048 ? 20 : 12) : 24) /*stsc*/ + 8 /*mdat*/;
  int bytesWritten = 0;

  for (int i = estimHeaderSize; i > 0; i -= STAT_HEADER_SIZE)
  {
    bytesWritten += _WRITE (m_fileHandle, m_staticHeader, __min (i, STAT_HEADER_SIZE));
  }
  m_mediaOffset = bytesWritten; // first frame will be written at this offset
  m_postLength -= __min (m_postLength, extraDelay);

  return bytesWritten;
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
  reset (frameLength, pregapLength, __min (USHRT_MAX, raPeriod), sampleRate);

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
  m_staticHeader[556] = (uint8_t) numChannels * 3; // = 6144 bits/channel

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
  push32BitValue (STSX_BSIZE + 4); // + 4count
  m_dynamicHeader.push_back (0x73); m_dynamicHeader.push_back (0x74);
  m_dynamicHeader.push_back (0x73); m_dynamicHeader.push_back (0x7A); // stsz
  push32BitValue (0);
  push32BitValue (0);
  push32BitValue (0);

  return 0; // correct operation
}

void BasicMP4Writer::reset (const unsigned frameLength, const unsigned pregapLength, const unsigned raPeriod, const unsigned sampleRate)
{
  m_ascSizeM5    = 0;
  m_frameCount   = 0;
  m_frameLength  = frameLength;
  m_mediaOffset  = 0;  // offset of first 'mdat' data byte serialized to file
  m_mediaSize    = 0; // total length of 'mdat' (access unit) payload in file
  m_preLength    = pregapLength;
  m_postLength   = sampleRate / 200u;
  m_rndAccPeriod = raPeriod;
  m_sampleRate   = sampleRate;
  m_dynamicHeader.clear ();
  m_rndAccOffsets.clear ();
#ifndef NO_PREROLL_DATA
  m_ipfCfgOffsets.clear ();
#endif

  if (m_fileHandle != -1) _SEEK (m_fileHandle, 0, 0 /*SEEK_SET*/);
}

#ifndef NO_PREROLL_DATA
int BasicMP4Writer::updateIPFs (const uint8_t* ascUcBuf, const uint32_t ascUcLength, const uint32_t ucOffset)
{
  const uint8_t bw = uint8_t (ascUcLength - ucOffset);
  int bytesWritten = 0, configsWritten = 0;

  if ((m_fileHandle == -1) || (ascUcBuf == nullptr) || (ascUcLength == 0) || (ascUcLength <= ucOffset))
  {
    return 1; // invalid file handle or IPF config parameter
  }

  // write updated UsacConfig() to AudioPreRoll() extensions
  for (uint32_t i = 0; i < (uint32_t) m_ipfCfgOffsets.size (); i++)
  {
    const uint32_t configOffset = m_ipfCfgOffsets.at (i);

    if (configOffset > 0) // this AU is an IPF
    {
      _SEEK (m_fileHandle, m_rndAccOffsets.at (i << 1) + m_mediaOffset + configOffset, 0 /*SEEK_SET*/);
      bytesWritten += _WRITE (m_fileHandle, &ascUcBuf[ucOffset], bw);
      configsWritten++;
    }
  }

  return (bytesWritten != configsWritten * bw ? 1 : 0);
}
#endif
