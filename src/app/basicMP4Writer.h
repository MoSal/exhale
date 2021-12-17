/* basicMP4Writer.h - header file for class with basic MPEG-4 file writing capability
 * written by C. R. Helmrich, last modified in 2021 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under the exhale Copyright License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2021 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#ifndef _BASIC_MP4_WRITER_H_
#define _BASIC_MP4_WRITER_H_

#include "exhaleAppPch.h"

// constant data sizes in bytes
#define STAT_HEADER_SIZE   576
#define STSX_BSIZE        0x10
#define UDTA_BSIZE        0x61 // udta: 0 to turn off!
#define ESDS_BSIZE  0x00, 0x36 // esds: 54 (+ m_ascSizeM5 later)
#define MP4A_BSIZE  0x00, 0x5A // mp4a: 36 + ESDS_BSIZE
#define STSD_BSIZE  0x00, 0x6A // mp4a: 16 + MP4A_BSIZE
#define STBL_BSIZE  0x00, 0x92 // stbl: 8 + 32 + STSD_BSIZE (+ rem later)
#define MINF_BSIZE  0x00, 0xCE // minf: 8 + 16 + 36 + STBL_BSIZE
#define MDIA_BSIZE  0x01, 0x1A // mdia: 8 + 32 + 36 + MINF_BSIZE
#define TRAK_BSIZE  0x01, 0xA2 // trak: 8 + 92 + 36 + MDIA_BSIZE
#define MOOV_BSIZE  0x02, 0x2E // moov: 8 +108 + 24 + TRAK_BSIZE

// basic MPEG-4 write-out class
class BasicMP4Writer
{
private:

  // member variables
  unsigned m_ascSizeM5;  // ASC + UsacConfig byte-size - 5
  int      m_fileHandle;
  unsigned m_frameCount;
  unsigned m_frameLength;
  uint32_t m_mediaOffset;  // offset of first mdat payload
  uint32_t m_mediaSize; // number of bytes of mdat content
  unsigned m_preLength;   // encoding look-ahead, pre-roll
  unsigned m_postLength; // decoding look-ahead, post-roll
  unsigned m_rndAccPeriod;  // random-access (RA) interval
  unsigned m_sampleRate;
  uint8_t  m_staticHeader[STAT_HEADER_SIZE]; // fixed-size
  std::vector <uint8_t> m_dynamicHeader; // variable-sized
  std::vector <uint32_t> m_rndAccOffsets; // random access
#ifndef NO_PREROLL_DATA
  std::vector <uint8_t> m_ipfCfgOffsets; // IPF UsacConfig
#endif

  // helper function
  void push32BitValue (const uint32_t value); // to header

public:

  // constructor
  BasicMP4Writer () { m_fileHandle = -1;  reset (0, 0, 0, 0); }
  // destructor
#ifdef NO_PREROLL_DATA
  ~BasicMP4Writer() { m_dynamicHeader.clear (); m_rndAccOffsets.clear (); }
#else
  ~BasicMP4Writer() { m_dynamicHeader.clear (); m_rndAccOffsets.clear (); m_ipfCfgOffsets.clear (); }
#endif
  // public functions
  int addFrameAU (const uint8_t* byteBuf, const uint32_t byteCount);
  int finishFile (const unsigned avgBitrate, const unsigned maxBitrate, const uint32_t audioLength,
                  const uint32_t modifTime = 0, const uint8_t* ascBuf = nullptr);
  unsigned getFrameCount () const { return m_frameCount; }
  int initHeader (const uint32_t audioLength, const unsigned extraDelay);
  unsigned open  (const int mp4FileHandle, const unsigned sampleRate,  const unsigned numChannels,
                  const unsigned bitDepth, const unsigned frameLength, const unsigned pregapLength,
                  const unsigned raPeriod, const uint8_t* ascBuf,      const unsigned ascSize,
                  const uint32_t creatTime = 0, const char vbrQuality = 0);
  void     reset (const unsigned frameLength, const unsigned pregapLength, const unsigned raPeriod, const unsigned sampleRate);
#ifndef NO_PREROLL_DATA
  int updateIPFs (const uint8_t* ascUcBuf, const uint32_t ascUcLength, const uint32_t ucOffset);
#endif
}; // BasicMP4Writer

#endif // _BASIC_MP4_WRITER_H_
