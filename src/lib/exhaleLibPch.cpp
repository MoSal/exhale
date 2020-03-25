/* exhaleLibPch.cpp - pre-compiled source file for classes of exhaleLib coding library
 * written by C. R. Helmrich, last modified in 2019 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under a Modified BSD-Style License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2020 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#include "exhaleLibPch.h"

// public bit-stream functions
void OutputStream::reset () // clear writer states and byte buffer
{
  heldBitChunk = 0;
  heldBitCount = 0;
  stream.clear ();
}

void OutputStream::write (const uint32_t bitChunk, const uint8_t bitCount)
{
  if (bitCount > 32) return; // only a maximum of 32 bits is writable at once

  const uint8_t totalBitCount   = bitCount + heldBitCount;
  const uint8_t totalByteCount  = totalBitCount >> 3;  // to be written
  const uint8_t newHeldBitCount = totalBitCount & 7; // not yet written
  const uint8_t newHeldBitChunk = (bitChunk << (8 - newHeldBitCount)) & UCHAR_MAX;

  if (totalByteCount == 0) // not enough bits to write, only update held bits
  {
    heldBitChunk |= newHeldBitChunk;
  }
  else // write bits
  {
    const uint32_t writtenChunk = (heldBitChunk << uint32_t ((bitCount - newHeldBitCount) & ~7)) | (bitChunk >> newHeldBitCount);
    switch (totalByteCount)
    {
      case 4: stream.push_back (writtenChunk >> 24);
      case 3: stream.push_back (writtenChunk >> 16);
      case 2: stream.push_back (writtenChunk >> 8);
      case 1: stream.push_back (writtenChunk);
    }
    heldBitChunk = newHeldBitChunk;
  }
  heldBitCount = newHeldBitCount;
}

// ISO/IEC 23003-3, Table 67
static const unsigned allowedSamplingRates[USAC_NUM_SAMPLE_RATES] = {
  96000, 88200, 64000, 48000, 44100, 32000, 24000, 22050, 16000, 12000, 11025,  8000, 7350, // AAC
  57600, 51200, 40000, 38400, 34150, 28800, 25600, 20000, 19200, 17075, 14400, 12800, 9600 // USAC
};

// public sampling rate functions
char toSamplingFrequencyIndex (const unsigned samplingRate)
{
  for (char i = 0; i < AAC_NUM_SAMPLE_RATES; i++)
  {
    if (samplingRate == allowedSamplingRates[(int) i]) // AAC rate
    {
      return i;
    }
#if !RESTRICT_TO_AAC
    if (samplingRate == allowedSamplingRates[i + AAC_NUM_SAMPLE_RATES] && (samplingRate % 19200) == 0) // Baseline USAC
    {
      return i + AAC_NUM_SAMPLE_RATES + 2;  // skip reserved entry
    }
#endif
  }
  return -1; // no index found
}

unsigned toSamplingRate (const char samplingFrequencyIndex)
{
#if RESTRICT_TO_AAC
  if ((samplingFrequencyIndex < 0) || (samplingFrequencyIndex >= AAC_NUM_SAMPLE_RATES))
#else
  if ((samplingFrequencyIndex < 0) || (samplingFrequencyIndex >= USAC_NUM_SAMPLE_RATES + 2))
#endif
  {
    return 0; // invalid index
  }
  return allowedSamplingRates[samplingFrequencyIndex > AAC_NUM_SAMPLE_RATES ? samplingFrequencyIndex - 2 : samplingFrequencyIndex];
}
