/* exhaleLibPch.cpp - pre-compiled source file for classes of exhaleLib coding library
 * written by C. R. Helmrich, last modified in 2020 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under the exhale Copyright License
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
  if (bitCount == 0) return; // nothing to do for length 0, max. length is 32

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

static int8_t quantizeSbrEnvelopeLevel (const uint64_t energy, const uint32_t divisor, const uint8_t noiseLevel)
{
  const double ener = (double) energy / (double) divisor;

  return (ener > 8192.0 ? int8_t (1.375 - 0.03125 * noiseLevel + 6.64385619 * log10 (ener)) - 26 : 0);
}

// public SBR related functions
int32_t getSbrEnvelopeAndNoise (int32_t* const sbrLevels, const uint8_t specFlat5b, const uint8_t tempFlat5b, const bool lr, const bool ind,
                                const uint8_t specFlatSte, const int32_t tmpValSte, const uint32_t frameSize, int32_t* sbrData)
{
  const uint64_t enValue[8] = {(uint64_t) sbrLevels[22] * (uint64_t) sbrLevels[22],
                               (uint64_t) sbrLevels[23] * (uint64_t) sbrLevels[23],
                               (uint64_t) sbrLevels[24] * (uint64_t) sbrLevels[24],
                               (uint64_t) sbrLevels[25] * (uint64_t) sbrLevels[25],
                               (uint64_t) sbrLevels[26] * (uint64_t) sbrLevels[26],
                               (uint64_t) sbrLevels[27] * (uint64_t) sbrLevels[27],
                               (uint64_t) sbrLevels[28] * (uint64_t) sbrLevels[28],
                               (uint64_t) sbrLevels[11] * (uint64_t) sbrLevels[11]};
  const uint64_t envTmp0[1] = { enValue[0] + enValue[1] + enValue[2] + enValue[3] +
                                enValue[4] + enValue[5] + enValue[6] + enValue[7]};
  const uint64_t envTmp1[2] = {(enValue[0] + enValue[1] + enValue[2] + enValue[3]) << 1,
                               (enValue[4] + enValue[5] + enValue[6] + enValue[7]) << 1};
  const uint64_t envTmp2[4] = {(enValue[0] + enValue[1]) << 2, (enValue[2] + enValue[3]) << 2,
                               (enValue[4] + enValue[5]) << 2, (enValue[6] + enValue[7]) << 2};
  const uint64_t envTmp3[8] = { enValue[0] << 3, enValue[1] << 3, enValue[2] << 3, enValue[3] << 3,
                                enValue[4] << 3, enValue[5] << 3, enValue[6] << 3, enValue[7] << 3};
  uint64_t errTmp[4] = {0, 0, 0, 0};
  uint64_t errBest;
  int32_t  tmpBest;
  uint8_t  t;

  for (t = 0; t < 8; t++) // get energy errors due to temporal merging
  {
    const int64_t ref = enValue[t] << 3;

    errTmp[0] += abs ((int64_t) envTmp0[t >> 3] - ref); // abs() since
    errTmp[1] += abs ((int64_t) envTmp1[t >> 2] - ref); // both values
    errTmp[2] += abs ((int64_t) envTmp2[t >> 1] - ref); // are already
    errTmp[3] += abs ((int64_t) envTmp3[t >> 0] - ref); // squares
  }
  errBest = errTmp[0];
  tmpBest = 0;

  for (t = 1; t < 3; t++) // find tmp value for minimal weighted error
  {
    if ((errTmp[t] << t) < errBest)
    {
      errBest = errTmp[t] << t;
      tmpBest = t;
    }
  }
  if ((errBest >> 3) > envTmp0[0]) tmpBest = (lr ? 2 : 3);

  if (tmpBest < tmpValSte) tmpBest = tmpValSte;

  /*Q*/if (tmpBest == 0)  // quantized envelopes for optimal tmp value
  {
    sbrData[0] = quantizeSbrEnvelopeLevel (envTmp0[0], frameSize, specFlat5b);
  }
  else if (tmpBest == 1)
  {
    sbrData[0] = quantizeSbrEnvelopeLevel (envTmp1[0], frameSize, specFlat5b);
    sbrData[1] = quantizeSbrEnvelopeLevel (envTmp1[1], frameSize, specFlat5b);
  }
  else if (tmpBest == 2)
  {
    sbrData[0] = quantizeSbrEnvelopeLevel (envTmp2[0], frameSize, specFlat5b);
    sbrData[1] = quantizeSbrEnvelopeLevel (envTmp2[1], frameSize, specFlat5b);
    sbrData[2] = quantizeSbrEnvelopeLevel (envTmp2[2], frameSize, specFlat5b);
    sbrData[3] = quantizeSbrEnvelopeLevel (envTmp2[3], frameSize, specFlat5b);
  }
  else // (tmpBest == 3)
  {
    sbrData[0] = quantizeSbrEnvelopeLevel (envTmp3[0], frameSize, specFlat5b);
    sbrData[1] = quantizeSbrEnvelopeLevel (envTmp3[1], frameSize, specFlat5b);
    sbrData[2] = quantizeSbrEnvelopeLevel (envTmp3[2], frameSize, specFlat5b);
    sbrData[3] = quantizeSbrEnvelopeLevel (envTmp3[3], frameSize, specFlat5b);
    sbrData[4] = quantizeSbrEnvelopeLevel (envTmp3[4], frameSize, specFlat5b);
    sbrData[5] = quantizeSbrEnvelopeLevel (envTmp3[5], frameSize, specFlat5b);
    sbrData[6] = quantizeSbrEnvelopeLevel (envTmp3[6], frameSize, specFlat5b);
    sbrData[7] = quantizeSbrEnvelopeLevel (envTmp3[7], frameSize, specFlat5b);
  }

  // quantized noise level for up to two temporal units, 30 = no noise
  t = __min (30, __max (specFlat5b, specFlatSte));
  sbrData[8] = ((int32_t) t << 13) | ((int32_t) t << 26);
#if ENABLE_INTERTES
  if ((t < 12) && (tempFlat5b > (lr ? 15 : 26)) && (tmpBest < 3))
  {
    sbrData[8] |= (1 << (1 << tmpBest)) - 1; // TODO: inter-TES mode = inv. filter mode?
  }
#endif
  memcpy (&sbrLevels[20], &sbrLevels[10] /*last*/, 10 * sizeof (int32_t)); // update the
  memcpy (&sbrLevels[10], sbrLevels /*& current*/, 10 * sizeof (int32_t)); // delay line

  tmpBest <<= 21; // config bits
  for (t = 0; t < (1 << (tmpBest >> 21)); t++)
  {
    const int32_t sbrEnvel = sbrData[t] & 127;
    const int32_t d = sbrLevels[30] - sbrEnvel; // two length-2 words!

    if ((t > 0 || !ind) && (d == 0 || d == 1))
    {
      tmpBest |= 1 << (12 + t); // delta-time coding flag for envelope
      sbrData[t] -= sbrEnvel;
      sbrData[t] |= d | (d << 7) | (d << 9) | (d << 11) | (d << 13) | (d << 15);
    }
    sbrLevels[30] = sbrEnvel;
  }
  for (t = 0; t < ((tmpBest >> 21) == 0 ? 1 : 2); t++)
  {
    const int32_t sbrNoise = (sbrData[8] >> (13 * (t + 1))) & 31;

    if ((t > 0 || !ind) && (sbrNoise == sbrLevels[31]))
    {
      tmpBest |= 1 << (4 + t); // and delta-time coding flag for noise
      sbrData[8] -= sbrNoise << (13 * (t + 1));
    }
    sbrLevels[31] = sbrNoise;
  }

  return tmpBest;
}

// public sampling rate functions
int8_t toSamplingFrequencyIndex (const unsigned samplingRate)
{
  for (int8_t i = 0; i < AAC_NUM_SAMPLE_RATES; i++)
  {
    if (samplingRate == allowedSamplingRates[i])  // (HE-)AAC rate
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

unsigned toSamplingRate (const int8_t samplingFrequencyIndex)
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
