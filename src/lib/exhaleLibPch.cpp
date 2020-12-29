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

// ISO/IEC 14496-3, Annex 4.A.6.1
static const uint8_t deltaHuffSbrF[14] = {0xFB, 0x7C, 0x3C, 0x1C, 0x0C, 0x05, 0x01, 0x00, 0x04, 0x0D, 0x1D, 0x3D, 0xFA, 0xFC};
static const uint8_t deltaHuffSbrT[14] = {0xFD, 0x7D, 0x3D, 0x1D, 0x0D, 0x05, 0x01, 0x00, 0x04, 0x0C, 0x1C, 0x3C, 0x7C, 0xFC};

// static SBR related functions
static int32_t getSbrDeltaBitCount (const int32_t delta, const bool dt)
{
  return (delta == 5 && !dt ? 8 : abs (delta) + 2 + (delta >> 31));
}

static int32_t getSbrDeltaHuffCode (const int32_t delta, const bool dt)
{
  const unsigned u = __max (-7, __min (6, delta)) + 7;

  return int32_t (dt ? deltaHuffSbrT[u] : deltaHuffSbrF[u]);
}

static int8_t getSbrQuantizedLevel (const double energy, const uint32_t divisor, const uint8_t noiseLevel)
{
  const double ener = energy / divisor;

  return (ener > 8192.0 ? int8_t (1.375 - 0.03125 * noiseLevel + 6.64385619 * log10 (ener)) - 26 : 0);
}

static int32_t packSbr3BandQuantLevels (const uint64_t energy, const uint32_t divisor, const uint8_t noiseLevel,
                                        const uint32_t ratioL, const uint32_t ratioM, const uint32_t ratioH)
{
  int8_t val[4] = {getSbrQuantizedLevel (ratioL * (double) energy, divisor * 6144, noiseLevel),
                   getSbrQuantizedLevel (ratioM * (double) energy, divisor * 6144, noiseLevel),
                   getSbrQuantizedLevel (ratioH * (double) energy, divisor * 6144, noiseLevel), 0};
  unsigned iMax = 0;

  val[3] = val[0]; // find last peak band having maximum value
  if (val[3] <= val[1]) { iMax = 1; val[3] = val[1]; }
  if (val[3] <= val[2]) { iMax = 2; }

  if (iMax > 0) // limit delta-value increases below peak band
  {
    if (val[iMax - 1] + 6 < val[iMax]) val[iMax - 1] = val[iMax] - 6;
    if ((iMax == 2) && (val[0] + 6 < val[1])) val[0] = val[1] - 6;
  }
  if (iMax < 2) // limit delta-value decreases above peak band
  {
    if (val[iMax + 1] + 7 < val[iMax]) val[iMax + 1] = val[iMax] - 7;
    if ((iMax == 0) && (val[2] + 7 < val[1])) val[2] = val[1] - 7;
  }
  val[1] = __max (val[1], __min (val[0], val[2]));

  return int32_t (val[0]) | (int32_t (val[1]) << 8) | (int32_t (val[2]) << 16);
}

static bool useSbr3BandDeltaTimeCoding (const int32_t dataCurr, const int32_t dataPrev)
{
  const int8_t curr[3] = {dataCurr & SCHAR_MAX, (dataCurr >> 8) & SCHAR_MAX, (dataCurr >> 16) & SCHAR_MAX};
  const int8_t prev[3] = {dataPrev & SCHAR_MAX, (dataPrev >> 8) & SCHAR_MAX, (dataPrev >> 16) & SCHAR_MAX};
  int32_t bDf = 7, bDt = 0;

  for (unsigned u = 0; u < 3; u++) // check delta range, count
  {
    if ((curr[u] > prev[u] + 6) || (curr[u] < prev[u] - 7)) return false;
    if (u > 0) bDf += getSbrDeltaBitCount (curr[u] - curr[u - 1], false);
    /*deltaT*/ bDt += getSbrDeltaBitCount (curr[u] - prev[u], true);
  }

  return (bDf > bDt);
}

// public SBR related functions
int32_t getSbrEnvelopeAndNoise (int32_t* const sbrLevels, const uint8_t specFlat5b, const uint8_t tempFlat5b, const bool lr, const bool ind,
                                const uint8_t specFlatSte, const int32_t tmpValSte, const uint32_t frameSize, int32_t* sbrData)
{
  const uint64_t enValue[8] = {square (sbrLevels[21]), square (sbrLevels[22]), square (sbrLevels[23]), square (sbrLevels[24]),
                               square (sbrLevels[25]), square (sbrLevels[26]), square (sbrLevels[27]), square (sbrLevels[10])};
  const uint64_t envTmp0[1] = { enValue[0] + enValue[1] + enValue[2] + enValue[3] +
                                enValue[4] + enValue[5] + enValue[6] + enValue[7]};
  const uint64_t envTmp1[2] = {(enValue[0] + enValue[1] + enValue[2] + enValue[3]) << 1,
                               (enValue[4] + enValue[5] + enValue[6] + enValue[7]) << 1};
  const uint64_t envTmp2[4] = {(enValue[0] + enValue[1]) << 2, (enValue[2] + enValue[3]) << 2,
                               (enValue[4] + enValue[5]) << 2, (enValue[6] + enValue[7]) << 2};
  const uint64_t envTmp3[8] = { enValue[0] << 3, enValue[1] << 3, enValue[2] << 3, enValue[3] << 3,
                                enValue[4] << 3, enValue[5] << 3, enValue[6] << 3, enValue[7] << 3};
  const uint32_t env3BandsL = sbrLevels[28] & USHRT_MAX;
  const uint32_t env3BandsM = sbrLevels[28] >> 16;
  const uint32_t env3BandsH = sbrLevels[29] & USHRT_MAX;
  uint64_t errTmp[4] = {0, 0, 0, 0};
  uint64_t errBest;
  int32_t  tmpBest, bitCount, diff;
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
    sbrData[0] = packSbr3BandQuantLevels (envTmp0[0], frameSize, specFlat5b, env3BandsL, env3BandsM, env3BandsH);
  }
  else if (tmpBest == 1)
  {
    for (t = 0; t < 2; t++)
    {
      sbrData[t] = packSbr3BandQuantLevels (envTmp1[t], frameSize, specFlat5b, env3BandsL, env3BandsM, env3BandsH);
    }
  }
  else if (tmpBest == 2)
  {
    for (t = 0; t < 4; t++)
    {
      sbrData[t] = packSbr3BandQuantLevels (envTmp2[t], frameSize, specFlat5b, env3BandsL, env3BandsM, env3BandsH);
    }
  }
  else // (tmpBest == 3)
  {
    for (t = 0; t < 8; t++)
    {
      sbrData[t] = packSbr3BandQuantLevels (envTmp3[t], frameSize, specFlat5b, env3BandsL, env3BandsM, env3BandsH);
    }
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
    const int32_t curr = sbrData[t], prev = sbrLevels[30];

    if ((t > 0 /*|| !ind*/) && useSbr3BandDeltaTimeCoding (curr, prev))
    {
      tmpBest |= 1 << (12 + t); // delta-time coding flag for envelope

      diff = (curr & SCHAR_MAX) - (prev & SCHAR_MAX);
      sbrData[t] = getSbrDeltaHuffCode (diff, true);  bitCount = 8; // see bitStreamWriter.cpp, ll. 186, 213
      diff = ((curr >> 8) & SCHAR_MAX) - ((prev >> 8) & SCHAR_MAX);
      sbrData[t] |= getSbrDeltaHuffCode (diff, true) << bitCount;  bitCount += getSbrDeltaBitCount (diff, true);
      diff = ((curr >> 16) & SCHAR_MAX) - ((prev >> 16) & SCHAR_MAX);
      sbrData[t] |= getSbrDeltaHuffCode (diff, true) << bitCount;  bitCount += getSbrDeltaBitCount (diff, true);
    }
    else // delta-frequency coding
    {
      diff = curr & SCHAR_MAX;
      sbrData[t] = diff;  bitCount = 8; // first envelope is PCM coded
      diff = ((curr >> 8) & SCHAR_MAX) - (curr & SCHAR_MAX);
      sbrData[t] |= getSbrDeltaHuffCode (diff, false) << bitCount;  bitCount += getSbrDeltaBitCount (diff, false);
      diff = ((curr >> 16) & SCHAR_MAX) - ((curr >> 8) & SCHAR_MAX);
      sbrData[t] |= getSbrDeltaHuffCode (diff, false) << bitCount;  bitCount += getSbrDeltaBitCount (diff, false);
    }
    sbrData[t] |= 1 << bitCount; // MSB delimiter for bitstream writer
    sbrLevels[30] = curr;
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
