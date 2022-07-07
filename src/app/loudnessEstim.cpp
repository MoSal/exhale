/* loudnessEstim.cpp - source file for class with ITU-R BS.1770-4 loudness level estimation
 * written by C. R. Helmrich, last modified in 2020 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under the exhale Copyright License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2021 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#include "exhaleAppPch.h"
#include "loudnessEstim.h"

#if LE_ACCURATE_CALC
static const int64_t kFilterCoeffs[4][8] = { // first 4: numerator (16->-32 bit), last 4: denominator (32 bit)
  { -974848000, 1329463296, -808124416, 185073664,  -946145519, 1253229559, -741406502, 165888314}, // <=32k
  {-1001717760, 1403125760, -874840064, 204996608,  -980574000, 1345222905, -822320842, 189236866}, // 44.1k
  {-1007550464, 1419378688, -889847808, 209584128,  -988032194, 1365543311, -840618073, 194671779}, // 48.0k
  {-1024000000, 1465647104, -933036032, 222953472, -1009281050, 1424244875, -894338686, 210939497}  // >=64k
};
#endif

// constructor
LoudnessEstimator::LoudnessEstimator (int32_t* const inputPcmData,         const unsigned bitDepth /*= 24*/,
                                      const unsigned sampleRate /*= 44k*/, const unsigned numChannels /*= 2*/)
{
#if LE_ACCURATE_CALC
  m_filterCoeffs  = kFilterCoeffs[sampleRate <= 44100 ? (sampleRate <= 32000 ? 0 : 1) : (sampleRate <= 48000 ? 2 : 3)];
  m_filterFactor  = (sampleRate < 48000 ? (48200 - sampleRate) >> 12 : (48000 - sampleRate) >> 14);
#else
  m_filterFactor  = 224 + (__min (SHRT_MAX, (int) sampleRate - 47616) >> 10);
#endif
  m_gbHopSize64   = (__min (163519, sampleRate) + 320) / 640; // 100 msec
  m_gbNormFactor  = (m_gbHopSize64 == 0 ? 0 : 1.0f / (4.0f * m_gbHopSize64));
  m_inputChannels = __min (8, numChannels);
  m_inputMaxValue = 1 << (__min (24, bitDepth) - 1);
  m_inputPcmData  = inputPcmData;

  reset ();
#if LE_ACCURATE_CALC
  for (unsigned ch = 0; ch < 8; ch++)
  {
    memset (m_filterMemI[ch], 0, 4 * sizeof (int32_t));
    memset (m_filterMemO[ch], 0, 4 * sizeof (int32_t));
  }
#else
  for (unsigned ch = 0; ch < 8; ch++) m_filterMemoryI[ch] = m_filterMemoryO[ch] = 0;
#endif
}

// public functions
uint32_t LoudnessEstimator::addNewPcmData (const unsigned samplesPerChannel)
{
  const unsigned frameSize64  = samplesPerChannel >> 6; // in units of 64
  const unsigned numSamples64 = 1 << 6; // sub-frame size (64, of course)
  const int32_t* chSig        = m_inputPcmData;
#if LE_ACCURATE_CALC
  const int64_t* filtI        = m_filterCoeffs;
  const int64_t* filtO        = &m_filterCoeffs[4];
#endif
  uint64_t* newQuarterPower   = m_powerValue[3];
  unsigned  ch, f, s;

  if ((chSig == nullptr) || (frameSize64 == 0))
  {
    return 1; // invalid sample pointer or frame size
  }

  // de-interleave and K-filter incoming audio samples in sub-frame units
  for (f = 0; f < frameSize64; f++) // sub-frame loop
  {
    for (s = 0; s < numSamples64; s++) // sample loop
    {
      for (ch = 0; ch < m_inputChannels; ch++)
      {
#if LE_ACCURATE_CALC
        // accurate K-filter according to ITU-R BS.1770-4, Annex 1 (2015)
        int32_t* const i = m_filterMemI[ch];
        int32_t* const o = m_filterMemO[ch];
        const int64_t pi = filtI[0] * i[0] + filtI[1] * i[1] + filtI[2] * i[2] + filtI[3] * i[3] -
                           filtO[0] * o[0] - filtO[1] * o[1] - filtO[2] * o[2] - filtO[3] * o[3];
        const int64_t to = (pi < 0 ? (1 << 28) - 1 : 0); // trunc. offset
        const int32_t xi = (*(chSig++)) * (1 << 2);
        const int32_t yi = xi + int32_t ((pi + (xi == 0 ? to : (1 << 27))) >> 28);
        const uint32_t a = abs (xi >> 2);

        i[3] = i[2];   i[2] = i[1];   i[1] = i[0];   i[0] = xi; // update
        o[3] = o[2];   o[2] = o[1];   o[1] = o[0];   o[0] = yi; // memory
#else
        // simplified K-filter, including 500-Hz high-pass pre-processing
        const int32_t xi = *(chSig++);
        const int32_t yi = xi - m_filterMemoryI[ch] + ((128 + m_filterFactor * m_filterMemoryO[ch]) >> 8);
        const uint32_t a = abs (xi);

        m_filterMemoryI[ch] = xi;
        m_filterMemoryO[ch] = yi;
#endif
        newQuarterPower[ch] += (int64_t) yi * (int64_t) yi;

        if (m_inputPeakValue < a) m_inputPeakValue = a; // get peak level
      }
    } // s

    if (++m_gbHopLength64 >= m_gbHopSize64) // completed 100-msec quarter
    {
      const float thrA = LE_THRESH_ABS * (float) m_inputMaxValue * (float) m_inputMaxValue;
      uint64_t zij, zj = 0;

      for (ch = 0; ch < m_inputChannels; ch++)  // sum 64-sample averages
      {
        zij = (m_powerValue[0][ch] + m_powerValue[1][ch] + m_powerValue[2][ch] + newQuarterPower[ch] + (1u << 5)) >> 6;
        zj += (ch > 2 ? (16u + 45 * zij) >> 5 : zij); // weighting by G_i
      }

      if (zj * m_gbNormFactor > thrA) // use sqrt (block RMS) if lj > -70
      {
        if (m_gbRmsValues.size () < INT_MAX) m_gbRmsValues.push_back (uint32_t (sqrt (zj * m_gbNormFactor) + 0.5f));
      }

      for (ch = 0; ch < m_inputChannels; ch++) // set up new gating block
      {
        m_powerValue[0][ch] = m_powerValue[1][ch];
        m_powerValue[1][ch] = m_powerValue[2][ch];
        m_powerValue[2][ch] = newQuarterPower[ch];
        newQuarterPower[ch] = 0;
      }
      m_gbHopLength64 = 0;
    }
  }

  return 0; // no error
}

uint32_t LoudnessEstimator::getStatistics (const bool includeWarmUp /*= false*/)
{
  const uint32_t numVectorValues = (uint32_t) m_gbRmsValues.size ();
  const uint32_t numWarmUpBlocks = includeWarmUp ? 0 : 3;
  const uint32_t numGatingBlocks = __max (numWarmUpBlocks, numVectorValues) - numWarmUpBlocks;
  const uint16_t maxValueDivisor = __max (1u, m_inputMaxValue >> 16);
  const uint16_t peakValue16Bits = __min (USHRT_MAX, (m_inputPeakValue + (maxValueDivisor >> 1)) / maxValueDivisor);
  uint32_t i, numBlocks = 0;
  float thrR, zg;

  if (numGatingBlocks == 0) return peakValue16Bits;  // no loudness stats

  const float normFac = 1.0f / numGatingBlocks; // prevents loop overflow

  // calculate arithmetic average of blocks satisfying absolute threshold
  for (zg = 0.0f, i = numWarmUpBlocks; i < numVectorValues; i++)
  {
    zg += normFac * (float) m_gbRmsValues.at (i) * (float) m_gbRmsValues.at (i);
  }
  if (zg < LE_THRESH_ABS) return peakValue16Bits;

  thrR = LE_THRESH_REL * zg; // find blocks satisfying relative threshold
  for (zg = 0.0f, i = numWarmUpBlocks; i < numVectorValues; i++)
  {
    const float p = (float) m_gbRmsValues.at (i) * (float) m_gbRmsValues.at (i);

    if (p > thrR) { zg += normFac * p;  numBlocks++; }
  }
  if (zg < LE_THRESH_ABS) return peakValue16Bits;

  zg = LE_LUFS_OFFSET + 10.0f * (float) log10 (zg / (normFac * numBlocks * (float) m_inputMaxValue * (float) m_inputMaxValue));
#if LE_ACCURATE_CALC
  zg -= m_filterFactor * 0.046875f; // for sample rates other than 48 kHz
#endif
  i  = __max (0, int32_t ((zg + 100.0f) * 512.0f + 0.5f)); // map to uint

  return (__min (USHRT_MAX, i) << 16) | peakValue16Bits; // L = i/512-100
}
