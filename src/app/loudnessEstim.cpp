/* loudnessEstim.cpp - source file for class with ITU-R BS.1770-4 loudness level estimation
 * written by C. R. Helmrich, last modified in 2020 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under a Modified BSD-Style License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2020 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#include "exhaleAppPch.h"
#include "loudnessEstim.h"

// constructor
LoudnessEstimator::LoudnessEstimator (int32_t* const inputPcmData,           const unsigned bitDepth /*= 24*/,
                                      const unsigned sampleRate /*= 44100*/, const unsigned numChannels /*= 2*/)
{
  m_filterFactor  = 224 + (__min (SHRT_MAX, (int) sampleRate - 47616) >> 10);
  m_gbHopSize64   = (__min (163519, sampleRate) + 320) / 640; // 100 msec
  m_gbNormFactor  = (m_gbHopSize64 == 0 ? 0 : 1.0f / (4.0f * m_gbHopSize64));
  m_inputChannels = __min (8, numChannels);
  m_inputMaxValue = 1 << (__min (24, bitDepth) - 1);
  m_inputPcmData  = inputPcmData;

  reset ();
  for (unsigned ch = 0; ch < 8; ch++) m_filterMemoryI[ch] = m_filterMemoryO[ch] = 0;
}

// public functions
uint32_t LoudnessEstimator::addNewPcmData (const unsigned samplesPerChannel)
{
  const unsigned frameSize64  = samplesPerChannel >> 6; // in units of 64
  const unsigned numSamples64 = 1 << 6; // sub-frame size (64, of course)
  const int32_t* chSig        = m_inputPcmData;
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
        // simplified K-filter, including 500-Hz high-pass pre-processing
        const int32_t xi = *(chSig++);
        const int32_t yi = xi - m_filterMemoryI[ch] + ((128 + m_filterFactor * m_filterMemoryO[ch]) >> 8);
        const uint32_t a = abs (xi);

        m_filterMemoryI[ch] = xi;
        m_filterMemoryO[ch] = yi;
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

  zg = LE_LUFS_OFFSET + 10.0f * log10 (zg / (normFac * numBlocks * (float) m_inputMaxValue * (float) m_inputMaxValue));
  i  = __max (0, int32_t ((zg + 100.0f) * 512.0f + 0.5f)); // map to uint

  return (__min (USHRT_MAX, i) << 16) | peakValue16Bits; // L = i/512-100
}
