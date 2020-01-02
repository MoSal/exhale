/* specAnalysis.cpp - source file for class providing spectral analysis of MCLT signals
 * written by C. R. Helmrich, last modified in 2019 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under a Modified BSD-Style License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2020 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#include "exhaleLibPch.h"
#include "specAnalysis.h"

// static helper function
static inline uint32_t packAvgSpecAnalysisStats (const uint64_t sumAvgBand, const uint64_t sumMaxBand,
                                                 const unsigned char predGain,
                                                 const uint16_t idxMaxSpec, const uint16_t idxLpStart)
{
  // temporal flatness, normalized for a value of 256 for a linear prediction gain of 1 (0 dB)
  const unsigned flatTemp = predGain;
  // spectral flatness, normalized for a value of 256 for steady low or mid-frequency sinusoid
  const int32_t  flatSpec = 256 - int (((sumAvgBand + SA_EPS) * 402) / (sumMaxBand + SA_EPS));

  return (flatTemp << 24) | (CLIP_UCHAR (flatSpec) << 16) | (__min (2047, idxMaxSpec) << 5) | __min (31, idxLpStart);
}

// constructor
SpecAnalyzer::SpecAnalyzer ()
{
  for (unsigned ch = 0; ch < USAC_MAX_NUM_CHANNELS; ch++)
  {
    m_bandwidthOff[ch] = 0;
    m_numAnaBands [ch] = 0;
    m_specAnaStats[ch] = 0;
    memset (m_parCorCoeffs[ch], 0, MAX_PREDICTION_ORDER * sizeof (short));
  }
  m_tnsPredictor = nullptr;
}

// public functions
unsigned SpecAnalyzer::getLinPredCoeffs (short parCorCoeffs[MAX_PREDICTION_ORDER], const unsigned channelIndex)  // returns best filter order
{
  unsigned bestOrder = MAX_PREDICTION_ORDER, predGainCurr, predGainPrev;

  if ((parCorCoeffs == nullptr) || (channelIndex >= USAC_MAX_NUM_CHANNELS))
  {
    return 0; // invalid arguments error
  }
  memcpy (parCorCoeffs, m_parCorCoeffs[channelIndex], MAX_PREDICTION_ORDER * sizeof (short));

  predGainCurr = (m_tnsPredGains[channelIndex] >> 24) & UCHAR_MAX;
  predGainPrev = (m_tnsPredGains[channelIndex] >> 16) & UCHAR_MAX;
  while ((bestOrder > 1) && (predGainPrev >= predGainCurr))  // find lowest-order gain maximum
  {
    bestOrder--;
    predGainCurr = predGainPrev;
    predGainPrev = (m_tnsPredGains[channelIndex] >> (8 * bestOrder - 16)) & UCHAR_MAX;
  }
  return ((bestOrder == 1) && (m_parCorCoeffs[channelIndex][0] == 0) ? 0 : bestOrder);
}

unsigned SpecAnalyzer::getMeanAbsValues (const int32_t* const mdctSignal, const int32_t* const mdstSignal, const unsigned nSamplesInFrame,
                                         const unsigned channelIndex, const uint16_t* const bandStartOffsets, const unsigned nBands,
                                         uint32_t* const meanBandValues)
{
  if ((mdctSignal == nullptr) || (bandStartOffsets == nullptr) || (meanBandValues == nullptr) || (channelIndex >= USAC_MAX_NUM_CHANNELS) ||
     (nSamplesInFrame > 2048) || (nSamplesInFrame < 2) || (nBands > nSamplesInFrame))
  {
    return 1; // invalid arguments error
  }

  if (mdstSignal != nullptr) // use complex-valued spectral data
  {
    for (unsigned b = 0; b < nBands; b++)
    {
      const unsigned bandOffset  = __min (nSamplesInFrame, bandStartOffsets[b]);
      const unsigned bandWidth   = __min (nSamplesInFrame, bandStartOffsets[b + 1]) - bandOffset;
      const unsigned anaBandIdx  = bandOffset >> SA_BW_SHIFT;

      if ((anaBandIdx < m_numAnaBands[channelIndex]) && (bandOffset == (anaBandIdx << SA_BW_SHIFT)) && ((bandWidth & (SA_BW - 1)) == 0))
      {
        const uint32_t* const anaAbsVal = &m_meanAbsValue[channelIndex][anaBandIdx];

        // data available from previous call to spectralAnalysis
        meanBandValues[b] = (bandWidth == SA_BW ? *anaAbsVal : uint32_t (((int64_t) anaAbsVal[0] + (int64_t) anaAbsVal[1] + 1) >> 1));
      }
      else // no previous data available, compute mean magnitude
      {
        const int32_t* const bMdct = &mdctSignal[bandOffset];
        const int32_t* const bMdst = &mdstSignal[bandOffset];
        uint64_t sumAbsVal = 0;

        for (int s = bandWidth - 1; s >= 0; s--)
        {
#if SA_EXACT_COMPLEX_ABS
          const double  complexSqr = (double) bMdct[s] * (double) bMdct[s] + (double) bMdst[s] * (double) bMdst[s];
          const unsigned absSample = unsigned (sqrt (complexSqr) + 0.5);
#else
          const unsigned absReal   = abs (bMdct[s]); // Richard Lyons, 1997; en.wikipedia.org/
          const unsigned absImag   = abs (bMdst[s]); // wiki/Alpha_max_plus_beta_min_algorithm
          const unsigned absSample = (absReal > absImag ? absReal + ((absImag * 3) >> 3) : absImag + ((absReal * 3) >> 3));
#endif
          sumAbsVal += absSample;
        }
        // average spectral sample magnitude across current band
        meanBandValues[b] = uint32_t ((sumAbsVal + (bandWidth >> 1)) / bandWidth);
      }
    } // for b
  }
  else // no imaginary part available, real-valued spectral data
  {
    for (unsigned b = 0; b < nBands; b++)
    {
      const unsigned bandOffset  = __min (nSamplesInFrame, bandStartOffsets[b]);
      const unsigned bandWidth   = __min (nSamplesInFrame, bandStartOffsets[b + 1]) - bandOffset;
      const int32_t* const bMdct = &mdctSignal[bandOffset];
      uint64_t sumAbsVal = 0;

      for (int s = bandWidth - 1; s >= 0; s--)
      {
        const unsigned absSample = abs (bMdct[s]);

        sumAbsVal += absSample;
      }
      // average spectral sample magnitude across frequency band
      meanBandValues[b] = uint32_t ((sumAbsVal + (bandWidth >> 1)) / bandWidth);
    } // for b
  }
  m_numAnaBands[channelIndex] = 0; // mark spectral data as used

  return 0; // no error
}

void SpecAnalyzer::getSpecAnalysisStats (uint32_t avgSpecAnaStats[USAC_MAX_NUM_CHANNELS], const unsigned nChannels)
{
  if ((avgSpecAnaStats == nullptr) || (nChannels > USAC_MAX_NUM_CHANNELS))
  {
    return;
  }
  memcpy (avgSpecAnaStats, m_specAnaStats, nChannels * sizeof (uint32_t));
}

void SpecAnalyzer::getSpectralBandwidth (uint16_t bandwidthOffset[USAC_MAX_NUM_CHANNELS], const unsigned nChannels)
{
  if ((bandwidthOffset == nullptr) || (nChannels > USAC_MAX_NUM_CHANNELS))
  {
    return;
  }
  memcpy (bandwidthOffset, m_bandwidthOff, nChannels * sizeof (uint16_t));
}

unsigned SpecAnalyzer::initLinPredictor (LinearPredictor* const linPredictor)
{
  if (linPredictor == nullptr)
  {
    return 1; // invalid arguments error
  }
  m_tnsPredictor = linPredictor;

  return 0; // no error
}

#if SA_OPT_WINDOW_GROUPING
unsigned SpecAnalyzer::optimizeGrouping (const unsigned channelIndex, const unsigned prefBandwidth, const unsigned prefGroupingIndex)
{
  const uint32_t* meanAbsValCurr = m_meanAbsValue[channelIndex];
  const uint32_t numAnaBandsInCh = m_numAnaBands [channelIndex];
  unsigned grpIdxCurr = prefGroupingIndex, maxBands, numBands;
  uint64_t energyCurrHF, energyPrefHF;
  uint32_t energyCurrLF, energyPrefLF;
  unsigned b;

  if ((prefBandwidth > 2048) || (grpIdxCurr == 0) || (grpIdxCurr >= 8) || (channelIndex >= USAC_MAX_NUM_CHANNELS) || (numAnaBandsInCh == 0))
  {
    return 8; // invalid arguments error, or pypassing
  }

  numBands = numAnaBandsInCh >> 3;
  maxBands = numAnaBandsInCh << SA_BW_SHIFT;  // available bandwidth, equal to nSamplesInFrame
  maxBands = (numBands * __min (maxBands, prefBandwidth) + (maxBands >> 1)) / maxBands;

  if (maxBands * numBands == 0) return 8; // low/no BW
  if (grpIdxCurr < 7) grpIdxCurr++; // after transient

  meanAbsValCurr += grpIdxCurr * numBands;
  grpIdxCurr++;
  energyPrefLF = meanAbsValCurr[0] >> 1; // - 6 dB
  energyPrefHF = 0;
  for (b = maxBands - 1; b > 0; b--)  // avoid LF band
  {
    energyPrefHF += meanAbsValCurr[b];
  }
  energyPrefHF >>= 1; // - 6 dB

  do // check whether HF or LF transient starts earlier than preferred grouping index suggests
  {
    meanAbsValCurr -= numBands;
    grpIdxCurr--;
    energyCurrLF = meanAbsValCurr[0];
    energyCurrHF = 0;
    for (b = maxBands - 1; b > 0; b--) // prev. window
    {
      energyCurrHF += meanAbsValCurr[b];
    }
  }
  while ((grpIdxCurr > 1) && (energyCurrHF >= energyPrefHF) && (energyCurrLF >= energyPrefLF));

  return __min (grpIdxCurr, prefGroupingIndex); // final optimized grouping index
}
#endif // SA_OPT_WINDOW_GROUPING

unsigned SpecAnalyzer::spectralAnalysis (const int32_t* const mdctSignals[USAC_MAX_NUM_CHANNELS],
                                         const int32_t* const mdstSignals[USAC_MAX_NUM_CHANNELS],
                                         const unsigned nChannels, const unsigned nSamplesInFrame, const unsigned samplingRate,
                                         const unsigned lfeChannelIndex /*= USAC_MAX_NUM_CHANNELS*/) // to skip an LFE channel
{
  const unsigned lpcStopBand16k = (samplingRate <= 32000 ? nSamplesInFrame : (32000 * nSamplesInFrame) / samplingRate) >> SA_BW_SHIFT;
  const unsigned thresholdSlope = (48000 + SA_EPS * samplingRate) / 96000;
  const unsigned thresholdStart = samplingRate >> 15;

  if ((mdctSignals == nullptr) || (nChannels > USAC_MAX_NUM_CHANNELS) || (lfeChannelIndex > USAC_MAX_NUM_CHANNELS) ||
      (nSamplesInFrame > 2048) || (nSamplesInFrame < 2) || (samplingRate < 7350) || (samplingRate > 96000))
  {
    return 1; // invalid arguments error
  }

  for (unsigned ch = 0; ch < nChannels; ch++)
  {
    const int32_t* const chMdct = mdctSignals[ch];
    const int32_t* const chMdst = (mdstSignals == nullptr ? nullptr : mdstSignals[ch]);
// --- get L1 norm and max value in each band
    uint16_t idxMaxSpec = 0;
    uint64_t sumAvgBand = 0;
    uint64_t sumMaxBand = 0;
    uint32_t valMaxSpec = 0;
    int b;

    if (ch == lfeChannelIndex) // no analysis
    {
      m_bandwidthOff[ch] = LFE_MAX;
      m_numAnaBands [ch] = 0;
      m_specAnaStats[ch] = 0; // flat/stationary frame
      continue;
    }

    m_bandwidthOff[ch] = 0;
    m_numAnaBands [ch] = nSamplesInFrame >> SA_BW_SHIFT;

    for (b = m_numAnaBands[ch] - 1; b >= 0; b--)
    {
      const uint16_t         offs = b << SA_BW_SHIFT; // start offset of current analysis band
      const int32_t* const  bMdct = &chMdct[offs];
      const int32_t* const  bMdst = (chMdst == nullptr ? nullptr : &chMdst[offs]);
      uint16_t  maxAbsIdx = 0;
      uint32_t  maxAbsVal = 0, tmp = UINT_MAX;
      uint64_t  sumAbsVal = 0;

      if (bMdst != nullptr) // complex-valued spectrum
      {
        for (int s = SA_BW - 1; s >= 0; s--)
        {
          // sum absolute values of complex signal, derive L1 norm, peak value, and peak index
#if SA_EXACT_COMPLEX_ABS
          const double  complexSqr = (double) bMdct[s] * (double) bMdct[s] + (double) bMdst[s] * (double) bMdst[s];
          const unsigned absSample = unsigned (sqrt (complexSqr) + 0.5);
#else
          const unsigned absReal   = abs (bMdct[s]); // Richard Lyons, 1997; en.wikipedia.org/
          const unsigned absImag   = abs (bMdst[s]); // wiki/Alpha_max_plus_beta_min_algorithm
          const unsigned absSample = (absReal > absImag ? absReal + ((absImag * 3) >> 3) : absImag + ((absReal * 3) >> 3));
#endif
          sumAbsVal += absSample;
          if (offs + s > 0) // exclude DC from max/min
          {
            if (maxAbsVal < absSample) // maximum data
            {
              maxAbsVal = absSample;
              maxAbsIdx = (uint16_t) s;
            }
            if (tmp/*min*/> absSample) // minimum data
            {
              tmp/*min*/= absSample;
            }
          } // b > 0
        }
      }
      else  // real-valued spectrum, no imaginary part
      {
        for (int s = SA_BW - 1; s >= 0; s--)
        {
          // obtain absolute values of real signal, derive L1 norm, peak value, and peak index
          const unsigned absSample = abs (bMdct[s]);

          sumAbsVal += absSample;
          if (offs + s > 0) // exclude DC from max/min
          {
            if (maxAbsVal < absSample) // maximum data
            {
              maxAbsVal = absSample;
              maxAbsIdx = (uint16_t) s;
            }
            if (tmp/*min*/> absSample) // minimum data
            {
              tmp/*min*/= absSample;
            }
          }
        }
      }
      // bandwidth detection
      if ((m_bandwidthOff[ch] == 0) && (maxAbsVal > __max (thresholdSlope * (thresholdStart + b), SA_EPS)))
      {
        m_bandwidthOff[ch] = __max (maxAbsIdx + 5/*guard*/, SA_BW) + offs;
        m_bandwidthOff[ch] = __min (m_bandwidthOff[ch], nSamplesInFrame);
      }
      // save mean magnitude
      tmp/*mean*/ = uint32_t ((sumAbsVal + (1 << (SA_BW_SHIFT - 1))) >> SA_BW_SHIFT);
      m_meanAbsValue[ch][b] = tmp;
      // spectral statistics
      if (b > 0)
      {
        sumAvgBand += tmp;
        sumMaxBand += maxAbsVal;
      }
      if (valMaxSpec < maxAbsVal)
      {
        valMaxSpec = maxAbsVal;
        idxMaxSpec = maxAbsIdx + offs;
      }
    } // for b

// --- spectral analysis statistics for frame
    b = 1;
    while (((unsigned) b + 1 < lpcStopBand16k) && ((uint64_t) m_meanAbsValue[ch][b] * (m_numAnaBands[ch] - 1) > sumAvgBand)) b++;
    b = __min (m_bandwidthOff[ch], b << SA_BW_SHIFT);

    // obtain prediction gain across spectrum
    m_tnsPredGains[ch] = m_tnsPredictor->calcParCorCoeffs (&chMdct[b], __min (m_bandwidthOff[ch], lpcStopBand16k << SA_BW_SHIFT) - b,
                                                           MAX_PREDICTION_ORDER, m_parCorCoeffs[ch]);
    m_specAnaStats[ch] = packAvgSpecAnalysisStats (sumAvgBand, sumMaxBand, m_tnsPredGains[ch] >> 24, idxMaxSpec, (unsigned) b >> SA_BW_SHIFT);
  } // for ch

  return 0; // no error
}
