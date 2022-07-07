/* specAnalysis.cpp - source file for class providing spectral analysis of MCLT signals
 * written by C. R. Helmrich, last modified in 2021 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under the exhale Copyright License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2021 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#include "exhaleLibPch.h"
#include "specAnalysis.h"

// static helper functions
static inline uint64_t complexAbs (const int32_t realPart, const int32_t imagPart)
{
#if SA_EXACT_COMPLEX_ABS
  const double complexSqr = (double) realPart * (double) realPart + (double) imagPart * (double) imagPart;

  return uint64_t (sqrt (complexSqr) + 0.5);
#else
  const uint64_t absReal  = abs (realPart); // Richard Lyons, 1997; en.wikipedia.org/
  const uint64_t absImag  = abs (imagPart); // wiki/Alpha_max_plus_beta_min_algorithm

  return (absReal > absImag ? absReal + ((absImag * 3) >> 3) : absImag + ((absReal * 3) >> 3));
#endif
}

static inline uint32_t packAvgSpecAnalysisStats (const uint64_t sumAvgBand, const uint64_t sumMaxBand,
                                                 const uint8_t  predGain,
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
    m_magnCorrPrev[ch] = 0;
    m_magnSpectra [ch] = nullptr;
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
  while ((predGainPrev >= predGainCurr) && --bestOrder > 1)  // find lowest-order gain maximum
  {
    predGainCurr = predGainPrev;
    predGainPrev = (m_tnsPredGains[channelIndex] >> (8 * bestOrder - 16)) & UCHAR_MAX;
  }
  return ((bestOrder == 1) && (m_parCorCoeffs[channelIndex][0] == 0) ? 0 : bestOrder);
}

unsigned SpecAnalyzer::getMeanAbsValues (const int32_t* const mdctSignal, const int32_t* const mdstSignal, const unsigned nSamplesInFrame,
                                         const unsigned channelIndex, const uint16_t* const bandStartOffsets, const unsigned nBands,
                                         uint32_t* const meanBandValues)
{
  if ((mdctSignal == nullptr) || (bandStartOffsets == nullptr) || (meanBandValues == nullptr) || (channelIndex > USAC_MAX_NUM_CHANNELS) ||
     (nSamplesInFrame > 2048) || (nSamplesInFrame < 2) || (nBands > nSamplesInFrame))
  {
    return 1; // invalid arguments error
  }

  if (mdstSignal != nullptr) // use complex-valued spectral data
  {
    for (unsigned b = 0; b < nBands; b++)
    {
      const unsigned bandOffset = __min (nSamplesInFrame, bandStartOffsets[b]);
      const unsigned bandWidth  = __min (nSamplesInFrame, bandStartOffsets[b + 1]) - bandOffset;
      const unsigned anaBandIdx = bandOffset >> SA_BW_SHIFT;

      if ((channelIndex < USAC_MAX_NUM_CHANNELS) && (anaBandIdx < m_numAnaBands[channelIndex]) &&
          (bandOffset == (anaBandIdx << SA_BW_SHIFT)) && ((bandWidth & (SA_BW - 1)) == 0))
      {
        const uint32_t* const anaAbsVal = &m_meanAbsValue[channelIndex][anaBandIdx];

        // data available from previous call to spectralAnalysis
        meanBandValues[b] = (bandWidth == SA_BW ? *anaAbsVal : uint32_t (((int64_t) anaAbsVal[0] + (int64_t) anaAbsVal[1] + 1) >> 1));
      }
      else // no previous data available, compute mean magnitude
      {
        const int32_t* bMdct = &mdctSignal[bandOffset];
        const int32_t* bMdst = &mdstSignal[bandOffset];
        uint64_t   sumAbsVal = 0;

        for (unsigned s = bandWidth; s > 0; s--) sumAbsVal += complexAbs (*(bMdct++), *(bMdst++));

        // average spectral sample magnitude across current band
        meanBandValues[b] = uint32_t ((sumAbsVal + (bandWidth >> 1)) / bandWidth);
      }
    } // for b
  }
  else // no imaginary part available, real-valued spectral data
  {
    int64_t prevMdct = mdctSignal[bandStartOffsets[0] + ((bandStartOffsets[0] > 0) && (channelIndex < USAC_MAX_NUM_CHANNELS) ? -1 : 1)];

    for (unsigned b = 0; b < nBands; b++)
    {
      const unsigned bandOffset = __min (nSamplesInFrame, bandStartOffsets[b]);
      const unsigned bandWidth  = __min (nSamplesInFrame, bandStartOffsets[b + 1]) - bandOffset;
      const int32_t* bMdct = &mdctSignal[bandOffset];
      const int32_t* bNext = &bMdct[1];
      uint64_t   sumAbsVal = (bandStartOffsets[b + 1] >= nSamplesInFrame ? abs (bMdct[bandWidth - 1]) : 0);

      for (int s = bandWidth - (bandStartOffsets[b + 1] >= nSamplesInFrame ? 1 : 0); s > 0; s--)
      {
        // based on S. Merdjani, L. Daudet, "Estimation of Frequency from MDCT-Encoded Files,"
        // DAFx-03, 2003, http://www.eecs.qmul.ac.uk/legacy/dafx03/proceedings/pdfs/dafx01.pdf
        sumAbsVal += complexAbs (*bMdct, int32_t ((*bNext - prevMdct) >> 1));
        bNext++;  prevMdct = *(bMdct++);
      }
      // average spectral sample magnitude across frequency band
      meanBandValues[b] = uint32_t ((sumAbsVal + (bandWidth >> 1)) / bandWidth);
    } // for b
  }
  if (channelIndex < USAC_MAX_NUM_CHANNELS) m_numAnaBands[channelIndex] = 0; // mark data used

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

unsigned SpecAnalyzer::initSigAnaMemory (LinearPredictor* const linPredictor, const unsigned nChannels, const unsigned maxTransfLength)
{
  if (linPredictor == nullptr)
  {
    return 1; // invalid arguments error
  }
  m_tnsPredictor = linPredictor;

  for (unsigned ch = 0; ch < nChannels; ch++)
  {
    if ((m_magnSpectra[ch] = (uint32_t*) malloc (maxTransfLength * sizeof (uint32_t))) == nullptr)
    {
      return 2; // mem. allocation error
    }
    memset (m_magnSpectra[ch], 0, maxTransfLength * sizeof (uint32_t));
  }
  return 0; // no error
}

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

unsigned SpecAnalyzer::spectralAnalysis (const int32_t* const mdctSignals[USAC_MAX_NUM_CHANNELS],
                                         const int32_t* const mdstSignals[USAC_MAX_NUM_CHANNELS],
                                         const unsigned nChannels, const unsigned nSamplesInFrame, const unsigned samplingRate,
                                         const unsigned lfeChannelIndex /*= USAC_MAX_NUM_CHANNELS*/) // to skip an LFE channel
{
  const uint64_t anaBwOffset = SA_BW >> 1;
  const unsigned lpcStopBand16k = (samplingRate <= 32000 ? nSamplesInFrame : (32000 * nSamplesInFrame) / samplingRate) >> SA_BW_SHIFT;
  const unsigned thresholdSlope = (48000 + SA_EPS * samplingRate) / 96000;
  const unsigned thresholdStart = samplingRate >> 15;

  if ((mdctSignals == nullptr) || (mdstSignals == nullptr) || (nChannels > USAC_MAX_NUM_CHANNELS) || (lfeChannelIndex > USAC_MAX_NUM_CHANNELS) ||
      (nSamplesInFrame > 2048) || (nSamplesInFrame <= 127) || (samplingRate < 7350) || (samplingRate > 96000))
  {
    return 1; // invalid arguments error
  }

  for (unsigned ch = 0; ch < nChannels; ch++)
  {
    const int32_t* const chMdct = mdctSignals[ch];
    const int32_t* const chMdst = mdstSignals[ch];
    uint32_t* const   chPrvMagn = m_magnSpectra[ch];
    const bool improvedSfmEstim = (chPrvMagn != nullptr);
    uint16_t currMC = 0, numMC = 0; // channel average

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
      const int32_t* const  bMdst = &chMdst[offs];
      uint32_t* const     prvMagn = (improvedSfmEstim ? &chPrvMagn[offs] : nullptr);
      uint16_t maxAbsIdx = 0;
      uint32_t maxAbsVal = 0, tmp = UINT_MAX;
      uint64_t sumAbsVal = 0;
      uint64_t sumAbsPrv = 0;
      uint64_t sumPrdCP  = 0, sumPrdCC = 0, sumPrdPP = 0;
      double ncp, dcc, dpp;

      for (int s = SA_BW - 1; s >= 0; s--)
      {
        // sum absolute values of complex spectrum, derive L1 norm, peak value, and peak index
        const uint64_t absSample = complexAbs (bMdct[s], bMdst[s]);

        if (improvedSfmEstim)   // correlation between current and previous magnitude spectrum
        {
          const uint64_t prvSample = prvMagn[s];

          sumPrdCP += (absSample * prvSample + anaBwOffset) >> SA_BW_SHIFT;
          sumPrdCC += (absSample * absSample + anaBwOffset) >> SA_BW_SHIFT;
          sumPrdPP += (prvSample * prvSample + anaBwOffset) >> SA_BW_SHIFT;

          sumAbsPrv += prvSample;
          prvMagn[s] = (uint32_t) absSample;
        }
        sumAbsVal += absSample;
        if (offs + s > 0) // exclude DC from max & min
        {
          if (maxAbsVal < absSample) // update maximum
          {
            maxAbsVal = (uint32_t) absSample;
            maxAbsIdx = (uint16_t) s;
          }
          if (tmp/*min*/> absSample) // update minimum
          {
            tmp/*min*/= (uint32_t) absSample;
          }
        }
      } // for s

      // bandwidth detection
      if ((m_bandwidthOff[ch] == 0) && (maxAbsVal > __max (thresholdSlope * (thresholdStart + b), SA_EPS)))
      {
        m_bandwidthOff[ch] = __max (maxAbsIdx + 5/*guard*/, SA_BW) + offs;
        m_bandwidthOff[ch] = __min (m_bandwidthOff[ch], nSamplesInFrame);
      }
      // save mean magnitude
      tmp/*mean*/ = uint32_t ((sumAbsVal + anaBwOffset) >> SA_BW_SHIFT);
      m_meanAbsValue[ch][b] = tmp;
      // spectral statistics
      if (improvedSfmEstim && (b > 0) && ((unsigned) b < lpcStopBand16k))
      {
        dcc = double (tmp);
        dpp = double ((sumAbsPrv + anaBwOffset) >> SA_BW_SHIFT);
        ncp = (sumPrdCP + dcc * dpp) * SA_BW - sumAbsVal * dpp - sumAbsPrv * dcc;
        dcc = (sumPrdCC + dcc * dcc) * SA_BW - sumAbsVal * dcc - sumAbsVal * dcc;
        dpp = (sumPrdPP + dpp * dpp) * SA_BW - sumAbsPrv * dpp - sumAbsPrv * dpp;
        sumPrdCP = uint64_t ((ncp <= 0.0) || (dcc * dpp <= 0.0) ? 0 : 0.5 + (256.0 * ncp * ncp) / (dcc * dpp));

        currMC += (uint16_t) __min (UCHAR_MAX, sumPrdCP); numMC++; // temporal correlation sum
      }
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
    if (samplingRate < 27713) sumAvgBand -= m_meanAbsValue[ch][b++];

    while (((unsigned) b + 1 < lpcStopBand16k) && ((uint64_t) m_meanAbsValue[ch][b] * (m_numAnaBands[ch] - 1) > sumAvgBand)) b++;
    b = __min (m_bandwidthOff[ch], b << SA_BW_SHIFT);
    if (samplingRate < 27713) sumAvgBand += m_meanAbsValue[ch][1];

    // obtain prediction gain across spectrum
    m_tnsPredGains[ch] = m_tnsPredictor->calcParCorCoeffs (&chMdct[b], __min (m_bandwidthOff[ch], lpcStopBand16k << SA_BW_SHIFT) - b,
                                                           MAX_PREDICTION_ORDER, m_parCorCoeffs[ch]);
    m_specAnaStats[ch] = packAvgSpecAnalysisStats (sumAvgBand, sumMaxBand, m_tnsPredGains[ch] >> 24, idxMaxSpec, (unsigned) b >> SA_BW_SHIFT);

    if (improvedSfmEstim)
    {
      if (numMC > 1) currMC = (currMC + (numMC >> 1)) / numMC;// smoothed temporal correlation
      valMaxSpec = (currMC + m_magnCorrPrev[ch] + 1) >> 1;
      m_magnCorrPrev[ch] = (uint8_t) currMC; // update

      if (valMaxSpec > ((m_specAnaStats[ch] >> 16) & UCHAR_MAX)) m_specAnaStats[ch] = (m_specAnaStats[ch] & 0xFF00FFFF) | (valMaxSpec << 16);
    }
  } // for ch

  return 0; // no error
}

int16_t SpecAnalyzer::stereoSigAnalysis (const int32_t* const mdctSignal1, const int32_t* const mdctSignal2,
                                         const int32_t* const mdstSignal1, const int32_t* const mdstSignal2,
                                         const unsigned nSamplesMax, const unsigned nSamplesInFrame, const bool shortTransforms,
                                         uint8_t* const stereoCorrValue) // per-band LR correlation
{
  const uint64_t anaBwOffset = SA_BW >> 1;
  const uint16_t numAnaBands = (shortTransforms ? nSamplesInFrame : nSamplesMax) >> SA_BW_SHIFT;
  const uint16_t numAnaModul = (shortTransforms ? numAnaBands >> 3 : numAnaBands + 1);
  int16_t b;

  if ((mdctSignal1 == nullptr) || (mdctSignal2 == nullptr) || (mdstSignal1 == nullptr) || (mdstSignal2 == nullptr) ||
      (nSamplesInFrame > 2048) || (nSamplesMax > 2048) || (numAnaBands == 0) || (numAnaModul == 0))
  {
    b = SHRT_MIN; // invalid arguments error
  }
  else
  {
    const uint16_t tempPreAnaPC = stereoCorrValue[0];
    uint16_t currPC = 0, numPC  = 0;
    uint64_t sumReM = 0, sumReS = 0;// mid-side RMS distribution

    for (b = numAnaBands - 1; b >= 0; b--)
    {
      const uint16_t anaBandModul = b % numAnaModul;  // to exclude first and last window band
      const uint16_t         offs = b << SA_BW_SHIFT; // start offset of current analysis band
      const int32_t* const lbMdct = &mdctSignal1[offs];
      const int32_t* const lbMdst = &mdstSignal1[offs];
      const int32_t* const rbMdct = &mdctSignal2[offs];
      const int32_t* const rbMdst = &mdstSignal2[offs];
      uint64_t sumMagnL = 0, sumMagnR = 0; // temporary RMS sums
      uint64_t sumPrdLR = 0, sumPrdLL = 0, sumPrdRR = 0;
      uint64_t sumRealL = 0, sumRealR = 0;
      uint64_t sumRealM = 0, sumRealS = 0, sumPrdMS; // mid-side
      double nlr, dll, drr;

      for (int s = SA_BW - 1; s >= 0; s--)
      {
        const uint64_t absMagnL = complexAbs (lbMdct[s], lbMdst[s]);
        const uint64_t absMagnR = complexAbs (rbMdct[s], rbMdst[s]);

        sumRealL += abs (lbMdct[s]);
        sumRealR += abs (rbMdct[s]);
        sumRealM += abs (lbMdct[s] + rbMdct[s]); // i.e., 2*mid,
        sumRealS += abs (lbMdct[s] - rbMdct[s]); // i.e., 2*side

        sumMagnL += absMagnL;
        sumMagnR += absMagnR;
        sumPrdLR += (absMagnL * absMagnR + anaBwOffset) >> SA_BW_SHIFT;
        sumPrdLL += (absMagnL * absMagnL + anaBwOffset) >> SA_BW_SHIFT;
        sumPrdRR += (absMagnR * absMagnR + anaBwOffset) >> SA_BW_SHIFT;
      } // for s

      sumRealL = (sumRealL + anaBwOffset) >> SA_BW_SHIFT; // avg
      sumRealR = (sumRealR + anaBwOffset) >> SA_BW_SHIFT;
      sumRealM = (sumRealM + anaBwOffset) >> SA_BW_SHIFT;
      sumRealS = (sumRealS + anaBwOffset) >> SA_BW_SHIFT;
      nlr = double (sumRealL * sumRealR) * 0.46875; // tuned for uncorrelated full-scale noise
      sumPrdMS = uint64_t (nlr > double (sumRealM * sumRealS) ? 256.0 : 0.5 + (512.0 * nlr) / __max (1.0, double (sumRealM * sumRealS)));

      dll = double ((sumMagnL + anaBwOffset) >> SA_BW_SHIFT);
      drr = double ((sumMagnR + anaBwOffset) >> SA_BW_SHIFT);
      nlr = (sumPrdLR + dll * drr) * SA_BW - sumMagnL * drr - sumMagnR * dll;
      dll = (sumPrdLL + dll * dll) * SA_BW - sumMagnL * dll - sumMagnL * dll;
      drr = (sumPrdRR + drr * drr) * SA_BW - sumMagnR * drr - sumMagnR * drr;
      sumPrdLR = uint64_t ((nlr <= 0.0) || (dll * drr <= 0.0) ? 0 : 0.5 + (256.0 * nlr * nlr) / (dll * drr));

      stereoCorrValue[b] = (uint8_t) __min (UCHAR_MAX, __max (sumPrdMS, sumPrdLR)); // in band

      if ((anaBandModul > 0) && (anaBandModul + 1 < numAnaModul)) // in frame (averaged below)
      {
        currPC += stereoCorrValue[b]; numPC++;
        sumReM += sumRealM;
        sumReS += sumRealS;
      }
    } // for b

    for (b = numAnaBands; b < int16_t (nSamplesInFrame >> SA_BW_SHIFT); b++)
    {
      stereoCorrValue[b] = UCHAR_MAX; // to allow joint-stereo coding at very high frequencies
    }

    if (numPC > 1) currPC = (currPC + (numPC >> 1)) / numPC; // frame's perceptual correlation

    if (currPC < tempPreAnaPC) currPC = (currPC + tempPreAnaPC + 1) >> 1;

    b = (int16_t) currPC * (sumReS * 2 > sumReM * 3 ? -1 : 1);  // negation implies side > mid
  }

  return b;
}
