/* specAnalysis.cpp - source file for class providing spectral analysis of MCLT signals
 * written by C. R. Helmrich, last modified in 2020 - see License.htm for legal notices
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
#if SA_IMPROVED_SFM_ESTIM
    m_magnCorrPrev[ch] = 0;
    m_magnSpectra [ch] = nullptr;
#endif
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
  if ((mdctSignal == nullptr) || (bandStartOffsets == nullptr) || (meanBandValues == nullptr) || (channelIndex > USAC_MAX_NUM_CHANNELS) ||
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

      if ((channelIndex < USAC_MAX_NUM_CHANNELS) && (anaBandIdx < m_numAnaBands[channelIndex]) &&
          (bandOffset == (anaBandIdx << SA_BW_SHIFT)) && ((bandWidth & (SA_BW - 1)) == 0))
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

          sumAbsVal += uint64_t (sqrt (complexSqr) + 0.5);
#else
          const uint64_t absReal   = abs (bMdct[s]); // Richard Lyons, 1997; en.wikipedia.org/
          const uint64_t absImag   = abs (bMdst[s]); // wiki/Alpha_max_plus_beta_min_algorithm

          sumAbsVal += (absReal > absImag ? absReal + ((absImag * 3) >> 3) : absImag + ((absReal * 3) >> 3));
#endif
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
      uint64_t sumAbsVal = (bandStartOffsets[b + 1] == nSamplesInFrame ? abs (bMdct[bandWidth - 1]) : 0);

      for (int s = bandWidth - (bandStartOffsets[b + 1] == nSamplesInFrame ? 2 : 1); s >= 0; s--)
      {
        // based on S. Merdjani, L. Daudet, "Estimation of Frequency from MDCT-Encoded Files,"
        // DAFx-03, 2003, http://www.eecs.qmul.ac.uk/legacy/dafx03/proceedings/pdfs/dafx01.pdf
        const uint64_t absReal   = abs (bMdct[s]); // Richard Lyons, 1997; see also code above
        const uint64_t absEstIm  = abs (bMdct[s + 1] - bMdct[s - 1]) >> 1; // s - 1 may be -1!

        sumAbsVal += (absReal > absEstIm ? absReal + ((absEstIm * 3) >> 3) : absEstIm + ((absReal * 3) >> 3));
      }
      if ((b == 0) && (bandWidth > 0)) // correct estimate at DC
      {
        const uint64_t absReal   = abs (bMdct[0]); // Richard Lyons, 1997; see also code above
        const uint64_t absEstIm  = abs (bMdct[1] - bMdct[-1]) >> 1; // if s - 1 was -1 earlier

        sumAbsVal -= (absReal > absEstIm ? absReal + ((absEstIm * 3) >> 3) : absEstIm + ((absReal * 3) >> 3));
        sumAbsVal += absReal;
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

unsigned SpecAnalyzer::initSigAnaMemory (LinearPredictor* const linPredictor, const unsigned nChannels, const unsigned maxTransfLength)
{
  if (linPredictor == nullptr)
  {
    return 1; // invalid arguments error
  }
  m_tnsPredictor = linPredictor;
#if SA_IMPROVED_SFM_ESTIM
  for (unsigned ch = 0; ch < nChannels; ch++)
  {
    if ((m_magnSpectra[ch] = (uint32_t*) malloc (maxTransfLength * sizeof (uint32_t))) == nullptr)
    {
      return 2; // mem. allocation error
    }
    memset (m_magnSpectra[ch], 0, maxTransfLength * sizeof (uint32_t));
  }
#endif
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
  const uint64_t anaBwOffset = SA_BW >> 1;
  const unsigned lpcStopBand16k = (samplingRate <= 32000 ? nSamplesInFrame : (32000 * nSamplesInFrame) / samplingRate) >> SA_BW_SHIFT;
  const unsigned thresholdSlope = (48000 + SA_EPS * samplingRate) / 96000;
  const unsigned thresholdStart = samplingRate >> 15;

  if ((mdctSignals == nullptr) || (mdstSignals == nullptr) || (nChannels > USAC_MAX_NUM_CHANNELS) || (lfeChannelIndex > USAC_MAX_NUM_CHANNELS) ||
      (nSamplesInFrame > 2048) || (nSamplesInFrame < 2) || (samplingRate < 7350) || (samplingRate > 96000))
  {
    return 1; // invalid arguments error
  }

  for (unsigned ch = 0; ch < nChannels; ch++)
  {
    const int32_t* const chMdct = mdctSignals[ch];
    const int32_t* const chMdst = mdstSignals[ch];
#if SA_IMPROVED_SFM_ESTIM
    uint32_t* const   chPrvMagn = m_magnSpectra[ch];
    const bool improvedSfmEstim = (chPrvMagn != nullptr);
    uint16_t currMC = 0, numMC = 0; // channel average
#endif
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
#if SA_IMPROVED_SFM_ESTIM
      uint32_t* const     prvMagn = (improvedSfmEstim ? &chPrvMagn[offs] : nullptr);
#endif
      uint16_t maxAbsIdx = 0;
      uint32_t maxAbsVal = 0, tmp = UINT_MAX;
      uint64_t sumAbsVal = 0;
#if SA_IMPROVED_SFM_ESTIM
      uint64_t sumAbsPrv = 0;
      uint64_t sumPrdCP  = 0, sumPrdCC = 0, sumPrdPP = 0;
      double ncp, dcc, dpp;
#endif
      for (int s = SA_BW - 1; s >= 0; s--)
      {
        // sum absolute values of complex spectrum, derive L1 norm, peak value, and peak index
#if SA_EXACT_COMPLEX_ABS
        const double  complexSqr = (double) bMdct[s] * (double) bMdct[s] + (double) bMdst[s] * (double) bMdst[s];
        const uint32_t absSample = uint32_t (sqrt (complexSqr) + 0.5);
#else
        const uint64_t absReal   = abs (bMdct[s]);   // Richard Lyons, 1997; en.wikipedia.org/
        const uint64_t absImag   = abs (bMdst[s]);   // wiki/Alpha_max_plus_beta_min_algorithm
        const uint32_t absSample = uint32_t (absReal > absImag ? absReal + ((absImag * 3) >> 3) : absImag + ((absReal * 3) >> 3));
#endif
#if SA_IMPROVED_SFM_ESTIM
        if (improvedSfmEstim)   // correlation between current and previous magnitude spectrum
        {
          const uint64_t prvSample = prvMagn[s];

          sumPrdCP += ((uint64_t) absSample * prvSample + anaBwOffset) >> SA_BW_SHIFT;
          sumPrdCC += ((uint64_t) absSample * absSample + anaBwOffset) >> SA_BW_SHIFT;
          sumPrdPP += ((uint64_t) prvSample * prvSample + anaBwOffset) >> SA_BW_SHIFT;

          sumAbsPrv += prvSample;
          prvMagn[s] = absSample;
        }
#endif
        sumAbsVal += absSample;
        if (offs + s > 0) // exclude DC from max & min
        {
          if (maxAbsVal < absSample) // update maximum
          {
            maxAbsVal = absSample;
            maxAbsIdx = (uint16_t) s;
          }
          if (tmp/*min*/> absSample) // update minimum
          {
            tmp/*min*/= absSample;
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
#if SA_IMPROVED_SFM_ESTIM
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
#endif
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
#if SA_IMPROVED_SFM_ESTIM
    if (improvedSfmEstim)
    {
      if (numMC > 1) currMC = (currMC + (numMC >> 1)) / numMC;// smoothed temporal correlation
      valMaxSpec = (currMC + m_magnCorrPrev[ch] + 1) >> 1;
      m_magnCorrPrev[ch] = (uint8_t) currMC; // update

      if (valMaxSpec > ((m_specAnaStats[ch] >> 16) & UCHAR_MAX)) m_specAnaStats[ch] = (m_specAnaStats[ch] & 0xFF00FFFF) | (valMaxSpec << 16);
    }
#endif
  } // for ch

  return 0; // no error
}

int16_t SpecAnalyzer::stereoSigAnalysis (const int32_t* const mdctSignal1, const int32_t* const mdctSignal2,
                                         const int32_t* const mdstSignal1, const int32_t* const mdstSignal2,
                                         const unsigned nSamplesMax, const unsigned nSamplesInFrame, const bool shortTransforms,
                                         uint8_t* const stereoCorrValue /*= nullptr*/) // per-band perceptual correlation data
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
    uint16_t currPC = 0, numPC = 0; // frame-average correlation
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
        const uint32_t absRealL  = abs (lbMdct[s]);
        const uint32_t absRealR  = abs (rbMdct[s]);
#if SA_EXACT_COMPLEX_ABS
        const double complexSqrL = (double) lbMdct[s] * (double) lbMdct[s] + (double) lbMdst[s] * (double) lbMdst[s];
        const uint32_t absMagnL  = uint32_t (sqrt (complexSqrL) + 0.5);
        const double complexSqrR = (double) rbMdct[s] * (double) rbMdct[s] + (double) rbMdst[s] * (double) rbMdst[s];
        const uint32_t absMagnR  = uint32_t (sqrt (complexSqrR) + 0.5);
#else
        const uint64_t absImagL  = abs (lbMdst[s]);  // Richard Lyons, 1997; en.wikipedia.org/
        const uint64_t absImagR  = abs (rbMdst[s]);  // wiki/Alpha_max_plus_beta_min_algorithm
        const uint64_t absMagnL  = (absRealL > absImagL ? absRealL + ((absImagL * 3) >> 3) : absImagL + ((absRealL * 3) >> 3));
        const uint64_t absMagnR  = (absRealR > absImagR ? absRealR + ((absImagR * 3) >> 3) : absImagR + ((absRealR * 3) >> 3));
#endif
        sumRealL += absRealL;
        sumRealR += absRealR;
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

    b = (int16_t) currPC * (sumReS * 2 > sumReM * 3 ? -1 : 1);  // negation implies side > mid
  }

  return b;
}
