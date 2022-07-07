/* tempAnalysis.cpp - source file for class providing temporal analysis of PCM signals
 * written by C. R. Helmrich, last modified in 2021 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under the exhale Copyright License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2021 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#include "exhaleLibPch.h"
#include "tempAnalysis.h"

static const int16_t lpfc12[65] = {  // 50% low-pass filter coefficients
  // 269-pt. sinc windowed by 0.409 * cos(0*pi.*t) - 0.5 * cos(2*pi.*t) + 0.091 * cos(4*pi.*t)
  17887, -27755, 16590, -11782, 9095, -7371, 6166, -5273, 4582, -4029, 3576, -3196, 2873,
  -2594, 2350, -2135, 1944, -1773, 1618, -1478, 1351, -1235, 1129, -1032, 942, -860, 784,
  -714, 650, -591, 536, -485, 439, -396, 357, -321, 287, -257, 229, -204, 181, -160, 141,
  -124, 108, -95, 82, -71, 61, -52, 44, -37, 31, -26, 21, -17, 14, -11, 8, -6, 5, -3, 2, -1, 1
};

static const int16_t lpfc34[128] = { // 25% low-pass filter coefficients
  // see also A. H. Nuttall, "Some Windows with Very Good Sidelobe Behavior," IEEE, Feb. 1981.
  3 /*<<16*/, 26221, -8914, 19626, 0, -11731, 13789, -8331, 0, 6431, -8148, 5212, 0, -4360,
  5688, -3728, 0, 3240, -4291, 2849, 0, -2529, 3378, -2260, 0, 2032, -2729, 1834, 0, -1662,
  2240, -1510, 0, 1375, -1856, 1253, 0, -1144, 1546, -1045, 0, 955, -1292, 873, 0, -798,
  1079, -729, 0, 666, -900, 608, 0, -555, 748, -505, 0, 459, -620, 418, 0, -379, 510, -343,
  0, 310, -417, 280, 0, -252, 338, -227, 0, 203, -272, 182, 0, -162, 216, -144, 0, 128, -170,
  113, 0, -100, 132, -88, 0, 77, -101, 67, 0, -58, 76, -50, 0, 43, -56, 37, 0, -31, 41, -26,
  0, 22, -28, 18, 0, -15, 19, -12, 0, 10, -12, 8, 0, -6, 7, -4, 0, 3, -4, 2, 0, -1, 2, -1
};

// static helper functions
static uint64_t updateAbsStats (const int32_t* const chSig, const int nSamples, unsigned* const maxAbsVal, int16_t* const maxAbsIdx)
{
  const int32_t* const chSigM1 = chSig - 1; // for first-order high-pass
  uint64_t sumAbs = 0;

  for (int s = nSamples - 1; s >= 0; s--)
  {
    // compute absolute values of high-pass signal, obtain L1 norm, peak value, and peak index
    const unsigned absSample = abs (chSig[s] - chSigM1[s]);

    sumAbs += absSample;
    if (*maxAbsVal < absSample)
    {
      *maxAbsVal = absSample;
      *maxAbsIdx = (int16_t) s;
    }
  }
  return sumAbs;
}

static uint64_t applyPitchPred (const int32_t* const chSig, const int nSamples, const int pitchLag, const int pitchSign = 1)
{
  const int32_t* const chSigM1 = chSig - 1; // for first-order high-pass
  const int32_t* const plSig   = chSig - pitchLag; // & pitch prediction
  const int32_t* const plSigM1 = plSig - 1;
  uint64_t sumAbs = 0;

  for (int s = nSamples - 1; s >= 0; s--)
  {
    // compute absolute values of pitch-predicted high-pass signal, obtain L1 norm, peak value
    sumAbs += abs (chSig[s] - chSigM1[s] - pitchSign * (plSig[s] - plSigM1[s]));
  }
  return sumAbs;
}

static inline uint32_t packAvgTempAnalysisStats (const uint64_t avgAbsHpL,  const uint64_t avgAbsHpR, const unsigned avgAbsHpP,
                                                 const uint64_t avgAbsPpLR, const unsigned maxAbsHpLR)
{
  // spectral flatness, normalized for a value of 256 for noise-like, spectrally flat waveform
  const unsigned flatSpec = 256 - int ((int64_t (avgAbsPpLR/*L+R sum*/ + TA_EPS) * 256) / (int64_t (avgAbsHpL + avgAbsHpR + TA_EPS)));
  // temporal flatness, normalized for a value of 256 for steady low or mid-frequency sinusoid
  const int32_t  flatTemp = 256 - int ((int64_t (avgAbsHpL + avgAbsHpR + TA_EPS) * 402) / (int64_t (maxAbsHpLR/*L+R sum*/ + TA_EPS)));
  // temporal stationarity, two sides, normalized for values of 256 for L1-stationary waveform
  const int32_t  statTmpL = 256 - int (((__min  (avgAbsHpP, avgAbsHpL) + TA_EPS) * 256) / ((__max  (avgAbsHpP, avgAbsHpL) + TA_EPS)));
  const int32_t  statTmpR = 256 - int (((__min  (avgAbsHpL, avgAbsHpR) + TA_EPS) * 256) / ((__max  (avgAbsHpL, avgAbsHpR) + TA_EPS)));

  return (CLIP_UCHAR (flatSpec) << 24) | (CLIP_UCHAR (flatTemp) << 16) | (CLIP_UCHAR (statTmpL) << 8) | CLIP_UCHAR (statTmpR);
}

static inline int16_t packTransLocWithPitchLag (const unsigned maxAbsValL, const unsigned maxAbsValR, const unsigned maxAbsValP,
                                                const int16_t  maxAbsIdxL, const int16_t  maxAbsIdxR, const int16_t  optPitchLag)
{
  if ((maxAbsValP * 5 < maxAbsValL * 2) || (maxAbsValL * 5 < maxAbsValR * 2)) // has transient
  {
    return (((maxAbsValR > maxAbsValL ? maxAbsIdxR : maxAbsIdxL) << 4) & 0xF800) | __min (2047, optPitchLag);
  }
  return -1 * optPitchLag; // has no transient
}

// constructor
TempAnalyzer::TempAnalyzer ()
{
  for (unsigned ch = 0; ch < USAC_MAX_NUM_CHANNELS; ch++)
  {
    m_avgAbsHpPrev[ch] = 0;
    m_maxAbsHpPrev[ch] = 0;
    m_maxHfLevPrev[ch] = 0;
    m_maxIdxHpPrev[ch] = 1;
    m_pitchLagPrev[ch] = 0;
    m_tempAnaStats[ch] = 0;
    m_transientLoc[ch] = -1;

    memset (m_filtSampPrev[ch], 0, 6 * sizeof (int64_t));
  }
}

// public functions
void TempAnalyzer::getTempAnalysisStats (uint32_t avgTempAnaStats[USAC_MAX_NUM_CHANNELS], const unsigned nChannels)
{
  if ((avgTempAnaStats == nullptr) || (nChannels > USAC_MAX_NUM_CHANNELS))
  {
    return;
  }
  memcpy (avgTempAnaStats, m_tempAnaStats, nChannels * sizeof (uint32_t));
}

void TempAnalyzer::getTransientAndPitch (int16_t transIdxAndPitch[USAC_MAX_NUM_CHANNELS], const unsigned nChannels)
{
  if ((transIdxAndPitch == nullptr) || (nChannels > USAC_MAX_NUM_CHANNELS))
  {
    return;
  }
  memcpy (transIdxAndPitch, m_transientLoc, nChannels * sizeof (int16_t));
}

uint8_t TempAnalyzer::stereoPreAnalysis (const int32_t* const timeSignals[2], const uint8_t specFlatness[2], const unsigned nSamplesInSig)
{
  const double   offsetSfmLR  = __max (0.0, ((double) specFlatness[0] + specFlatness[1] - 256.0) * 0.5);
  const int32_t* const sigL   = timeSignals[0] + (nSamplesInSig >> 1);
  const int32_t* const sigLM1 = sigL - 1;
  const int32_t* const sigR   = timeSignals[1] + (nSamplesInSig >> 1);
  const int32_t* const sigRM1 = sigR - 1;
  int64_t hpNextL = sigL[nSamplesInSig] - sigLM1[nSamplesInSig];
  int64_t hpNextR = sigR[nSamplesInSig] - sigRM1[nSamplesInSig];
  int64_t sumSqrL = hpNextL * hpNextL, sumSqrR = hpNextR * hpNextR;
  int64_t sumPC00 = (hpNextL * hpNextR) >> 1, sumPC01 = 0, sumPC10 = 0;
  double d;

  for (int s = nSamplesInSig - 1; s >= 0; s--)
  {
    // compute correlation between high-pass channel signals with and without 1 smp time delay
    const int64_t hpL = sigL[s] - sigLM1[s];
    const int64_t hpR = sigR[s] - sigRM1[s];

    sumSqrL += hpL * hpL;
    sumSqrR += hpR * hpR;
    sumPC00 += hpL * hpR;
    sumPC01 += hpL * hpNextR;
    sumPC10 += hpR * hpNextL;

    hpNextL = hpL;
    hpNextR = hpR;
  }

  if (sumSqrL < nSamplesInSig || sumSqrR < nSamplesInSig) return 0; // stop on low-level input

  sumPC00 = abs (sumPC00);
  sumPC01 = abs (sumPC01);
  sumPC10 = abs (sumPC10);

  d = 256.0 * __max (sumPC00, __max (sumPC01, sumPC10)); // max. corr. regardless of the delay

  return (uint8_t) __max (0.0, d / sqrt ((double) sumSqrL * sumSqrR) - offsetSfmLR);
}

unsigned TempAnalyzer::temporalAnalysis (const int32_t* const timeSignals[USAC_MAX_NUM_CHANNELS], const unsigned nChannels,
                                         const int nSamplesInFrame, const unsigned lookaheadOffset, const uint8_t sbrShift,
                                         int32_t* const lrCoreTimeSignals[USAC_MAX_NUM_CHANNELS] /*= nullptr*/, // if using SBR
                                         const unsigned lfeChannelIndex /*= USAC_MAX_NUM_CHANNELS*/)  // to skip an LFE channel
{
  const bool applyResampler = (sbrShift > 0 && lrCoreTimeSignals != nullptr);
  const int halfFrameOffset = nSamplesInFrame >> 1;
  const int resamplerOffset = (int) lookaheadOffset - 128;

  if ((timeSignals == nullptr) || (nChannels > USAC_MAX_NUM_CHANNELS) || (lfeChannelIndex > USAC_MAX_NUM_CHANNELS) || (sbrShift > 1) ||
      (nSamplesInFrame > 2048) || (nSamplesInFrame <= 128 * sbrShift) || (lookaheadOffset > 4096) || (lookaheadOffset <= 256u * sbrShift))
  {
    return 1;
  }

  for (unsigned ch = 0; ch < nChannels; ch++)
  {
    const int32_t* const chSig   = &timeSignals[ch][lookaheadOffset];
    const int32_t* const chSigM1 = chSig - 1; // for first-order high-pass
    const int32_t* const chSigPH = chSig + halfFrameOffset;
// --- get L1 norm and pitch lag of both sides
    uint64_t sumAbsValL = 0,  sumAbsValR = 0;
    unsigned maxAbsValL = 0,  maxAbsValR = 0;
    int32_t  maxHfrLevL = 0,  maxHfrLevR = 0;
    int16_t  maxAbsIdxL = 0,  maxAbsIdxR = 0;
    int      splitPtL   = 0;
    int      splitPtC   = halfFrameOffset;
    int      splitPtR   = nSamplesInFrame;
    uint64_t ue[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0}; // sub-fr. unit energies
    unsigned uL0 = abs (chSig[splitPtL    ] - chSigM1[splitPtL    ]);
    unsigned uL1 = abs (chSig[splitPtC - 1] - chSigM1[splitPtC - 1]);
    unsigned uR0 = abs (chSig[splitPtC    ] - chSigM1[splitPtC    ]);
    unsigned uR1 = abs (chSig[splitPtR - 1] - chSigM1[splitPtR - 1]);
    unsigned u; // temporary value - register?

    if (applyResampler && lrCoreTimeSignals[ch] != nullptr) // downsampler
    {
      /*LF*/int32_t* lrSig = &lrCoreTimeSignals[ch][resamplerOffset >> sbrShift];
      const int32_t* hrSig = &timeSignals[ch][resamplerOffset];
      int64_t* const rPrev = m_filtSampPrev[ch];
      uint64_t     subSumL = 0, subSumM = 0, subSumH = 0;

      for (int i = nSamplesInFrame >> sbrShift; i > 0; i--, lrSig++, hrSig += 2)
      {
        int64_t r  = ((int64_t) hrSig[0] * (1 << 17)) + (hrSig[-1] + (int64_t) hrSig[1]) * -2*SHRT_MIN;
        int16_t s;

        for (u = 65, s = 129; u > 0; s -= 2) r += (hrSig[-s] + (int64_t) hrSig[s]) * lpfc12[--u];

        *lrSig = int32_t ((r + (1 << 17)) >> 18); // low-pass at half rate
        if (*lrSig < -8388608) *lrSig = -8388608;
        else
        if (*lrSig >  8388607) *lrSig =  8388607;

        if ((i & 1) != 0) // compute quarter-rate mid-frequency SBR signal
        {
          r  = ((3 * (int64_t) hrSig[0]) * (1 << 16)) - (hrSig[-1] + (int64_t) hrSig[1]) * SHRT_MIN - r;
          r += (hrSig[-2] + (int64_t) hrSig[2]) * SHRT_MIN;

          for (s = 127; s > 0; s--/*u = s*/) r += (hrSig[-s] + (int64_t) hrSig[s]) * lpfc34[s];

          r = (r + (1 << 17)) >> 18; // SBR env. band-pass at quarter rate
          ue[i >> 7] += square (r);

          // calculate 3 SBR subband envelope energies (low, mid and high)
          subSumL += square ((6 * rPrev[2] + 5 * (rPrev[1] + rPrev[3]) + 3 * (rPrev[0] + rPrev[4]) + (r + rPrev[5]) + 8) >> 4);
          subSumM += square ((2 * rPrev[2] - (rPrev[0] + rPrev[4]) + 2) >> 2);
          subSumH += square ((6 * rPrev[2] - 5 * (rPrev[1] + rPrev[3]) + 3 * (rPrev[0] + rPrev[4]) - (r + rPrev[5]) + 8) >> 4);

          rPrev[5] = rPrev[4];  rPrev[4] = rPrev[3];  rPrev[3] = rPrev[2];
          rPrev[2] = rPrev[1];  rPrev[1] = rPrev[0];  rPrev[0] = r;
        }
      }

      if (ch != lfeChannelIndex) // calculate overall and unit-wise levels
      {
        const unsigned numUnits = nSamplesInFrame >> (sbrShift + 7);
        int32_t* const hfrLevel = &lrCoreTimeSignals[ch][(resamplerOffset + nSamplesInFrame) >> sbrShift];

        for (u = numUnits; u > 0;  )
        {
          ue[8] += ue[--u];
          hfrLevel[numUnits - 1 - u] = int32_t (0.5 + sqrt ((double) ue[u]));
        }

        if (ue[8] < 1) ue[8] = 1;  // low, mid, high subband energy ratios
        hfrLevel[numUnits]   = int32_t (0.5 + __min (USHRT_MAX, (21845.3 * subSumL) / ue[8]));
        hfrLevel[numUnits]  |= int32_t (0.5 + __min ( SHRT_MAX, (21845.3 * subSumM) / ue[8])) << 16;
        hfrLevel[numUnits+1] = int32_t (0.5 + __min (USHRT_MAX, (21845.3 * subSumH) / ue[8]));

        for (u = numUnits >> 1; u > 0;  ) // stabilize transient detection
        {
          u--;
          if (maxHfrLevL < hfrLevel[u]) /* update max. */ maxHfrLevL = hfrLevel[u];
          if (maxHfrLevR < hfrLevel[u + (numUnits >> 1)]) maxHfrLevR = hfrLevel[u + (numUnits >> 1)];
        }
      }
    }

    if (ch == lfeChannelIndex)  // no analysis
    {
      m_tempAnaStats[ch] = 0; // flat/stationary frame
      m_transientLoc[ch] = -1;
      continue;
    }

    do // find last sample of left-side region
    {
      sumAbsValL += (u = uL1);
      splitPtC--;
    }
    while ((splitPtC > /*start +*/1) && (uL1 = abs (chSig[splitPtC - 1] - chSigM1[splitPtC - 1])) < u);

    do // find first sample of left-side range
    {
      sumAbsValL += (u = uL0);
      splitPtL++;
    }
    while ((splitPtL < splitPtC - 1) && (uL0 = abs (chSig[splitPtL] - chSigM1[splitPtL])) < u);

    sumAbsValL += updateAbsStats (&chSig[splitPtL], splitPtC - splitPtL, &maxAbsValL, &maxAbsIdxL);
    maxAbsIdxL += splitPtL; // left-side stats
    if ((maxAbsIdxL == 1) && (maxAbsValL <= u))
    {
      maxAbsValL = u;
      maxAbsIdxL--;
    }

    splitPtC = halfFrameOffset;

    do // find last sample of right-side region
    {
      sumAbsValR += (u = uR1);
      splitPtR--;
    }
    while ((splitPtR > splitPtC + 1) && (uR1 = abs (chSig[splitPtR - 1] - chSigM1[splitPtR - 1])) < u);

    do // find first sample of right-side range
    {
      sumAbsValR += (u = uR0);
      splitPtC++;
    }
    while ((splitPtC < splitPtR - 1) && (uR0 = abs (chSig[splitPtC] - chSigM1[splitPtC])) < u);

    sumAbsValR += updateAbsStats (&chSig[splitPtC], splitPtR - splitPtC, &maxAbsValR, &maxAbsIdxR);
    maxAbsIdxR += splitPtC; // right-side stats
    if ((maxAbsIdxR == halfFrameOffset + 1) && (maxAbsValR <= u))
    {
      maxAbsValR = u;
      maxAbsIdxR--;
    }

// --- find best pitch lags minimizing L1 norms
    if (sumAbsValL == 0 && sumAbsValR == 0)
    {
      m_tempAnaStats[ch] = 0; // flat/stationary frame
      m_transientLoc[ch] = -1;
      // re-init stats history for this channel
      m_avgAbsHpPrev[ch] = 0;
      m_maxAbsHpPrev[ch] = 0;
      m_maxIdxHpPrev[ch] = 1;
      m_pitchLagPrev[ch] = 0;
    }
    else // nonzero signal in the current frame
    {
      const int maxAbsIdxP = __max ((int) m_maxIdxHpPrev[ch] - nSamplesInFrame, 1 - (int) lookaheadOffset);
      uint64_t   sumAbsHpL = sumAbsValL,  sumAbsHpR = sumAbsValR; // after high-pass filter
      uint64_t   sumAbsPpL = sumAbsValL,  sumAbsPpR = sumAbsValR; // after pitch prediction
      int pLag,  pLagBestR = 0,  pSgn;

      // test left-side pitch lag on this frame
      pLag = __min (maxAbsIdxL - maxAbsIdxP, (int) lookaheadOffset - 1);
      pSgn = (((chSig[maxAbsIdxL] - chSigM1[maxAbsIdxL] > 0) && (chSig[maxAbsIdxP] - chSigM1[maxAbsIdxP] < 0)) ||
              ((chSig[maxAbsIdxL] - chSigM1[maxAbsIdxL] < 0) && (chSig[maxAbsIdxP] - chSigM1[maxAbsIdxP] > 0)) ? -1 : 1);
      if ((sumAbsValL = applyPitchPred (chSig, halfFrameOffset, pLag, pSgn)) < sumAbsPpL)
      {
        sumAbsPpL = sumAbsValL; // left side
      }
      if ((sumAbsValR = applyPitchPred (chSigPH, halfFrameOffset, pLag, pSgn)) < sumAbsPpR)
      {
        sumAbsPpR = sumAbsValR; // right side
        pLagBestR = pLag;
      }
      // test right-side pitch lag on the frame
      pLag = __min (maxAbsIdxR - maxAbsIdxL, (int) lookaheadOffset - 1);
      pSgn = (((chSig[maxAbsIdxR] - chSigM1[maxAbsIdxR] > 0) && (chSig[maxAbsIdxL] - chSigM1[maxAbsIdxL] < 0)) ||
              ((chSig[maxAbsIdxR] - chSigM1[maxAbsIdxR] < 0) && (chSig[maxAbsIdxL] - chSigM1[maxAbsIdxL] > 0)) ? -1 : 1);
      if ((sumAbsValL = applyPitchPred (chSig, halfFrameOffset, pLag, pSgn)) < sumAbsPpL)
      {
        sumAbsPpL = sumAbsValL; // left side
      }
      if ((sumAbsValR = applyPitchPred (chSigPH, halfFrameOffset, pLag, pSgn)) < sumAbsPpR)
      {
        sumAbsPpR = sumAbsValR; // right side
        pLagBestR = pLag;
      }
      // try previous frame's lag on this frame
      pLag = (m_pitchLagPrev[ch] > 0 ? (int) m_pitchLagPrev[ch] : __min (halfFrameOffset, (int) lookaheadOffset - 1));
      pSgn = (((chSig[maxAbsIdxL] - chSigM1[maxAbsIdxL] > 0) && (chSig[maxAbsIdxL-pLag] - chSigM1[maxAbsIdxL-pLag] < 0)) ||
              ((chSig[maxAbsIdxL] - chSigM1[maxAbsIdxL] < 0) && (chSig[maxAbsIdxL-pLag] - chSigM1[maxAbsIdxL-pLag] > 0)) ? -1 : 1);
      if ((sumAbsValL = applyPitchPred (chSig, halfFrameOffset, pLag, pSgn)) < sumAbsPpL)
      {
        sumAbsPpL = sumAbsValL; // left side
      }
      if ((sumAbsValR = applyPitchPred (chSigPH, halfFrameOffset, pLag, pSgn)) < sumAbsPpR)
      {
        sumAbsPpR = sumAbsValR; // right side
        pLagBestR = pLag;
      }
      if (pLagBestR >= halfFrameOffset) // half
      {
        pLag = pLagBestR >> 1;
        pSgn = (((chSig[maxAbsIdxR] - chSigM1[maxAbsIdxR] > 0) && (chSig[maxAbsIdxR-pLag] - chSigM1[maxAbsIdxR-pLag] < 0)) ||
                ((chSig[maxAbsIdxR] - chSigM1[maxAbsIdxR] < 0) && (chSig[maxAbsIdxR-pLag] - chSigM1[maxAbsIdxR-pLag] > 0)) ? -1 : 1);
        if ((sumAbsValL = applyPitchPred (chSig, halfFrameOffset, pLag, pSgn)) < sumAbsPpL)
        {
          sumAbsPpL = sumAbsValL; // left side
        }
        if ((sumAbsValR = applyPitchPred (chSigPH, halfFrameOffset, pLag, pSgn)) < sumAbsPpR)
        {
          sumAbsPpR = sumAbsValR; // right side
          pLagBestR = pLag;
        }
      }

      // convert L1 norms into average values
      sumAbsHpL = (sumAbsHpL + unsigned (halfFrameOffset >> 1)) / unsigned (halfFrameOffset);
      sumAbsHpR = (sumAbsHpR + unsigned (halfFrameOffset >> 1)) / unsigned (halfFrameOffset);
      sumAbsPpL = (sumAbsPpL + unsigned (halfFrameOffset >> 1)) / unsigned (halfFrameOffset);
      sumAbsPpR = (sumAbsPpR + unsigned (halfFrameOffset >> 1)) / unsigned (halfFrameOffset);
// --- temporal analysis statistics for frame
      m_tempAnaStats[ch] = packAvgTempAnalysisStats (sumAbsHpL,  sumAbsHpR,  m_avgAbsHpPrev[ch],
                                                     sumAbsPpL + sumAbsPpR,  maxAbsValL + maxAbsValR);
      u = maxAbsValR;
      if ((m_maxHfLevPrev[ch] < (maxHfrLevL >> 4)) || (maxHfrLevL < (maxHfrLevR >> 4))) // HF
      {
        maxAbsValL = maxHfrLevL;
        maxAbsValR = maxHfrLevR;
        m_maxAbsHpPrev[ch] = m_maxHfLevPrev[ch];
      }
      else
      {
        memset (ue, 0, 8 * sizeof (uint64_t));
        for (u = nSamplesInFrame - 1; u > 0; u--) ue[u >> 8] += abs (chSig[u] - chSigM1[u]);

        sumAbsValL = ue[0];
        sumAbsValR = (uint64_t) maxAbsValL + (uint64_t) maxAbsValR;
        for (u = (nSamplesInFrame >> 8) - 1; u > 0; u--) sumAbsValL = __min (sumAbsValL, ue[u]);

        u = maxAbsValR;
        if (sumAbsValL < sumAbsValR * (1u + (nSamplesInFrame >> 10)) && m_maxAbsHpPrev[ch] > TA_EPS) m_maxAbsHpPrev[ch] = TA_EPS;
      }
      m_transientLoc[ch] = packTransLocWithPitchLag (maxAbsValL, maxAbsValR, m_maxAbsHpPrev[ch],
                                                     maxAbsIdxL, maxAbsIdxR, __max (1, pLagBestR));
      // update stats history for this channel
      m_avgAbsHpPrev[ch] = (unsigned) sumAbsHpR;
      m_maxAbsHpPrev[ch] = u;
      m_maxIdxHpPrev[ch] = (unsigned) maxAbsIdxR;
      m_pitchLagPrev[ch] = (unsigned) pLagBestR;
    } // if sumAbsValL == 0 && sumAbsValR == 0

    if (applyResampler) m_maxHfLevPrev[ch] = maxHfrLevR;
  } // ch

  return 0; // no error
}
