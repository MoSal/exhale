/* tempAnalysis.cpp - source file for class providing temporal analysis of PCM signals
 * written by C. R. Helmrich, last modified in 2020 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under the exhale Copyright License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2020 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#include "exhaleLibPch.h"
#include "tempAnalysis.h"

static const int16_t lffc2x[65] = { // low-frequency filter coefficients
  // 269-pt. sinc windowed by 0.409 * cos(0*pi.*t) - 0.5 * cos(2*pi.*t) + 0.091 * cos(4*pi.*t)
  17887, -27755, 16590, -11782, 9095, -7371, 6166, -5273, 4582, -4029, 3576, -3196, 2873,
  -2594, 2350, -2135, 1944, -1773, 1618, -1478, 1351, -1235, 1129, -1032, 942, -860, 784,
  -714, 650, -591, 536, -485, 439, -396, 357, -321, 287, -257, 229, -204, 181, -160, 141,
  -124, 108, -95, 82, -71, 61, -52, 44, -37, 31, -26, 21, -17, 14, -11, 8, -6, 5, -3, 2, -1, 1
};

// static helper functions
static unsigned updateAbsStats (const int32_t* const chSig, const int nSamples, unsigned* const maxAbsVal, int16_t* const maxAbsIdx)
{
  const int32_t* const chSigM1 = chSig - 1; // for first-order high-pass
  unsigned sumAbs = 0;

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

static unsigned applyPitchPred (const int32_t* const chSig, const int nSamples, const int pitchLag, const int pitchSign = 1)
{
  const int32_t* const chSigM1 = chSig - 1; // for first-order high-pass
  const int32_t* const plSig   = chSig - pitchLag; // & pitch prediction
  const int32_t* const plSigM1 = plSig - 1;
  unsigned sumAbs = 0;

  for (int s = nSamples - 1; s >= 0; s--)
  {
    // compute absolute values of pitch-predicted high-pass signal, obtain L1 norm, peak value
    sumAbs += abs (chSig[s] - chSigM1[s] - pitchSign * (plSig[s] - plSigM1[s]));
  }
  return sumAbs;
}

static inline uint32_t packAvgTempAnalysisStats (const unsigned avgAbsHpL,  const unsigned avgAbsHpR, const unsigned avgAbsHpP,
                                                 const unsigned avgAbsPpLR, const unsigned maxAbsHpLR)
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
    m_maxIdxHpPrev[ch] = 1;
    m_pitchLagPrev[ch] = 0;
    m_tempAnaStats[ch] = 0;
    m_transientLoc[ch] = -1;
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

unsigned TempAnalyzer::temporalAnalysis (const int32_t* const timeSignals[USAC_MAX_NUM_CHANNELS], const unsigned nChannels,
                                         const int nSamplesInFrame, const unsigned lookaheadOffset, const uint8_t sbrShift,
                                         int32_t* const lrCoreTimeSignals[USAC_MAX_NUM_CHANNELS] /*= nullptr*/, // if using SBR
                                         const unsigned lfeChannelIndex /*= USAC_MAX_NUM_CHANNELS*/)  // to skip an LFE channel
{
  const bool applyResampler = (sbrShift > 0 && lrCoreTimeSignals != nullptr);
  const int halfFrameOffset = nSamplesInFrame >> 1;
  const int resamplerOffset = (int) lookaheadOffset - 128;

  if ((timeSignals == nullptr) || (nChannels > USAC_MAX_NUM_CHANNELS) || (lfeChannelIndex > USAC_MAX_NUM_CHANNELS) || (sbrShift > 1) ||
      (nSamplesInFrame > 2048) || (nSamplesInFrame < 2) || (lookaheadOffset > 4096) || (lookaheadOffset <= 256u * sbrShift))
  {
    return 1;
  }

  for (unsigned ch = 0; ch < nChannels; ch++)
  {
    const int32_t* const chSig   = &timeSignals[ch][lookaheadOffset];
    const int32_t* const chSigM1 = chSig - 1; // for first-order high-pass
    const int32_t* const chSigPH = chSig + halfFrameOffset;
// --- get L1 norm and pitch lag of both sides
    unsigned sumAbsValL = 0,  sumAbsValR = 0;
    unsigned maxAbsValL = 0,  maxAbsValR = 0;
    int16_t  maxAbsIdxL = 0,  maxAbsIdxR = 0;
    int      splitPtL   = 0;
    int      splitPtC   = halfFrameOffset;
    int      splitPtR   = nSamplesInFrame;
    unsigned uL0 = abs (chSig[splitPtL    ] - chSigM1[splitPtL    ]);
    unsigned uL1 = abs (chSig[splitPtC - 1] - chSigM1[splitPtC - 1]);
    unsigned uR0 = abs (chSig[splitPtC    ] - chSigM1[splitPtC    ]);
    unsigned uR1 = abs (chSig[splitPtR - 1] - chSigM1[splitPtR - 1]);
    unsigned u; // temporary value - register?

    if (applyResampler && lrCoreTimeSignals[ch] != nullptr) // downsampler
    {
      /*LF*/int32_t* lrSig = &lrCoreTimeSignals[ch][resamplerOffset >> sbrShift]; // low-rate,
      const int32_t* hrSig = &timeSignals[ch][resamplerOffset]; // high-rate input time signal

      for (int i = nSamplesInFrame >> sbrShift; i > 0; i--, lrSig++, hrSig += 2)
      {
        int64_t r  = ((int64_t) hrSig[0] << 17) + (hrSig[-1] + (int64_t) hrSig[1]) * -2*SHRT_MIN;
        int16_t s;

        for (u = 65, s = 129; u > 0; s -= 2) r += (hrSig[-s] + (int64_t) hrSig[s]) * lffc2x[--u];

        *lrSig = int32_t ((r + (1 << 17)) >> 18); // low-pass and low-rate
// TODO: bandpass
        if (*lrSig < -8388608) *lrSig = -8388608;
        else
        if (*lrSig >  8388607) *lrSig =  8388607;
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
      m_maxAbsHpPrev[ch] = 0; // maxAbsValR
      m_maxIdxHpPrev[ch] = 1; // maxAbsIdxR
      m_pitchLagPrev[ch] = 0; // pLagBestR
    }
    else // nonzero signal in the current frame
    {
      const int maxAbsIdxP = __max ((int) m_maxIdxHpPrev[ch] - nSamplesInFrame, 1 - (int) lookaheadOffset);
      unsigned   sumAbsHpL = sumAbsValL,  sumAbsHpR = sumAbsValR; // after high-pass filter
      unsigned   sumAbsPpL = sumAbsValL,  sumAbsPpR = sumAbsValR; // after pitch prediction
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
      m_transientLoc[ch] = packTransLocWithPitchLag (maxAbsValL, maxAbsValR, m_maxAbsHpPrev[ch],
                                                     maxAbsIdxL, maxAbsIdxR, __max (1, pLagBestR));
      // update stats history for this channel
      m_avgAbsHpPrev[ch] = sumAbsHpR;
      m_maxAbsHpPrev[ch] = maxAbsValR;
      m_maxIdxHpPrev[ch] = (unsigned) maxAbsIdxR;
      m_pitchLagPrev[ch] = (unsigned) pLagBestR;
    } // if sumAbsValL == 0 && sumAbsValR == 0
  } // for ch

  return 0; // no error
}
