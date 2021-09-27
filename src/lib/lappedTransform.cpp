/* lappedTransform.cpp - source file for class providing time-frequency transformation
 * written by C. R. Helmrich, last modified in 2019 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under the exhale Copyright License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2021 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#include "exhaleLibPch.h"
#include "lappedTransform.h"
#ifdef _MSC_VER
# include <intrin.h> // _BitScanReverse

# pragma intrinsic (_BitScanReverse)
#endif

// static helper functions
static short* createPermutTable (const short tableSize)
{
  const short lOver2 = tableSize >> 1;
  short* permutTable = nullptr;
  short  i = 0;

  if ((permutTable = (short*) malloc (tableSize * sizeof (short))) == nullptr)
  {
    return nullptr; // allocation error
  }

  permutTable[0] = 0;
  for (short s = 1; s < tableSize; s++)
  {
    short l = lOver2;

    while (i >= l)
    {
      i -= l;
      l >>= 1;
    }
    permutTable[s] = (i += l);
  }

  return permutTable;
}

static inline int shortIntLog2 (uint16_t s)
{
#ifdef _MSC_VER
  unsigned long l;

  _BitScanReverse (&l, s);

  return (int) l;
#else
  // fast base-2 integer logarithm by Todd Lehman (Stack Overflow), July 2014
  int i = 0;
# define S(k) if (s >= (1u << k)) { i += k; s >>= k; }
  S (8); S (4); S (2); S (1);
# undef S
  return i;
#endif
}

// private helper functions
void LappedTransform::applyHalfSizeFFT (int32_t* const iR/*eal*/, int32_t* const iI/*mag*/, const bool shortTransform) // works in-place
{
  // int32 FFT version based on http://paulbourke.net/miscellaneous/dft, 1993
  const int    l = (shortTransform ? m_transfLengthS : m_transfLengthL) >> 1;
  const short* p = (shortTransform ? m_fftPermutS : m_fftPermutL); // look-up
  int l2 = 1, l3 = m_transfLengthL >> 1;

  if (iR == nullptr)
  {
    return; // null-pointer input error
  }

  // sort input with permutation look-up table
  if (iI != nullptr)
  {
    for (int i = l - 1; i >= 0; i--)
    {
      const int j = p[i];

      if (j > i) // swap input data at i and j
      {
        const int32_t iRTmp = iR[i]; // use re-
        const int32_t iITmp = iI[i]; // gisters

        iR[i] = iR[j];
        iI[i] = iI[j];
        iR[j] = iRTmp;
        iI[j] = iITmp;
      }
    }
  }
  else // real-valued input, no imaginary part
  {
    for (int i = l - 1; i >= 0; i--)
    {
      const int j = p[i];

      if (j > i) // swap real input at i and j
      {
        const int32_t iRTmp = iR[i];

        iR[i] = iR[j];
        iR[j] = iRTmp;
      }
      iI[i] = 0;  // zero imaginary input data
    }
  }

  // get length-l Fast Fourier Transform (FFT)
  for (int k = shortIntLog2 ((uint16_t) l) - 1; k >= 0; k--)
  {
    const int l1 = l2;
    l2 <<= 1;
    l3 >>= 1;

    for (int j = l1 - 1; j >= 0; j--)
    {
      const int       jTl3 = j * l3;
      const int64_t cosjl3 = m_fftHalfCos[jTl3]; // cos/sin
      const int64_t sinjl3 = m_fftHalfSin[jTl3]; // look-up

      for (int i = j; i < l; i += l2)
      {
        const int     iPl1 = i + l1;
        const int32_t rotR = int32_t ((cosjl3 * iR[iPl1] + sinjl3 * iI[iPl1] + LUT_OFFSET) >> LUT_SHIFT); // clockwise
        const int32_t rotI = int32_t ((cosjl3 * iI[iPl1] - sinjl3 * iR[iPl1] + LUT_OFFSET) >> LUT_SHIFT); // rotation

        iR[iPl1] = iR[i] + rotR;  iR[i] -= rotR;
        iI[iPl1] = iI[i] + rotI;  iI[i] -= rotI;
      }
    }
  }
}

void LappedTransform::windowAndFoldInL (const int32_t* inputL, const bool shortTransform, const bool kbdWindowL, const bool lowOverlapL,
                                        const bool mdstKernel, int32_t* const output)
{
  const unsigned ws = (kbdWindowL ? 1 : 0); // shape
  const int32_t* wl = (lowOverlapL ? m_timeWindowS[ws] : m_timeWindowL[ws]);
  const int Mo2     = (shortTransform ? m_transfLengthS : m_transfLengthL) >> 1;
  const int Mm1     = Mo2 * 2 - 1;
  const int Mo2m1   = Mo2 - 1;
  const int Mo2mO   = (lowOverlapL ? Mo2 - (m_transfLengthS >> 1) : 0);
  const int Mm1mO   = Mm1 - Mo2mO; // overlap offset
  int n;

  if (mdstKernel) // time-reversal and TDA sign flip
  {
    for (n = Mo2m1; n >= Mo2mO; n--) // windowed pt.
    {
      const int64_t i64 = (int64_t) inputL[Mm1 - n] * wl[Mm1mO - n] + (int64_t) inputL[n] * wl[n - Mo2mO];

      output[Mo2m1 - n] = int32_t ((i64 + WIN_OFFSET) >> WIN_SHIFT);
    }
    for (/*Mo2mO-1*/; n >= 0; n--) // unwindowed pt.
    {
      output[Mo2m1 - n] = (inputL[Mm1 - n] + 2) >> 2;
    }
  }
  else // MDCT kernel, no time-reversal or sign flip
  {
    for (n = Mo2m1; n >= Mo2mO; n--) // windowed pt.
    {
      const int64_t i64 = (int64_t) inputL[Mm1 - n] * wl[Mm1mO - n] - (int64_t) inputL[n] * wl[n - Mo2mO];

      output[Mo2 + n]   = int32_t ((i64 + WIN_OFFSET) >> WIN_SHIFT);
    }
    for (/*Mo2mO-1*/; n >= 0; n--) // unwindowed pt.
    {
      output[Mo2 + n]   = (inputL[Mm1 - n] + 2) >> 2;
    }
  }
}

void LappedTransform::windowAndFoldInR (const int32_t* inputR, const bool shortTransform, const bool kbdWindowR, const bool lowOverlapR,
                                        const bool mdstKernel, int32_t* const output)
{
  const unsigned ws = (kbdWindowR ? 1 : 0); // shape
  const int32_t* wr = (lowOverlapR ? m_timeWindowS[ws] : m_timeWindowL[ws]);
  const int Mo2     = (shortTransform ? m_transfLengthS : m_transfLengthL) >> 1;
  const int Mm1     = Mo2 * 2 - 1;
  const int Mo2m1   = Mo2 - 1;
  const int Mo2mO   = (lowOverlapR ? Mo2 - (m_transfLengthS >> 1) : 0);
  const int Mm1mO   = Mm1 - Mo2mO; // overlap offset
  int n;

  if (mdstKernel) // time-reversal and TDA sign flip
  {
    for (n = Mo2m1; n >= Mo2mO; n--) // windowed pt.
    {
      const int64_t i64 = (int64_t) inputR[n] * wr[Mm1mO - n] - (int64_t) inputR[Mm1 - n] * wr[n - Mo2mO];

      output[Mo2 + n]   = int32_t ((i64 + WIN_OFFSET) >> WIN_SHIFT);
    }
    for (/*Mo2mO-1*/; n >= 0; n--) // unwindowed pt.
    {
      output[Mo2 + n]   = (inputR[n] + 2) >> 2;
    }
  }
  else // MDCT kernel, no time-reversal or sign flip
  {
    for (n = Mo2m1; n >= Mo2mO; n--) // windowed pt.
    {
      const int64_t i64 = (int64_t) inputR[n] * wr[Mm1mO - n] + (int64_t) inputR[Mm1 - n] * wr[n - Mo2mO];

      output[Mo2m1 - n] = int32_t ((i64 + WIN_OFFSET) >> WIN_SHIFT);
    }
    for (/*Mo2mO-1*/; n >= 0; n--) // unwindowed pt.
    {
      output[Mo2m1 - n] = (inputR[n] + 2) >> 2;
    }
  }
}

// constructor
LappedTransform::LappedTransform ()
{
  // initialize all helper buffers
  m_dctRotCosL = nullptr;
  m_dctRotCosS = nullptr;
  m_dctRotSinL = nullptr;
  m_dctRotSinS = nullptr;
  m_fftHalfCos = nullptr;
  m_fftHalfSin = nullptr;
  m_fftPermutL = nullptr;
  m_fftPermutS = nullptr;
  m_tempIntBuf = nullptr;

  // initialize all window buffers
  for (short s = 0; s < 2; s++)
  {
    m_timeWindowL[s] = nullptr;
    m_timeWindowS[s] = nullptr;
  }
  m_transfLengthL = 0;
  m_transfLengthS = 0;
}

// destructor
LappedTransform::~LappedTransform ()
{
  // free allocated helper buffers
  MFREE (m_dctRotCosL);
  MFREE (m_dctRotCosS);
  MFREE (m_dctRotSinL);
  MFREE (m_dctRotSinS);
  MFREE (m_fftHalfCos);
  MFREE (m_fftHalfSin);
  MFREE (m_fftPermutL);
  MFREE (m_fftPermutS);
  m_tempIntBuf = nullptr;
}

// public functions
unsigned LappedTransform::applyNegDCT4 (int32_t* const signal, const bool shortTransform) // works in-place
{
  // int32 negative-DCT-IV version based on http://www.ee.columbia.edu/~marios/mdct/mdct_giraffe.html, 2003
  // NOTE: amplifies short-transform results (8 times shorter than long-transform results) by a factor of 8
  const int lm1   = (shortTransform ? m_transfLengthS : m_transfLengthL) - 1;
  const int lm1o2 = lm1 >> 1;
  const int32_t* rotatCos = (shortTransform ? m_dctRotCosS : m_dctRotCosL);
  const int32_t* rotatSin = (shortTransform ? m_dctRotSinS : m_dctRotSinL);
  const int64_t rotOffset = (shortTransform ? LUT_OFFSET>>3 : LUT_OFFSET);
  const int64_t  rotShift = (shortTransform ? LUT_SHIFT - 3 : LUT_SHIFT);
  int32_t* const tempReal = m_tempIntBuf;
  int32_t* const tempImag = &m_tempIntBuf[lm1o2 + 1];
  int i, i2;

  if (signal == nullptr)
  {
    return 1; // null-pointer input error
  }

  // resort, separate, pre-twiddle signal
  for (i = lm1o2, i2 = lm1 - 1; i >= 0; i--, i2 -= 2)
  {
    const int64_t c = rotatCos[i];
    const int64_t s = rotatSin[i];
    const int64_t e = signal[i2/*even*/];
    const int64_t o = signal[lm1 - i2];

    tempReal[i] = int32_t ((rotOffset + e * c - o * s) >> rotShift);
    tempImag[i] = int32_t ((rotOffset + o * c + e * s) >> rotShift);
  }

  applyHalfSizeFFT (tempReal, tempImag, shortTransform);

  // post-twiddle, combine, resort output
  for (i = lm1o2, i2 = lm1 - 1; i >= 0; i--, i2 -= 2)
  {
    const int64_t c = rotatCos[i];
    const int64_t s = rotatSin[i];
    const int64_t e = tempReal[i];
    const int64_t o = tempImag[i];

    signal[i2/*even*/] = int32_t ((LUT_OFFSET + o * s - e * c) >> LUT_SHIFT);
    signal[lm1 - i2]   = int32_t ((LUT_OFFSET + e * s + o * c) >> LUT_SHIFT);
  }

  return 0; // no error
}

unsigned LappedTransform::applyMCLT (const int32_t* timeSig, const bool eightTransforms, bool kbdWindowL, const bool kbdWindowR,
                                     const bool lowOverlapL, const bool lowOverlapR, int32_t* const outMdct, int32_t* const outMdst)
{
  if ((timeSig == nullptr) || (outMdct == nullptr) || (outMdst == nullptr))
  {
    return 1; // invalid arguments error
  }

  if (eightTransforms)  // short windows
  {
    const int32_t* tSigS = &timeSig[(m_transfLengthL - m_transfLengthS) >> 1];
    int32_t*    outMdctS = outMdct;
    int32_t*    outMdstS = outMdst;

    for (unsigned w = 0; w < 8; w++)
    {
      windowAndFoldInL (tSigS /*window half 1*/, true, kbdWindowL, lowOverlapL, false, outMdctS);
      windowAndFoldInR (&tSigS[m_transfLengthS], true, kbdWindowR, lowOverlapR, false, outMdctS);
      windowAndFoldInL (tSigS /*window half 1*/, true, kbdWindowL, lowOverlapL, true, outMdstS);
      windowAndFoldInR (&tSigS[m_transfLengthS], true, kbdWindowR, lowOverlapR, true, outMdstS);
      // MDCT via DCT-IV
      applyNegDCT4 (outMdctS, true);
      // MDST via DCT-IV
      applyNegDCT4 (outMdstS, true);
#if 1 // not needed here
      for (int i = m_transfLengthS - 2; i >= 0; i -= 2)
      {
        outMdstS[i] *= -1; // DCT to DST
      }
#endif
      kbdWindowL = kbdWindowR; // only first window uses last frame's shape
      tSigS     += m_transfLengthS;
      outMdctS  += m_transfLengthS;
      outMdstS  += m_transfLengthS;
    }
  }
  else  // 1 long window
  {
    windowAndFoldInL (timeSig /*window half 1*/, false, kbdWindowL, lowOverlapL, false, outMdct);
    windowAndFoldInR (&timeSig[m_transfLengthL], false, kbdWindowR, lowOverlapR, false, outMdct);
    windowAndFoldInL (timeSig /*window half 1*/, false, kbdWindowL, lowOverlapL, true, outMdst);
    windowAndFoldInR (&timeSig[m_transfLengthL], false, kbdWindowR, lowOverlapR, true, outMdst);
    // MDCT using DCT-IV
    applyNegDCT4 (outMdct, false);
    // MDST using DCT-IV
    applyNegDCT4 (outMdst, false);
#if 1 // not needed here
    for (int i = m_transfLengthL - 2; i >= 0; i -= 2)
    {
      outMdst[i] *= -1; // DCT to DST
    }
#endif
  }

  return 0; // no error
}

unsigned LappedTransform::initConstants (int32_t* const tempIntBuf, int32_t* const timeWindowL[2], int32_t* const timeWindowS[2],
                                         const unsigned maxTransfLength)
{
  const short  halfLength = short (maxTransfLength >> 1);
  const short  sixtLength = short (maxTransfLength >> 4);
  const double dNormL     = 3.141592653589793 / (2.0 * halfLength);
  const double dNormS     = 3.141592653589793 / (2.0 * sixtLength);
  const double dNormL4    = dNormL * 4.0;
  short s;

  if ((tempIntBuf == nullptr) || (timeWindowL == nullptr) || (timeWindowS == nullptr) ||
      (maxTransfLength < 128) || (maxTransfLength > 8192) || (maxTransfLength & (maxTransfLength - 1)))
  {
    return 1; // invalid arguments error
  }

  m_transfLengthL = 2 * halfLength;
  m_transfLengthS = 2 * sixtLength;

  if ((m_dctRotCosL = (int32_t*) malloc (halfLength * sizeof (int32_t))) == nullptr ||
      (m_dctRotCosS = (int32_t*) malloc (sixtLength * sizeof (int32_t))) == nullptr ||
      (m_dctRotSinL = (int32_t*) malloc (halfLength * sizeof (int32_t))) == nullptr ||
      (m_dctRotSinS = (int32_t*) malloc (sixtLength * sizeof (int32_t))) == nullptr ||
      (m_fftHalfCos = (int32_t*) malloc ((halfLength >> 1) * sizeof (int32_t))) == nullptr ||
      (m_fftHalfSin = (int32_t*) malloc ((halfLength >> 1) * sizeof (int32_t))) == nullptr ||
      (m_fftPermutL = createPermutTable (halfLength)) == nullptr ||
      (m_fftPermutS = createPermutTable (sixtLength)) == nullptr)
  {
    return 2; // memory allocation error
  }

  // obtain cosine and sine coefficients
  for (s = 0; s < halfLength; s++)
  {
    m_dctRotCosL[s] = int32_t (cos (dNormL * (s + 0.125)) * (INT_MAX + 1.0) + 0.5);
    m_dctRotSinL[s] = int32_t (sin (dNormL * (s + 0.125)) * INT_MIN - 0.5);
  }
  for (s = 0; s < sixtLength; s++)
  {
    m_dctRotCosS[s] = int32_t (cos (dNormS * (s + 0.125)) * (INT_MAX + 1.0) + 0.5);
    m_dctRotSinS[s] = int32_t (sin (dNormS * (s + 0.125)) * INT_MIN - 0.5);
  }

  for (s = 0; s < m_transfLengthS; s++)
  {
    m_fftHalfSin[s] = int32_t (sin (dNormL4 * s) * INT_MIN - 0.5);
    m_fftHalfCos[m_transfLengthS + s] = -m_fftHalfSin[s];
  }
  // complete missing entries by copying
  m_fftHalfSin[s] = INT_MIN;
  m_fftHalfCos[0] = INT_MIN;
  for (s = 1; s < m_transfLengthS; s++)
  {
    m_fftHalfSin[m_transfLengthS + s] = m_fftHalfSin[m_transfLengthS - s];
    m_fftHalfCos[m_transfLengthS - s] = m_fftHalfSin[s];
  }

  // adopt helper/window buffer pointers
  m_tempIntBuf = tempIntBuf;
  for (s = 0; s < 2; s++)
  {
    m_timeWindowL[s] = timeWindowL[s];
    m_timeWindowS[s] = timeWindowS[s];
  }

  return 0; // no error
}
