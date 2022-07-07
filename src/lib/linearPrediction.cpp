/* linearPrediction.cpp - source file for class providing linear prediction capability
 * written by C. R. Helmrich, last modified in 2021 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under the exhale Copyright License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2021 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#include "exhaleLibPch.h"
#include "linearPrediction.h"

// reconstructed 3-bit TNS coefficients
static const short tnsQuantCoeff3[9 /*2^3+1*/] = { // = round (2^11 * sin (x * pi / (x < 0 ? 9 : 7)))
  -2017, -1774, -1316, -700,  0,  889,  1601,  1997,  1997
};

// reconstructed 4-bit TNS coefficients
static const short tnsQuantCoeff4[17/*2^4+1*/] = { // = round (2^11 * sin (x * pi / (x < 0 ? 17 : 15)))
  -2039, -1970, -1833, -1634, -1380, -1078, -740, -376,  0,  426,  833,  1204,  1522,  1774,  1948,  2037,  2037
};

static const short* tnsQuantCoeff[2/*coefRes*/] = {tnsQuantCoeff3, tnsQuantCoeff4};

// ISO/IEC 14496-3, Sec. 4.6.9.3, 3-bit
static const int8_t tnsQuantIndex3[SCHAR_MAX + 1] = { // = round (asin (x / 64) * (x < 0 ? 9 : 7) / pi)
  -4, -4, -4, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2,
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
   1,  1,  1,  1,  1,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  3,  3,  3,  3,  3,  3,  3
};

// ISO/IEC 14496-3, Sec. 4.6.9.3, 4-bit
static const int8_t tnsQuantIndex4[SCHAR_MAX + 1] = { // = round (asin (x / 64) * (x < 0 ? 17 : 15) / pi)
  -8, -7, -7, -7, -6, -6, -6, -6, -6, -5, -5, -5, -5, -5, -5, -5, -4, -4, -4, -4, -4, -4, -4, -4, -4, -3, -3, -3, -3, -3, -3, -3,
  -3, -3, -3, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  3,
   3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  4,  4,  4,  4,  4,  4,  4,  4,  4,  5,  5,  5,  5,  5,  5,  5,  6,  6,  6,  6,  7,  7
};

static const int8_t* tnsQuantIndex[2/*coefRes*/] = {tnsQuantIndex3, tnsQuantIndex4};

// static helper functions
static int quantizeParCorCoeffs (const short* const parCorCoeffs, const uint16_t nCoeffs, const short bitDepth, int8_t* const quantCoeffs,
                                 const bool lowRes)
{
  const short   bitShift = bitDepth - 7;
  const unsigned  tabIdx = (lowRes ? 0 : 1);
  const int8_t tabOffset = 4 << tabIdx;
  const short*  coeffTab = tnsQuantCoeff[tabIdx];
  const int8_t* indexTab = tnsQuantIndex[tabIdx];
  int dist0, dist1, distTotal = 0;

  for (uint16_t s = 0; s < nCoeffs; s++)
  {
    const short   coeff = (bitShift < 0 ? parCorCoeffs[s] << -bitShift : parCorCoeffs[s] >> bitShift);
    const int8_t coeff1 = indexTab[coeff + 1 + (SCHAR_MAX >> 1)];
    const int8_t coeff0 = (coeff1 <= -tabOffset ? coeff1 : coeff1 - 1);

    dist0 = (int) coeffTab[coeff0 + tabOffset] - parCorCoeffs[s];
    dist0 *= dist0;
    dist1 = (int) coeffTab[coeff1 + tabOffset] - parCorCoeffs[s];
    dist1 *= dist1;
    if (dist0 < dist1) // use down-rounded coeff
    {
      quantCoeffs[s] = coeff0;
      distTotal += dist0;
    }
    else // dist0 >= dist1, use up-rounded coeff
    {
      quantCoeffs[s] = ((dist0 == dist1) && (abs (coeff0) < abs (coeff1)) ? coeff0 : coeff1);
      distTotal += dist1;
    }
  } // for s

  return distTotal;  // total quantization error
}

// constructor
LinearPredictor::LinearPredictor ()
{
  // initialize member variables
  memset (m_tempBuf, 0, 2 * MAX_PREDICTION_ORDER * sizeof (int64_t));
}

// public functions
uint32_t LinearPredictor::calcParCorCoeffs (const int32_t* const anaSignal, const uint16_t nAnaSamples, const uint16_t nCoeffs,
                                            short* const parCorCoeffs)  // returns 256 - 256 / prediction gain per filter order, or 0
{
  // int64 version of algorithm in Figure 1 of J. LeRoux and C. Gueguen, "A Fixed Point Computation
  // of Partial Correlation Coefficients," IEEE Trans. ASSP, vol. 27, no. 3, pp. 257-259, June 1977.
  int64_t* const acf = m_tempBuf; // correlation
  uint32_t pg[MAX_PREDICTION_ORDER] = {0, 0, 0, 0};
  int64_t  pgDen, pgOff; // for prediction gains
  short s = nAnaSamples - 1;

  if ((anaSignal == nullptr) || (parCorCoeffs == nullptr) || (nCoeffs == 0) || (nCoeffs > MAX_PREDICTION_ORDER) || (nAnaSamples <= nCoeffs))
  {
    return 0;
  }

  if (nCoeffs >= 2) // high predictor order, 2-4
  {
    int64_t* const EN = &m_tempBuf[0];
    int64_t* const EP = &m_tempBuf[MAX_PREDICTION_ORDER];
    int64_t sampleHO = anaSignal[s];

    // calculate autocorrelation function values
    acf[0]  = sampleHO * sampleHO;
    sampleHO = anaSignal[--s];
    acf[0] += sampleHO * sampleHO;
    acf[1]  = sampleHO * anaSignal[s + 1];
    sampleHO = anaSignal[--s];
    acf[0] += sampleHO * sampleHO;
    acf[1] += sampleHO * anaSignal[s + 1];
    acf[2]  = sampleHO * anaSignal[s + 2];

    if (nCoeffs == 2)
    {
      for (s--; s >= 0; s--)
      {
        const int64_t sample = anaSignal[s];

        acf[0] += sample * sample;
        acf[1] += sample * anaSignal[s + 1];
        acf[2] += sample * anaSignal[s + 2];
      }
    }
    else // order 3-4
    {
      sampleHO = anaSignal[--s];
      acf[0] += sampleHO * sampleHO;
      acf[1] += sampleHO * anaSignal[s + 1];
      acf[2] += sampleHO * anaSignal[s + 2];
      acf[3]  = sampleHO * anaSignal[s + 3];

      if (nCoeffs == 3)
      {
        for (s--; s >= 0; s--)
        {
          const int64_t sample = anaSignal[s];

          acf[0] += sample * sample;
          acf[1] += sample * anaSignal[s + 1];
          acf[2] += sample * anaSignal[s + 2];
          acf[3] += sample * anaSignal[s + 3];
        }
      }
      else // order 4
      {
        acf[4] = 0;
        for (s--; s >= 0; s--)
        {
          const int64_t sample = anaSignal[s];

          acf[0] += sample * sample;
          acf[1] += sample * anaSignal[s + 1];
          acf[2] += sample * anaSignal[s + 2];
          acf[3] += sample * anaSignal[s + 3];
          acf[4] += sample * anaSignal[s + 4];
        }
      }
    }

    // reduce correlation value range to <32 bit
    acf[0] = (acf[0] - INT_MIN/*eps*/) >> 31;

    EP[3] = (acf[4] - (INT_MIN >> 1)) >> 31;
    EN[3] = (acf[3] - (INT_MIN >> 1)) >> 31;
    EP[2] = EN[3];
    EN[2] = (acf[2] - (INT_MIN >> 1)) >> 31;
    EP[1] = EN[2];
    EN[1] = (acf[1] - (INT_MIN >> 1)) >> 31;
    EP[0] = EN[1];  // finish EP buffer creation
    EN[0] = acf[0]; // finish EN buffer creation

    pgOff = acf[0];
    pgDen = pgOff << 1;

    for (s = 0; s < (short) nCoeffs; s++)
    {
      uint16_t p;
      // ParCor coefficient Ks & prediction gain
      int64_t Ks = (EP[s] * -512) / (EN[0] + LP_EPS);

      Ks = CLIP_PM (Ks, 511);  // enforce |Ks|<1
      parCorCoeffs[s] = (short) Ks;

      for (p = s; p < nCoeffs; p++)
      {
        sampleHO  =  EN[p - s]; // use register?
        EN[p - s] = (EN[p - s] * (1 << 9)) + EP[p] * Ks;
        EP[p]     = (EP[p] * (1 << 9))  + sampleHO * Ks;
      }
      if (s > 0 && EN[0] < 0)  // EN wrap-around
      {
        pg[s] = pg[s - 1]; // "invent" some gain
      }
      else
      {
        sampleHO = (1 << (s * 9)) >> 1;// offset
        pg[s] = (pgOff <= LP_EPS ? 0 : (1 << 8) - uint32_t ((((EN[0] + sampleHO) >> (s * 9)) + pgOff) / pgDen));
      }
    } // for s
  }
  else  // nCoeffs == 1, minimum predictor order
  {
    int64_t sampleMO = anaSignal[s];

    acf[0] = sampleMO * sampleMO;
    acf[1] = 0;
    for (s--; s >= 0; s--)
    {
      const int64_t sample = anaSignal[s];

      acf[0] += sample * sample;
      acf[1] += sample * anaSignal[s + 1];
    }

    // reduce correlation value range to <32 bit
    acf[0] = (acf[0] - INT_MIN/*eps*/) >> 31;
    acf[1] = (acf[1] - (INT_MIN >> 1)) >> 31;

    pgOff = acf[0];
    pgDen = pgOff << 1;

    // 1st-order coefficient K & prediction gain
    sampleMO = (acf[1] * -512) / (acf[0] + LP_EPS);

    sampleMO = CLIP_PM (sampleMO, 511); // |K|<1
    parCorCoeffs[0] = (short) sampleMO;

    sampleMO = (acf[0] << 9) + acf[1] * sampleMO;
    pg[0] = (pgOff <= LP_EPS ? 0 : (1 << 8) - uint32_t ((sampleMO + pgOff) / pgDen));
  }

  return (CLIP_UCHAR (pg[3]) << 24) | (CLIP_UCHAR (pg[2]) << 16) | (CLIP_UCHAR (pg[1]) << 8) | CLIP_UCHAR (pg[0]);
}

uint8_t LinearPredictor::calcOptTnsCoeffs (short* const parCorCoeffs, int8_t* const quantCoeffs, bool* const lowCoeffRes,
                                           const uint16_t maxOrder, const uint8_t predGain, const uint8_t tonality /*= 0*/,
                                           const uint16_t parCorCoeffBitDepth /*= 10*/) // returns optimized filter order for TNS
{
  const short bitShift = LP_DEPTH - (short) parCorCoeffBitDepth;
  short       shortBuf[MAX_PREDICTION_ORDER];
  uint16_t s, order = __min (maxOrder, MAX_PREDICTION_ORDER);
  int d, i;

  if ((parCorCoeffs == nullptr) || (quantCoeffs == nullptr) || (maxOrder == 0) || (maxOrder > MAX_PREDICTION_ORDER) || (parCorCoeffBitDepth < 2) || (bitShift < 0))
  {
    if (quantCoeffs) memset (quantCoeffs, 0, order * sizeof (int8_t));

    return 0;   // invalid input arguments error
  }

  // determine direct-form filter damping factor
  parCorCoeffs[0] *= 1 << bitShift;
  i = abs (parCorCoeffs[0]);
  for (s = 1; s < order; s++)
  {
    parCorCoeffs[s] *= 1 << bitShift; // scale coeff
    i = __max (i, abs (parCorCoeffs[s]));
  }
  for (/*s*/; s < MAX_PREDICTION_ORDER; s++)
  {
    parCorCoeffs[s] = 0;  // similarParCorCoeffs
  }

  if (predGain < 41 + (tonality >> 3)) // 1.5 dB
  {
    memset (quantCoeffs, 0, order * sizeof (int8_t));

    return 0;  // LPC prediction gain is too low
  }

  d = (3 << (LP_DEPTH - 1)) >> 2;

  if (i > d) // apply direct-form filter damping
  {
    i = ((i >> 1) + d * (1 << LP_SHIFT)) / i;
    d = i;
    if (parCorToLpCoeffs (parCorCoeffs, order, shortBuf, LP_DEPTH) > 0)
    {
      return 0;  // coefficient conversion error
    }
    for (s = 0; s < order; s++)
    {
      shortBuf[s] = short ((LP_OFFSET + d * shortBuf[s]) >> LP_SHIFT);
      d = (LP_OFFSET + d * i) >> LP_SHIFT;
    }
    // verify order and go back to ParCor domain
    if ((s > 0) && (shortBuf[s - 1] == 0))
    {
      order--;
    }
    if (lpToParCorCoeffs (shortBuf, order, parCorCoeffs, LP_DEPTH) > 0)
    {
      return 0;  // coefficient conversion error
    }
  }

  // high-res quantizer, obtain coeff distortion
  d = quantizeParCorCoeffs (parCorCoeffs, order, LP_DEPTH, quantCoeffs, false);

  if ((lowCoeffRes != nullptr) && (quantizeParCorCoeffs (parCorCoeffs, order, LP_DEPTH, (int8_t* const) m_tempBuf, true) < d))
  {
    // low-res quantizer yields lower distortion
    *lowCoeffRes = true;
    memcpy (quantCoeffs, m_tempBuf, order * sizeof (int8_t));
  }

  for (; order > 0; order--) // return opt order
  {
    if (quantCoeffs[order - 1] != 0) return (uint8_t) order;
  }
  return 0;
}

unsigned LinearPredictor::lpToParCorCoeffs (/*mod!*/short* const lpCoeffs, const uint16_t nCoeffs, short* const parCorCoeffs,
                                            const uint16_t parCorCoeffBitDepth /*= 10*/)
{
  int* const intBuf = (int* const) m_tempBuf;
  const int  shift  = parCorCoeffBitDepth - 1;
  const int  offset = 1 << (shift - 1);

  if ((lpCoeffs == nullptr) || (parCorCoeffs == nullptr) || (nCoeffs == 0) || (nCoeffs > MAX_PREDICTION_ORDER) || (parCorCoeffBitDepth < 2))
  {
    return 1;  // error
  }

  for (uint16_t p, s = nCoeffs - 1; s > 0; s--)
  {
    const int i = (parCorCoeffs[s] = lpCoeffs[s]);
    const int d = (1 << shift) - ((offset + i * i) >> shift);
    const int o = d >> 1;

    if (d <= 0) return s; // invalid coefficient

    for (p = 0; p < s; p++)
    {
      intBuf[p] = lpCoeffs[s - 1 - p];
    }
    for (p = 0; p < s; p++)
    {
      lpCoeffs[p] = short ((o + ((int) lpCoeffs[p] * (1 << shift)) - intBuf[p] * i) / d);
    }
  }
  parCorCoeffs[0] = lpCoeffs[0];

  return 0; // no error
}

unsigned LinearPredictor::parCorToLpCoeffs (const short* const parCorCoeffs, const uint16_t nCoeffs, short* const lpCoeffs,
                                            const uint16_t parCorCoeffBitDepth /*= 10*/)
{
  int* const intBuf = (int* const) m_tempBuf;
  const int  shift  = parCorCoeffBitDepth - 1;
  const int  offset = 1 << (shift - 1);

  if ((parCorCoeffs == nullptr) || (lpCoeffs == nullptr) || (nCoeffs == 0) || (nCoeffs > MAX_PREDICTION_ORDER) || (parCorCoeffBitDepth < 2))
  {
    return 1;  // error
  }

  lpCoeffs[0] = parCorCoeffs[0];
  for (uint16_t p, s = 1; s < nCoeffs; s++)
  {
    const int i = (lpCoeffs[s] = parCorCoeffs[s]);

    if (abs (i) > (1 << shift)) return s; // > 1

    for (p = 0; p < s; p++)
    {
      intBuf[p] = lpCoeffs[s - 1 - p];
    }
    for (p = 0; p < s; p++)
    {
      lpCoeffs[p] += short ((offset + intBuf[p] * i) >> shift);
    }
  }
  return 0; // no error
}

unsigned LinearPredictor::quantTnsToLpCoeffs (const int8_t* const quantTnsCoeffs, const uint16_t nCoeffs, const bool lowCoeffRes,
                                              short* const parCorCoeffs, short* const lpCoeffs)
{
  const unsigned  tabIdx = (lowCoeffRes ? 0 : 1);
  const int8_t tabOffset = 4 << tabIdx;
  const short*  coeffTab = tnsQuantCoeff[tabIdx];

  if ((quantTnsCoeffs == nullptr) || (parCorCoeffs == nullptr) || (lpCoeffs == nullptr) || (nCoeffs == 0) || (nCoeffs > MAX_PREDICTION_ORDER))
  {
    return 1;  // error
  }

  for (uint16_t s = 0; s < nCoeffs; s++)
  {
    parCorCoeffs[s] = coeffTab[CLIP_PM (quantTnsCoeffs[s], tabOffset) + tabOffset];
  }

  return parCorToLpCoeffs (parCorCoeffs, nCoeffs, lpCoeffs, LP_DEPTH);
}

bool LinearPredictor::similarParCorCoeffs (const short* const parCorCoeffs1, const short* const parCorCoeffs2, const uint16_t nCoeffs,
                                           const uint16_t parCorCoeffBitDepth /*= 10*/)
{
  unsigned sumAbsDiff = 0;

  if ((parCorCoeffs1 == nullptr) || (parCorCoeffs2 == nullptr) || (nCoeffs == 0) || (nCoeffs > MAX_PREDICTION_ORDER) || (parCorCoeffBitDepth < 2))
  {
    return false; // error
  }

  for (uint16_t s = 0; s < nCoeffs; s++)
  {
    sumAbsDiff += abs (parCorCoeffs1[s] - parCorCoeffs2[s]);
  }

  return (sumAbsDiff + 12u * nCoeffs < ((4u * nCoeffs) << (parCorCoeffBitDepth >> 1)));
}
