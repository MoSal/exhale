/* quantization.cpp - source file for class with nonuniform quantization functionality
 * written by C. R. Helmrich, last modified in 2020 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under a Modified BSD-Style License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2020 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#include "exhaleLibPch.h"
#include "quantization.h"

#define EC_TRAIN (0 && EC_TRELLIS_OPT_CODING) // for RDOC testing

// static helper functions
static inline short getBitCount (EntropyCoder& entrCoder, const int sfIndex, const int sfIndexPred,
                                 const uint8_t groupLength, const uint8_t* coeffQuant,
                                 const uint16_t coeffOffset, const uint16_t numCoeffs)
{
  unsigned bitCount = (sfIndex != UCHAR_MAX && sfIndexPred == UCHAR_MAX ? 8 : entrCoder.indexGetBitCount (sfIndex - sfIndexPred));

  if (groupLength == 1) // include arithmetic coding in bit count
  {
#if EC_TRELLIS_OPT_CODING
    const unsigned bitsStart = (entrCoder.arithGetCtxState () >> 17) & 31;
#endif
    bitCount += entrCoder.arithCodeSigMagn (coeffQuant, coeffOffset, numCoeffs);
#if EC_TRELLIS_OPT_CODING
    bitCount += (entrCoder.arithGetCtxState () >> 17) & 31;
    bitCount -= __min (bitsStart, bitCount); // +new-old m_acBits
#endif
  }

  return (short) __min (SHRT_MAX, bitCount); // exclude sign bits
}

#if EC_TRELLIS_OPT_CODING && !EC_TRAIN
static inline double getLagrangeValue (const uint16_t rateIndex) // RD optimization constant
{
  return (95.0 + rateIndex * rateIndex) * 0.0009765625; // / 1024
}
#endif

// private helper functions
double SfbQuantizer::getQuantDist (const unsigned* const coeffMagn, const uint8_t scaleFactor,
                                   const uint8_t* const coeffQuant, const uint16_t numCoeffs)
{
  const double stepSizeDiv = m_lutSfNorm[scaleFactor];
  double dDist = 0.0;

  for (int i = numCoeffs - 1; i >= 0; i--)
  {
    const double d = m_lutXExp43[coeffQuant[i]] - coeffMagn[i] * stepSizeDiv;

    dDist += d * d; // TODO: do this in fixed-point and with SIMD
  }

  // consider quantization step-size in calculation of distortion
  return dDist * m_lut2ExpX4[scaleFactor] * m_lut2ExpX4[scaleFactor];
}

uint8_t SfbQuantizer::quantizeMagnSfb (const unsigned* const coeffMagn, const uint8_t scaleFactor,
                                      /*mod*/uint8_t* const coeffQuant, const uint16_t numCoeffs,
#if EC_TRELLIS_OPT_CODING
                                       EntropyCoder* const arithmCoder, const uint16_t coeffOffset,
#endif
                                       short* const sigMaxQ /*= nullptr*/, short* const sigNumQ /*= nullptr*/)
{
  const double stepSizeDiv = m_lutSfNorm[scaleFactor];
  double  dNum = 0.0, dDen = 0.0;
  short sf, maxQ = 0, numQ = 0;

  for (int i = numCoeffs - 1; i >= 0; i--)
  {
    const double normalizedMagn = (double) coeffMagn[i] * stepSizeDiv;
#if SFB_QUANT_FAST_POW
    short q;

    if (normalizedMagn < 28.5)  // fast approximate pow (d, 0.75)
    {
      // based on code from: N. N. Schraudolph, "A Fast, Compact Approximation of the Expo-
      // nential Function," Neural Comput., vol. 11, pp. 853-862, 1998 and M. Ankerl, 2007,
      // https://martin.ankerl.com/2007/10/04/optimized-pow-approximation-for-java-and-c-c/
      union { double d; int32_t i[2]; } u = { normalizedMagn };

      u.i[1] = int32_t (0.75 * (u.i[1] - 1072632447) + 1072632447.0);
      u.i[0] = 0;
      q = short (u.d + (u.d < 1.0 ? 0.3822484 : 0.2734375));
    }
    else
    {
      q = short (SFB_QUANT_OFFSET + pow (normalizedMagn, 0.75));
    }
#else
    short q = short (SFB_QUANT_OFFSET + pow (normalizedMagn, 0.75));
#endif

    if (q > 0)
    {
      if (q >= SCHAR_MAX)
      {
        if (maxQ < q)
        {
          maxQ = q; // find maximum quantized magnitude in vector
        }
        q = SCHAR_MAX;
      }
      else
      {
        const double diffRoundD = m_lutXExp43[q    ] - normalizedMagn;
        const double diffRoundU = m_lutXExp43[q + 1] - normalizedMagn;

        if (diffRoundU * diffRoundU < diffRoundD * diffRoundD)
        {
          q++; // round-up gives lower distortion than round-down
        }
      }
      if (maxQ < q)
      {
        maxQ = q;
      }
      numQ++;
      dNum += m_lutXExp43[q] * normalizedMagn;
      dDen += m_lutXExp43[q] * m_lutXExp43[q];
    }
#if SFB_QUANT_PERCEPT_OPT
    else // q == 0, assume perceptual transparency for code below
    {
      dNum += normalizedMagn * normalizedMagn;
      dDen += normalizedMagn * normalizedMagn;
    }
#endif
    coeffQuant[i] = (uint8_t) q;
  }

  if (sigMaxQ) *sigMaxQ = maxQ; // max. quantized value magnitude
  if (sigNumQ) *sigNumQ = numQ; // nonzero coeff. count (L0 norm)

  sf = scaleFactor;
  // compute least-squares optimal modifier added to scale factor
  if (dNum > SF_THRESH_POS * dDen) sf++;
  else
  if (dNum < SF_THRESH_NEG * dDen) sf--;

#if EC_TRELLIS_OPT_CODING
  if (arithmCoder && (sf > 0) && (maxQ <= SCHAR_MAX)) // use RDOC
  {
    EntropyCoder& entrCoder = *arithmCoder;
#if EC_TRAIN
    const uint32_t codStart = entrCoder.arithGetCodState ();
    const uint32_t ctxStart = entrCoder.arithGetCtxState ();
    uint32_t bitCount = entrCoder.arithCodeSigMagn (&coeffQuant[-((int) coeffOffset)], coeffOffset, numCoeffs);

    bitCount += (entrCoder.arithGetCtxState () >> 17) & 31; // refinement: +new-old m_acBits
    bitCount -= __min ((ctxStart >> 17) & 31, bitCount);
    bitCount += (uint32_t) numQ;  // add sign bits for completion

    entrCoder.arithSetCodState (codStart);  // back to last state
    entrCoder.arithSetCtxState (ctxStart);
#else
    uint32_t bitCount = (uint32_t) numQ;
#endif
    if ((bitCount = quantizeMagnRDOC (entrCoder, (uint8_t) sf, bitCount, coeffOffset, coeffMagn, numCoeffs, coeffQuant)) > 0)
    {
      numQ = bitCount & SHRT_MAX;

      if ((numQ > 0) && (sf < m_maxSfIndex)) // nonzero-quantized
      {
        const double magnNormDiv = m_lutSfNorm[sf];

        dNum = dDen = 0.0;
        for (int i = numCoeffs - 1; i >= 0; i--)
        {
          const double normalizedMagn = (double) coeffMagn[i] * magnNormDiv;
          const uint8_t q = coeffQuant[i];

          if (q > 0)
          {
            dNum += m_lutXExp43[q] * normalizedMagn;
            dDen += m_lutXExp43[q] * m_lutXExp43[q];
          }
#if SFB_QUANT_PERCEPT_OPT
          else   // assume perceptual transparency for code below
          {
            dNum += normalizedMagn * normalizedMagn;
            dDen += normalizedMagn * normalizedMagn;
          }
#endif
        }

        // re-compute least-squares optimal scale factor modifier
        if (dNum > SF_THRESH_POS * dDen) sf++;
#if !SFB_QUANT_PERCEPT_OPT
        else
        if (dNum < SF_THRESH_NEG * dDen) sf--; // reduces SFB RMS
#endif
      } // if nonzero

      if (sigMaxQ) *sigMaxQ = (numQ > 0 ? maxQ : 0); // a new max
      if (sigNumQ) *sigNumQ = numQ; // a new nonzero coeff. count
    }
  }
#endif // EC_TRELLIS_OPT_CODING

#if SFB_QUANT_PERCEPT_OPT
  if ((numQ > 0) && (sf > 0 && sf <= scaleFactor)) // recover RMS
  {
    const double magnNormDiv = m_lutSfNorm[sf];

    dNum = 0.0;  // dDen has normalized energy after quantization
    for (int i = numCoeffs - 1; i >= 0; i--)
    {
      const double normalizedMagn = (double) coeffMagn[i] * magnNormDiv;

      dNum += normalizedMagn * normalizedMagn; // TODO: use SIMD?
    }

    if (dNum > SF_THRESH_POS * SF_THRESH_POS * dDen) sf++;
  }
#endif
  return (uint8_t) __max (0, sf); // optimized scale factor index
}

#if EC_TRELLIS_OPT_CODING
uint32_t SfbQuantizer::quantizeMagnRDOC (EntropyCoder& entropyCoder, const uint8_t optimalSf, const unsigned targetBitCount,
                                         const uint16_t coeffOffset, const unsigned* const coeffMagn, // initial MDCT magnitudes
                                         const uint16_t numCoeffs, uint8_t* const quantCoeffs) // returns updated SFB statistics
{
  // numTuples: num of trellis stages. Based on: A. Aggarwal, S. L. Regunathan, and K. Rose,
  // "Trellis-Based Optimization of MPEG-4 Advanced Audio Coding," in Proc. IEEE Workshop on
  // Speech Coding, pp. 142-144, Sep. 2000. Modified for arithmetic instead of Huffman coder
  const uint32_t  codStart = entropyCoder.arithGetCodState ();
  const uint32_t  ctxStart = entropyCoder.arithGetCtxState (); // before call to getBitCount
  const double stepSizeDiv = m_lutSfNorm[optimalSf];
  const uint16_t numStates = 4; // 4 reduction types: [0, 0], [0, -1], [-1, 0], and [-1, -1]
  const uint16_t numTuples = numCoeffs >> 1;
  uint8_t* const quantRate = &m_coeffTemp[((unsigned) m_maxSize8M1 + 1) << 3];
  uint32_t prevCodState[4] = {0, 0, 0, 0};
  uint32_t prevCtxState[4] = {0, 0, 0, 0};
  double   prevVtrbCost[4] = {0, 0, 0, 0};
  uint32_t tempCodState[4] = {0, 0, 0, 0};
  uint32_t tempCtxState[4] = {0, 0, 0, 0};
  double   tempVtrbCost[4] = {0, 0, 0, 0};
  double   quantDist[32][4];   // TODO: dynamic memory allocation
  uint8_t* const optimalIs = (uint8_t* const) (quantDist[32-1]);
  uint8_t  tempQuant[4], numQ; // for tuple/SFB sign bit counting
  unsigned tempBitCount, tuple, is;
  int ds;
#if EC_TRAIN
  double refSfbDist = 0.0, tempSfbDist = 0.0;
#else
  const double lambda = getLagrangeValue (m_rateIndex);
#endif

  if ((coeffMagn == nullptr) || (quantCoeffs == nullptr) || (optimalSf > m_maxSfIndex) || (numTuples == 0) || (numTuples > 32) ||
      (targetBitCount == 0)  || (targetBitCount > SHRT_MAX))
  {
    return 0; // invalid input error
  }

  // save third-last tuple value, required due to an insufficiency of arithGet/SetCtxState()
  if (coeffOffset > 5) tempQuant[3] = entropyCoder.arithGetTuplePtr ()[(coeffOffset >> 1) - 3];

  for (tuple = 0; tuple < numTuples; tuple++) // tuple-wise non-weighted distortion and rate
  {
    const uint16_t  tupleStart = tuple << 1;
    const uint16_t tupleOffset = coeffOffset + tupleStart;
    const double   normalMagnA = (double) coeffMagn[tupleStart    ] * stepSizeDiv;
    const double   normalMagnB = (double) coeffMagn[tupleStart + 1] * stepSizeDiv;
    uint8_t  coeffQuantA = quantCoeffs[tupleStart];
    uint8_t  coeffQuantB = quantCoeffs[tupleStart + 1];

    for (is = 0; is < numStates; is++)  // populate tuple trellis
    {
      uint8_t* const mag = (is != 0 ? tempQuant : quantCoeffs) - (int) tupleOffset; // see arithCodeSigMagn()
      uint8_t*  currRate = &quantRate[(is + tuple * numStates) * numStates];
      double diffA, diffB;

      if (is != 0) // test reduction of quantized MDCT magnitudes
      {
        const uint8_t redA = is >> 1;
        const uint8_t redB = is &  1;

        if ((redA > 0 && coeffQuantA != 1) || (redB > 0 && coeffQuantB != 1))  // avoid path
        {
          tempCodState[is] = tempCodState[0];
          tempCtxState[is] = tempCtxState[0];
          for (ds = numStates - 1; ds >= 0; ds--)
          {
            currRate[ds] = UCHAR_MAX;
          }
          continue;
        }
        tempQuant[0] = (coeffQuantA -= redA);
        tempQuant[1] = (coeffQuantB -= redB);
      }
      diffA = m_lutXExp43[coeffQuantA] - normalMagnA;
      diffB = m_lutXExp43[coeffQuantB] - normalMagnB;
      quantDist[tuple][is] = diffA * diffA + diffB * diffB;

      numQ  = (coeffQuantA > 0 ? 1 : 0) + (coeffQuantB > 0 ? 1 : 0);

      if (tuple == 0) // first tuple, with tupleStart == sfbStart
      {
        entropyCoder.arithSetCodState (codStart); // start of SFB
        entropyCoder.arithSetCtxState (ctxStart);
        tempBitCount = entropyCoder.arithCodeSigMagn (mag, tupleOffset, 2);

        tempBitCount += (entropyCoder.arithGetCtxState () >> 17) & 31;  // +new-old m_acBits
        tempBitCount -= __min ((ctxStart >> 17) & 31, tempBitCount);
        tempBitCount += numQ; // add sign bits to finish estimate

        for (ds = numStates - 1; ds >= 0; ds--)
        {
          currRate[ds] = (uint8_t) __min (UCHAR_MAX, tempBitCount);
        }
      }
      else // tuple > 0, rate depends on decisions for last tuple
      {
        for (ds = numStates - 1; ds >= 0; ds--)
        {
          if (quantRate[(ds + (tuple-1) * numStates) * numStates] >= UCHAR_MAX)// avoid path
          {
            currRate[ds] = UCHAR_MAX;

            continue;
          }

          entropyCoder.arithSetCodState (prevCodState[ds]);
          entropyCoder.arithSetCtxState (prevCtxState[ds], tupleOffset);
          tempBitCount = entropyCoder.arithCodeSigMagn (mag, tupleOffset, 2);

          tempBitCount += (entropyCoder.arithGetCtxState () >> 17) & 31;// +new-old m_acBits
          tempBitCount -= __min ((prevCtxState[ds] >> 17) & 31, tempBitCount);
          tempBitCount += numQ; // + sign bits to finish estimate

          currRate[ds] = (uint8_t) __min (UCHAR_MAX, tempBitCount);
        }
      }
      // statistically best place to save states is after ds == 0
      tempCodState[is] = entropyCoder.arithGetCodState ();
      tempCtxState[is] = entropyCoder.arithGetCtxState ();
    } // for is
#if EC_TRAIN
    refSfbDist += quantDist[tuple][0];
#endif
    memcpy (prevCodState, tempCodState, numStates * sizeof (uint32_t));
    memcpy (prevCtxState, tempCtxState, numStates * sizeof (uint32_t));
  } // for tuple

  entropyCoder.arithSetCodState (codStart); // back to last state
  entropyCoder.arithSetCtxState (ctxStart, coeffOffset);
  // restore third-last tuple value, see insufficiency note above
  if (coeffOffset > 5) entropyCoder.arithGetTuplePtr ()[(coeffOffset >> 1) - 3] = tempQuant[3];

#if EC_TRAIN
  tempBitCount = targetBitCount + 1; // Viterbi search for minimum distortion at target rate
  for (double lambda = 0.015625; (lambda <= 0.375) && (tempBitCount > targetBitCount); lambda += 0.0078125)
#endif
  {
    double* const  prevCost = prevVtrbCost;
#if !EC_TRAIN
    uint8_t* const prevPath = (uint8_t*) quantDist;// backtracker
#endif
    double   costMinIs = (double) UINT_MAX;
    unsigned pathMinIs = 0;
#if EC_TRAIN
    uint8_t prevPath[16*4];
    tempSfbDist = 0.0;
#endif

    for (is = 0; is < numStates; is++) // initialize minimum path
    {
      const uint8_t  currRate = quantRate[is * numStates];

      prevCost[is] = (currRate >= UCHAR_MAX ? (double) UINT_MAX : lambda * currRate + quantDist[0][is]);
      prevPath[is] = 0;
    }

    for (tuple = 1; tuple < numTuples; tuple++) // find min. path
    {
      double* const  currCost = tempVtrbCost;
      uint8_t* const currPath = &prevPath[tuple * numStates];

      for (is = 0; is < numStates; is++)  // tuple's minimum path
      {
        uint8_t* currRate = &quantRate[(is + tuple * numStates) * numStates];
        double  costMinDs = (double) UINT_MAX;
        uint8_t pathMinDs = 0;

        for (ds = numStates - 1; ds >= 0; ds--)    // transitions
        {
          const double costCurr = (currRate[ds] >= UCHAR_MAX ? (double) UINT_MAX : prevCost[ds] + lambda * currRate[ds]);

          if (costMinDs > costCurr)
          {
            costMinDs = costCurr;
            pathMinDs = (uint8_t) ds;
          }
        }
        if (costMinDs < UINT_MAX) costMinDs += quantDist[tuple][is];

        currCost[is] = costMinDs;
        currPath[is] = pathMinDs;
      } // for is

      memcpy (prevCost, currCost, numStates * sizeof (double)); // TODO: avoid memcpy, use pointer swapping instead for speed
    } // for tuple

    for (tempBitCount = is = 0; is < numStates; is++)  // minimum
    {
      if (costMinIs > prevCost[is])
      {
        costMinIs = prevCost[is];
        pathMinIs = is;
      }
    }

    for (tuple--; tuple > 0; tuple--)  // min-cost rate and types
    {
      const uint8_t* currPath = &prevPath[tuple * numStates];
      const uint8_t pathMinDs = currPath[pathMinIs];

      optimalIs[tuple] = (uint8_t) pathMinIs;
      tempBitCount    += quantRate[pathMinDs + (pathMinIs + tuple * numStates) * numStates];
#if EC_TRAIN
      tempSfbDist += quantDist[tuple][pathMinIs];
#endif
      pathMinIs = pathMinDs;
    }
    optimalIs[0]  = (uint8_t) pathMinIs;
    tempBitCount += quantRate[pathMinIs * numStates];
#if EC_TRAIN
    tempSfbDist += quantDist[0][pathMinIs];
#endif
  } // Viterbi search

#if EC_TRAIN
  if ((tempSfbDist <= refSfbDist) || (tempBitCount <= targetBitCount))
#endif
  {
#if !EC_TRAIN
    numQ = 0;
#endif
    for (tuple = 0; tuple < numTuples; tuple++) // re-quantize SFB with R/D optimal rounding
    {
      const uint16_t tupleStart = tuple << 1;
      const uint8_t  tupIs = optimalIs[tuple];
      uint8_t& coeffQuantA = quantCoeffs[tupleStart];
      uint8_t& coeffQuantB = quantCoeffs[tupleStart + 1];

      if (tupIs != 0) // optimal red of quantized MDCT magnitudes
      {
        coeffQuantA -= (tupIs >> 1);
        coeffQuantB -= (tupIs &  1);
      }
#if !EC_TRAIN
      numQ += (coeffQuantA > 0 ? 1 : 0) + (coeffQuantB > 0 ? 1 : 0);
#endif
    } // for tuple

#if EC_TRAIN
    return tempBitCount;
#else
    return (1u << 15) | numQ; // final stats: OK flag | sign bits
#endif
  }

  return targetBitCount;
}
#endif // EC_TRELLIS_OPT_CODING

// constructor
SfbQuantizer::SfbQuantizer ()
{
  // initialize all helper buffers
  m_coeffMagn = nullptr;
#if EC_TRELLIS_OPT_CODING
  m_coeffTemp = nullptr;
#endif
  m_lut2ExpX4 = nullptr;
  m_lutSfNorm = nullptr;
  m_lutXExp43 = nullptr;

  m_maxSfIndex = 0;
#if EC_TRELLIS_OPT_CODING
  m_numCStates = 0;

  for (unsigned b = 0; b < 52; b++)
  {
    m_quantDist[b] = nullptr;
    m_quantInSf[b] = nullptr;
    m_quantRate[b] = nullptr;
  }
#endif
}

// destructor
SfbQuantizer::~SfbQuantizer ()
{
  // free allocated helper buffers
  MFREE (m_coeffMagn);
#if EC_TRELLIS_OPT_CODING
  MFREE (m_coeffTemp);
#endif
  MFREE (m_lut2ExpX4);
  MFREE (m_lutSfNorm);
  MFREE (m_lutXExp43);
#if EC_TRELLIS_OPT_CODING

  for (unsigned b = 0; b < 52; b++)
  {
    MFREE (m_quantDist[b]);
    MFREE (m_quantInSf[b]);
    MFREE (m_quantRate[b]);
  }
#endif
}

// public functions
unsigned SfbQuantizer::initQuantMemory (const unsigned maxTransfLength,
#if EC_TRELLIS_OPT_CODING
                                        const uint8_t numSwb, const uint8_t bitRateMode, const unsigned samplingRate,
#endif
                                        const uint8_t maxScaleFacIndex /*= SCHAR_MAX*/)
{
  const unsigned numScaleFactors = (unsigned) maxScaleFacIndex + 1;
#if EC_TRELLIS_OPT_CODING
  const uint8_t numTrellisStates = (samplingRate < 44100 ? 8 - samplingRate / 16000 : 5) - __min (2, (bitRateMode + 2) >> 2); // states per SFB
  const uint8_t numSquaredStates = numTrellisStates * numTrellisStates;
  const uint16_t quantRateLength = (samplingRate < 28800 || samplingRate >= 57600 ? 512 : 256); // quantizeMagnRDOC()
#endif
  unsigned x;

  if ((maxTransfLength < 128) || (maxTransfLength > 2048) || (maxTransfLength & 7) || (maxScaleFacIndex == 0) || (maxScaleFacIndex > SCHAR_MAX))
  {
    return 1; // invalid arguments error
  }

  m_maxSfIndex = maxScaleFacIndex;

  if ((m_coeffMagn = (unsigned*) malloc (maxTransfLength * sizeof (unsigned))) == nullptr ||
#if EC_TRELLIS_OPT_CODING
      (m_coeffTemp = (uint8_t* ) malloc (maxTransfLength + quantRateLength  )) == nullptr ||
#endif
      (m_lut2ExpX4 = (double*  ) malloc (numScaleFactors * sizeof (double  ))) == nullptr ||
      (m_lutSfNorm = (double*  ) malloc (numScaleFactors * sizeof (double  ))) == nullptr ||
      (m_lutXExp43 = (double*  ) malloc ((SCHAR_MAX + 1) * sizeof (double  ))) == nullptr)
  {
    return 2; // memory allocation error
  }
#if EC_TRELLIS_OPT_CODING
  m_maxSize8M1 = (maxTransfLength >> 3) - 1;
  m_numCStates = numTrellisStates;
  m_rateIndex  = bitRateMode;

  for (x = 0; x < __min (52u, numSwb); x++)
  {
    if ((m_quantDist[x] = (double*  ) malloc (numTrellisStates * sizeof (double  ))) == nullptr ||
        (m_quantInSf[x] = (uint8_t* ) malloc (numTrellisStates * sizeof (uint8_t ))) == nullptr ||
        (m_quantRate[x] = (uint16_t*) malloc (numSquaredStates * sizeof (uint16_t))) == nullptr)
    {
      return 2;
    }
  }
#else
  memset (m_coeffTemp, 0, sizeof (m_coeffTemp));
#endif
  // calculate scale factor gain 2^(x/4)
  for (x = 0; x < numScaleFactors; x++)
  {
    m_lut2ExpX4[x] = pow (2.0, (double) x / 4.0);
    m_lutSfNorm[x] = 1.0 / m_lut2ExpX4[x];
  }
  // calculate dequantized coeff x^(4/3)
  for (x = 0; x < (SCHAR_MAX + 1); x++)
  {
    m_lutXExp43[x] = pow ((double) x, 4.0 / 3.0);
  }

  return 0; // no error
}

uint8_t SfbQuantizer::quantizeSpecSfb (EntropyCoder& entropyCoder, const int32_t* const inputCoeffs, const uint8_t grpLength,
                                       const uint16_t* const grpOffsets, uint32_t* const grpStats,  // quant./coding statistics
                                       const unsigned sfb, const uint8_t sfIndex, const uint8_t sfIndexPred /*= UCHAR_MAX*/,
                                       uint8_t* const quantCoeffs /*= nullptr*/) // returns the RD optimized scale factor index
{
#if EC_TRELLIS_OPT_CODING
  EntropyCoder* const entrCoder = (grpLength == 1 ? &entropyCoder : nullptr);
#endif
  uint8_t sfBest = sfIndex;

  if ((inputCoeffs == nullptr) || (grpOffsets == nullptr) || (sfb >= 52) || (sfIndex > m_maxSfIndex))
  {
    return UCHAR_MAX; // invalid input error
  }

#if EC_TRELLIS_OPT_CODING
  if (grpLength == 1) // references for RDOC
  {
    m_quantDist[sfb][1] = -1.0;
    m_quantInSf[sfb][1] = sfIndex;
    m_quantRate[sfb][1] = 0; // for sgn bits
    m_quantRate[sfb][0] = entropyCoder.arithGetCtxState () & USHRT_MAX; // ref start context
  }
#endif
  if ((sfIndex == 0) || (sfIndexPred <= m_maxSfIndex && sfIndex + INDEX_OFFSET < sfIndexPred))
  {
    const uint16_t   grpStart = grpOffsets[0];
    const uint16_t   sfbStart = grpOffsets[sfb];
    const uint16_t   sfbWidth = grpOffsets[sfb + 1] - sfbStart;
    uint32_t* const coeffMagn = &m_coeffMagn[sfbStart];

    for (int i = sfbWidth - 1; i >= 0; i--) // back up magnitudes. TODO: use SIMD for speed?
    {
      coeffMagn[i] = abs (inputCoeffs[sfbStart + i]);
    }

    if (quantCoeffs)
    {
      memset (&quantCoeffs[sfbStart], 0, sfbWidth * sizeof (uint8_t)); // SFB output zeroing
      if (grpStats) // approximate bit count
      {
        grpStats[sfb] = getBitCount (entropyCoder, 0, 0, grpLength, &quantCoeffs[grpStart], sfbStart - grpStart, sfbWidth);
      }
    }
    return sfIndexPred - (sfIndex == 0 ? 0 : INDEX_OFFSET); // save delta bits if applicable
  }
  else // nonzero sf, optimized quantization
  {
    const uint16_t   grpStart = grpOffsets[0];
    const uint16_t   sfbStart = grpOffsets[sfb];
    const uint16_t   sfbWidth = grpOffsets[sfb + 1] - sfbStart;
    const uint16_t   cpyWidth = sfbWidth * sizeof (uint8_t);
    uint32_t* const coeffMagn = &m_coeffMagn[sfbStart];
    uint32_t codStart = 0, ctxStart = 0;
    uint32_t codFinal = 0, ctxFinal = 0;
    double   distBest = 0, distCurr = 0;
    short    maxQBest = 0, maxQCurr = 0;
    short    numQBest = 0, numQCurr = 0;
#if EC_TRELLIS_OPT_CODING
    bool rdOptimQuant = (grpLength != 1);
#else
    bool rdOptimQuant = true;
#endif
    uint8_t* ptrBest  = &m_coeffTemp[0];
    uint8_t* ptrCurr  = &m_coeffTemp[100];
    uint8_t  sfCurr   = sfIndex;

    for (int i = sfbWidth - 1; i >= 0; i--) // back up magnitudes. TODO: use SIMD for speed?
    {
      coeffMagn[i] = abs (inputCoeffs[sfbStart + i]);
    }

// --- determine default quantization result using range limited scale factor as a reference
    sfBest = quantizeMagnSfb (coeffMagn, sfCurr, ptrBest, sfbWidth,
#if EC_TRELLIS_OPT_CODING
                              entrCoder, sfbStart - grpStart,
#endif
                              &maxQBest, &numQBest);

    if (maxQBest > SCHAR_MAX) // limit SNR via scale factor index
    {
      for (uint8_t c = 0; (c < 2) && (maxQBest > SCHAR_MAX); c++)  // very rarely done twice
      {
        sfCurr += getScaleFacOffset (pow ((double) maxQBest, 4.0 / 3.0) * 0.001566492688) + c; // / m_lutXExp43[SCHAR_MAX]
        sfBest = quantizeMagnSfb (coeffMagn, sfCurr, ptrBest, sfbWidth,
#if EC_TRELLIS_OPT_CODING
                                  entrCoder, sfbStart - grpStart,
#endif
                                  &maxQBest, &numQBest);
      }
      rdOptimQuant = false;
    }
    else if ((sfBest < sfCurr) && (sfBest != sfIndexPred)) // re-optimize above quantization
    {
      sfBest = quantizeMagnSfb (coeffMagn, --sfCurr, ptrBest, sfbWidth,
#if EC_TRELLIS_OPT_CODING
                                entrCoder, sfbStart - grpStart,
#endif
                                &maxQBest, &numQBest);

      rdOptimQuant &= (maxQBest <= SCHAR_MAX);
    }

#if EC_TRELLIS_OPT_CODING
    if (grpLength == 1) // ref masking level
    {
      m_quantInSf[sfb][1] = __min (sfCurr, m_maxSfIndex);
    }
#endif
    if (maxQBest == 0) // SFB was quantized to zero - zero output
    {
      if (quantCoeffs)
      {
        memset (&quantCoeffs[sfbStart], 0, cpyWidth);
        if (grpStats) // estimated bit count
        {
          grpStats[sfb] = getBitCount (entropyCoder, 0, 0, grpLength, &quantCoeffs[grpStart], sfbStart - grpStart, sfbWidth);
        }
      }
      return sfIndexPred; // repeat scale factor, save delta bits
    }

// --- check whether optimized quantization and coding results in lower rate-distortion cost
    distBest = getQuantDist (coeffMagn, sfBest, ptrBest, sfbWidth);

#if EC_TRELLIS_OPT_CODING
    if (grpLength == 1) // ref band-wise NMR
    {
      const double refSfbNmrDiv = m_lutSfNorm[m_quantInSf[sfb][1]];

      m_quantDist[sfb][1] = distBest * refSfbNmrDiv * refSfbNmrDiv;
      m_quantRate[sfb][1] = numQBest; // sgn
    }
#endif
    if (quantCoeffs)
    {
      memcpy (&quantCoeffs[sfbStart], ptrBest, cpyWidth);

      codStart = entropyCoder.arithGetCodState (); // start state
      ctxStart = entropyCoder.arithGetCtxState ();
      numQBest += getBitCount (entropyCoder, sfBest, sfIndexPred, grpLength, &quantCoeffs[grpStart], sfbStart - grpStart, sfbWidth);
      codFinal = entropyCoder.arithGetCodState (); // final state
      ctxFinal = entropyCoder.arithGetCtxState ();
    }
    rdOptimQuant &= (distBest > 0.0);

    if ((sfBest < sfCurr) && (sfBest != sfIndexPred) && rdOptimQuant) // R/D re-optimization
    {
#if EC_TRELLIS_OPT_CODING
      const double refSfbNmrDiv = m_lutSfNorm[sfCurr];
      const double lambda       = getLagrangeValue (m_rateIndex);
#endif
      sfCurr = quantizeMagnSfb (coeffMagn, sfCurr - 1, ptrCurr, sfbWidth,
#if EC_TRELLIS_OPT_CODING
                                entrCoder, sfbStart - grpStart,
#endif
                                &maxQCurr, &numQCurr);

      distCurr = getQuantDist (coeffMagn, sfCurr, ptrCurr, sfbWidth);
      if (quantCoeffs)
      {
        memcpy (&quantCoeffs[sfbStart], ptrCurr, cpyWidth);

        entropyCoder.arithSetCodState (codStart);  // reset state
        entropyCoder.arithSetCtxState (ctxStart);
        numQCurr += getBitCount (entropyCoder, sfCurr, sfIndexPred, grpLength, &quantCoeffs[grpStart], sfbStart - grpStart, sfbWidth);
      }

      // rate-distortion decision, using empirical Lagrange value
#if EC_TRELLIS_OPT_CODING
      if (distCurr * refSfbNmrDiv * refSfbNmrDiv + lambda * numQCurr < distBest * refSfbNmrDiv * refSfbNmrDiv + lambda * numQBest)
#else
      if ((maxQCurr <= maxQBest) && (numQCurr <= numQBest + (distCurr >= distBest ? -1 : short (0.5 + distBest / __max (1.0, distCurr)))))
#endif
      {
        maxQBest = maxQCurr;
        numQBest = numQCurr;
        sfBest   = sfCurr;
      }
      else if (quantCoeffs) // discard result, recover best trial
      {
        memcpy (&quantCoeffs[sfbStart], ptrBest, cpyWidth);

        entropyCoder.arithSetCodState (codFinal);  // reset state
        entropyCoder.arithSetCtxState (ctxFinal);
      }
    }

    if (grpStats)
    {
      grpStats[sfb] = ((uint32_t) maxQBest << 16) | numQBest; // max magnitude and bit count
    }
  } // if sfIndex == 0

  return __min (sfBest, m_maxSfIndex);
}

#if EC_TRELLIS_OPT_CODING
unsigned SfbQuantizer::quantizeSpecRDOC (EntropyCoder& entropyCoder, uint8_t* const optimalSf, const unsigned targetBitCount,
                                         const uint16_t* const grpOffsets, uint32_t* const grpStats,  // quant./coding statistics
                                         const unsigned numSfb, uint8_t* const quantCoeffs)  // returns RD optimization bit count
{
  // numSfb: number of trellis stages. Based on: A. Aggarwal, S. L. Regunathan, and K. Rose,
  // "Trellis-Based Optimization of MPEG-4 Advanced Audio Coding," in Proc. IEEE Workshop on
  // Speech Coding, pp. 142-144, Sep. 2000. Modified for arithmetic instead of Huffman coder
  const uint32_t codStart = USHRT_MAX << 16;
  const uint32_t ctxStart = m_quantRate[0][0]; // start context before call to quantizeSfb()
  const uint32_t codFinal = entropyCoder.arithGetCodState ();
  const uint32_t ctxFinal = entropyCoder.arithGetCtxState (); // after call to quantizeSfb()
  const uint16_t grpStart = grpOffsets[0];
  uint8_t* const inScaleFac = &m_coeffTemp[((unsigned) m_maxSize8M1 - 6) << 3];
  uint32_t  prevCodState[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  uint32_t  prevCtxState[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  uint8_t   prevScaleFac[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  double    prevVtrbCost[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  uint32_t  tempCodState[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  uint32_t  tempCtxState[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  uint8_t   tempScaleFac[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  double    tempVtrbCost[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  unsigned  tempBitCount, sfb, is;
  int ds;
#if EC_TRAIN
  double refGrpDist = 0.0, tempGrpDist = 0.0;
#else
  const double lambda = getLagrangeValue (m_rateIndex);
#endif

  if ((optimalSf == nullptr) || (quantCoeffs == nullptr) || (grpOffsets == nullptr) || (numSfb == 0) || (numSfb > 52) ||
      (targetBitCount == 0)  || (targetBitCount > SHRT_MAX))
  {
    return 0; // invalid input error
  }

  for (sfb = 0; sfb < numSfb; sfb++) // SFB-wise scale factor, weighted distortion, and rate
  {
    const uint8_t       refSf = m_quantInSf[sfb][1];
    const uint16_t    refNumQ = m_quantRate[sfb][1];
    const double refQuantDist = m_quantDist[sfb][1];
    const double refQuantNorm = m_lutSfNorm[refSf] * m_lutSfNorm[refSf];
    const uint16_t   sfbStart = grpOffsets[sfb];
    const uint16_t   sfbWidth = grpOffsets[sfb + 1] - sfbStart;
    const uint32_t* coeffMagn = &m_coeffMagn[sfbStart];
    uint8_t* const  tempQuant = &m_coeffTemp[sfbStart - grpStart];
    bool maxSnrReached = false;

    if (refQuantDist < 0.0) memset (tempQuant, 0, sfbWidth * sizeof (uint8_t));
#if EC_TRAIN
    else refGrpDist += refQuantDist;
#endif
    if (grpStats)
    {
      grpStats[sfb] = (grpStats[sfb] & (USHRT_MAX << 16)) | refNumQ; // keep magn, sign bits
    }

    for (is = 0; is < m_numCStates; is++) // populate SFB trellis
    {
      const uint8_t* mag = (is != 1 ? m_coeffTemp /*= tempQuant[grpStart - sfbStart]*/ : &quantCoeffs[grpStart]);
      double&   currDist = m_quantDist[sfb][is];
      uint16_t* currRate = &m_quantRate[sfb][is * m_numCStates];
      uint8_t     sfBest = optimalSf[sfb]; // optimal scalefactor
      short maxQCurr = 0, numQCurr = 0; // for sign bits counting

      if (refQuantDist < 0.0) // -1.0 means SFB is zero-quantized
      {
        currDist = -1.0;
        m_quantInSf[sfb][is] = refSf;
      }
      else if (is != 1) // quantization & distortion not computed
      {
        const uint8_t sfCurr = __max (0, __min (m_maxSfIndex, refSf + 1 - (int) is));

        currDist = -1.0;
        if ((sfCurr == 0) || maxSnrReached)
        {
          maxSnrReached = true;
        }
        else // sfCurr > 0 && sfCurr <= m_maxSfIndex, re-quantize
        {
          sfBest = quantizeMagnSfb (coeffMagn, sfCurr, tempQuant, sfbWidth,
                                    &entropyCoder, sfbStart - grpStart,
                                    &maxQCurr, &numQCurr);

          if (maxQCurr > SCHAR_MAX)
          {
            maxSnrReached = true; numQCurr = 0;
          }
          else
          {
            currDist = getQuantDist (coeffMagn, sfBest, tempQuant, sfbWidth) * refQuantNorm;
          }
        }
        if (currDist < 0.0) memset (tempQuant, 0, sfbWidth * sizeof (uint8_t));
        m_quantInSf[sfb][is] = sfCurr; // store initial scale fac
      }
      else // is == 1, quant. & dist. computed with quantizeSfb()
      {
        numQCurr = refNumQ;
      }

      if (sfb == 0) // first SFB, having sfbStart - grpStart == 0
      {
        entropyCoder.arithSetCodState (codStart);  // group start
        entropyCoder.arithSetCtxState (ctxStart, 0);
        tempBitCount = (maxSnrReached ? USHRT_MAX : numQCurr + getBitCount (entropyCoder, sfBest, UCHAR_MAX, 1, mag, 0, sfbWidth));

        for (ds = m_numCStates - 1; ds >= 0; ds--)
        {
          currRate[ds] = (uint16_t) tempBitCount;
        }
        tempCodState[is] = entropyCoder.arithGetCodState ();
        tempCtxState[is] = entropyCoder.arithGetCtxState ();
      }
      else // sfb > 0, rate depends on decisions in preceding SFB
      {
        for (ds = m_numCStates - 1; ds >= 0; ds--)
        {
          const uint16_t prevRate = m_quantRate[sfb - 1][ds * m_numCStates];

          entropyCoder.arithSetCodState (prevCodState[ds]);
          entropyCoder.arithSetCtxState (prevCtxState[ds], sfbStart - grpStart);
          tempBitCount = (maxSnrReached || (prevRate >= USHRT_MAX) ? USHRT_MAX : numQCurr + getBitCount (entropyCoder,
                           (numQCurr == 0 ? prevScaleFac[ds] : sfBest), prevScaleFac[ds], 1, mag, sfbStart - grpStart, sfbWidth));
          currRate[ds] = (uint16_t) tempBitCount;

          if (ds == 1) // statistically best place to save states
          {
            tempCodState[is] = entropyCoder.arithGetCodState ();
            tempCtxState[is] = entropyCoder.arithGetCtxState ();
          }
        }
      }
      tempScaleFac[is] = sfBest; // optimized factor for next SFB
    } // for is

    memcpy (prevCodState, tempCodState, m_numCStates * sizeof (uint32_t));
    memcpy (prevCtxState, tempCtxState, m_numCStates * sizeof (uint32_t));
    memcpy (prevScaleFac, tempScaleFac, m_numCStates * sizeof (uint8_t ));
  } // for sfb

  entropyCoder.arithSetCodState (codFinal); // back to last state
  entropyCoder.arithSetCtxState (ctxFinal, grpOffsets[numSfb] - grpStart);

#if EC_TRAIN
  tempBitCount = targetBitCount + 1; // Viterbi search for minimum distortion at target rate
  for (double lambda = 0.015625; (lambda <= 0.375) && (tempBitCount > targetBitCount); lambda += 0.0078125)
#endif
  {
    double* const  prevCost = prevVtrbCost;
    uint8_t* const prevPath = m_coeffTemp; // trellis backtracker
    double   costMinIs = (double) UINT_MAX;
    unsigned pathMinIs = 1;
#if EC_TRAIN
    tempGrpDist = 0.0;
#endif

    for (is = 0; is < m_numCStates; is++) // initial minimum path
    {
      const uint16_t currRate = m_quantRate[0][is * m_numCStates];

      prevCost[is] = (currRate >= USHRT_MAX ? (double) UINT_MAX : lambda * currRate + __max (0.0, m_quantDist[0][is]));
      prevPath[is] = 0;
    }

    for (sfb = 1; sfb < numSfb; sfb++) // search for minimum path
    {
      double* const  currCost = tempVtrbCost;
      uint8_t* const currPath = &prevPath[sfb * m_numCStates];

      for (is = 0; is < m_numCStates; is++) // SFB's minimum path
      {
        uint16_t* currRate = &m_quantRate[sfb][is * m_numCStates];
        double   costMinDs = (double) UINT_MAX;
        uint8_t  pathMinDs = 1;

        for (ds = m_numCStates - 1; ds >= 0; ds--) // transitions
        {
          const double costCurr = (currRate[ds] >= USHRT_MAX ? (double) UINT_MAX : prevCost[ds] + lambda * currRate[ds]);

          if (costMinDs > costCurr)
          {
            costMinDs = costCurr;
            pathMinDs = (uint8_t) ds;
          }
        }
        if (costMinDs < UINT_MAX) costMinDs += __max (0.0, m_quantDist[sfb][is]);

        currCost[is] = costMinDs;
        currPath[is] = pathMinDs;
      } // for is

      memcpy (prevCost, currCost, m_numCStates * sizeof (double)); // TODO: avoid memcpy, use pointer swapping instead for speed
    } // for sfb

    for (sfb--, is = 0; is < m_numCStates; is++) // group minimum
    {
      if (costMinIs > prevCost[is])
      {
        costMinIs = prevCost[is];
        pathMinIs = is;
      }
    }

    for (tempBitCount = 0; sfb > 0; sfb--) // min-cost group rate
    {
      const uint8_t* currPath = &prevPath[sfb * m_numCStates];
      const uint8_t pathMinDs = currPath[pathMinIs];

      inScaleFac[sfb] = (m_quantDist[sfb][pathMinIs] < 0.0 ? UCHAR_MAX : m_quantInSf[sfb][pathMinIs]);
      tempBitCount   +=  m_quantRate[sfb][pathMinDs + pathMinIs * m_numCStates];
#if EC_TRAIN
      tempGrpDist += __max (0.0, m_quantDist[sfb][pathMinIs]);
#endif
      pathMinIs = pathMinDs;
    }
    inScaleFac[0] = (m_quantDist[0][pathMinIs] < 0.0 ? UCHAR_MAX : m_quantInSf[0][pathMinIs]);
    tempBitCount +=  m_quantRate[0][pathMinIs * m_numCStates];
#if EC_TRAIN
    tempGrpDist += __max (0.0, m_quantDist[0][pathMinIs]);
#endif
  } // Viterbi search

#if EC_TRAIN
  if ((tempGrpDist <= refGrpDist) || (tempBitCount <= targetBitCount))
#endif
  {
    uint8_t sfIndexPred = UCHAR_MAX;

    if (grpStats)
    {
      entropyCoder.arithSetCodState (codStart);// set group start
      entropyCoder.arithSetCtxState (ctxStart, 0);

      tempBitCount = 0;
    }
    for (sfb = 0; sfb < numSfb; sfb++) // re-quantize spectrum with R/D optimized parameters
    {
      const uint16_t sfbStart = grpOffsets[sfb];
      const uint16_t sfbWidth = grpOffsets[sfb + 1] - sfbStart;

      if ((inScaleFac[sfb] == UCHAR_MAX) || (sfIndexPred <= m_maxSfIndex && inScaleFac[sfb] + INDEX_OFFSET < sfIndexPred))
      {
        memset (&quantCoeffs[sfbStart], 0, sfbWidth * sizeof (uint8_t));  // zero SFB output

        optimalSf[sfb] = sfIndexPred - (inScaleFac[sfb] == UCHAR_MAX ? 0 : INDEX_OFFSET);
      }
      else if (inScaleFac[sfb] != m_quantInSf[sfb][1]) // speedup
      {
        short maxQBest = 0, numQBest = 0;

        optimalSf[sfb] = quantizeMagnSfb (&m_coeffMagn[sfbStart], inScaleFac[sfb], &quantCoeffs[sfbStart], sfbWidth,
                                          &entropyCoder, sfbStart - grpStart,
                                          &maxQBest, &numQBest);

        if (maxQBest == 0) optimalSf[sfb] = sfIndexPred; // empty
        if (grpStats)
        {
          grpStats[sfb] = ((uint32_t) maxQBest << 16) | numQBest; // max magn. and sign bits
        }
      }

      if (grpStats) // complete statistics with per-SFB bit count
      {
        grpStats[sfb] += getBitCount (entropyCoder, optimalSf[sfb], sfIndexPred, 1, &quantCoeffs[grpStart], sfbStart - grpStart, sfbWidth);
        tempBitCount  += grpStats[sfb] & USHRT_MAX;
      }

      if ((sfb > 0) && (optimalSf[sfb] < UCHAR_MAX) && (sfIndexPred == UCHAR_MAX))
      {
        memset (optimalSf, optimalSf[sfb], sfb * sizeof (uint8_t)); // back-propagate factor
      }
      sfIndexPred = optimalSf[sfb];
    } // for sfb

    return tempBitCount + (grpStats ? 2 : 0); // last coding bits
  }

  return targetBitCount;
}
#endif // EC_TRELLIS_OPT_CODING
