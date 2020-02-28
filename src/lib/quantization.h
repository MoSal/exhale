/* quantization.h - header file for class with nonuniform quantization functionality
 * written by C. R. Helmrich, last modified in 2020 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under a Modified BSD-Style License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2020 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#ifndef _QUANTIZATION_H_
#define _QUANTIZATION_H_

#include "exhaleLibPch.h"
#include "entropyCoding.h"

// constants, experimental macros
#define FOUR_LOG102   13.28771238 // 4 / log10 (2)
#define SF_QUANT_OFFSET 0.4783662 // for scale fac
#define SFB_QUANT_FAST_POW      1 // faster pow ()
#define SFB_QUANT_PERCEPT_OPT   1 // psych. quant.
#if SFB_QUANT_FAST_POW
#define SFB_QUANT_OFFSET 0.496094 // 13 - 29^(3/4)
#else
#define SFB_QUANT_OFFSET 0.405396 // 1 - 0.5^(3/4)
#endif

// class for xHE-AAC quantization
class SfbQuantizer
{
private:

  // member variables
  unsigned* m_coeffMagn; // temp memory
#if EC_TRELLIS_OPT_CODING
  uint8_t   m_coeffTemp[1024]; // TODO?
#else
  uint8_t   m_coeffTemp[200]; // 40 * 5 - NOTE: increase this when maximum grpLength > 5
#endif
  double*   m_lut2ExpX4; // for 2^(X/4)
  double*   m_lutSfNorm; // 1 / 2^(X/4)
  double*   m_lutXExp43; // for X^(4/3)
  uint8_t   m_maxSfIndex; // 1,..., 127
#if EC_TRELLIS_OPT_CODING
  uint8_t   m_numCStates; // states/SFB
  uint16_t  m_rateIndex; // lambda mode
  // trellis memory, max. 8 KB @ num_swb=51
  double*   m_quantDist[52]; // quantizing distortion
  uint8_t*  m_quantInSf[52]; // initial scale factors
  uint16_t* m_quantRate[52]; // MDCT and SF bit count
#endif

  // helper functions
  double    getQuantDist (const unsigned* const coeffMagn, const uint8_t scaleFactor,
                          const uint8_t* const coeffQuant, const uint16_t numCoeffs);
  uint8_t   quantizeMagnSfb (const unsigned* const coeffMagn, const uint8_t scaleFactor,
                            /*mod*/uint8_t* const coeffQuant, const uint16_t numCoeffs,
#if EC_TRELLIS_OPT_CODING
                             EntropyCoder* const arithmCoder, const uint16_t coeffOffset,
#endif
                             short* const sigMaxQ = nullptr,  short* const sigNumQ = nullptr);
#if EC_TRELLIS_OPT_CODING
  uint32_t quantizeMagnRDOC (EntropyCoder& entropyCoder, const uint8_t optimalSf, const unsigned targetBitCount,
                             const uint16_t coeffOffset, const unsigned* const coeffMagn,  // initial MDCT magnitudes
                             const uint16_t numCoeffs, uint8_t* const quantCoeffs); // returns updated SFB statistics
#endif
public:

  // constructor
  SfbQuantizer ();
  // destructor
  ~SfbQuantizer ();
  // public functions
  unsigned* getCoeffMagnPtr ()                      const { return m_coeffMagn; }
  double*   getSfNormTabPtr ()                      const { return m_lutSfNorm; }
  uint8_t getScaleFacOffset (const double absValue) const { return uint8_t (SF_QUANT_OFFSET + FOUR_LOG102 * log10 (__max (1.0, absValue))); }
  unsigned  initQuantMemory (const unsigned maxTransfLength,
#if EC_TRELLIS_OPT_CODING
                             const uint8_t numSwb, const uint8_t bitRateMode,
#endif
                             const uint8_t maxScaleFacIndex = SCHAR_MAX);
  uint8_t   quantizeSpecSfb (EntropyCoder& entropyCoder, const int32_t* const inputCoeffs, const uint8_t grpLength,
                             const uint16_t* const grpOffsets, uint32_t* const grpStats,  // quant./coding statistics
                             const unsigned sfb, const uint8_t sfIndex, const uint8_t sfIndexPred = UCHAR_MAX,
                             uint8_t* const quantCoeffs = nullptr); // returns index of the RD optimized scale factor
#if EC_TRELLIS_OPT_CODING
  unsigned quantizeSpecRDOC (EntropyCoder& entropyCoder, uint8_t* const optimalSf, const unsigned targetBitCount,
                             const uint16_t* const grpOffsets, uint32_t* const grpStats,  // quant./coding statistics
                             const unsigned numSfb, uint8_t* const quantCoeffs); // returns RD optimization bit count
#endif
}; // SfbQuantizer

#endif // _QUANTIZATION_H_
