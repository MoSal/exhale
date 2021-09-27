/* lappedTransform.h - header file for class providing time-frequency transformation
 * written by C. R. Helmrich, last modified in 2019 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under the exhale Copyright License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2021 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#ifndef _LAPPED_TRANSFORM_H_
#define _LAPPED_TRANSFORM_H_

#include "exhaleLibPch.h"

// constants, experimental macros
#define LUT_OFFSET      (1 << 30)
#define LUT_SHIFT              31
#define WIN_OFFSET      (1 << 24)
#define WIN_SHIFT              25

// time-frequency transform class
class LappedTransform
{
private:

  // member variables
  int32_t* m_dctRotCosL;
  int32_t* m_dctRotCosS;
  int32_t* m_dctRotSinL;
  int32_t* m_dctRotSinS;
  int32_t* m_fftHalfCos;
  int32_t* m_fftHalfSin;
  short*   m_fftPermutL;
  short*   m_fftPermutS;
  int32_t* m_tempIntBuf;     // pointer to temporary helper buffer
  int32_t* m_timeWindowL[2]; // pointer to two long window halves
  int32_t* m_timeWindowS[2]; // pointer to two short window halves
  short    m_transfLengthL;
  short    m_transfLengthS;

  // helper functions
  void applyHalfSizeFFT (int32_t* const iR/*eal*/, int32_t* const iI/*mag*/, const bool shortTransform);
  void windowAndFoldInL (const int32_t* inputL, const bool shortTransform, const bool kbdWindowL, const bool lowOverlapL,
                         const bool mdstKernel, int32_t* const output);
  void windowAndFoldInR (const int32_t* inputR, const bool shortTransform, const bool kbdWindowR, const bool lowOverlapR,
                         const bool mdstKernel, int32_t* const output);

public:

  // constructor
  LappedTransform ();
  // destructor
  ~LappedTransform ();
  // public functions
  unsigned applyNegDCT4  (int32_t* const signal,  const bool shortTransform);
  unsigned applyMCLT     (const int32_t* timeSig, const bool eightTransforms, bool kbdWindowL, const bool kbdWindowR,
                          const bool lowOverlapL, const bool lowOverlapR, int32_t* const outMdct, int32_t* const outMdst);
  unsigned initConstants (int32_t* const tempIntBuf, int32_t* const timeWindowL[2], int32_t* const timeWindowS[2],
                          const unsigned maxTransfLength);
}; // LappedTransform

#endif // _LAPPED_TRANSFORM_H_
