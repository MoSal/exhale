/* linearPrediction.h - header file for class providing linear prediction capability
 * written by C. R. Helmrich, last modified in 2020 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under the exhale Copyright License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2021 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#ifndef _LINEAR_PREDICTION_H_
#define _LINEAR_PREDICTION_H_

#include "exhaleLibPch.h"

// constants, experimental macros
#define LP_EPS                  1
#define LP_SHIFT               15
#define LP_DEPTH               (1 + LP_SHIFT - MAX_PREDICTION_ORDER)
#define LP_OFFSET              (1 << (LP_SHIFT - 1))

// linear predictive filter class
class LinearPredictor
{
private:

  // temporary buffer
  int64_t m_tempBuf[2 * MAX_PREDICTION_ORDER];

public:

  // constructor
  LinearPredictor ();
  // destructor
  ~LinearPredictor () { }
  // public functions
  uint32_t calcParCorCoeffs (const int32_t* const anaSignal, const uint16_t nAnaSamples, const uint16_t nCoeffs,
                             short* const parCorCoeffs); // returns 256 - 256 / prediction gain per filter order, or 0
  uint8_t  calcOptTnsCoeffs (short* const parCorCoeffs, int8_t* const quantCoeffs, bool* const lowCoeffRes,
                             const uint16_t maxOrder, const uint8_t predGain, const uint8_t tonality = 0,
                             const uint16_t parCorCoeffBitDepth = 10); // returns optimized filter order for TNS
  unsigned lpToParCorCoeffs (short* const lpCoeffs, const uint16_t nCoeffs, short* const parCorCoeffs,
                             const uint16_t parCorCoeffBitDepth = 10);
  unsigned parCorToLpCoeffs (const short* const parCorCoeffs,  const uint16_t nCoeffs, short* const lpCoeffs,
                             const uint16_t parCorCoeffBitDepth = 10);
  unsigned quantTnsToLpCoeffs(const int8_t* const quantCoeffs, const uint16_t nCoeffs, const bool lowCoeffRes,
                             short* const parCorCoeffs, short* const lpCoeffs);
  bool  similarParCorCoeffs (const short* const parCorCoeffs1, const short* const parCorCoeffs2, const uint16_t nCoeffs,
                             const uint16_t parCorCoeffBitDepth = 10);
}; // LinearPredictor

#endif // _LINEAR_PREDICTION_H_
