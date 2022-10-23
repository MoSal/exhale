/* entropyCoding.h - header file for class with lossless entropy coding capability
 * written by C. R. Helmrich, last modified in 2022 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under the exhale Copyright License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2021 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#ifndef _ENTROPY_CODING_H_
#define _ENTROPY_CODING_H_

#include "exhaleLibPch.h"

// constants, experimental macro
#define ARITH_ESCAPE          16
#define ARITH_SIZE           742
#define INDEX_OFFSET          60
#define INDEX_SIZE           121
#define EC_TRELLIS_OPT_CODING  1

// lossless entropy coding class
class EntropyCoder
{
private:

  // member variables
  uint8_t* m_qcCurr;         // curr. window's quantized context q[1]
  uint8_t* m_qcPrev;         // prev. window's quantized context q[0]

  uint16_t m_acBits;         // bits_to_follow in arith_encode, 0..31
  uint16_t m_acHigh;         // high in arith_encode as in Annex B.25
  uint16_t m_acLow;          // low in arith_encode, as in Annex B.25
  uint16_t m_acSize;         // context window size (N/4 in Scl. 7.4)
  uint32_t m_csCurr;         // context state, see initWindowCoding()
  unsigned m_maxTupleLength; // maximum half-transform length (<4096)
  bool     m_shortTrafoCurr; // used to derive N in Scl. 7.4 and B.25
  bool     m_shortTrafoPrev; // used to derive previous_N in Scl. 7.4

  // helper functions
  unsigned arithCodeSymbol (const uint16_t symbol, const uint16_t* table, OutputStream* const stream = nullptr);
  unsigned arithGetContext (const unsigned ctx, const unsigned idx);
  unsigned arithMapContext (const bool arithResetFlag);
#if EC_TRELLIS_OPT_CODING
  void     arithSetContext (const unsigned newCtxState, const uint16_t sigEnd);
#endif

public:

  // constructor
  EntropyCoder ();
  // destructor
  ~EntropyCoder ();
  // public functions
  unsigned arithCodeSigMagn (const uint8_t* const magn, const uint16_t sigOffset, const uint16_t sigLength,
                             const bool arithFinish = false, OutputStream* const stream = nullptr);
#if EC_TRELLIS_OPT_CODING
  unsigned arithCodeSigTest (const uint8_t* const magn, const uint16_t sigOffset, const uint16_t sigLength); // +-m_acBits
  unsigned arithCodeTupTest (const uint8_t* const magn, const uint16_t sigOffset); // for sigLength of 2 - also +-m_acBits
#endif
  unsigned arithGetCodState () const                     { return ((unsigned) m_acHigh << 16) | (unsigned) m_acLow; }
  unsigned arithGetCtxState () const                     { return m_csCurr; }
  unsigned arithGetResetBit (const uint8_t* const magn, const uint16_t sigOffset, const uint16_t sigLength);
  char*    arithGetTuplePtr () const                     { return (char*) m_qcCurr; }
  void     arithResetMemory () { memset (m_qcPrev, 0, (m_maxTupleLength + 1) * sizeof (uint8_t)); m_acBits = 0; }
  void     arithSetCodState (const unsigned newCodState) { m_acHigh = newCodState >> 16; m_acLow = newCodState & USHRT_MAX; }
#if EC_TRELLIS_OPT_CODING
  void     arithSetCtxState (const unsigned newCtxState, const uint16_t sigOffset = 0) { arithSetContext (newCtxState, sigOffset >> 1); }
#else
  void     arithSetCtxState (const unsigned newCtxState) { m_csCurr = newCtxState; }
#endif
  unsigned indexGetBitCount (const int scaleFactorDelta) const;
  unsigned indexGetHuffCode (const int scaleFactorDelta) const;

  unsigned initCodingMemory (const unsigned maxTransfLength);
  unsigned initWindowCoding (const bool forceArithReset, const bool shortWin = false);

  bool     getIsShortWindow () const                     { return m_shortTrafoCurr; }
  void     setIsShortWindow (const bool shortWin)        { m_shortTrafoCurr = shortWin; }
}; // EntropyCoder

#endif // _ENTROPY_CODING_H_
