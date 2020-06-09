/* specGapFilling.h - header file for class with spectral gap filling coding methods
 * written by C. R. Helmrich, last modified in 2020 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under a Modified BSD-Style License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2020 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#ifndef _SPEC_GAP_FILLING_H_
#define _SPEC_GAP_FILLING_H_

#include "exhaleLibPch.h"
#include "quantization.h"

// constants, experimental macro
#define SGF_LIMIT 2*INDEX_OFFSET
#define SGF_OPT_SHORT_WIN_CALC 1
#define SGF_SF_PEAK_SMOOTHING  1

// MDCT-domain gap filling class
class SpecGapFiller
{
private:

  // member variables
  uint16_t  m_1stGapFillSfb;
  int16_t   m_1stNonZeroSfb[NUM_WINDOW_GROUPS];

public:

  // constructor
  SpecGapFiller ();
  // destructor
  ~SpecGapFiller () { }
  // public functions
  uint16_t  getFirstGapFillSfb () const { return m_1stGapFillSfb; }
  uint8_t   getSpecGapFillParams (const SfbQuantizer& sfbQuantizer, const uint8_t* const quantMagn,
                                  const uint8_t numSwbShort, SfbGroupData& grpData /*modified*/,
                                  const unsigned nSamplesInFrame, const bool saveRate, const uint8_t specFlat);
}; // SpecGapFiller

#endif // _SPEC_GAP_FILLING_H_
