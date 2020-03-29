/* bitAllocation.h - header file for class needed for psychoacoustic bit-allocation
 * written by C. R. Helmrich, last modified in 2020 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under a Modified BSD-Style License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2020 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#ifndef _BIT_ALLOCATION_H_
#define _BIT_ALLOCATION_H_

#include "exhaleLibPch.h"

// constants, experimental macros
#define BA_EPS                  1
#define BA_INTER_CHAN_SIM_MASK  0  // 5.1 cross-channel simultaneous masking

// class for audio bit-allocation
class BitAllocator
{
private:

  // member variables
  uint32_t m_avgStepSize[USAC_MAX_NUM_CHANNELS];
  uint8_t  m_avgSpecFlat[USAC_MAX_NUM_CHANNELS];
  uint8_t  m_avgTempFlat[USAC_MAX_NUM_CHANNELS];

public:

  // constructor
  BitAllocator ();
  // destructor
  ~BitAllocator () { }
  // public functions
  void getChAverageSpecFlat (uint8_t meanSpecFlatInCh[USAC_MAX_NUM_CHANNELS], const unsigned nChannels);
  void getChAverageTempFlat (uint8_t meanTempFlatInCh[USAC_MAX_NUM_CHANNELS], const unsigned nChannels);
  uint8_t       getScaleFac (const uint32_t sfbStepSize, const int32_t* const sfbSignal, const uint8_t sfbWidth,
                             const uint32_t sfbRmsValue);
  unsigned initSfbStepSizes (const SfbGroupData* const groupData[USAC_MAX_NUM_CHANNELS], const uint8_t numSwbShort,
                             const uint32_t specAnaStats[USAC_MAX_NUM_CHANNELS],
                             const uint32_t tempAnaStats[USAC_MAX_NUM_CHANNELS],
                             const unsigned nChannels, const unsigned samplingRate, uint32_t* const sfbStepSizes,
                             const unsigned lfeChannelIndex = USAC_MAX_NUM_CHANNELS, const bool tnsDisabled = false);
}; // BitAllocator

#endif // _BIT_ALLOCATION_H_
