/* exhaleAppPch.cpp - pre-compiled source file for source code of exhale application
 * written by C. R. Helmrich, last modified in 2019 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under a Modified BSD-Style License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2020 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#include "exhaleAppPch.h"

// ISO/IEC 23003-3 USAC Table 67
static const unsigned supportedSamplingRates[26] = {
  96000, 88200, 64000, 48000, 44100, 32000, 24000, 22050, 16000, 12000, 11025,  8000, 7350, // AAC
  57600, 51200, 40000, 38400, 34150, 28800, 25600, 20000, 19200, 17075, 14400, 12800, 9600 // USAC
};

// public sampling rate function
bool isSamplingRateSupported (const unsigned samplingRate)
{
  for (unsigned i = 0; i < 26; i++)
  {
    if (samplingRate == supportedSamplingRates[i]) return true;
  }
  return false; // not supported
}
