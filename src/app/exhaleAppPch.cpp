/* exhaleAppPch.cpp - pre-compiled source file for source code of exhale application
 * written by C. R. Helmrich, last modified in 2021 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under the exhale Copyright License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2021 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#include "exhaleAppPch.h"

// ISO/IEC 23003-3 USAC Table 67
static const unsigned supportedSamplingRates[16] = {
  96000, 88200, 64000, 48000, 44100, 32000, 24000, 22050, 16000, 12000, 11025, 8000, 7350, // AAC
  57600, 38400, 19200 // BL USAC
};

// public extrapolation function
void eaExtrapolate (int32_t* const pcmBuffer, const uint16_t pcmOffset, // start/end of PCM fades
                    const uint16_t frameSize, const uint16_t numChannels, const bool fadeIn /*= false*/)
{
  const int32_t delta = (fadeIn ? -1 : 1) * numChannels;
  const uint16_t size = (fadeIn ? pcmOffset : frameSize - pcmOffset);

  if ((pcmOffset == 0 && fadeIn) || (pcmOffset >= frameSize) || !pcmBuffer) return;

  for (uint16_t ch = 0; ch < numChannels; ch++)
  {
    int32_t* chPcmBuf = pcmBuffer + ch + (pcmOffset - (fadeIn ? 0 : 1)) * numChannels;
    int32_t  result32 = (pcmOffset == 0 ? 0 : *chPcmBuf * (1 << 8)); // input is known to be 24-bit PCM
    const int32_t s32 = result32 / size;

    for (uint16_t i = size; i > 0; i--) *(chPcmBuf += delta) = ((result32 -= s32) + 128) >> 8;
  }
}

// public sampling rate function
bool isSamplingRateSupported (const unsigned samplingRate)
{
  for (uint8_t i = 0; i < 16; i++)
  {
    if (samplingRate == supportedSamplingRates[i]) return true;
  }
  return false; // not supported
}
