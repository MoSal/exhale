/* exhaleApp.cpp - source file with main() routine for exhale application executable
 * written by C. R. Helmrich, last modified in 2023 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under the exhale Copyright License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2021 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#include "exhaleAppPch.h"
#include "basicMP4Writer.h"
#include "basicWavReader.h"
#include "loudnessEstim.h"
// #define USE_EXHALELIB_DLL (defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64))
#if USE_EXHALELIB_DLL
#include "exhaleDecl.h"
#else
#include "../lib/exhaleEnc.h"
#endif
#include "version.h"

#include <iostream>
#include <sys/stat.h>
#include <fcntl.h>
#include <time.h>
#if defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64)
#include <direct.h>
#include <windows.h>
#ifdef __MINGW32__
#include <share.h>
#endif

#define EXHALE_APP_WIN
#if defined (_MSC_VER) || defined (__INTEL_COMPILER) || defined (__MINGW32__) // || defined (__GNUC__)
#define EXHALE_APP_WCHAR
#ifndef _O_U16TEXT
# define _O_U16TEXT   0x20000
#endif
#define _ERROR1(fmt)      fwprintf_s (stderr, L##fmt)
#define _ERROR2(fmt, dat) fwprintf_s (stderr, L##fmt, dat)
#define _GETCWD _wgetcwd
#define _SOPENS _wsopen_s
#define _STRLEN  wcslen
#else
#define _ERROR1(fmt)      fprintf_s (stderr, fmt)
#define _ERROR2(fmt, dat) fprintf_s (stderr, fmt, dat)
#define _GETCWD  _getcwd
#define _SOPENS  _sopen_s
#define _STRLEN  strlen
#endif
#define EXHALE_TEXT_BLUE  (FOREGROUND_INTENSITY | FOREGROUND_BLUE | FOREGROUND_GREEN)
#define EXHALE_TEXT_PINK  (FOREGROUND_INTENSITY | FOREGROUND_BLUE | FOREGROUND_RED)
#else // Linux, MacOS, Unix
#define _ERROR1(fmt)      fprintf_s (stderr, fmt)
#define _ERROR2(fmt, dat) fprintf_s (stderr, fmt, dat)
#define _GETCWD  getcwd
#define _STRLEN  strlen

#if 0  // change this to "#if 1" to avoid garbage text in some terminals
#define EXHALE_TEXT_INIT  ""
#define EXHALE_TEXT_BLUE  ""
#define EXHALE_TEXT_PINK  ""
#else
#define EXHALE_TEXT_INIT  "\x1b[0m"
#define EXHALE_TEXT_BLUE  "\x1b[36m"
#define EXHALE_TEXT_PINK  "\x1b[35m"
#endif
#endif

// constants, experimental macros
#if LE_ACCURATE_CALC
#define EA_LOUD_INIT  16384u  // bsSamplePeakLevel = 0 & methodValue = 0
#else
#define EA_LOUD_INIT  16399u  // bsSamplePeakLevel = 0 & methodValue = 0
#endif
#define EA_LOUD_NORM -42.25f  // -100 + 57.75 of ISO 23003-4, Table A.48
#define EA_PEAK_NORM -96.33f  // 20 * log10(2^-16), 16-bit normalization
#define EA_PEAK_MIN   0.262f  // 20 * log10() + EA_PEAK_NORM = -108 dbFS
#define EA_USE_WORK_DIR    1  // 1: use working instead of app directory
#define ENABLE_STDOUT_LOAS 0  // 1: experimental LOAS packed pipe output
#define FULL_FRM_LOOKAHEAD   // on: encoder delay = zero or frame length

static const int16_t usfc2x[32] = { // 2x upsampling filter coefficients
  (83359-65536), -27563, 16273, -11344, 8541, -6708, 5403, -4419, 3647, -3025, 2514, -2088, 1730,
  -1428, 1173, -957, 775, -622, 494, -388, 300, -230, 172, -127, 91, -63, 43, -27, 16, -9, 4, -1
};

static const int16_t rsfc3x[128] = {// 3x resampling filter coefficients
  21846, 6711, (36099-32768), 0, -18000, -14370, 0, 10208, 8901, 0, -7062, -6389, 0, 5347, 4934, 0, -4258,
  -3977, 0, 3499, 3294, 0, -2937, -2780, 0, 2501, 2376, 0, -2151, -2050, 0, 1864, 1779, 0, -1623, -1551,
  0, 1417, 1355, 0, -1240, -1187, 0, 1086, 1040, 0, -952, -910, 0, 833, 797, 0, -728, -696, 0, 635, 607,
  0, -553, -528, 0, 480, 457, 0, -415, -395, 0, 358, 340, 0, -307, -291, 0, 262, 248, 0, -223, -211, 0,
  188, 177, 0, -158, -149, 0, 132, 124, 0, -109, -102, 0, 90, 84, 0, -73, -68, 0, 59, 55, 0, -47, -43,
  0, 37, 34, 0, -29, -26, 0, 22, 20, 0, -16, -15, 0, 12, 10, 0, -8, -7, 0, 5, 5, 0, -3, -3, 0, 2
};

static bool eaInitUpsampler2x (int32_t** upsampleBuffer, const uint16_t bitRateMode, const uint16_t sampleRate,
                               const uint16_t frameSize, const uint16_t numChannels)
{
  const uint16_t inLength = frameSize >> 1;
  const uint16_t chLength = inLength + (32 << 1);
  const bool useUpsampler = (frameSize > (32 << 1) && bitRateMode * 3675 > sampleRate);

  if (useUpsampler)
  {
    if ((*upsampleBuffer = (int32_t*) malloc (chLength * numChannels * sizeof (int32_t))) == nullptr) return false;

    for (uint16_t ch = 0; ch < numChannels; ch++)
    {
      memset (*upsampleBuffer + inLength + chLength * ch, 0, (chLength - inLength) * sizeof (int32_t));
    }
  }
  return useUpsampler;
}

static bool eaInitDownsampler (int32_t** resampleBuffer, const uint16_t bitRateMode, const uint16_t sampleRate,
                               const uint16_t frameSize, const uint16_t numChannels)
{
  const uint16_t inLength = (frameSize * 3u) >> 1;
  const uint16_t chLength = inLength + (frameSize >> 3);
  const bool useResampler = (frameSize >= 512 && bitRateMode <= 1 && sampleRate == 48000);

  if (useResampler)
  {
    if ((*resampleBuffer = (int32_t*) malloc (chLength * numChannels * sizeof (int32_t))) == nullptr) return false;

    for (uint16_t ch = 0; ch < numChannels; ch++)
    {
      memset (*resampleBuffer + inLength + chLength * ch, 0, (chLength - inLength) * sizeof (int32_t));
    }
  }
  return useResampler;
}

static void eaApplyUpsampler2x (int32_t* const pcmBuffer, int32_t* const upsampleBuffer,
                                const uint16_t frameSize, const uint16_t numChannels, const bool firstFrame = false)
{
  const int16_t lookahead = 32;
  const uint16_t inLength = (frameSize >> 1) + (firstFrame ? lookahead : 0);
  const uint16_t chLength = (frameSize >> 1) + (lookahead << 1);
  uint16_t ch;

  for (ch = 0; ch < numChannels; ch++) // step 1: add deinterleaved input samples to resampling buffer
  {
    int32_t* chPcmBuf = &pcmBuffer[ch];
    int32_t* chUpsBuf = &upsampleBuffer[chLength * ch];

    if (firstFrame) // construct leading sample values via extrapolation
    {
      for (int8_t i = 0; i < 32; i++) chUpsBuf[i] = (*chPcmBuf * i + (32 >> 1)) >> 5;
    }
    else
    {
      memcpy (chUpsBuf, &chUpsBuf[inLength], (chLength - inLength) * sizeof (int32_t)); // update memory
    }
    chUpsBuf += chLength - inLength;

    for (uint16_t i = inLength; i > 0; i--, chPcmBuf += numChannels, chUpsBuf++)
    {
      *chUpsBuf = *chPcmBuf; // deinterleave, store in resampling buffer
    }
  }

  for (ch = 0; ch < numChannels; ch++) // step 2: upsample, reinterleave, and save to PCM input buffer
  {
    /*in*/int32_t* chPcmBuf = &pcmBuffer[ch];
    const int32_t* chUpsBuf = &upsampleBuffer[chLength * ch + lookahead];

    for (uint16_t i = frameSize >> 1; i > 0; i--, chUpsBuf++)
    {
      int64_t r = (chUpsBuf[0] + (int64_t) chUpsBuf[1]) * (usfc2x[0] - 2 * SHRT_MIN);

      for (int16_t c = lookahead - 1; c > 0; c--)
      {
        r += (chUpsBuf[-c] + (int64_t) chUpsBuf[c + 1]) * usfc2x[c];
      }
      *chPcmBuf = *chUpsBuf; // no filtering necessary, just copy sample
      if (*chPcmBuf < MIN_VALUE_AUDIO24) *chPcmBuf = MIN_VALUE_AUDIO24;
      else
      if (*chPcmBuf > MAX_VALUE_AUDIO24) *chPcmBuf = MAX_VALUE_AUDIO24;
      chPcmBuf += numChannels;

      *chPcmBuf = int32_t ((r - 2 * SHRT_MIN) >> 17);  // interp. sample
      if (*chPcmBuf < MIN_VALUE_AUDIO24) *chPcmBuf = MIN_VALUE_AUDIO24;
      else
      if (*chPcmBuf > MAX_VALUE_AUDIO24) *chPcmBuf = MAX_VALUE_AUDIO24;
      chPcmBuf += numChannels;
    }
  }
}

static void eaApplyDownsampler (int32_t* const pcmBuffer, int32_t* const resampleBuffer,
                                const uint16_t frameSize, const uint16_t numChannels, const bool firstFrame = false)
{
  const int16_t lookahead = frameSize >> 4;
  const uint16_t inLength = ((frameSize * 3u) >> 1) + (firstFrame ? lookahead : 0);
  const uint16_t chLength = ((frameSize * 3u) >> 1) + (lookahead << 1);
  uint16_t ch;

  for (ch = 0; ch < numChannels; ch++) // step 1: add deinterleaved input samples to resampling buffer
  {
    int32_t* chPcmBuf = &pcmBuffer[ch];
    int32_t* chResBuf = &resampleBuffer[chLength * ch];

    if (firstFrame) // construct leading sample values via extrapolation
    {
      memset (chResBuf, 0, (lookahead - 32) * sizeof (int32_t));
      for (int8_t i = 0; i < 32; i++) chResBuf[lookahead - 32 + i] = (*chPcmBuf * i + (32 >> 1)) >> 5;
    }
    else
    {
      memcpy (chResBuf, &chResBuf[inLength], (chLength - inLength) * sizeof (int32_t)); // update memory
    }
    chResBuf += chLength - inLength;

    for (uint16_t i = inLength; i > 0; i--, chPcmBuf += numChannels, chResBuf++)
    {
      *chResBuf = *chPcmBuf; // deinterleave, store in resampling buffer
    }
  }

  for (ch = 0; ch < numChannels; ch++) // step 2: resample, reinterleave, and save to PCM input buffer
  {
    /*in*/int32_t* chPcmBuf = &pcmBuffer[ch];
    const int32_t* chResBuf = &resampleBuffer[chLength * ch + lookahead];

    for (uint16_t i = frameSize >> 1; i > 0; i--, chResBuf += 3)
    {
      int64_t r1 = (int64_t) chResBuf[0] * (rsfc3x[0] - 2 * SHRT_MIN) - (chResBuf[-1] + (int64_t) chResBuf[1]) * SHRT_MIN +
                   (int64_t) chResBuf[-lookahead] + (int64_t) chResBuf[lookahead];
      int64_t r2 = (chResBuf[1] + (int64_t) chResBuf[2]) * (rsfc3x[1] - 2 * SHRT_MIN);

      for (int16_t c = lookahead - 1; c > 0; c--)
      {
        r1 += (chResBuf[-c] + (int64_t) chResBuf[c]) * rsfc3x[c << 1];
        r2 += (chResBuf[1 - c] + (int64_t) chResBuf[c + 2]) * rsfc3x[(c << 1) + 1];
      }
      *chPcmBuf = int32_t ((r1 - 2 * SHRT_MIN) >> 17); // lowpass sample
      if (*chPcmBuf < MIN_VALUE_AUDIO24) *chPcmBuf = MIN_VALUE_AUDIO24;
      else
      if (*chPcmBuf > MAX_VALUE_AUDIO24) *chPcmBuf = MAX_VALUE_AUDIO24;
      chPcmBuf += numChannels;

      *chPcmBuf = int32_t ((r2 - 2 * SHRT_MIN) >> 17); // interp. sample
      if (*chPcmBuf < MIN_VALUE_AUDIO24) *chPcmBuf = MIN_VALUE_AUDIO24;
      else
      if (*chPcmBuf > MAX_VALUE_AUDIO24) *chPcmBuf = MAX_VALUE_AUDIO24;
      chPcmBuf += numChannels;
    }
  }
}

static uint16_t eaApplyLevelNorm (int32_t* /*o*/ pcmBuffer, uint16_t* oldLoudness, const uint16_t currLoudness,
                                  const uint16_t frameSize, const uint16_t numChannels)
{
  const int64_t gainOld = __min (MAX_VALUE_AUDIO24, int64_t (0.5f + MAX_VALUE_AUDIO24 * pow (10.0f, (77.0f - *oldLoudness / 512.f) / 20.0f)));
  const int64_t gainNew = __min (MAX_VALUE_AUDIO24, int64_t (0.5f + MAX_VALUE_AUDIO24 * pow (10.0f, (77.0f - currLoudness / 512.f) / 20.0f)));
  const int64_t gRShift = 32 + __min (3, frameSize >> 10);
  const int64_t gOffset = (int64_t) 1 << (gRShift - 1);
  uint16_t ch, i;

  if ((frameSize & (frameSize - 1)) || !pcmBuffer) return 0;

  for (i = frameSize; i > 0; i--)
  {
    const int64_t gainI = gainOld * i + gainNew * (frameSize - i);

    for (ch = numChannels; ch > 0; ch--, pcmBuffer++)
    {
      *pcmBuffer = int32_t ((gOffset + *pcmBuffer * gainI) >> gRShift);
    }
  }
  return (*oldLoudness = currLoudness);
}

#if ENABLE_STDOUT_LOAS
static uint16_t eaInitLoasHeader (uint8_t* const loasHeader, // sets up LATM/LOAS header, returns payload offset
                                  const uint8_t* const ascUcBuf, const uint32_t ascUcSize)
{
  if (ascUcSize == 0) return 0;

  loasHeader[0] = 0x56; // 11-bit sync. word
  loasHeader[1] = 0xE0; // AudioSyncStream()
  // audioMuxLengthBytes will be placed here
  loasHeader[3] = 0x20; // AudioMuxElement()
  loasHeader[4] = 0x00;
  memcpy (loasHeader + 5, ascUcBuf, ascUcSize); // 3 trailing bits
  loasHeader[4 + ascUcSize] &= 0xE0;
  loasHeader[4 + ascUcSize] |= 0x3; // first 5 bits after ASC + UC
  loasHeader[5 + ascUcSize] = 0xFC;

  return uint16_t (6 + ascUcSize);
}

static uint32_t eaWriteLoasFrame (const int fileHandle, uint8_t* const loasHeader, const uint16_t payloadOffset,
                                  const uint8_t* const auBuffer, const uint32_t auSize)
{
  const uint32_t audioMuxLengthBytes = payloadOffset + 1 + auSize / 255 + auSize;
  uint32_t size = payloadOffset, tmp = auSize;

  if (audioMuxLengthBytes > 8191) return 0;

  loasHeader[1] = uint8_t (audioMuxLengthBytes >> 8) | 0xE0; // set AudioSyncStream()
  loasHeader[2] = uint8_t (audioMuxLengthBytes & UCHAR_MAX);

  loasHeader[size++] = (uint8_t) __min (UCHAR_MAX, tmp); // write PayloadLengthInfo()
  while (tmp >= UCHAR_MAX)
  {
    tmp -= UCHAR_MAX;
    loasHeader[size++] = (uint8_t) __min (UCHAR_MAX, tmp);
  }
    if ( (uint32_t) _WRITE (fileHandle, loasHeader, size) != size) return 0;
  return (uint32_t) _WRITE (fileHandle, auBuffer, auSize); // then write PayloadMux()
}
#endif // ENABLE_STDOUT_LOAS

// main routine
#ifdef EXHALE_APP_WCHAR
extern "C" int wmain (const int argc, wchar_t* argv[])
#else
int main (const int argc, char* argv[])
#endif
{
  if (argc <= 0) return argc; // for safety

  const bool readStdin = (argc == 3 || argc == 5);
  BasicWavReader wavReader;
  int32_t* inPcmData = nullptr;  // 24-bit WAVE audio input buffer
  int32_t* inPcmRsmp = nullptr;  // temporary buffer for resampler
  uint8_t* outAuData = nullptr;  // access unit (AU) output buffer
  int   inFileHandle = -1, outFileHandle = -1;
  uint32_t loudStats = EA_LOUD_INIT;  // valid empty loudness data
#ifdef EXHALE_APP_WCHAR
# if EA_USE_WORK_DIR
  wchar_t*      currPath = nullptr;
# endif
  const wchar_t* exePath = argv[0];
#else
# if EA_USE_WORK_DIR
  char*         currPath = nullptr;
# endif
  const char*    exePath = argv[0];
#endif
  uint16_t i, exePathEnd = 0, loudMemory = 0; // for LUFS leveling
  uint16_t compatibleExtensionFlag = 0; // 0: disabled, 1: enabled
  uint16_t coreSbrFrameLengthIndex = 1; // 0: 768, 1: 1024 samples
  uint16_t variableCoreBitRateMode = 3; // 0: lowest... 9: highest
  uint8_t  zeroDelayForSbrEncoding = (argc >= 5 && (argv[2][0] == 's' || argv[2][0] == 'S') && argv[2][1] == 0 ? 1 : 0);
#if ENABLE_STDOUT_LOAS
  const bool writeStdout = (zeroDelayForSbrEncoding != 0 && argv[1][0] >= 'a' && argv[argc - 1][0] == '-' && argv[argc - 1][1] == 0);
  uint16_t loasMuxOffset = 0;
  uint8_t loasHeader[64] = {0};
#endif
  bool  enableLufsLevel = (argc >= 5 && (argv[2][0] == 'l' || argv[2][0] == 'L') && argv[2][1] == 0);
#ifdef EXHALE_APP_WIN
  const HANDLE hConsole = GetStdHandle (STD_OUTPUT_HANDLE);
  CONSOLE_SCREEN_BUFFER_INFO csbi;
#endif

  for (i = 0; (exePath[i] != 0) && (i < USHRT_MAX); i++)
  {
#ifdef EXHALE_APP_WIN
    if (exePath[i] == '\\') exePathEnd = i + 1;
#else
    if (exePath[i] == '/' ) exePathEnd = i + 1;
#endif
  }
#ifdef EXHALE_APP_WCHAR
  _setmode (_fileno (stderr), _O_U16TEXT);

  const wchar_t* const exeFileName = exePath + exePathEnd;
#else
  const char* const exeFileName = exePath + exePathEnd;
#endif
  if ((exeFileName[0] == 0) || (i == USHRT_MAX))
  {
    _ERROR1 (" ERROR reading executable name or path: the string is invalid!\n\n");

    return 32768;  // bad executable string
  }

  // print program header with compile info in plain text if we pass -V
#ifdef EXHALE_APP_WCHAR
  if ((argc > 1) && (wcscmp (argv[1], L"-V") == 0 || wcscmp (argv[1], L"-v") == 0))
#else
  if ((argc > 1) && (strcmp (argv[1], "-V") == 0 || strcmp (argv[1], "-v") == 0))
#endif
  {
#if defined (__arm__) || defined (__aarch64__) || defined (__arm64__)
    fprintf_s (stdout, "exhale %s.%s%s (ARM",
#elif defined (_WIN64) || defined (WIN64) || defined (_LP64) || defined (__LP64__) || defined (__x86_64) || defined (__x86_64__)
    fprintf_s (stdout, "exhale %s.%s%s (x64",
#else // 32-bit OS
    fprintf_s (stdout, "exhale %s.%s%s (x86",
#endif
               EXHALELIB_VERSION_MAJOR, EXHALELIB_VERSION_MINOR, EXHALELIB_VERSION_BUGFIX);
#ifdef EXHALE_APP_WCHAR
    if (wcscmp (argv[1], L"-V") == 0)
#else
    if (strcmp (argv[1], "-V") == 0)
#endif
    {
      char fts[] = __TIMESTAMP__; // append month and year of file time
      fts[7] = 0;
#ifdef EXHALE_APP_WCHAR
      fprintf_s (stdout, ", Unicode");
#endif
      fprintf_s (stdout, ", %s %s)\n", &fts[4], &fts[sizeof (fts) - 5]);
    }
    else fprintf_s (stdout, ")\n");

    return 0;
  }

  // print program header with compile info
#if ENABLE_STDOUT_LOAS
  if (writeStdout)
  {
# ifdef EXHALE_APP_WCHAR
    wchar_t dateStr[12] = {0};

    mbstowcs_s (nullptr, dateStr, 12, __DATE__, _TRUNCATE);
    _ERROR1 ("\n  ----------------------------------------------------------------------\n");
    _ERROR2 (" | exhale (stdout mode, built on %s) - written by C.R.Helmrich |\n", dateStr);
# else
    _ERROR1 ("\n  ----------------------------------------------------------------------\n");
    _ERROR2 (" | exhale (stdout mode, built on %s) - written by C.R.Helmrich |\n", __DATE__);
# endif
    _ERROR1 ("  ----------------------------------------------------------------------\n\n");
  }
  else
  {
#endif
  fprintf_s (stdout, "\n  ---------------------------------------------------------------------\n");
  fprintf_s (stdout, " | ");
#ifdef EXHALE_APP_WIN
  GetConsoleScreenBufferInfo (hConsole, &csbi); // save the text color
  SetConsoleTextAttribute (hConsole, EXHALE_TEXT_PINK); fprintf_s (stdout, "exhale");
  SetConsoleTextAttribute (hConsole, csbi.wAttributes); fprintf_s (stdout, " - ");
  SetConsoleTextAttribute (hConsole, EXHALE_TEXT_PINK); fprintf_s (stdout, "e");
  SetConsoleTextAttribute (hConsole, csbi.wAttributes); fprintf_s (stdout, "codis e");
  SetConsoleTextAttribute (hConsole, EXHALE_TEXT_PINK); fprintf_s (stdout, "x");
  SetConsoleTextAttribute (hConsole, csbi.wAttributes); fprintf_s (stdout, "tended ");
  SetConsoleTextAttribute (hConsole, EXHALE_TEXT_PINK); fprintf_s (stdout, "h");
  SetConsoleTextAttribute (hConsole, csbi.wAttributes); fprintf_s (stdout, "igh-efficiency ");
  SetConsoleTextAttribute (hConsole, EXHALE_TEXT_PINK); fprintf_s (stdout, "a");
  SetConsoleTextAttribute (hConsole, csbi.wAttributes); fprintf_s (stdout, "nd ");
  SetConsoleTextAttribute (hConsole, EXHALE_TEXT_PINK); fprintf_s (stdout, "l");
  SetConsoleTextAttribute (hConsole, csbi.wAttributes); fprintf_s (stdout, "ow-complexity ");
  SetConsoleTextAttribute (hConsole, EXHALE_TEXT_PINK); fprintf_s (stdout, "e");
  SetConsoleTextAttribute (hConsole, csbi.wAttributes); fprintf_s (stdout, "ncoder |\n");
#else
  fprintf_s (stdout, EXHALE_TEXT_PINK "exhale");
  fprintf_s (stdout, EXHALE_TEXT_INIT " - ");
  fprintf_s (stdout, EXHALE_TEXT_PINK "e");
  fprintf_s (stdout, EXHALE_TEXT_INIT "codis e");
  fprintf_s (stdout, EXHALE_TEXT_PINK "x");
  fprintf_s (stdout, EXHALE_TEXT_INIT "tended ");
  fprintf_s (stdout, EXHALE_TEXT_PINK "h");
  fprintf_s (stdout, EXHALE_TEXT_INIT "igh-efficiency ");
  fprintf_s (stdout, EXHALE_TEXT_PINK "a");
  fprintf_s (stdout, EXHALE_TEXT_INIT "nd ");
  fprintf_s (stdout, EXHALE_TEXT_PINK "l");
  fprintf_s (stdout, EXHALE_TEXT_INIT "ow-complexity ");
  fprintf_s (stdout, EXHALE_TEXT_PINK "e");
  fprintf_s (stdout, EXHALE_TEXT_INIT "ncoder |\n");
#endif
  fprintf_s (stdout, " |                                                                     |\n");
#if defined (__arm__) || defined (__aarch64__) || defined (__arm64__)
  fprintf_s (stdout, " | version %s.%s%s (ARM, built on %s) - written by C.R.Helmrich |\n",
#elif defined (_WIN64) || defined (WIN64) || defined (_LP64) || defined (__LP64__) || defined (__x86_64) || defined (__x86_64__)
  fprintf_s (stdout, " | version %s.%s%s (x64, built on %s) - written by C.R.Helmrich |\n",
#else // 32-bit OS
  fprintf_s (stdout, " | version %s.%s%s (x86, built on %s) - written by C.R.Helmrich |\n",
#endif
             EXHALELIB_VERSION_MAJOR, EXHALELIB_VERSION_MINOR, EXHALELIB_VERSION_BUGFIX, __DATE__);
  fprintf_s (stdout, "  ---------------------------------------------------------------------\n\n");
#if ENABLE_STDOUT_LOAS
  }
#endif

  // check arg. list, print usage if needed
  if ((argc < 3) || (argc > 6) || (argc > 1 && argv[1][1] != 0))
  {
    fprintf_s (stdout, " Copyright 2018-2023 C.R.Helmrich, project ecodis. See License.htm for details.\n\n");

    fprintf_s (stdout, " This software is made available under the exhale Copyright License and comes\n");
    fprintf_s (stdout, " with ABSOLUTELY NO WARRANTY. This software may be subject to other third-party\n");
    fprintf_s (stdout, " rights, including patent rights. No such rights are granted under this License.\n\n");
#ifdef EXHALE_APP_WIN
    SetConsoleTextAttribute (hConsole, EXHALE_TEXT_BLUE); fprintf_s (stdout, " Usage:\t");
    SetConsoleTextAttribute (hConsole, csbi.wAttributes);
#else
    fprintf_s (stdout, EXHALE_TEXT_BLUE " Usage:\t" EXHALE_TEXT_INIT);
#endif
#ifdef EXHALE_APP_WCHAR
    fwprintf_s (stdout, L"%s preset [inputWaveFile.wav] outputMP4File.m4a\n\n where\n\n", exeFileName);
#else
    fprintf_s (stdout, "%s preset [inputWaveFile.wav] outputMP4File.m4a\n\n where\n\n", exeFileName);
#endif
#ifdef EXHALE_APP_WIN
    fprintf_s (stdout, " preset\t=  # (0-9)  low-complexity ISO/MPEG-D Extended HE-AAC at 16ú#+48 kbit/s\n");
    fprintf_s (stdout, " \t     (a-g)  low-complexity Extended HE-AAC using eSBR at 12ú#+36 kbit/s\n");
#else
    fprintf_s (stdout, " preset\t=  # (0-9)  low-complexity ISO/MPEG-D Extended HE-AAC at 16*#+48 kbit/s\n");
    fprintf_s (stdout, " \t     (a-g)  low-complexity Extended HE-AAC using eSBR at 12*#+36 kbit/s\n");
#endif
    fprintf_s (stdout, "\n inputWaveFile.wav  lossless WAVE audio input, read from stdin if not specified\n\n");
    fprintf_s (stdout, " outputMP4File.m4a  encoded MPEG-4 bit-stream, extension should be .m4a or .mp4\n\n\n");
#ifdef EXHALE_APP_WIN
    SetConsoleTextAttribute (hConsole, EXHALE_TEXT_BLUE); fprintf_s (stdout, " Notes:\t");
    SetConsoleTextAttribute (hConsole, csbi.wAttributes);
#else
    fprintf_s (stdout, EXHALE_TEXT_BLUE " Notes:\t" EXHALE_TEXT_INIT);
#endif
    fprintf_s (stdout, "The above bit-rates are for stereo and change for mono or multichannel.\n");
#if !EA_USE_WORK_DIR
    if (exePathEnd > 0)
    {
# ifdef EXHALE_APP_WIN
#  ifdef EXHALE_APP_WCHAR
      fwprintf_s (stdout, L" \tUse filename prefix .\\ for the current directory if this executable was\n\tcalled with a path (call: %s).\n", exePath);
#  else
      fprintf_s (stdout, " \tUse filename prefix .\\ for the current directory if this executable was\n\tcalled with a path (call: %s).\n", exePath);
#  endif
# else
      fprintf_s (stdout, " \tUse filename prefix ./ for the current directory if this executable was\n\tcalled with a path (call: %s).\n", exePath);
# endif
    }
#endif
    return 0;  // no arguments, which is OK
  }

  // check preset mode, derive coder config
  if ((*argv[1] >= '0' && *argv[1] <= '9') || (*argv[1] >= 'a' && *argv[1] <= 'g'))
  {
    i = (uint16_t) argv[1][0];
    compatibleExtensionFlag = (i & 0x40) >> 6;
    coreSbrFrameLengthIndex = (i > 0x60 ? 5 : (i & 0x20) >> 5);
    variableCoreBitRateMode = (i & 0x0F) - (i >> 6);
  }
  else if (*argv[1] == '#') // default mode
  {
#if ENABLE_STDOUT_LOAS
    if (writeStdout)
    {
      _ERROR2 (" Default preset is specified, encoding to low-complexity xHE-AAC, preset mode %d\n\n", variableCoreBitRateMode);
    }
    else
#endif
    fprintf_s (stdout, " Default preset is specified, encoding to low-complexity xHE-AAC, preset mode %d\n\n", variableCoreBitRateMode);
  }
  else
  {
    _ERROR2 (" ERROR reading preset mode: character %s is not supported! Use 0-9 or a-g.\n\n", argv[1]);

    return 16384; // preset isn't supported
  }

  const bool enableSbrCoding = (coreSbrFrameLengthIndex >= 3); // SBR coding flag
  const unsigned frameLength = (3 + coreSbrFrameLengthIndex) << 8; // dec. output
  const unsigned startLength = (frameLength * 25) >> 4; // encoder PCM look-ahead

  if (readStdin) // configure stdin
  {
#ifdef EXHALE_APP_WIN
    inFileHandle = _fileno (stdin);
    if (_setmode (inFileHandle, _O_RDONLY | _O_BINARY) == -1)
    {
      _ERROR1 (" ERROR while trying to set stdin to binary mode! Has stdin been closed?\n\n");
      inFileHandle = -1;

      goto mainFinish; // stdin setup error
    }
#else
    inFileHandle = fileno (stdin);
#endif
  }
  else // argc = 4, open input file
  {
#ifdef EXHALE_APP_WCHAR
    const wchar_t* inFileName = argv[argc - 2];
#else
    const char* inFileName = argv[argc - 2];
#endif
    uint16_t    inPathEnd  = 0;

    for (i = 0; (inFileName[i] != 0) && (i < USHRT_MAX); i++)
    {
#ifdef EXHALE_APP_WIN
      if (inFileName[i] == '\\') inPathEnd = i + 1;
#else
      if (inFileName[i] == '/' ) inPathEnd = i + 1;
#endif
    }
    if ((inFileName[0] == 0) || (i == USHRT_MAX))
    {
      _ERROR1 (" ERROR reading input file name or path: the string is invalid!\n\n");

      goto mainFinish;  // bad input string
    }

    if (inPathEnd == 0) // name has no path
    {
#if EA_USE_WORK_DIR
# ifdef __linux__
      if ((currPath == nullptr) && (currPath = _GETCWD (NULL, 0)) != nullptr)
# else
      if ((currPath == nullptr) && (currPath = _GETCWD (NULL, 1)) != nullptr)
# endif
      {
        exePath = currPath;
        exePathEnd = (uint16_t) __min (USHRT_MAX - 1, _STRLEN (currPath));
# ifdef EXHALE_APP_WIN
        if (currPath[exePathEnd - 1] != '\\') currPath[exePathEnd++] = '\\';
# else
        if (currPath[exePathEnd - 1] != '/' ) currPath[exePathEnd++] = '/';
# endif
      }
#endif
#ifdef EXHALE_APP_WCHAR
      inFileName = (const wchar_t*) malloc ((exePathEnd + i + 1) * sizeof (wchar_t));  // 0-terminated
      memcpy ((void*) inFileName, exePath, exePathEnd * sizeof (wchar_t));  // prepend executable path
      memcpy ((void*)(inFileName + exePathEnd), argv[argc - 2], (i + 1) * sizeof (wchar_t));// to name
#else
      inFileName = (const char*) malloc ((exePathEnd + i + 1) * sizeof (char)); // 0-terminated string
      memcpy ((void*) inFileName, exePath, exePathEnd * sizeof (char)); // prepend executable path ...
      memcpy ((void*)(inFileName + exePathEnd), argv[argc - 2], (i + 1) * sizeof (char)); //...to name
#endif
    }

#ifdef EXHALE_APP_WIN
    if (_SOPENS (&inFileHandle, inFileName, _O_RDONLY | _O_SEQUENTIAL | _O_BINARY, _SH_DENYWR, _S_IREAD) != 0)
#else
    if ((inFileHandle = ::open (inFileName, O_RDONLY, 0666)) == -1)
#endif
    {
      _ERROR2 (" ERROR while trying to open input file %s! Does it already exist?\n\n", inFileName);
      inFileHandle = -1;
      if (inPathEnd == 0) free ((void*) inFileName);

      goto mainFinish;  // input file error
    }
    if (inPathEnd == 0) free ((void*) inFileName);
  }

#ifdef EXHALE_APP_WIN
  if ((wavReader.open (inFileHandle, startLength, readStdin ? LLONG_MAX : _filelengthi64 (inFileHandle)) != 0) ||
#else
  if ((wavReader.open (inFileHandle, startLength, readStdin ? LLONG_MAX : lseek (inFileHandle, 0, 2 /*SEEK_END*/)) != 0) ||
#endif
      (wavReader.getSampleRate () >= 1000 && wavReader.getSampleRate () < 22050 && enableSbrCoding) || (wavReader.getNumChannels () >= (enableSbrCoding ? 3u : 7u)))
  {
    _ERROR1 (" ERROR while trying to open WAVE file: invalid or unsupported audio format!\n\n");

    if (wavReader.getSampleRate () >= 1000 && wavReader.getSampleRate () < 22050 && enableSbrCoding)
    {
      _ERROR2 (" The sampling rate is %d kHz but encoding using eSBR requires at least 22 kHz.\n\n", wavReader.getSampleRate () / 1000);
    }
    if (wavReader.getNumChannels () >= 3 && enableSbrCoding)
    {
      _ERROR2 (" The channel count is %d but exhale can't encode multichannel audio with eSBR.\n\n", wavReader.getNumChannels ());
    }
    i = 8192; // return value

    goto mainFinish; // audio format invalid
  }
#if ENABLE_STDOUT_LOAS
  else if (writeStdout)  // configure stdout
  {
    if (wavReader.getNumChannels () != 2)
    {
      _ERROR1 (" ERROR during encoding! Input audio must be stereo for encoding to stdout!\n\n");
      i = 8192; // return value

      goto mainFinish; // stdout audio error
    }

    fflush (stdout);
# ifdef EXHALE_APP_WIN
    outFileHandle = _fileno (stdout);
    if (_setmode (outFileHandle, _O_BINARY) == -1)
    {
      _ERROR1 (" ERROR while trying to set stdout to binary mode! Has stdout been closed?\n\n");
      outFileHandle = -1;

      goto mainFinish; // stdout setup error
    }
# else
    outFileHandle = fileno (stdout);
# endif

    if ((wavReader.getSampleRate () > 32100 + (unsigned) variableCoreBitRateMode * 12000 + (variableCoreBitRateMode >> 2) * 3900) &&
        (variableCoreBitRateMode > 1 || wavReader.getSampleRate () != 48000) && !enableSbrCoding)
    {
      i = (variableCoreBitRateMode > 4 ? 96 : __min (64, 32 + variableCoreBitRateMode * 12));
# ifdef EXHALE_APP_WCHAR
      fwprintf_s (stderr, L" ERROR during encoding! Input sample rate must be <=%d kHz for preset mode %d!\n\n", i, variableCoreBitRateMode);
# else
      fprintf_s (stderr, " ERROR during encoding! Input sample rate must be <=%d kHz for preset mode %d!\n\n", i, variableCoreBitRateMode);
# endif
      i = 4096; // return value

      goto mainFinish; // ask for resampling
    }
    if ((wavReader.getSampleRate () > 32000) && !enableSbrCoding && (variableCoreBitRateMode <= 1))
    {
      if (wavReader.getSampleRate () == 48000)
      {
        _ERROR2 (" NOTE: Downsampling the input audio from 48 kHz to 32 kHz with preset mode %d\n\n", variableCoreBitRateMode);
      }
      else
      {
        _ERROR2 (" WARNING: The input sampling rate should be 32 kHz or less for preset mode %d!\n\n", variableCoreBitRateMode);
      }
    }
  }
#endif
  else // WAVE OK, open output file
  {
#ifdef EXHALE_APP_WCHAR
    const wchar_t* outFileName = argv[argc - 1];
#else
    const char* outFileName = argv[argc - 1];
#endif
    uint16_t    outPathEnd  = readStdin ? 1 : 0; // no path prepends when the input is read from stdin

    for (i = 0; (outFileName[i] != 0) && (i < USHRT_MAX); i++)
    {
#ifdef EXHALE_APP_WIN
      if (outFileName[i] == '\\') outPathEnd = i + 1;
#else
      if (outFileName[i] == '/' ) outPathEnd = i + 1;
#endif
    }
    if ((outFileName[0] == 0) || (i == USHRT_MAX))
    {
      _ERROR1 (" ERROR reading output file name or path: the string is invalid!\n\n");

      goto mainFinish;  // bad output string
    }

    if ((wavReader.getSampleRate () > 32100 + (unsigned) variableCoreBitRateMode * 12000 + (variableCoreBitRateMode >> 2) * 3900) &&
        (variableCoreBitRateMode > 1 || wavReader.getSampleRate () != 48000) && !enableSbrCoding)
    {
      i = (variableCoreBitRateMode > 4 ? 96 : __min (64, 32 + variableCoreBitRateMode * 12));
#ifdef EXHALE_APP_WCHAR
      fwprintf_s (stderr, L" ERROR during encoding! Input sample rate must be <=%d kHz for preset mode %d!\n\n", i, variableCoreBitRateMode);
#else
      fprintf_s (stderr, " ERROR during encoding! Input sample rate must be <=%d kHz for preset mode %d!\n\n", i, variableCoreBitRateMode);
#endif
      i = 4096; // return value

      goto mainFinish; // ask for resampling
    }
    if ((wavReader.getSampleRate () > 32000) && !enableSbrCoding && (variableCoreBitRateMode <= 1))
    {
      if (wavReader.getSampleRate () == 48000)
      {
        fprintf_s (stdout, " NOTE: Downsampling the input audio from 48 kHz to 32 kHz with preset mode %d\n\n", variableCoreBitRateMode);
      }
      else
      {
        _ERROR2 (" WARNING: The input sampling rate should be 32 kHz or less for preset mode %d!\n\n", variableCoreBitRateMode);
      }
    }

    if (outPathEnd == 0) // name has no path
    {
#if EA_USE_WORK_DIR
# ifdef __linux__
      if ((currPath != exePath) && (currPath = _GETCWD (NULL, 0)) != nullptr)
# else
      if ((currPath != exePath) && (currPath = _GETCWD (NULL, 1)) != nullptr)
# endif
      {
        exePath = currPath;
        exePathEnd = (uint16_t) __min (USHRT_MAX - 1, _STRLEN (currPath));
# ifdef EXHALE_APP_WIN
        if (currPath[exePathEnd - 1] != '\\') currPath[exePathEnd++] = '\\';
# else
        if (currPath[exePathEnd - 1] != '/' ) currPath[exePathEnd++] = '/';
# endif
      }
#endif
#ifdef EXHALE_APP_WCHAR
      outFileName = (const wchar_t*) malloc ((exePathEnd + i + 1) * sizeof (wchar_t));  // 0-terminated
      memcpy ((void*) outFileName, exePath, exePathEnd * sizeof (wchar_t));  // prepend executable path
      memcpy ((void*)(outFileName + exePathEnd), argv[argc - 1], (i + 1) * sizeof (wchar_t));// to name
#else
      outFileName = (const char*) malloc ((exePathEnd + i + 1) * sizeof (char)); // 0-terminated string
      memcpy ((void*) outFileName, exePath, exePathEnd * sizeof (char)); // prepend executable path ...
      memcpy ((void*)(outFileName + exePathEnd), argv[argc - 1], (i + 1) * sizeof (char)); //...to name
#endif
    }

    i = (readStdin ? O_RDWR : O_WRONLY);
#ifdef EXHALE_APP_WIN
    if (_SOPENS (&outFileHandle, outFileName, i | _O_SEQUENTIAL | _O_CREAT | _O_EXCL | _O_BINARY, _SH_DENYRD, _S_IWRITE) != 0)
#else
    if ((outFileHandle = ::open (outFileName, i | O_CREAT | O_EXCL, 0666)) == -1)
#endif
    {
      _ERROR2 (" ERROR while trying to open output file %s! Does it already exist?\n\n", outFileName);
      outFileHandle = -1;
      if (outPathEnd == 0) free ((void*) outFileName);

      goto mainFinish;  // output file error
    }
    if (outPathEnd == 0) free ((void*) outFileName);
  }

  // enforce executable specific constraints
  i = __min (USHRT_MAX, wavReader.getSampleRate ());
  if (i == 22050 && enableSbrCoding) loudStats |= 0xC0000000; // fix dflt_freq_scale in SBR header for 22050-Hz input

  if ((wavReader.getNumChannels () > 3 || enableSbrCoding) && (i == 57600 || i == 38400 || i == 28800 || i == 19200)) // BL USAC
  {
#ifdef EXHALE_APP_WCHAR
    fwprintf_s (stderr, L" ERROR: exhale does not support %d-channel coding with %d Hz sampling rate.\n\n", wavReader.getNumChannels (), i);
#else
    fprintf_s (stderr, " ERROR: exhale does not support %d-channel coding with %d Hz sampling rate.\n\n", wavReader.getNumChannels (), i);
#endif
    goto mainFinish; // encoder config error
  }
  else
  {
    const unsigned numChannels = wavReader.getNumChannels ();
    const unsigned inSampDepth = wavReader.getBitDepth ();
    const unsigned sbrEncDelay = (enableSbrCoding ? 962 : 0);
    const bool enableUpsampler = eaInitUpsampler2x (&inPcmRsmp, variableCoreBitRateMode, i, frameLength, numChannels);
    const bool enableResampler = (enableSbrCoding ? false : // no 3:2 downsampling required when encoding in SBR mode
                                 eaInitDownsampler (&inPcmRsmp, variableCoreBitRateMode, i, frameLength, numChannels));
    const uint16_t firstLength = uint16_t (enableUpsampler ? (frameLength >> 1) + 32 : (enableResampler ? startLength : frameLength));
    const unsigned inFrameSize = (enableResampler ? startLength : frameLength) * sizeof (int32_t); // max buffer size
    const unsigned resampRatio = (enableResampler ? 3 : 1); // for resampling ratio
    const unsigned resampShift = (enableResampler || enableUpsampler ? 1 : 0);
#ifdef FULL_FRM_LOOKAHEAD
    const uint16_t inPadLength = uint16_t ((((frameLength << 1) - startLength) * resampRatio) >> resampShift) +
                                 (enableSbrCoding ? firstLength - (sbrEncDelay >> resampShift) - (enableUpsampler ? 32 : 0) : 0);
#endif
    const int64_t expectLength = (wavReader.getDataBytesLeft () << resampShift) / int64_t ((numChannels * inSampDepth * resampRatio) >> 3);

    if (enableUpsampler) // notify by printf
    {
#if ENABLE_STDOUT_LOAS
      if (writeStdout) // relocate to stderr
      {
# ifdef EXHALE_APP_WCHAR
        fwprintf_s (stderr, L" NOTE: Upsampling the input audio from %d kHz to %d kHz with preset mode %d\n\n",
# else
        fprintf_s (stderr, " NOTE: Upsampling the input audio from %d kHz to %d kHz with preset mode %d\n\n",
# endif
                    i / 1000, i / 500, variableCoreBitRateMode);
      }
      else
#endif
      fprintf_s (stdout, " NOTE: Upsampling the input audio from %d kHz to %d kHz with preset mode %d\n\n", i / 1000, i / 500, variableCoreBitRateMode);
    }

    if (variableCoreBitRateMode == 0)
    {
      if (enableSbrCoding)
        _ERROR1 (" WARNING: The usage of preset mode a is not recommended since the audio quality\n resulting from preset a does not reflect the full capabilities of the Extended\n HE-AAC standard. Therefore, use the lowest bit-rate modes only when necessary!\n\n");
      else
        _ERROR1 (" WARNING: The usage of preset mode 0 is not recommended since the audio quality\n resulting from preset 0 does not reflect the full capabilities of the Extended\n HE-AAC standard. Therefore, use the lowest bit-rate modes only when necessary!\n\n");
    }

    // allocate dynamic frame memory buffers
    inPcmData = (int32_t*) malloc (inFrameSize * numChannels); // max frame in size
#ifdef NO_PREROLL_DATA
    outAuData = (uint8_t*) malloc ((6144 >> 3) * numChannels); // max frame AU size
#else
    outAuData = (uint8_t*) malloc ((9984 >> 3) * numChannels); // max frame AU size
#endif
    if ((inPcmData == nullptr) || (outAuData == nullptr))
    {
      _ERROR1 (" ERROR while trying to allocate dynamic memory! Not enough free RAM available!\n\n");
      i = 2048; // return value

      goto mainFinish; // memory alloc error
    }

#ifdef FULL_FRM_LOOKAHEAD
    memset (inPcmData, 0, inPadLength * numChannels * sizeof (int32_t)); // padding

    if (inPadLength + wavReader.read (inPcmData + inPadLength * numChannels, firstLength - inPadLength) != firstLength)
#else
    if (wavReader.read (inPcmData, firstLength) != firstLength) // full first frame
#endif
    {
      _ERROR1 (" ERROR while trying to encode input audio data! The audio stream is too short!\n\n");
      i = 1024; // return value

      goto mainFinish; // audio is too short
    }
    else // start coding loop, show progress
    {
      const unsigned sampleRate  = (wavReader.getSampleRate () << resampShift) / resampRatio;
      const bool userIndepPeriod = (argc >= 5 && argv[3][0] > '0' && argv[3][0] <= '9' && argv[3][1] >= '0' && argv[3][1] <= '9' && argv[3][2] == 0);
      const unsigned indepPeriod = (userIndepPeriod ? 10 * (argv[3][0] - 48) + (argv[3][1] - 48) : (sampleRate < 48000 ? sampleRate - 320u : 50u << 10u) / frameLength);
#if ENABLE_STDOUT_LOAS
      const unsigned mod3Percent = (writeStdout ? 0 : unsigned ((expectLength * (3 + (coreSbrFrameLengthIndex & 3))) >> 17));
      uint32_t byteCount = 0, bw = (numChannels < 7 ? loudStats | (writeStdout ? 0x4A0022CB /*-23 LUFS*/ : 0) : 0);
#else
      const unsigned mod3Percent = unsigned ((expectLength * (3 + (coreSbrFrameLengthIndex & 3))) >> 17);
      uint32_t byteCount = 0, bw = (numChannels < 7 ? loudStats : 0);
#endif
      uint32_t br, bwMax = 0, bwTmp = 0; // br will hold bytes read and/or bit-rate
      uint32_t headerRes = 0;
      // initialize LoudnessEstimator object
      LoudnessEstimator loudnessEst (inPcmData, 24 /*bit*/, sampleRate, numChannels);
      // open & prepare ExhaleEncoder object
#if USE_EXHALELIB_DLL
      ExhaleEncAPI&  exhaleEnc = *exhaleCreate (inPcmData, outAuData, sampleRate, numChannels, frameLength, indepPeriod, variableCoreBitRateMode
#else
      ExhaleEncoder  exhaleEnc (inPcmData, outAuData, sampleRate, numChannels, frameLength, indepPeriod, variableCoreBitRateMode
#endif
                                + (enableUpsampler && (variableCoreBitRateMode < 9) ? 1 : 0)
#if !RESTRICT_TO_AAC
                              , !(argc >= 5 && (argv[2][0] == 'n' || argv[2][0] == 'N') && argv[2][1] == 0), compatibleExtensionFlag > 0
#endif
                                );
      BasicMP4Writer mp4Writer; // .m4a file

      // init encoder, generate UsacConfig()
      memset (outAuData, 0, 108 * sizeof (uint8_t));  // max. allowed ASC + UC size
#ifdef FULL_FRM_LOOKAHEAD
      if (!enableSbrCoding) zeroDelayForSbrEncoding = 0;

      if (zeroDelayForSbrEncoding)
      {
        // resample PCM priming if necessary
        if (enableUpsampler) eaApplyUpsampler2x (inPcmData, inPcmRsmp, frameLength, numChannels, true);
        else
        if (enableResampler) eaApplyDownsampler (inPcmData, inPcmRsmp, frameLength, numChannels, true);
      }
      // extrapolate samples in padding region of first frame since exhaleLib can't
      // take over this job when inPadLength > 0. Improves gapless playback.
      else if (inPadLength > 0) eaExtrapolate (inPcmData, inPadLength, frameLength, numChannels, true); // fade-in

      // signal 1-frame skip and PCM priming
      outAuData[0] = 1 | zeroDelayForSbrEncoding * (uint8_t) __min (254, (firstLength - inPadLength) << (resampShift + 1));
#endif
      i = exhaleEnc.initEncoder (outAuData, &bw); // bw stores actual ASC + UC size
#ifdef FULL_FRM_LOOKAHEAD
      if ((i == 0) && (zeroDelayForSbrEncoding))
      {
        if (wavReader.read (inPcmData, (frameLength * resampRatio) >> resampShift) == 0) i = 1; // discard priming
      }
#endif

#if ENABLE_STDOUT_LOAS
      if ((i == 0) && writeStdout)
      {
        br = 0; // init frame count & header
        if ((loasMuxOffset = eaInitLoasHeader (loasHeader, outAuData, bw)) == 0) i = 1;
        else enableLufsLevel = true;
      }
      else
#endif
      if ((i |= mp4Writer.open (outFileHandle, sampleRate, numChannels, inSampDepth, frameLength,
#ifdef FULL_FRM_LOOKAHEAD
                                (frameLength << (enableSbrCoding && !zeroDelayForSbrEncoding ? 1 : 0))
#else
                                startLength + sbrEncDelay
#endif
#ifndef NO_PREROLL_DATA
                                - frameLength
#endif
                              , indepPeriod, outAuData, bw, (time (nullptr) + 2082844800) & UINT_MAX, (char) variableCoreBitRateMode)) != 0)
      {
        _ERROR2 (" ERROR while trying to initialize exhale encoder: error value %d was returned!\n\n", i);
        i <<= 2; // return value
#if USE_EXHALELIB_DLL
        exhaleDelete (&exhaleEnc);
#endif
        goto mainFinish; // coder init error
      }

      if (*argv[1] != '#') // user-def. mode
      {
#if ENABLE_STDOUT_LOAS
        if (writeStdout)  // print to stderr
        {
# ifdef EXHALE_APP_WCHAR
          fwprintf_s (stderr, L" Encoding %d-kHz %d-channel %d-bit WAVE to low-complexity xHE-AAC at %d kbit/s\n\n", sampleRate / 1000,
# else
          fprintf_s (stderr, " Encoding %d-kHz %d-channel %d-bit WAVE to low-complexity xHE-AAC at %d kbit/s\n\n", sampleRate / 1000,
# endif
                      numChannels, inSampDepth, __min (5, numChannels) * (((24 + variableCoreBitRateMode * 8) * (enableSbrCoding ? 3 : 4)) >> 2));
        }
        else
#endif
        fprintf_s (stdout, " Encoding %d-kHz %d-channel %d-bit WAVE to low-complexity xHE-AAC at %d kbit/s\n\n", sampleRate / 1000,
                   numChannels, inSampDepth, __min (5, numChannels) * (((24 + variableCoreBitRateMode * 8) * (enableSbrCoding ? 3 : 4)) >> 2));
      }
      if (!readStdin && (mod3Percent > 0))
      {
#ifdef EXHALE_APP_WIN
        SetConsoleTextAttribute (hConsole, EXHALE_TEXT_BLUE);
        fprintf_s (stdout, " Progress: ");
        SetConsoleTextAttribute (hConsole, csbi.wAttributes); // initial text color
        fprintf_s (stdout, "-");  fflush (stdout);
#else
        fprintf_s (stdout, EXHALE_TEXT_BLUE " Progress: " EXHALE_TEXT_INIT "-"); fflush (stdout);
#endif
      }

#if ENABLE_STDOUT_LOAS
      if (!readStdin && !writeStdout) // reserve space required for MP4 file header
#else
      if (!readStdin) // reserve space for MP4 file header. TODO: nasty, avoid this
#endif
      {
        if ((headerRes = (uint32_t) mp4Writer.initHeader (uint32_t (__min (UINT_MAX - startLength, expectLength)), sbrEncDelay >> 2)) < 666)
        {
          _ERROR2 ("\n ERROR while trying to write MPEG-4 bit-stream header: stopped after %d bytes!\n\n", headerRes);
          i = 3; // return value
#if USE_EXHALELIB_DLL
          exhaleDelete (&exhaleEnc);
#endif
          goto mainFinish; // writeout error
        }
      }
#if EE_MORE_MSE
      std::cout << "\n" << "headerSizeBytes " << (headerRes - 34) << "\n";
#endif
      i = 1; // for progress bar

      // resample initial frame if necessary
      if (enableUpsampler) eaApplyUpsampler2x (inPcmData, inPcmRsmp, frameLength, numChannels, !zeroDelayForSbrEncoding);
      else
      if (enableResampler) eaApplyDownsampler (inPcmData, inPcmRsmp, frameLength, numChannels, !zeroDelayForSbrEncoding);

      // initial frame, encode look-ahead AU
      if ((bw = exhaleEnc.encodeLookahead ()) < 3)
      {
        _ERROR2 ("\n ERROR while trying to create first audio frame: error value %d was returned!\n\n", bw);
        i = 2; // return value
#if USE_EXHALELIB_DLL
        exhaleDelete (&exhaleEnc);
#endif
        goto mainFinish; // coder-time error
      }
#ifdef FULL_FRM_LOOKAHEAD
      if (zeroDelayForSbrEncoding)
      {
        if (loudnessEst.addNewPcmData (frameLength))
        {
# if USE_EXHALELIB_DLL
          exhaleDelete (&exhaleEnc);
# endif
          goto mainFinish; // loudness error
        }
      }
      wavReader.read (inPcmData, (frameLength * resampRatio) >> resampShift); // discard the initial look-ahead AU

      // resample leading frame if necessary
      if (enableUpsampler) eaApplyUpsampler2x (inPcmData, inPcmRsmp, frameLength, numChannels);
      else
      if (enableResampler) eaApplyDownsampler (inPcmData, inPcmRsmp, frameLength, numChannels);

      // leading frame, actual look-ahead AU
      if ((bw = exhaleEnc.encodeFrame ()) < 3)
      {
        _ERROR2 ("\n ERROR while trying to create first audio frame: error value %d was returned!\n\n", bw);
        i = 2; // return value
# if USE_EXHALELIB_DLL
        exhaleDelete (&exhaleEnc);
# endif
        goto mainFinish; // coder-time error
      }
#endif
      bwTmp = bw;
#ifdef NO_PREROLL_DATA
      if (bwMax < bw) bwMax = bw;

      // write first AU, add frame to header
      if ((mp4Writer.addFrameAU (outAuData, bw) != (int) bw) || loudnessEst.addNewPcmData (frameLength))
      {
# if USE_EXHALELIB_DLL
        exhaleDelete (&exhaleEnc);
# endif
        goto mainFinish;   // writeout error
      }
      byteCount += bw;
#else
      if (loudnessEst.addNewPcmData (frameLength))
      {
# if USE_EXHALELIB_DLL
        exhaleDelete (&exhaleEnc);
# endif
        goto mainFinish; // estimation error
      }
#endif

      while (wavReader.read (inPcmData, (frameLength * resampRatio) >> resampShift) > 0) // read a new audio frame
      {
        // resample audio frame if necessary
        if (enableUpsampler) eaApplyUpsampler2x (inPcmData, inPcmRsmp, frameLength, numChannels);
        else
        if (enableResampler) eaApplyDownsampler (inPcmData, inPcmRsmp, frameLength, numChannels);

        // frame coding loop, encode next AU
        loudnessEst.addNewPcmData (frameLength);
        if (enableLufsLevel) eaApplyLevelNorm (inPcmData, &loudMemory, loudnessEst.getStatistics () >> 16, frameLength, numChannels);

        if ((bw = exhaleEnc.encodeFrame ()) < 3)
        {
          _ERROR2 ("\n ERROR while trying to create audio frame: error value %d was returned!\n\n", bw);
          i = 2; // return value
#if USE_EXHALELIB_DLL
          exhaleDelete (&exhaleEnc);
#endif
          goto mainFinish; // encoding error
        }
        bwTmp = (enableSbrCoding ? bw : (bwTmp + bw) >> 1u);
        if (bwMax < bwTmp) bwMax = bwTmp;
        bwTmp = bw;

        // write new AU, add frame to header
#if ENABLE_STDOUT_LOAS
        if (writeStdout)
        {
          if (eaWriteLoasFrame (outFileHandle, loasHeader, loasMuxOffset, outAuData, bw) != bw)
          {
# if USE_EXHALELIB_DLL
            exhaleDelete (&exhaleEnc);
# endif
            goto mainFinish; // writer error
          }
          if (br < UINT_MAX) br++;
        }
        else
#endif
        if (mp4Writer.addFrameAU (outAuData, bw) != (int) bw)
        {
#if USE_EXHALELIB_DLL
          exhaleDelete (&exhaleEnc);
#endif
          goto mainFinish; // writeout error
        }
        byteCount += bw;

        if (!readStdin && (mod3Percent > 0) && !(mp4Writer.getFrameCount () % mod3Percent))
        {
          if ((i++) < (enableSbrCoding ? 17 : 34))
          {
            fprintf_s (stdout, "-");  fflush (stdout);
          }
        }
      } // frame loop

      // resample the last frame if necessary
      if (enableUpsampler) eaApplyUpsampler2x (inPcmData, inPcmRsmp, frameLength, numChannels);
      else
      if (enableResampler) eaApplyDownsampler (inPcmData, inPcmRsmp, frameLength, numChannels);

      // end of coding loop, encode final AU
      loudnessEst.addNewPcmData (frameLength);
      if (enableLufsLevel) eaApplyLevelNorm (inPcmData, &loudMemory, loudnessEst.getStatistics () >> 16, frameLength, numChannels);

      if ((bw = exhaleEnc.encodeFrame ()) < 3)
      {
        _ERROR2 ("\n ERROR while trying to create audio frame: error value %d was returned!\n\n", bw);
        i = 2; // return value
#if USE_EXHALELIB_DLL
        exhaleDelete (&exhaleEnc);
#endif
        goto mainFinish; // coder-time error
      }
      bwTmp = (enableSbrCoding ? bw : (bwTmp + bw) >> 1u);
      if (bwMax < bwTmp) bwMax = bwTmp;
      bwTmp = bw;

      // write final AU, add frame to header
#if ENABLE_STDOUT_LOAS
      if (writeStdout)
      {
        if (eaWriteLoasFrame (outFileHandle, loasHeader, loasMuxOffset, outAuData, bw) != bw)
        {
# if USE_EXHALELIB_DLL
          exhaleDelete (&exhaleEnc);
# endif
          goto mainFinish; // writeout error
        }
        if (br < UINT_MAX) br++;
      }
      else
#endif
      if (mp4Writer.addFrameAU (outAuData, bw) != (int) bw)
      {
#if USE_EXHALELIB_DLL
        exhaleDelete (&exhaleEnc);
#endif
        goto mainFinish;   // writeout error
      }
      byteCount += bw;

      const int64_t actualLength = (wavReader.getDataBytesRead () << resampShift) / int64_t ((numChannels * inSampDepth * resampRatio) >> 3);
      const int64_t inFileLength = wavReader.getDataBytesRead () / int64_t ((numChannels * inSampDepth) >> 3);
      const unsigned inFrmLength = (frameLength * resampRatio) >> resampShift;
      const unsigned resampDelay = (enableUpsampler ? 32 : (enableResampler ? 64 : 0));
#ifdef FULL_FRM_LOOKAHEAD
      const unsigned flushLength = (inFileLength - resampDelay + inPadLength) % inFrmLength;
#else
      const unsigned flushLength = (inFileLength - resampDelay) % inFrmLength;
#endif
      if ((flushLength + ((startLength * resampRatio) >> resampShift) - inFrmLength + resampDelay
        - (resampDelay >> 6)/*rnd*/+ wavReader.getSampleRate () / 200 > inFrmLength
        - ((2 + sbrEncDelay * 3) >> 2)) || (flushLength == 0))  // flush last frame
      {
        memset (inPcmData, 0, inFrameSize * numChannels);

        // resample flush frame if necessary
        if (enableUpsampler) eaApplyUpsampler2x (inPcmData, inPcmRsmp, frameLength, numChannels);
        else
        if (enableResampler) eaApplyDownsampler (inPcmData, inPcmRsmp, frameLength, numChannels);

        // flush remaining audio into new AU
        // no loudnessEst.addNewPcmData call
        if (enableLufsLevel) eaApplyLevelNorm (inPcmData, &loudMemory, loudnessEst.getStatistics () >> 16, frameLength, numChannels);

        if ((bw = exhaleEnc.encodeFrame ()) < 3)
        {
          _ERROR2 ("\n ERROR while trying to create last audio frame: error value %d was returned!\n\n", bw);
          i = 2; // return value
#if USE_EXHALELIB_DLL
          exhaleDelete (&exhaleEnc);
#endif
          goto mainFinish; // encoding error
        }
        bwTmp = (enableSbrCoding ? bw : (bwTmp + bw) >> 1u);
        if (bwMax < bwTmp) bwMax = bwTmp;
        bwTmp = bw;

        // the flush AU, add frame to header
#if ENABLE_STDOUT_LOAS
        if (writeStdout)
        {
          if (eaWriteLoasFrame (outFileHandle, loasHeader, loasMuxOffset, outAuData, bw) != bw)
          {
# if USE_EXHALELIB_DLL
            exhaleDelete (&exhaleEnc);
# endif
            goto mainFinish; // writer error
          }
          if (br < UINT_MAX) br++;
        }
        else
#endif
        if (mp4Writer.addFrameAU (outAuData, bw) != (int) bw)
        {
#if USE_EXHALELIB_DLL
          exhaleDelete (&exhaleEnc);
#endif
          goto mainFinish; // writeout error
        }
        byteCount += bw;
      } // trailing frame

#if ENABLE_STDOUT_LOAS
      if (readStdin && !writeStdout) // reserve space necessary for MP4 file header
#else
      if (readStdin) // reserve space for MP4 file header (is there an easier way?)
#endif
      {
        int64_t pos = _SEEK (outFileHandle, 0, 1 /*SEEK_CUR*/);

        if ((headerRes = (uint32_t) mp4Writer.initHeader (uint32_t (__min (UINT_MAX - startLength, actualLength)), sbrEncDelay >> 2)) < 666)
        {
          _ERROR2 ("\n ERROR while trying to write MPEG-4 bit-stream header: stopped after %d bytes!\n\n", headerRes);
          i = 3; // return value
#if USE_EXHALELIB_DLL
          exhaleDelete (&exhaleEnc);
#endif
          goto mainFinish; // writeout error
        }
        // move AU data forward to make room for actual MP4 header at start of file
        br = inFrameSize * numChannels;
        while ((pos -= br) > 0) // move loop
        {
          _SEEK (outFileHandle, pos, 0 /*SEEK_SET*/);
          bw = _READ (outFileHandle, inPcmData, br);
          _SEEK (outFileHandle, pos + headerRes, 0 /*SEEK_SET*/);
          bw = _WRITE(outFileHandle, inPcmData, br);
        }
        if ((br = (uint32_t) __max (0, pos + br)) > 0) // remainder of data to move
        {
          _SEEK (outFileHandle, 0, 0 /*SEEK_SET*/);
          bw = _READ (outFileHandle, inPcmData, br);
          _SEEK (outFileHandle, headerRes, 0 /*SEEK_SET*/);
          bw = _WRITE(outFileHandle, inPcmData, br);
        }
      }
      i = 0; // no errors

      // loudness and sample peak of program
      bw = loudStats;
      loudStats = loudnessEst.getStatistics ();
#if ENABLE_STDOUT_LOAS
      if ((numChannels < 7) && !writeStdout)
#else
      if (numChannels < 7)
#endif
      {
        // quantize for loudnessInfo() reset
        const uint32_t qLoud = uint32_t (enableLufsLevel ? 139/* -23 LUFS */ : 4.0f * __max (0.0f, (loudStats >> 16) / 512.f + EA_LOUD_NORM) + 0.5f);
        const uint32_t qPeak = uint32_t (32.0f * (20.0f - 20.0f * log10 (__max (EA_PEAK_MIN, float (loudStats & USHRT_MAX))) - EA_PEAK_NORM) + 0.5f);
        // NOTE: In case of enableLufsLevel, the input peak is also the approximate
        // target peak - just as the -23.0 LUFS is the approximate target loudness.

        // recreate ASC + UC + loudness data
        bw |= (qPeak << 18) | (qLoud << 6) | 11; // measurementSystem & reliability
        memset (outAuData, 0, 108 * sizeof (uint8_t)); // max allowed ASC + UC size
        i = exhaleEnc.initEncoder (outAuData, &bw); // with finished loudnessInfo()
#ifndef NO_PREROLL_DATA
        if (i == 0)
        {
          i = __min (USHRT_MAX, wavReader.getSampleRate ());
          i = (uint16_t) mp4Writer.updateIPFs (outAuData, bw, (i == 57600 || i == 38400 || i == 28800 || i == 19200 /*BL USAC*/ ? 6 : 3));
        }
#endif
      }
      // mean & max. bit-rate of encoded AUs
#if ENABLE_STDOUT_LOAS
      if (writeStdout)
      {
        br = uint32_t (((actualLength >> 1) + 8 * (byteCount + loasMuxOffset * (int64_t) br) * sampleRate) / actualLength);
        bw = 0; // print encoding statistics

        _ERROR2 (" Done, actual average %.1f kbit/s,", (float) br * 0.001f);
        _ERROR2 (" program loudness in %.2f => out -23 LUFS\n\n", __max (3u, loudStats >> 16) / 512.f - 100.0f);
      }
      else
      {
#endif
      br = uint32_t (((actualLength >> 1) + 8 * (byteCount + 4 * (int64_t) mp4Writer.getFrameCount ()) * sampleRate) / actualLength);
      bw = uint32_t (((frameLength  >> 1) + 8 * (bwMax + 4u /* max. 2048 AU size + stsz as bit-rate */) * sampleRate) / frameLength);
      bw = mp4Writer.finishFile (br, bw, uint32_t (__min (UINT_MAX - startLength, actualLength)), (time (nullptr) + 2082844800) & UINT_MAX,
                                 (i == 0) && (numChannels < 7) ? outAuData : nullptr);
      // print out collected file statistics
      if (enableSbrCoding)
      {
        fprintf_s (stdout, " Done, actual average incl. SBR data %.2f kbit/s\n\n", (float) br * 0.001f);
      }
      else
      {
        fprintf_s (stdout, " Done, actual average %.1f kbit/s\n\n", (float) br * 0.001f);
      }
      if (numChannels < 7)
      {
        fprintf_s (stdout, " Input statistics:  File loudness %.2f LUFS,\tsample peak level %.2f dBFS\n\n",
                   __max (3u, loudStats >> 16) / 512.f - 100.0f, 20.0f * log10 (__max (EA_PEAK_MIN, float (loudStats & USHRT_MAX))) + EA_PEAK_NORM);
      }
#if ENABLE_STDOUT_LOAS
      } // writeStdout
#endif

      if (!readStdin && (actualLength != expectLength || bw != headerRes))
      {
        if (actualLength != expectLength)
#ifdef EXHALE_APP_WCHAR
        fwprintf_s (stderr, L" WARNING: %lld sample frames read but %lld sample frames expected!\n", (long long) actualLength, (long long) expectLength);
#else
        fprintf_s (stderr, " WARNING: %lld sample frames read but %lld sample frames expected!\n", (long long) actualLength, (long long) expectLength);
#endif
        if (bw != headerRes) _ERROR1 (" WARNING: The encoded MPEG-4 bit-stream is likely to be unreadable!\n");
        _ERROR1 ("\n");
      }
#if USE_EXHALELIB_DLL
      exhaleDelete (&exhaleEnc);
#endif
    } // end coding loop and stats print-out
  }

mainFinish:

  // free all dynamic memory
  MFREE (inPcmData);
  MFREE (inPcmRsmp);
#if EA_USE_WORK_DIR
  MFREE (currPath);
#endif
  MFREE (outAuData);

  // close input file
  if (inFileHandle != -1)
  {
    if (_CLOSE (inFileHandle) != 0)
    {
      if (readStdin) // stdin
      {
        _ERROR1 (" ERROR while trying to close stdin stream! Has it already been closed?\n\n");
      }
      else  // argc = 4, file
      {
        _ERROR2 (" ERROR while trying to close input file %s! Does it still exist?\n\n", argv[argc - 2]);
      }
    }
    inFileHandle = 0;
  }
  // close output file
  if (outFileHandle != -1)
  {
#if ENABLE_STDOUT_LOAS
    if (!writeStdout)
#endif
    if (_CLOSE (outFileHandle) != 0)
    {
      _ERROR2 (" ERROR while trying to close output file %s! Does it still exist?\n\n", argv[argc - 1]);
    }
    outFileHandle = 0;
  }

  return (inFileHandle | outFileHandle | i);
}
