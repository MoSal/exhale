/* exhaleApp.cpp - source file with main() routine for exhale application executable
 * written by C. R. Helmrich, last modified in 2020 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under a Modified BSD-Style License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2020 Christian R. Helmrich, project ecodis. All rights reserved.
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
#include <windows.h>
#ifdef __MINGW32__
#include <share.h>
#endif

#define EXHALE_TEXT_BLUE  (FOREGROUND_INTENSITY | FOREGROUND_BLUE | FOREGROUND_GREEN)
#define EXHALE_TEXT_PINK  (FOREGROUND_INTENSITY | FOREGROUND_BLUE | FOREGROUND_RED)
#else // Linux, MacOS, Unix
#define EXHALE_TEXT_INIT  "\x1b[0m"
#define EXHALE_TEXT_BLUE  "\x1b[36m"
#define EXHALE_TEXT_PINK  "\x1b[35m"
#endif

// constants, experimental macros
#define EA_LOUD_INIT  16399u  // bsSamplePeakLevel = 0 & methodValue = 0
#define EA_LOUD_NORM -42.25f  // -100 + 57.75 of ISO 23003-4, Table A.48
#define EA_PEAK_NORM -96.33f  // 20 * log10(2^-16), 16-bit normalization
#define EA_PEAK_MIN   0.262f  // 20 * log10() + EA_PEAK_NORM = -108 dbFS
#define IGNORE_WAV_LENGTH  0  // 1: ignore input size indicators (nasty)
#define XHE_AAC_LOW_DELAY  0  // 1: allow encoding with 768 frame length

// main routine
int main (const int argc, char* argv[])
{
  if (argc <= 0) return argc; // for safety

  const bool readStdin = (argc == 3);
  BasicWavReader wavReader;
  int32_t* inPcmData = nullptr;  // 24-bit WAVE audio input buffer
  uint8_t* outAuData = nullptr;  // access unit (AU) output buffer
  int   inFileHandle = -1, outFileHandle = -1;
  uint32_t loudStats = EA_LOUD_INIT;  // valid empty loudness data
  uint16_t i, exePathEnd = 0;
  uint16_t compatibleExtensionFlag = 0; // 0: disabled, 1: enabled
  uint16_t coreSbrFrameLengthIndex = 1; // 0: 768, 1: 1024 samples
  uint16_t variableCoreBitRateMode = 3; // 0: lowest... 9: highest
#if defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64)
  const HANDLE hConsole = GetStdHandle (STD_OUTPUT_HANDLE);
  CONSOLE_SCREEN_BUFFER_INFO csbi;
#endif

  for (i = 0; (argv[0][i] != 0) && (i < USHRT_MAX); i++)
  {
#if defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64)
    if (argv[0][i] == '\\') exePathEnd = i + 1;
#else // Linux, MacOS, Unix
    if (argv[0][i] == '/' ) exePathEnd = i + 1;
#endif
  }
  const char* const exeFileName = argv[0] + exePathEnd;

  if ((exeFileName[0] == 0) || (i == USHRT_MAX))
  {
    fprintf_s (stderr, " ERROR reading executable name or path: the string is invalid!\n\n");

    return 32768;  // bad executable string
  }

  // print program header with compile info in plain text if we pass -V
  if ((argc > 1) && (strcmp (argv[1], "-V") == 0 || strcmp (argv[1], "-v") == 0))
  {
#if defined (_WIN64) || defined (WIN64) || defined (_LP64) || defined (__LP64__) || defined (__x86_64) || defined (__x86_64__)
    fprintf_s (stdout, "exhale %s.%s%s (x64)\n",
#else // 32-bit OS
    fprintf_s (stdout, "exhale %s.%s%s (x86)\n",
#endif
               EXHALELIB_VERSION_MAJOR, EXHALELIB_VERSION_MINOR, EXHALELIB_VERSION_BUGFIX);
    return 0;
  }

  // print program header with compile info
  fprintf_s (stdout, "\n  ---------------------------------------------------------------------\n");
  fprintf_s (stdout, " | ");
#if defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64)
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
#else // Linux, MacOS, Unix
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
#if defined (_WIN64) || defined (WIN64) || defined (_LP64) || defined (__LP64__) || defined (__x86_64) || defined (__x86_64__)
  fprintf_s (stdout, " | version %s.%s%s (x64, built on %s) - written by C.R.Helmrich |\n",
#else // 32-bit OS
  fprintf_s (stdout, " | version %s.%s%s (x86, built on %s) - written by C.R.Helmrich |\n",
#endif
             EXHALELIB_VERSION_MAJOR, EXHALELIB_VERSION_MINOR, EXHALELIB_VERSION_BUGFIX, __DATE__);
  fprintf_s (stdout, "  ---------------------------------------------------------------------\n\n");

  // check arg. list, print usage if needed
  if ((argc < 3) || (argc > 4) || (argc > 1 && argv[1][1] != 0))
  {
    fprintf_s (stdout, " Copyright 2018-2020 C.R.Helmrich, project ecodis. See License.htm for details.\n\n");

    fprintf_s (stdout, " This software is being made available under a Modified BSD License and comes\n");
    fprintf_s (stdout, " with ABSOLUTELY NO WARRANTY. This software may be subject to other third-party\n");
    fprintf_s (stdout, " rights, including patent rights. No such rights are granted under this License.\n\n");
#if defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64)
    SetConsoleTextAttribute (hConsole, EXHALE_TEXT_BLUE); fprintf_s (stdout, " Usage:\t");
    SetConsoleTextAttribute (hConsole, csbi.wAttributes);
#else
    fprintf_s (stdout, EXHALE_TEXT_BLUE " Usage:\t" EXHALE_TEXT_INIT);
#endif
    fprintf_s (stdout, "%s preset [inputWaveFile.wav] outputMP4File.m4a\n\n where\n\n", exeFileName);
#if defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64)
    fprintf_s (stdout, " preset\t=  # (1-9)  low-complexity standard compliant xHE-AAC at 16ú#+48 kbit/s\n");
# if XHE_AAC_LOW_DELAY
//  fprintf_s (stdout, " \t     (a-i)  low-complexity compatible xHE-AAC with BE at 16ú#+48 kbit/s\n");
    fprintf_s (stdout, " \t     (A-I)  41ms low-delay compatible xHE-AAC with BE at 16ú#+48 kbit/s\n");
# endif
#else
    fprintf_s (stdout, " preset\t=  # (1-9)  low-complexity standard compliant xHE-AAC at 16*#+48 kbit/s\n");
# if XHE_AAC_LOW_DELAY
//  fprintf_s (stdout, " \t     (a-i)  low-complexity compatible xHE-AAC with BE at 16*#+48 kbit/s\n");
    fprintf_s (stdout, " \t     (A-I)  41ms low-delay compatible xHE-AAC with BE at 16*#+48 kbit/s\n");
# endif
#endif
    fprintf_s (stdout, "\n inputWaveFile.wav  lossless WAVE audio input, read from stdin if not specified\n\n");
    fprintf_s (stdout, " outputMP4File.m4a  encoded MPEG-4 bit-stream, extension should be .m4a or .mp4\n\n\n");
#if defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64)
    SetConsoleTextAttribute (hConsole, EXHALE_TEXT_BLUE); fprintf_s (stdout, " Notes:\t");
    SetConsoleTextAttribute (hConsole, csbi.wAttributes);
#else
    fprintf_s (stdout, EXHALE_TEXT_BLUE " Notes:\t" EXHALE_TEXT_INIT);
#endif
    fprintf_s (stdout, "The above bit-rates are for stereo and change for mono or multichannel.\n");
    if (exePathEnd > 0)
    {
#if defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64)
      fprintf_s (stdout, " \tUse filename prefix .\\ for the current directory if this executable was\n\tcalled with a path (call: %s).\n", argv[0]);
#else // Linux, MacOS, Unix
      fprintf_s (stdout, " \tUse filename prefix ./ for the current directory if this executable was\n\tcalled with a path (call: %s).\n", argv[0]);
#endif
    }
    return 0;  // no arguments, which is OK
  }

  // check preset mode, derive coder config
#if XHE_AAC_LOW_DELAY
  if ((*argv[1] >= '1' && *argv[1] <= '9') || (*argv[1] >= 'a' && *argv[1] <= 'i') || (*argv[1] >= 'A' && *argv[1] <= 'I'))
#else
  if ((*argv[1] >= '1' && *argv[1] <= '9') || (*argv[1] >= 'a' && *argv[1] <= 'i'))
#endif
  {
    i = (uint16_t) argv[1][0];
    compatibleExtensionFlag = (i & 0x40) >> 6;
    coreSbrFrameLengthIndex = (i & 0x20) >> 5;
    variableCoreBitRateMode = (i & 0x0F);
  }
  else if (*argv[1] == '#') // default mode
  {
    fprintf_s (stdout, " Default preset is specified, encoding to low-complexity xHE-AAC, preset mode %d\n\n", variableCoreBitRateMode);
  }
  else
  {
#if XHE_AAC_LOW_DELAY
    fprintf_s (stderr, " ERROR reading preset mode: character %s is not supported! Use 1-9 or A-I.\n\n", argv[1]);
#else
    fprintf_s (stderr, " ERROR reading preset mode: character %s is not supported! Please use 1-9.\n\n", argv[1]);
#endif
    return 16384; // preset isn't supported
  }

  const unsigned frameLength = (3 + coreSbrFrameLengthIndex) << 8;

  if (readStdin) // configure stdin
  {
#if defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64)
    inFileHandle = _fileno (stdin);
    if (_setmode (inFileHandle, _O_RDONLY | _O_BINARY) == -1)
    {
      fprintf_s (stderr, " ERROR while trying to set stdin to binary mode! Has stdin been closed?\n\n");
      inFileHandle = -1;

      goto mainFinish; // stdin setup error
    }
#else // Linux, MacOS, Unix
    inFileHandle = fileno (stdin);
#endif
  }
  else // argc = 4, open input file
  {
    const char* inFileName = argv[2];
    uint16_t    inPathEnd  = 0;

    for (i = 0; (inFileName[i] != 0) && (i < USHRT_MAX); i++)
    {
#if defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64)
      if (inFileName[i] == '\\') inPathEnd = i + 1;
#else // Linux, MacOS, Unix
      if (inFileName[i] == '/' ) inPathEnd = i + 1;
#endif
    }
    if ((inFileName[0] == 0) || (i == USHRT_MAX))
    {
      fprintf_s (stderr, " ERROR reading input file name or path: the string is invalid!\n\n");

      goto mainFinish;  // bad input string
    }

    if (inPathEnd == 0) // name has no path
    {
      inFileName = (const char*) malloc ((exePathEnd + i + 1) * sizeof (char)); // 0-terminated string
      memcpy ((void*) inFileName, argv[0], exePathEnd * sizeof (char)); // prepend executable path ...
      memcpy ((void*)(inFileName + exePathEnd), argv[2], (i + 1) * sizeof (char)); // ... to file name
    }

#if defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64)
    if (_sopen_s (&inFileHandle, inFileName, _O_RDONLY | _O_SEQUENTIAL | _O_BINARY, _SH_DENYWR, _S_IREAD) != 0)
#else // Linux, MacOS, Unix
    if ((inFileHandle = ::open (inFileName, O_RDONLY, 0666)) == -1)
#endif
    {
      fprintf_s (stderr, " ERROR while trying to open input file %s! Does it already exist?\n\n", inFileName);
      inFileHandle = -1;
      if (inPathEnd == 0) free ((void*) inFileName);

      goto mainFinish;  // input file error
    }
    if (inPathEnd == 0) free ((void*) inFileName);
  }

#if defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64)
  if ((wavReader.open (inFileHandle, frameLength, readStdin ? LLONG_MAX : _filelengthi64 (inFileHandle)) != 0) ||
#else // Linux, MacOS, Unix
  if ((wavReader.open (inFileHandle, frameLength, readStdin ? LLONG_MAX : lseek (inFileHandle, 0, 2 /*SEEK_END*/)) != 0) ||
#endif
      (wavReader.getNumChannels () >= 7))
  {
    fprintf_s (stderr, " ERROR while trying to open WAVE file: invalid or unsupported audio format!\n\n");
    i = 8192; // return value

    goto mainFinish; // audio format invalid
  }
  else // WAVE OK, open output file
  {
    const char* outFileName = argv[argc - 1];
    uint16_t    outPathEnd  = readStdin ? 1 : 0; // no path prepends when the input is read from stdin

    for (i = 0; (outFileName[i] != 0) && (i < USHRT_MAX); i++)
    {
#if defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64)
      if (outFileName[i] == '\\') outPathEnd = i + 1;
#else // Linux, MacOS, Unix
      if (outFileName[i] == '/' ) outPathEnd = i + 1;
#endif
    }
    if ((outFileName[0] == 0) || (i == USHRT_MAX))
    {
      fprintf_s (stderr, " ERROR reading output file name or path: the string is invalid!\n\n");

      goto mainFinish;  // bad output string
    }

    if ((variableCoreBitRateMode < 2) && (wavReader.getSampleRate () > 32000))
    {
      fprintf_s (stderr, " ERROR during encoding! Input sample rate must be <=32 kHz for preset mode %d!\n\n", variableCoreBitRateMode);
      i = 4096; // return value

      goto mainFinish; // resample to 32 kHz
    }
    if ((variableCoreBitRateMode < 4) && (wavReader.getSampleRate () > 48000))
    {
      fprintf_s (stderr, " ERROR during encoding! Input sample rate must be <=48 kHz for preset mode %d!\n\n", variableCoreBitRateMode);
      i = 4096; // return value

      goto mainFinish; // resample to 44 kHz
    }

    if (outPathEnd == 0) // name has no path
    {
      outFileName = (const char*) malloc ((exePathEnd + i + 1) * sizeof (char)); // 0-terminated string
      memcpy ((void*) outFileName, argv[0], exePathEnd * sizeof (char)); // prepend executable path ...
      memcpy ((void*)(outFileName + exePathEnd), argv[argc - 1], (i + 1) * sizeof (char)); //...to name
    }

    i = (readStdin ? O_RDWR : O_WRONLY);
#if defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64)
    if (_sopen_s (&outFileHandle, outFileName, i | _O_SEQUENTIAL | _O_CREAT | _O_EXCL | _O_BINARY, _SH_DENYRD, _S_IWRITE) != 0)
#else // Linux, MacOS, Unix
    if ((outFileHandle = ::open (outFileName, i | O_CREAT | O_EXCL, 0666)) == -1)
#endif
    {
      fprintf_s (stderr, " ERROR while trying to open output file %s! Does it already exist?\n\n", outFileName);
      outFileHandle = -1;
      if (outPathEnd == 0) free ((void*) outFileName);

      goto mainFinish;  // output file error
    }
    if (outPathEnd == 0) free ((void*) outFileName);
  }

  // enforce executable specific constraints
  i = __min (USHRT_MAX, wavReader.getSampleRate ());
  if ((wavReader.getNumChannels () > 3) && (i == 57600 || i == 51200 || i == 40000 || i == 38400 || i == 34150 ||
      i == 28800 || i == 25600 || i == 20000 || i == 19200 || i == 17075 || i == 14400 || i == 12800 || i == 9600))
  {
    fprintf_s (stderr, " ERROR: exhale does not support %d-channel coding with %d Hz sampling rate.\n\n", wavReader.getNumChannels (), i);

    goto mainFinish; // encoder config error
  }
  else
  {
    const unsigned startLength = (frameLength * 25) >> 4; // encoder PCM look-ahead
    const unsigned numChannels = wavReader.getNumChannels ();
    const unsigned inFrameSize = frameLength * sizeof (int32_t);
    const unsigned inSampDepth = wavReader.getBitDepth ();
    const int64_t expectLength = wavReader.getDataBytesLeft () / int64_t (numChannels * inSampDepth >> 3);

    // allocate dynamic frame memory buffers
    inPcmData = (int32_t*) malloc (inFrameSize * numChannels); // max frame in size
    outAuData = (uint8_t*) malloc ((6144 >> 3) * numChannels); // max frame AU size
    if ((inPcmData == nullptr) || (outAuData == nullptr))
    {
      fprintf_s (stderr, " ERROR while trying to allocate dynamic memory! Not enough free RAM available!\n\n");
      i = 2048; // return value

      goto mainFinish; // memory alloc error
    }

    if (wavReader.read (inPcmData, frameLength) != frameLength) // full first frame
    {
      fprintf_s (stderr, " ERROR while trying to encode input audio data! The audio stream is too short!\n\n");
      i = 1024; // return value

      goto mainFinish; // audio is too short
    }
    else // start coding loop, show progress
    {
      const unsigned sampleRate  = wavReader.getSampleRate ();
      const unsigned indepPeriod = (sampleRate < 48000 ? sampleRate / frameLength : 45 /*for 50-Hz video, use 50 for 60-Hz video*/);
      const unsigned mod3Percent = unsigned ((expectLength * (3 + coreSbrFrameLengthIndex)) >> 17);
      uint32_t byteCount = 0, bw = (numChannels < 7 ? loudStats : 0);
      uint32_t br, bwMax = 0; // br will be used to hold bytes read and/or bit-rate
      uint32_t headerRes = 0;
      // initialize LoudnessEstimator object
      LoudnessEstimator loudnessEst (inPcmData, 24 /*bit*/, sampleRate, numChannels);
      // open & prepare ExhaleEncoder object
#if USE_EXHALELIB_DLL
      ExhaleEncAPI&  exhaleEnc = *exhaleCreate (inPcmData, outAuData, sampleRate, numChannels, frameLength, indepPeriod, variableCoreBitRateMode +
#else
      ExhaleEncoder  exhaleEnc (inPcmData, outAuData, sampleRate, numChannels, frameLength, indepPeriod, variableCoreBitRateMode
#endif
#if !RESTRICT_TO_AAC
                              , true /*noise filling*/, compatibleExtensionFlag > 0
#endif
                                );
      BasicMP4Writer mp4Writer; // .m4a file

      // init encoder, generate UsacConfig()
      memset (outAuData, 0, 108 * sizeof (uint8_t));  // max. allowed ASC + UC size
      i = exhaleEnc.initEncoder (outAuData, &bw); // bw stores actual ASC + UC size

      if ((i |= mp4Writer.open (outFileHandle, sampleRate, numChannels, inSampDepth, frameLength, startLength,
                                indepPeriod, outAuData, bw, time (nullptr) & UINT_MAX, (char) variableCoreBitRateMode)) != 0)
      {
        fprintf_s (stderr, " ERROR while trying to initialize xHE-AAC encoder: error value %d was returned!\n\n", i);
        i <<= 2; // return value
#if USE_EXHALELIB_DLL
        exhaleDelete (&exhaleEnc);
#endif
        goto mainFinish; // coder init error
      }

      if (*argv[1] != '#') // user-def. mode
      {
        fprintf_s (stdout, " Encoding %d-kHz %d-channel %d-bit WAVE to low-complexity xHE-AAC at %d kbit/s\n\n",
                   sampleRate / 1000, numChannels, inSampDepth, __min (4, numChannels) * (24 + variableCoreBitRateMode * 8));
      }
      if (!readStdin && (mod3Percent > 0))
      {
#if defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64)
        SetConsoleTextAttribute (hConsole, EXHALE_TEXT_BLUE);
        fprintf_s (stdout, " Progress: ");
        SetConsoleTextAttribute (hConsole, csbi.wAttributes); // initial text color
        fprintf_s (stdout, "-");  fflush (stdout);
#else
        fprintf_s (stdout, EXHALE_TEXT_BLUE " Progress: " EXHALE_TEXT_INIT "-"); fflush (stdout);
#endif
      }
#if !IGNORE_WAV_LENGTH
      if (!readStdin) // reserve space for MP4 file header. TODO: nasty, avoid this
      {
        if ((headerRes = (uint32_t) mp4Writer.initHeader (uint32_t (__min (UINT_MAX - startLength, expectLength)))) < 666)
        {
          fprintf_s (stderr, "\n ERROR while trying to write MPEG-4 bit-stream header: stopped after %d bytes!\n\n", headerRes);
          i = 3; // return value
# if USE_EXHALELIB_DLL
          exhaleDelete (&exhaleEnc);
# endif
          goto mainFinish; // writeout error
        }
      }
#endif
      i = 1; // for progress bar

      // initial frame, encode look-ahead AU
      if ((bw = exhaleEnc.encodeLookahead ()) < 3)
      {
        fprintf_s (stderr, "\n ERROR while trying to create first xHE-AAC frame: error value %d was returned!\n\n", bw);
        i = 2; // return value
#if USE_EXHALELIB_DLL
        exhaleDelete (&exhaleEnc);
#endif
        goto mainFinish; // coder-time error
      }
      if (bwMax < bw) bwMax = bw;
      // write first AU, add frame to header
      if ((mp4Writer.addFrameAU (outAuData, bw) != (int) bw) || loudnessEst.addNewPcmData (frameLength))
      {
#if USE_EXHALELIB_DLL
        exhaleDelete (&exhaleEnc);
#endif
        goto mainFinish;   // writeout error
      }
      byteCount += bw;

      while (wavReader.read (inPcmData, frameLength) > 0) // read a new audio frame
      {
        // frame coding loop, encode next AU
        if ((bw = exhaleEnc.encodeFrame ()) < 3)
        {
          fprintf_s (stderr, "\n ERROR while trying to create xHE-AAC frame: error value %d was returned!\n\n", bw);
          i = 2; // return value
#if USE_EXHALELIB_DLL
          exhaleDelete (&exhaleEnc);
#endif
          goto mainFinish; // encoding error
        }
        if (bwMax < bw) bwMax = bw;
        // write new AU, add frame to header
        if ((mp4Writer.addFrameAU (outAuData, bw) != (int) bw) || loudnessEst.addNewPcmData (frameLength))
        {
#if USE_EXHALELIB_DLL
          exhaleDelete (&exhaleEnc);
#endif
          goto mainFinish; // writeout error
        }
        byteCount += bw;

        if (!readStdin && (mod3Percent > 0) && !(mp4Writer.getFrameCount () % mod3Percent))
        {
          if ((i++) < 34) // for short files
          {
            fprintf_s (stdout, "-");  fflush (stdout);
          }
        }
      } // frame loop

      // end of coding loop, encode final AU
      if ((bw = exhaleEnc.encodeFrame ()) < 3)
      {
        fprintf_s (stderr, "\n ERROR while trying to create xHE-AAC frame: error value %d was returned!\n\n", bw);
        i = 2; // return value
#if USE_EXHALELIB_DLL
        exhaleDelete (&exhaleEnc);
#endif
        goto mainFinish; // coder-time error
      }
      if (bwMax < bw) bwMax = bw;
      // write final AU, add frame to header
      if ((mp4Writer.addFrameAU (outAuData, bw) != (int) bw) || loudnessEst.addNewPcmData (frameLength))
      {
#if USE_EXHALELIB_DLL
        exhaleDelete (&exhaleEnc);
#endif
        goto mainFinish;   // writeout error
      }
      byteCount += bw;

      const int64_t actualLength = wavReader.getDataBytesRead () / int64_t (numChannels * inSampDepth >> 3);

      if (((actualLength + startLength) % frameLength) > 0) // flush trailing audio
      {
        memset (inPcmData, 0, inFrameSize * numChannels);

        // flush remaining audio into new AU
        if ((bw = exhaleEnc.encodeFrame ()) < 3)
        {
          fprintf_s (stderr, "\n ERROR while trying to create last xHE-AAC frame: error value %d was returned!\n\n", bw);
          i = 2; // return value
#if USE_EXHALELIB_DLL
          exhaleDelete (&exhaleEnc);
#endif
          goto mainFinish; // encoding error
        }
        if (bwMax < bw) bwMax = bw;
        // the flush AU, add frame to header
        if (mp4Writer.addFrameAU (outAuData, bw) != (int) bw) // no loudness update
        {
#if USE_EXHALELIB_DLL
          exhaleDelete (&exhaleEnc);
#endif
          goto mainFinish; // writeout error
        }
        byteCount += bw;
      } // trailing frame

#if !IGNORE_WAV_LENGTH
      if (readStdin) // reserve space for MP4 file header (is there an easier way?)
#endif
      {
        int64_t pos = _SEEK (outFileHandle, 0, 1 /*SEEK_CUR*/);

        if ((headerRes = (uint32_t) mp4Writer.initHeader (uint32_t (__min (UINT_MAX - startLength, actualLength)))) < 666)
        {
          fprintf_s (stderr, "\n ERROR while trying to write MPEG-4 bit-stream header: stopped after %d bytes!\n\n", headerRes);
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
      loudStats = loudnessEst.getStatistics ();
      if (numChannels < 7)
      {
        // quantize for loudnessInfo() reset
        const uint32_t qLoud = uint32_t (4.0f * __max (0.0f, (loudStats >> 16) / 512.f + EA_LOUD_NORM) + 0.5f);
        const uint32_t qPeak = uint32_t (32.0f * (20.0f - 20.0f * log10 (__max (EA_PEAK_MIN, float (loudStats & USHRT_MAX))) - EA_PEAK_NORM) + 0.5f);

        // recreate ASC + UC + loudness data
        bw = EA_LOUD_INIT | (qPeak << 18) | (qLoud << 6); // measurementSystem is 3
        memset (outAuData, 0, 108 * sizeof (uint8_t)); // max allowed ASC + UC size
        i = exhaleEnc.initEncoder (outAuData, &bw); // with finished loudnessInfo()
      }
      // mean & max. bit-rate of encoded AUs
      br = uint32_t (((actualLength >> 1) + 8 * (byteCount + 4 * (int64_t) mp4Writer.getFrameCount ()) * sampleRate) / actualLength);
      bw = uint32_t (((frameLength  >> 1) + 8 * (bwMax + 4u /* maximum AU size + stsz as a bit-rate */) * sampleRate) / frameLength);
      bw = mp4Writer.finishFile (br, bw, uint32_t (__min (UINT_MAX - startLength, actualLength)), time (nullptr) & UINT_MAX,
                                 (i == 0) && (numChannels < 7) ? outAuData : nullptr);
      // print out collected file statistics
      fprintf_s (stdout, " Done, actual average %.1f kbit/s\n\n", (float) br * 0.001f);
      if (numChannels < 7)
      {
        fprintf_s (stdout, " Input statistics: Mobile loudness %.2f LUFS,\tsample peak level %.2f dBFS\n\n",
                   __max (3u, loudStats >> 16) / 512.f - 100.0f, 20.0f * log10 (__max (EA_PEAK_MIN, float (loudStats & USHRT_MAX))) + EA_PEAK_NORM);
      }

      if (!readStdin && (actualLength != expectLength || bw != headerRes))
      {
        fprintf_s (stderr, " WARNING: %lld sample frames read but %lld sample frames expected!\n", (long long) actualLength, (long long) expectLength);
        if (bw != headerRes) fprintf_s (stderr, "          The encoded MPEG-4 bit-stream is likely to be unreadable!\n");
        fprintf_s (stderr, "\n");
      }
#if USE_EXHALELIB_DLL
      exhaleDelete (&exhaleEnc);
#endif
    } // end coding loop and stats print-out
  }

mainFinish:

  // free all dynamic memory
  if (inPcmData != nullptr)
  {
    free ((void*) inPcmData);
    inPcmData = nullptr;
  }
  if (outAuData != nullptr)
  {
    free ((void*) outAuData);
    outAuData = nullptr;
  }
  // close input file
  if (inFileHandle != -1)
  {
    if (_CLOSE (inFileHandle) != 0)
    {
      if (readStdin) // stdin
      {
        fprintf_s (stderr, " ERROR while trying to close stdin stream! Has it already been closed?\n\n");
      }
      else  // argc = 4, file
      {
        fprintf_s (stderr, " ERROR while trying to close input file %s! Does it still exist?\n\n", argv[2]);
      }
    }
    inFileHandle = 0;
  }
  // close output file
  if (outFileHandle != -1)
  {
    if (_CLOSE (outFileHandle) != 0)
    {
      fprintf_s (stderr, " ERROR while trying to close output file %s! Does it still exist?\n\n", argv[argc - 1]);
    }
    outFileHandle = 0;
  }

  return (inFileHandle | outFileHandle | i);
}
