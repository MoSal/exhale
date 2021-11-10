/* exhaleAppPch.h - pre-compiled header file for source code of exhale application
 * written by C. R. Helmrich, last modified in 2021 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under the exhale Copyright License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2021 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#ifndef _EXHALE_APP_PCH_H_
#define _EXHALE_APP_PCH_H_

#include <limits.h> // for .._MAX, .._MIN
#include <math.h>   // for log, pow, sqrt
#include <stdint.h> // for (u)int8_t, (u)int16_t, (u)int32_t, (u)int64_t
#include <stdlib.h> // for abs, div, calloc, malloc, free, (__)max, (__)min, (s)rand
#include <string.h> // for memcpy, memset
#include <vector>   // for std::vector <>
#if defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64)
# include <io.h>

# define _CLOSE _close
# define _READ  _read
# define _SEEK  _lseeki64
# define _WRITE _write
#else // Linux, MacOS, Unix
# include <unistd.h>

# define _CLOSE ::close
# define _READ  ::read
# define _SEEK  ::lseek
# define _WRITE ::write
#endif

#ifndef __max
# define __max(a, b)           ((a) > (b) ? (a) : (b))
#endif
#ifndef __min
# define __min(a, b)           ((a) < (b) ? (a) : (b))
#endif
#if !defined (fprintf_s) && !defined (__MINGW32__)
# define fprintf_s             fprintf
#endif
#if !defined (fwprintf_s) && !defined (__MINGW32__)
# define fwprintf_s            fwprintf
#endif
#ifndef MFREE
# define MFREE(x)              if (x != nullptr) { free ((void*) x); x = nullptr; }
#endif

// public extrapolation function
void eaExtrapolate (int32_t* const pcmBuffer, const uint16_t pcmOffset,
                    const uint16_t frameSize, const uint16_t numChannels, const bool fadeIn = false);

// public sampling rate function
bool isSamplingRateSupported (const unsigned samplingRate);

#endif // _EXHALE_APP_PCH_H_
