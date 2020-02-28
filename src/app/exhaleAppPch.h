/* exhaleAppPch.h - pre-compiled header file for source code of exhale application
 * written by C. R. Helmrich, last modified in 2019 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under a Modified BSD-Style License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2020 Christian R. Helmrich, project ecodis. All rights reserved.
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
#ifndef fprintf_s
# define fprintf_s             fprintf
#endif

// public sampling rate function
bool isSamplingRateSupported (const unsigned samplingRate);

#endif // _EXHALE_APP_PCH_H_
