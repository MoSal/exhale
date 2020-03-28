/* exhaleLibPch.h - pre-compiled header file for classes of exhaleLib coding library
 * written by C. R. Helmrich, last modified in 2020 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under a Modified BSD-Style License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2020 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#ifndef _EXHALE_LIB_PCH_H_
#define _EXHALE_LIB_PCH_H_

#include <limits.h> // for .._MAX, .._MIN
#include <math.h>   // for pow, sin, sqrt
#include <stdint.h> // for (u)int8_t, (u)int16_t, (u)int32_t, (u)int64_t
#include <stdlib.h> // for abs, div, calloc, malloc, free, (__)max, (__)min, (s)rand
#include <string.h> // for memcpy, memset
#include <vector>   // for std::vector <>

// constants, experimental macros
#define AAC_NUM_SAMPLE_RATES   13
#define MAX_PREDICTION_ORDER    4
#define MAX_NUM_SWB_LFE         6
#define MAX_NUM_SWB_LONG       51
#define MAX_NUM_SWB_SHORT      15
#define MIN_NUM_SWB_SHORT      12
#define NUM_WINDOW_GROUPS       4 // must be between 4 and 8
#define USAC_MAX_NUM_CHANNELS   8
#define USAC_MAX_NUM_ELCONFIGS 13
#define USAC_MAX_NUM_ELEMENTS   5
#define USAC_NUM_FREQ_TABLES    6
#define USAC_NUM_SAMPLE_RATES  (2 * AAC_NUM_SAMPLE_RATES)

#define RESTRICT_TO_AAC         0 // allow only AAC tool-set

#if RESTRICT_TO_AAC
# define LFE_MAX               12
#else
# define LFE_MAX               24
#endif

#ifndef __max
# define __max(a, b)           ((a) > (b) ? (a) : (b))
#endif
#ifndef __min
# define __min(a, b)           ((a) < (b) ? (a) : (b))
#endif
#ifndef CLIP_PM
# define CLIP_PM(x, clip)      (__max (-clip, __min (+clip, x)))
#endif
#ifndef CLIP_UCHAR
# define CLIP_UCHAR(x)         (__max (0, __min (UCHAR_MAX, x)))
#endif

#ifndef MFREE
# define MFREE(x)              if (x != nullptr) { free ((void*) x); x = nullptr; }
#endif

// usacElementType[el] definition
typedef enum ELEM_TYPE : int8_t
{
  ID_EL_UNDEF = -1,
  ID_USAC_SCE = 0, // single-channel element (CCI_1_CH)
  ID_USAC_CPE = 1, // channel-pair element (CCI_2_CH)
  ID_USAC_LFE = 2, // low-frequency effects element (CCI_1_CH) with ONLY_LONG_SEQUENCE, no TNS
  ID_USAC_EXT = 3  // extension element (not to be used here!)
} ELEM_TYPE;

// window_sequence defines for ICS
typedef enum USAC_WSEQ : uint8_t
{
  ONLY_LONG   = 0, // one symmetric long window  .´`.
  LONG_START  = 1, // asymmet. long-start window .´|_
  EIGHT_SHORT = 2, // 8 symmetric short windows  _MM_
  LONG_STOP   = 3, // asymmet. long-stop window  _|`.
  STOP_START  = 4  // symmet. stop-start window  _||_
} USAC_WSEQ;

// window_shape definition for ICS
typedef enum USAC_WSHP : uint8_t
{
  WINDOW_SINE = 0, // half-sine window = sin(pi*((0:L-1)'+0.5)/L)
  WINDOW_KBD  = 1  // Kaiser-Bessel derived (KBD) window by Dolby
} USAC_WSHP;

// ics_info(): channel data struct
struct IcsInfo
{
  uint8_t   maxSfb;         // max_sfb(1)
  uint8_t   windowGrouping; // scale_factor_grouping (index)
  USAC_WSEQ windowSequence; // window_sequence
  USAC_WSHP windowShape;    // window_shape
};

// tns_data(): channel data struct
struct TnsData
{
  int8_t    coeff[3][MAX_PREDICTION_ORDER];
  int16_t   coeffParCor[MAX_PREDICTION_ORDER];
  bool      coeffResLow;  // means coef_res[w]=0
  uint8_t   filteredWindow; // filtered window w
  bool      filterDownward[3]; // direction[f]=1
  uint8_t   filterLength[3]; // filter length[f]
  uint8_t   filterOrder[3];   // filter order[f]
  uint8_t   numFilters;   // n_filt for window w
};

// scale factor group. data struct
struct SfbGroupData
{
  uint16_t  numWindowGroups; // 1 | NUM_WINDOW_GROUPS, num_window_groups
  uint16_t  sfbOffsets[1+MAX_NUM_SWB_SHORT * NUM_WINDOW_GROUPS];
  uint32_t  sfbRmsValues[MAX_NUM_SWB_SHORT * NUM_WINDOW_GROUPS];
  uint8_t   scaleFactors[MAX_NUM_SWB_SHORT * NUM_WINDOW_GROUPS]; // sf[]
  uint8_t   sfbsPerGroup; // max_sfb(1) duplicate needed by BitAllocator
  uint8_t   windowGroupLength[NUM_WINDOW_GROUPS]; // window_group_length
};

// UsacCoreCoderData() data struct
struct CoreCoderData
{
  bool          commonMaxSfb;  // common_max_sfb in StereoCoreToolInfo()
  bool          commonTnsData;     // common_tns in StereoCoreToolInfo()
  bool          commonWindow;   // common_window in StereoCoreToolInfo()
  ELEM_TYPE     elementType;   // usacElementType in UsacDecoderConfig()
  SfbGroupData  groupingData[2];  // window grouping and SFB offset data
  IcsInfo       icsInfoCurr[2];   // current ics_info() for each channel
  IcsInfo       icsInfoPrev[2];  // previous ics_info() for each channel
#if !RESTRICT_TO_AAC
  uint8_t       specFillData[2];  // noise filling data for each channel
#endif
  uint8_t       stereoConfig;  // cplx_pred_data() config: pred_dir etc.
  uint8_t       stereoData[MAX_NUM_SWB_SHORT * NUM_WINDOW_GROUPS];
  uint8_t       stereoMode;   // ms_mask_present in StereoCoreToolInfo()
  bool          tnsActive;    // tns_active flag in StereoCoreToolInfo()
  TnsData       tnsData[2];       // current tns_data() for each channel
  bool          tnsOnLeftRight;     // tns_on_lr in StereoCoreToolInfo()
};

// bit-stream encoding data struct
struct OutputStream
{
  uint8_t heldBitChunk; // bits not yet flushed to buffer
  uint8_t heldBitCount; // number of bits not yet flushed
  std::vector <uint8_t> stream; // FIFO bit-stream buffer
  // constructor
  OutputStream () { reset (); }
  // destructor
  ~OutputStream() { stream.clear (); }
  // public functions
  void reset ();   // clear writer states and byte buffer
  void write (const uint32_t bitChunk, const uint8_t bitCount);
}; // OutputStream

// fast calculation of sqrt (256 - x): (4 + eightTimesSqrt256Minus[x]) >> 3, for 0 <= x <= 255
const uint8_t eightTimesSqrt256Minus[256] = {
  128, 128, 127, 127, 127, 127, 126, 126, 126, 126, 125, 125, 125, 125, 124, 124, 124, 124, 123, 123, 123, 123, 122, 122, 122, 122,
  121, 121, 121, 121, 120, 120, 120, 119, 119, 119, 119, 118, 118, 118, 118, 117, 117, 117, 116, 116, 116, 116, 115, 115, 115, 115,
  114, 114, 114, 113, 113, 113, 113, 112, 112, 112, 111, 111, 111, 111, 110, 110, 110, 109, 109, 109, 109, 108, 108, 108, 107, 107,
  107, 106, 106, 106, 106, 105, 105, 105, 104, 104, 104, 103, 103, 103, 102, 102, 102, 102, 101, 101, 101, 100, 100, 100,  99,  99,
   99,  98,  98,  98,  97,  97,  97,  96,  96,  96,  95,  95,  95,  94,  94,  94,  93,  93,  93,  92,  92,  92,  91,  91,  91,  90,
   90,  89,  89,  89,  88,  88,  88,  87,  87,  87,  86,  86,  85,  85,  85,  84,  84,  84,  83,  83,  82,  82,  82,  81,  81,  80,
   80,  80,  79,  79,  78,  78,  78,  77,  77,  76,  76,  75,  75,  75,  74,  74,  73,  73,  72,  72,  72,  71,  71,  70,  70,  69,
   69,  68,  68,  67,  67,  66,  66,  65,  65,  64,  64,  63,  63,  62,  62,  61,  61,  60,  60,  59,  59,  58,  58,  57,  57,  56,
   55,  55,  54,  54,  53,  52,  52,  51,  51,  50,  49,  49,  48,  47,  47,  46,  45,  45,  44,  43,  42,  42,  41,  40,  39,  38,
   38,  37,  36,  35,  34,  33,  32,  31,  30,  29,  28,  27,  25,  24,  23,  21,  20,  18,  16,  14,  11,   8
};

// fast calculation of x / den: (x * oneTwentyEightOver[den]) >> 7, accurate for 0 <= x <= 162
const uint8_t oneTwentyEightOver[14] = {0, 128, 64, 43, 32, 26, 22, 19, 16, 15, 13, 12, 11, 10};

// public sampling rate functions
char toSamplingFrequencyIndex (const unsigned samplingRate);
unsigned toSamplingRate (const char samplingFrequencyIndex);

#endif // _EXHALE_LIB_PCH_H_
