/* exhaleDecl.h - header file with declarations for exhale DLL ex-/import under Windows
 * written by C. R. Helmrich, last modified in 2020 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under the exhale Copyright License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2021 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#ifndef _EXHALE_DECL_H_
#define _EXHALE_DECL_H_

#include <stdint.h> /* for (u)int8_t, (u)int16_t, (u)int32_t, (u)int64_t */
#if !defined(_MSC_VER) || _MSC_VER >= 1800
#include <stdbool.h>
#endif

#if defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64)
# ifdef EXHALE_DYN_LINK
#  define EXHALE_DECL __declspec (dllexport)
# else
#  define EXHALE_DECL __declspec (dllimport)
# endif
#else
# define EXHALE_DECL
#endif

#ifdef __cplusplus
struct ExhaleEncAPI
{
  /* initializer */
  virtual unsigned initEncoder (unsigned char* const audioConfigBuffer, uint32_t* const audioConfigBytes = nullptr) = 0;
  /* lookahead encoder */
  virtual unsigned encodeLookahead () = 0;
  /* frame encoder */
  virtual unsigned encodeFrame () = 0;
  /* destructor */
  virtual ~ExhaleEncAPI () { }
};

extern "C"
{
#else /* C, not C++ */
struct ExhaleEncAPI; /* opaque type */
typedef struct ExhaleEncAPI ExhaleEncAPI;
#endif

/* C constructor */
EXHALE_DECL ExhaleEncAPI* exhaleCreate (int32_t* const, unsigned char* const, const unsigned, const unsigned,
                                        const unsigned, const unsigned, const unsigned, const bool, const bool);
/* C destructor */
EXHALE_DECL unsigned exhaleDelete (ExhaleEncAPI*);

/* C initializer */
EXHALE_DECL unsigned exhaleInitEncoder (ExhaleEncAPI*, unsigned char* const, uint32_t* const);

/* C lookahead encoder */
EXHALE_DECL unsigned exhaleEncodeLookahead (ExhaleEncAPI*);

/* C frame encoder */
EXHALE_DECL unsigned exhaleEncodeFrame (ExhaleEncAPI*);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* _EXHALE_DECL_H_ */
