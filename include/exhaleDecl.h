/* exhaleDecl.h - header file with declarations for exhale DLL ex-/import under Windows
 * written by C. R. Helmrich, last modified in 2020 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under a Modified BSD-Style License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2020 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#ifndef _EXHALE_DECL_H_
#define _EXHALE_DECL_H_

#include "../src/lib/exhaleEnc.h"

// DLL constructor
extern "C" EXHALE_DECL ExhaleEncoder* exhaleCreate (int32_t* const, unsigned char* const, const unsigned, const unsigned,
                                                    const unsigned, const unsigned, const unsigned, const bool, const bool);

// DLL destructor
extern "C" EXHALE_DECL unsigned exhaleDelete (ExhaleEncoder*);

// DLL initializer
extern "C" EXHALE_DECL unsigned exhaleInitEncoder (ExhaleEncoder*, unsigned char* const, uint32_t* const);

// DLL lookahead encoder
extern "C" EXHALE_DECL unsigned exhaleEncodeLookahead (ExhaleEncoder*);

// DLL frame encoder
extern "C" EXHALE_DECL unsigned exhaleEncodeFrame (ExhaleEncoder*);

#endif // _EXHALE_DECL_H_
