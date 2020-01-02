/* exhaleDecl.h - header file with declarations for exhale DLL export under Windows
 * written by C. R. Helmrich, last modified in 2019 - see License.htm for legal notices
 *
 * The copyright in this software is being made available under a Modified BSD-Style License
 * and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 * party rights, including patent rights. No such rights are granted under this License.
 *
 * Copyright (c) 2018-2020 Christian R. Helmrich, project ecodis. All rights reserved.
 */

#ifndef _EXHALE_DECL_H_
#define _EXHALE_DECL_H_

#if defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64)
# ifdef EXHALE_DYN_LINK
#  define EXHALE_DECL __declspec (dllexport)
# else
#  define EXHALE_DECL
# endif
#endif

#endif // _EXHALE_DECL_H_
