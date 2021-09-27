## makefile - code library make-file for compiling exhale on Linux and MacOS platforms
 # written by C. R. Helmrich, last modified in 2020 - see License.htm for legal notices
 #
 # The copyright in this software is being made available under the exhale Copyright License
 # and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 # party rights, including patent rights. No such rights are granted under this License.
 #
 # Copyright (c) 2018-2021 Christian R. Helmrich, project ecodis. All rights reserved.
 ##

# define as source code library
CONFIG    = LIBRARY

# source and output directories
DIR_BIN   = ../../bin
DIR_OBJ   = ../../build
DIR_INC   = ../../include
DIR_LIB   = ../../lib
DIR_SRC   = ../../src/lib

# build with large file support
DEFS      = -DMSYS_LINUX -DMSYS_UNIX_LARGEFILE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64

# name of product / binary file
PRD_NAME  = exhale

# name of temporary object file
OBJS      = \
	$(DIR_OBJ)/bitAllocation.o \
	$(DIR_OBJ)/bitStreamWriter.o \
	$(DIR_OBJ)/entropyCoding.o \
	$(DIR_OBJ)/exhaleEnc.o \
	$(DIR_OBJ)/exhaleLibPch.o \
	$(DIR_OBJ)/lappedTransform.o \
	$(DIR_OBJ)/linearPrediction.o \
	$(DIR_OBJ)/quantization.o \
	$(DIR_OBJ)/specAnalysis.o \
	$(DIR_OBJ)/specGapFilling.o \
	$(DIR_OBJ)/stereoProcessing.o \
	$(DIR_OBJ)/tempAnalysis.o \

# define libraries to link with
LIBS      = -lpthread
DYN_LIBS  = -ldl
STAT_LIBS =

DYN_DEBUG_LIBS    =
DYN_RELEASE_LIBS  =
STAT_DEBUG_LIBS   =
STAT_RELEASE_LIBS =

# include common makefile.base
include ../makefile.base
