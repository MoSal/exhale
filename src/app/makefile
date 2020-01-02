## makefile - application make-file for compiling exhale on Linux and MacOS platforms
 # written by C. R. Helmrich, last modified 2019 - see License.txt for legal notices
 #
 # The copyright in this software is being made available under a Modified BSD License
 # and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 # party rights, including patent rights. No such rights are granted under this License.
 #
 # Copyright (c) 2018-2019 Christian R. Helmrich, project ecodis. All rights reserved.
 ##

# define as console application
CONFIG    = CONSOLE

# source and output directories
DIR_BIN   = ../../bin
DIR_OBJ   = ../../build
DIR_INC   = ../../include
DIR_LIB   = ../../lib
DIR_SRC   = ../../src/app

# build with large file support
DEFS      = -DMSYS_LINUX -DMSYS_UNIX_LARGEFILE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64

# name of product / binary file
PRD_NAME  = Exhale

# name of temporary object file
OBJS      = \
	$(DIR_OBJ)/basicMP4Writer.o \
	$(DIR_OBJ)/basicWavReader.o \
	$(DIR_OBJ)/exhaleApp.o \
	$(DIR_OBJ)/exhaleAppPch.o \

# define libraries to link with
LIBS      = -ldl
DYN_LIBS  =
STAT_LIBS = -lpthread

DYN_DEBUG_LIBS        = -lExhaleDynd
DYN_DEBUG_PREREQS     = $(DIR_LIB)/libExhaleDynd.a
STAT_DEBUG_LIBS       = -lExhaled
STAT_DEBUG_PREREQS    = $(DIR_LIB)/libExhaled.a

DYN_RELEASE_LIBS      = -lExhale
DYN_RELEASE_PREREQS   = $(DIR_LIB)/libExhale.a
STAT_RELEASE_LIBS     = -lExhale
STAT_RELEASE_PREREQS  = $(DIR_LIB)/libExhale.a

# include common makefile.base
include ../makefile.base
