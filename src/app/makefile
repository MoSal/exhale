## makefile - application make-file for compiling exhale on Linux and MacOS platforms
 # written by C. R. Helmrich, last modified in 2020 - see License.htm for legal notices
 #
 # The copyright in this software is being made available under the exhale Copyright License
 # and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 # party rights, including patent rights. No such rights are granted under this License.
 #
 # Copyright (c) 2018-2021 Christian R. Helmrich, project ecodis. All rights reserved.
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
PRD_NAME  = exhale

# name of temporary object file
OBJS      = \
	$(DIR_OBJ)/basicMP4Writer.o \
	$(DIR_OBJ)/basicWavReader.o \
	$(DIR_OBJ)/exhaleApp.o \
	$(DIR_OBJ)/exhaleAppPch.o \
	$(DIR_OBJ)/loudnessEstim.o \

# define libraries to link with
LIBS      = -ldl
DYN_LIBS  =
STAT_LIBS = -lpthread

DYN_DEBUG_LIBS        = -l$(PRD_NAME)Dynd
DYN_DEBUG_PREREQS     = $(DIR_LIB)/lib$(PRD_NAME)Dynd.a
STAT_DEBUG_LIBS       = -l$(PRD_NAME)d
STAT_DEBUG_PREREQS    = $(DIR_LIB)/lib$(PRD_NAME)d.a

DYN_RELEASE_LIBS      = -l$(PRD_NAME)
DYN_RELEASE_PREREQS   = $(DIR_LIB)/lib$(PRD_NAME).a
STAT_RELEASE_LIBS     = -l$(PRD_NAME)
STAT_RELEASE_PREREQS  = $(DIR_LIB)/lib$(PRD_NAME).a

# include common makefile.base
include ../makefile.base
