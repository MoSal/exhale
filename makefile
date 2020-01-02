## makefile - master user make-file for compiling exhale on Linux and MacOS platforms
 # written by C. R. Helmrich, last modified 2019 - see License.txt for legal notices
 #
 # The copyright in this software is being made available under a Modified BSD License
 # and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 # party rights, including patent rights. No such rights are granted under this License.
 #
 # Copyright (c) 2018-2019 Christian R. Helmrich, project ecodis. All rights reserved.
 ##

## BUILD32=1: compile for 32-bit platforms, BUILD32=0: compile for 64-bit platforms
BUILD32?=0

export BUILD32

all:
	$(MAKE) -C src/lib  MM32=$(BUILD32)
	$(MAKE) -C src/app  MM32=$(BUILD32)

debug:
	$(MAKE) -C src/lib  debug MM32=$(BUILD32)
	$(MAKE) -C src/app  debug MM32=$(BUILD32)

release:
	$(MAKE) -C src/lib  release MM32=$(BUILD32)
	$(MAKE) -C src/app  release MM32=$(BUILD32)

clean:
	$(MAKE) -C src/lib  clean MM32=$(BUILD32)
	$(MAKE) -C src/app  clean MM32=$(BUILD32)
