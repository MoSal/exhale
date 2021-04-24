## makefile - master user make-file for compiling exhale on Linux and MacOS platforms
 # written by C. R. Helmrich, last modified in 2021 - see License.htm for legal notices
 # Universal 2 support for macOS added by Christopher Snowhill on hydrogenaud.io in 2021
 #
 # The copyright in this software is being made available under the exhale Copyright License
 # and comes with ABSOLUTELY NO WARRANTY. This software may be subject to other third-
 # party rights, including patent rights. No such rights are granted under this License.
 #
 # Copyright (c) 2018-2021 Christian R. Helmrich, project ecodis. All rights reserved.
 ##

## BUILD32=1: compile for 32-bit platforms, BUILD32=0: compile for 64-bit platforms
BUILD32?=0

## UNIVERSAL2=1: compile for both x86_64 and arm64 on macOS, UNIVERSAL2=0: compile for native architecture on all platforms
UNIVERSAL2?=0

export BUILD32
export UNIVERSAL2

all:
	$(MAKE) -C src/lib  MM32=$(BUILD32) UNIVERSAL2=$(UNIVERSAL2)
	$(MAKE) -C src/app  MM32=$(BUILD32) UNIVERSAL2=$(UNIVERSAL2)

debug:
	$(MAKE) -C src/lib  debug MM32=$(BUILD32) UNIVERSAL2=$(UNIVERSAL2)
	$(MAKE) -C src/app  debug MM32=$(BUILD32) UNIVERSAL2=$(UNIVERSAL2)

release:
	$(MAKE) -C src/lib  release MM32=$(BUILD32) UNIVERSAL2=$(UNIVERSAL2)
	$(MAKE) -C src/app  release MM32=$(BUILD32) UNIVERSAL2=$(UNIVERSAL2)

clean:
	$(MAKE) -C src/lib  clean MM32=$(BUILD32) UNIVERSAL2=$(UNIVERSAL2)
	$(MAKE) -C src/app  clean MM32=$(BUILD32) UNIVERSAL2=$(UNIVERSAL2)
