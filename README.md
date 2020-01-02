exhale
======

exhale, which is an acronym for "Ecodis eXtended High-efficiency And
Low-complexity Encoder", is a lightweight library and application to
encode uncompressed WAVE-format audio files into MPEG-4 format files
complying with the ISO/IEC 23003-3 (MPEG-D) Unified Speech and Audio
Coding (USAC, also known as Extended High-Efficiency AAC) standard.

exhale currently makes use of all frequency-domain (FD) coding tools
in the scalefactor based MDCT processing path, except for predictive
joint stereo, which is already being worked on. It objective is high
quality mono, stereo, and multichannel coding at medium and high bit
rates, so the parametric USAC coding tools (ACELP, TCX, Enhanced SBR
and MPEG Surround with Unified Stereo coding) won't be integrated. 

____________________________________________________________________


Copyright
---------

(c) 2020 Christian R. Helmrich, project ecodis. All rights reserved.


License
-------

exhale is being made available under an open-source license which is
similar to the 3-clause BSD license but modified to address specific
lar aspects dictated by the nature and output of this application.

The license text and release notes for the current version 1.0.0 can
be found in the "include" subdirectory of the exhale distribution.


Compilation
-----------

This section describes how to compile the exhale source code into an
executable application under Linux and Microsoft Windows. The binary
application files will show up in a newly created `bin` subdirectory
of the exhale distribution directory and/or a subdirectory thereof.

### Linux and MacOS (GNU Compiler Collection, gcc):

In a terminal, change to the exhale distribution directory and enter

`
make release
`

to build a release-mode executable with the default (usually 64-bit)
configuration. A 32-bit debug-mode executable can be built by typing

`
make debug BUILD32=1
`

### Microsoft Windows (Visual Studio 2012 and later):

Doubleclick the exhale_vs2012.sln file to open the project in Visual
Studio. Once it's loaded, rightclick on `exhaleApp` in the "Solution
Explorer" window on the right-hand side, then select `Set as StartUp
Project`. Now simply press `F7` to build the solution in debug mode.

For fastest encoding speed, please select `Release` and `x64` before
building the solution. This will create a release-mode 64-bit binary.


Usage
-----

This section describes how to run the exhale application either from
the command-line or using a third-party software providing WAVE data
to exhale's standard input pipe (stdin), such as foobar2000.

### Standalone (command-line):

In a terminal, change to exhale's `bin` subdirectory and then enter

`./exhale` (on Linux and MacOS) or `exhaleApp.exe` (on Windows)

to print out usage information. As an example, the following command

`exhaleApp.exe 5 C:\Music\Input.wav C:\Music\Output.m4a`

converts file Input.wav to file Output.m4a at roughly 128 kbit/s (if
the input file is two-channel stereo) and in xHE-AAC audio format.

### Third-party stdin (foobar2000):

will follow soon
