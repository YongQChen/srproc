# srproc
srproc is utility program that PROCesses next-gen multiplexed 
Sequencing Reads that have been aligned to reference sequences
and reports read counts per reference sequences

Sequencing read data are expected to have the following fields:
1-8:	SequencerName	flowCell	lane	tile	x	y	tag	read
9-15:	sequence	quality	ref	pos	direction	map	mapQuality

See the files in testData for examples

# build
The code currently build under Microsoft Visual Studio 2019 with 
additional dependency on Dirent which can be downloaded from
[GitHub](https://github.com/tronkko/dirent/releases)
Add dirent.h file in your include directory.

The code may need to be modified to work on Linux/UNIX environment.
A makefile is provided, but has not been tested.

# Copying
srproc may be freely distributed under the MIT license.  See the LICENSE
file for details.
