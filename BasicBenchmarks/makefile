# makefile for doing some benchmarks

# Compiler
CC=icc

# Compilation Flags
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings	
CFLAGS = -g -Wall -openmp -vec-report=3 -O3

BMTOOLS = ../BenchmarkTools

INCLUDES = $(BMTOOLS)


all: bmarks bmarks-mic

bmarks: BasicBenchmarks.cc BasicBenchmarks.h $(BMTOOLS)/BenchmarkTools.cc $(BMTOOLS)/BenchmarkTools.h
	$(CC) $(CFLAGS) BasicBenchmarks.cc $(BMTOOLS)/BenchmarkTools.cc -o bmarks -I$(INCLUDES)

bmarks-mic: BasicBenchmarks.cc BasicBenchmarks.h $(BMTOOLS)/BenchmarkTools.cc $(BMTOOLS)/BenchmarkTools.h
	$(CC) $(CFLAGS) -mmic BasicBenchmarks.cc $(BMTOOLS)/BenchmarkTools.cc -o bmarks-mic -I$(INCLUDES)
	scp bmarks-mic root@mic0:

run: bmarks bmarks-mic
	./bmarks
	ssh root@mic0 ./bmarks-mic


clean:
	rm -f bmarks

clean-mic:
	rm -f bmakrs-mic

clean-all:
	rm -f bmarks bmarks-mic	