# makefile for doing some benchmarks

# Compiler
CC=icc

# Compilation Flags
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings	
CFLAGS = -g -Wall -openmp -vec-report=3 -O3

BMTOOLS = ../BenchmarkTools

INCLUDES = $(BMTOOLS)


main: BasicBenchmarks.cc BasicBenchmarks.h $(BMTOOLS)/BenchmarkTools.cc $(BMTOOLS)/BenchmarkTools.h
	$(CC) $(CFLAGS) BasicBenchmarks.cc $(BMTOOLS)/BenchmarkTools.cc -o main -I$(INCLUDES)


clean:
	rm -f main	