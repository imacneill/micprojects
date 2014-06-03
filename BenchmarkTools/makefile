# makefile for testing the benchmark tools

# Compiler
CC=icc

# Compilation Flags
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings	
CFLAGS = -g -Wall


test: test.cc BenchmarkTools.cc BenchmarkTools.h
	$(CC) $(CFLAGS1) test.cc BenchmarkTools.cc -o test


clean:
	rm -f test	