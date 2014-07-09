# makefile for doing some benchmarks

# Compiler
CC=icpc

# Compilation Flags
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings	
CFLAGS = -g -Wall -openmp -vec-report=3 -O3

BMTOOLS = ../BenchmarkTools

INCLUDES = $(BMTOOLS)

all: flops flops-mic

flops: flopstest.cc $(BMTOOLS)/BenchmarkTools.cc $(BMTOOLS)/BenchmarkTools.h
	$(CC) $(CFLAGS) flopstest.cc $(BMTOOLS)/BenchmarkTools.cc -o flops -I$(INCLUDES) -D__ISHOST__
	#$(CC) $(CFLAGS) flopstest.cc $(BMTOOLS)/BenchmarkTools.cc -I$(INCLUDES) -D__ISHOST__

flops-mic: flopstest.cc $(BMTOOLS)/BenchmarkTools.cc $(BMTOOLS)/BenchmarkTools.h
	$(CC) $(CFLAGS) -mmic flopstest.cc $(BMTOOLS)/BenchmarkTools.cc -o flops-mic -I$(INCLUDES) -D__ISMIC__
	#$(CC) $(CFLAGS) -mmic flopstest.cc $(BMTOOLS)/BenchmarkTools.cc -I$(INCLUDES) -D__ISMIC__
	scp flops-mic root@mic0:

run: flops flops-mic
	./flops
	ssh root@mic0 ./flops-mic


clean:
	rm -f flops

clean-mic:
	rm -f flops-mic

clean-all:
	rm -f flops flops-mic	