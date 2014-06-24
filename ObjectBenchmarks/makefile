# makefile for doing some benchmarks

# Compiler
CC=icpc

# Compilation Flags
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings	
CFLAGS = -g -Wall -openmp -vec-report=5 -O3 -std=gnu++0x

BMTOOLS = ../BenchmarkTools

INCLUDES = $(BMTOOLS)



all: obj-bmarks obj-bmarks-mic

obj-bmarks: ObjectBenchmarks.cc ObjectBenchmarks.h Seeds.h Seeds.cc $(BMTOOLS)/BenchmarkTools.cc $(BMTOOLS)/BenchmarkTools.h
	$(CC) $(CFLAGS) ObjectBenchmarks.cc ObjectBenchmarks.h Seeds.h Seeds.cc $(BMTOOLS)/BenchmarkTools.cc $(BMTOOLS)/BenchmarkTools.h -o obj-bmarks -I$(INCLUDES)

obj-bmarks-mic: ObjectBenchmarks.cc ObjectBenchmarks.h Seeds.h Seeds.cc $(BMTOOLS)/BenchmarkTools.cc $(BMTOOLS)/BenchmarkTools.h
	$(CC) $(CFLAGS) -mmic ObjectBenchmarks.cc ObjectBenchmarks.h Seeds.h Seeds.cc $(BMTOOLS)/BenchmarkTools.cc $(BMTOOLS)/BenchmarkTools.h -o obj-bmarks-mic -I$(INCLUDES)
	scp obj-bmarks-mic root@mic0:

run: obj-bmarks
	./obj-bmarks

run-mic: obj-bmarks-mic
	ssh root@mic0 ./obj-bmarks-mic

run-all: obj-bmarks obj-bmarks-mic
	./obj-bmarks
	ssh root@mic0 ./obj-bmarks-mic


clean:
	rm -f obj-bmarks

clean-mic:
	rm -f obj-bmakrs-mic

clean-all:
	rm -f obj-bmarks obj-bmarks-mic	