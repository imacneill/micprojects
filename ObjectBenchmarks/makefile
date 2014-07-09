# makefile for doing some benchmarks

# Compiler
# use icpc because icc has an annoying habit of not knowing when a header file is c++ or not
CC=icpc

VECHOST = -mavx
VECMIC = -mmic
# Compilation Flags
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings	
#  -std=c++11 -std=gnu++11 or -std=gnu++0x should add c++ 11 support
CFLAGS = -g -Wall -openmp -vec-report=3 -O3 -std=gnu++0x

BMTOOLS = ../BenchmarkTools/BenchmarkTools.cc ../BenchmarkTools/BenchmarkTools.h

INCLUDES = ../BenchmarkTools

ROOT = `root-config --libs --cflags`


all: obj-bmarks obj-bmarks-mic

% : %.C
	g++ -g -Wall -std=gnu++0x $(ROOT) -o $@ $<

%.o : %.cc *.h
	$(CC) $(CFLAGS) $(VECHOST) -I$(INCLUDES) -c -o $@ $<

%.om : %.cc *.h
	$(CC) $(CFLAGS) $(VECMIC) -I$(INCLUDES) -c -o $@ $<

obj-bmarks: ObjectBenchmarks.o Seeds.o $(BMTOOLS)
	$(CC) $(CFLAGS) $(VECHOST) -o $@ $^

obj-bmarks-mic: ObjectBenchmarks.om Seeds.om $(BMTOOLS)
	$(CC) $(CFLAGS) $(VECMIC) -o $@ $^
	scp $@ mic0:

run: obj-bmarks
	./obj-bmarks

run-mic: obj-bmarks-mic
	ssh mic0 mkdir -p output
	ssh mic0 ./obj-bmarks-mic
	scp mic0:~/output/* output_mic/

run-all: obj-bmarks obj-bmarks-mic
	./obj-bmarks
	ssh mic0 mkdir -p output
	ssh mic0 ./obj-bmarks-mic
	scp mic0:~/output/* output_mic/


clean:
	rm -f obj-bmarks

clean-mic:
	rm -f obj-bmakrs-mic

clean-all:
	rm -f obj-bmarks obj-bmarks-mic	*.o *.om