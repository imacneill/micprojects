# makefile for doing some benchmarks

# Compiler
CC=icc

# Compilation Flags
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings	
CFLAGS = -g -Wall -openmp -vec-report=3 -O3



all: main_test main_test-mic

main_test: main.cc main.h
	$(CC) $(CFLAGS) main.cc -o main_test

main_test-mic: main.cc main.h
	$(CC) $(CFLAGS) -mmic main.cc -o main_test-mic
	scp main_test-mic root@mic0:	

run: main_test main_test-mic
	./main_test
	ssh root@mic0 ./main_test-mic


clean:
	rm -f main_test

clean-mic:
	rm -f main_test-mic

clean-all:
	rm -f main_test main_test-mic	