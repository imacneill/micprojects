# makefile for SeedAlgorithm Code

# Compiler
CC=icc

# Compilation Flags
CFLAGS1=-openmp -vec-report=3 -O3 -I ../tools
CFLAGS2=-openmp -mmic -vec-report=3 -O3 -I ../tools

all:
	$(CC) $(CFLAGS1) ArraySeedAlgorithm.C  ../tools/ReadSeedFile.C -o ../output/arrayTest
	$(CC) $(CFLAGS1) StructSeedAlgorithm.C ../tools/ReadSeedFile.C -o ../output/structTest
	$(CC) $(CFLAGS1) AOSOA.C               ../tools/ReadSeedFile.C -o ../output/aosoaTest	
	$(CC) $(CFLAGS2) ArraySeedAlgorithm.C  ../tools/ReadSeedFile.C -o ../output/arrayTest_mmic
	$(CC) $(CFLAGS2) StructSeedAlgorithm.C ../tools/ReadSeedFile.C -o ../output/structTest_mmic
	$(CC) $(CFLAGS2) AOSOA.C               ../tools/ReadSeedFile.C -o ../output/aosoaTest_mmic

array:
	$(CC) $(CFLAGS1) ArraySeedAlgorithm.C  ../tools/ReadSeedFile.C -o ../output/arrayTest
	$(CC) $(CFLAGS2) ArraySeedAlgorithm.C  ../tools/ReadSeedFile.C -o ../output/arrayTest_mmic

struct:
	$(CC) $(CFLAGS1) StructSeedAlgorithm.C ../tools/ReadSeedFile.C -o ../output/structTest	
	$(CC) $(CFLAGS2) StructSeedAlgorithm.C ../tools/ReadSeedFile.C -o ../output/structTest_mmic	

aosoa:
	$(CC) $(CFLAGS1) AOSOA.C               ../tools/ReadSeedFile.C -o ../output/aosoaTest	
	$(CC) $(CFLAGS2) AOSOA.C               ../tools/ReadSeedFile.C -o ../output/aosoaTest_mmic	


clean:
	rm -f ../output/arrayTest
	rm -f ../output/structTest
	rm -f ../output/aosoaTest
	rm -f ../output/arrayTest_mmic
	rm -f ../output/structTest_mmic
	rm -f ../output/aosoaTest_mmic