#  Seed Algorithm Code

## Motivation
There are three structures of code within SeedAlgorithm.  

-  Array-based code
-  Struct-based code
-  Array of Structures of Arrays (AOSOA) - based code

This is designed to test the parallel processing speeds under different data architectures.  


## Running the Code 

1.  In the src/ directory,  
    *$  make*
2.  Move to the output directory,   
    *$  cd ../output/*  
3.  Transfer executables to the mic0   
    *$  scp 'filename' mic0:~/bin/*  
4.  Move to mic0   
    *$  ssh mic0*   
    *%  cd bin/*
5.  Set environment variables and execute file**   
    *%   export OMP_NUM_THREADS=122*   
    *%   ./'filename'*  


 **This relies on the data file, seedDataFile.txt, being in data folder.  Executable looks for data in, '../data/seedDataFile.txt'