// Standard C
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <sys/time.h>
#include <stdlib.h>

// Standard C++
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cmath>

// User
#include "Seed.h"
#include "ReadSeedFile.h"


/////////////////////////////////////
//             Methods             //
/////////////////////////////////////

// Returns the current wall clock time
double dtime()
{
    double tseconds     = 0.0;
    struct timeval mytime;
    gettimeofday(&mytime,(struct timezone*)0);
    tseconds = (double)(mytime.tv_sec + mytime.tv_usec*1.0e-6);
    return( tseconds );
}

// Extend the length of the vector by a multiplicative factor
std::vector<seed_t> ExtendVector(std::vector<seed_t> & init_vec, int factor)
{
    std::vector<seed_t> ext_vec;

    for(int i = 0; i < factor; i++)
    {
        ext_vec.insert(ext_vec.end(), init_vec.begin(), init_vec.end());
    }
    
    return ext_vec;
}    

// Define Struct of Small Arrays
struct smallArrays_t
{
    float hit0_radius[16];
    float hit1_radius[16];
    float hit2_radius[16];
    float hit0_z[16];
    float hit1_z[16];
    float hit2_z[16];
};

/* To Do
1.  Create one SOA every 16 seeds
2.  Fill properly
3.  Figure out dependencies warnings
*/


// Returns an AOSOA of small array-size N
smallArrays_t * MakeAOSOA(const std::vector<seed_t> & seed_vector, const int N)
{
    // Declare the array pointer for the AOSOA
    smallArrays_t * AOSOA;

    // Determine the size of the AOSOA
    int seed_entries   = seed_vector.size();
    int num_structs    = seed_entries/N;
    int residual_seeds = seed_entries % N;

    std::cout<<"seed_entries: " << seed_entries << std::endl;
    std::cout<<"num_structs: " << num_structs << std::endl;
    std::cout<<"residual_seeds: " << residual_seeds << std::endl;

    // If the seed_entries are not evenly divisible, increase size of AOSOA to hold remainder
    if(residual_seeds  != 0)
        num_structs++;

    // Allocate memory for AOSOA
    AOSOA = (smallArrays_t *) _mm_malloc(sizeof(smallArrays_t) * (num_structs), 64);
    int k = -1;  // Index for the AOSOA

    // Iterate over the seed_vector, add smallArrays_t structs into AOSOA
    //for(int i = 0; i < seed_entries; i++)
    for(int i = 0; i < 100; i++)
    {
        // Index for the small arrays
        int j = i % N;        
        
        std::cout<<"i, j, k: " << i <<", " << j << ", " << k << std::endl;

        // Create an instance of the Structure of Arrays (SOA)
        std::cout<<"New SOA created"<<std::endl;
        smallArrays_t SOA;  

        // Every 16 elements, do the following procedures
        if(j == 0) {k++;}    
        
        // Fill the SOA from the seed_vector
        SOA.hit0_radius[j] = seed_vector[i].hit0_radius;
        SOA.hit1_radius[j] = seed_vector[i].hit1_radius;
        SOA.hit2_radius[j] = seed_vector[i].hit2_radius;
        SOA.hit0_z[j]      = seed_vector[i].hit0_z;
        SOA.hit1_z[j]      = seed_vector[i].hit1_z;
        SOA.hit2_z[j]      = seed_vector[i].hit2_z;

        // Add the SOA to the AOSOA after all 16 array elements have been filled        
        if(j == 15)
        {
            std::cout<<"Add SOA to AOSOA at position " << k <<std::endl;
            AOSOA[k] = SOA;
        }    
            
    }    

    //for(int i = 0; i < num_structs)

}


// Times the Residuals Calculation Using an Array Of Structs Of Arrays (AOSOA)
/*  Accessing memory consecutively is the fastest way to access memory on Intel Phi
    This improves cache efficiency, reduces the number of Translational Lookaside Buffer (TLB)
    misses, etc.
    AOSOA allows efficienct vectorization, while not overloading the TLB, nor spreading accesses
    across many pages. 
*/
void ComputeAOSOAResidual(smallArrays_t *, bool square = false)
{
    return;
}

/*
void ComputeAOSOAResidual(const std::vector<seed_t> & seed_vector, bool square = false, const int N)
{
    if(seed_vector.empty() == true)
    {
        std::cout<< "Seed vector empty" << std::endl;
        return;    
    }        
    
    double tstart = 0, tstop = 0, ttime = 0;
    int i = 0, num_threads = 0, chunk = 0;
    double gflops_res   = 0.0;
    double gflops_sqr   = 0.0;

    int seed_entries = seed_vector.size();

    // Determine the divisibility of the seed vector 
    int num_structs    = seed_entries/N;
    int residual_seeds = seed_entries % N;

    if(residual_seeds != 0) 
        std::cout<<"Non-even division of seeds into structs"<<std::endl;

    // Define array pointers
    float * hit0_z_array;
    float * hit1_z_array;
    float * hit2_z_array;
    float * hit0_radius_array;
    float * hit1_radius_array;
    float * hit2_radius_array;
    float * residual_array;

    // Define the array to hold the structs of small arrays
    smallArrays_t * small_array_struct_array;

    // Allocate memory 
    hit0_z_array      = (float*) _mm_malloc(sizeof(float)*seed_entries, 64); 
    hit1_z_array      = (float*) _mm_malloc(sizeof(float)*seed_entries, 64); 
    hit2_z_array      = (float*) _mm_malloc(sizeof(float)*seed_entries, 64); 
    hit0_radius_array = (float*) _mm_malloc(sizeof(float)*seed_entries, 64); 
    hit1_radius_array = (float*) _mm_malloc(sizeof(float)*seed_entries, 64); 
    hit2_radius_array = (float*) _mm_malloc(sizeof(float)*seed_entries, 64);  
    residual_array    = (float*) _mm_malloc(sizeof(float)*seed_entries, 64); 

    // Allocate memory for num_structs + 1, and have one partially filled struct 
    small_array_struct_array = (smallArrays_t*) _mm_malloc(sizeof(smallArrays_t)*(num_structs + 1), 64);

    #pragma omp parallel
    #pragma omp master
    num_threads = omp_get_num_threads();
    chunk = seed_entries/num_threads;

    std::cout<<"Number of threads = " << num_threads << std::endl;
    
    #pragma omp parallel for simd private(i)
    for(int i = 0; i < seed_entries; i++)
    {    
        hit0_radius_array[i] = seed_vector[i].hit0_radius; 
        hit1_radius_array[i] = seed_vector[i].hit1_radius; 
        hit2_radius_array[i] = seed_vector[i].hit2_radius; 
        hit0_z_array[i]      = seed_vector[i].hit0_z; 
        hit1_z_array[i]      = seed_vector[i].hit1_z; 
        hit2_z_array[i]      = seed_vector[i].hit2_z; 

        
        if(i % N == 0)
        {
            // If the index is divisible by N, begin filling a new a new struct of arrays
            smallArrays_t array_struct;
            array_struct.hit0_radius[i] = seed_vector[i].hit0_radius;
            array_struct.hit1_radius[i] = seed_vector[i].hit1_radius;
            array_struct.hit2_radius[i] = seed_vector[i].hit2_radius;
            array_struct.hit0_z[i]      = seed_vector[i].hit0_z;
            array_struct.hit1_z[i]      = seed_vector[i].hit1_z;
            array_struct.hit2_z[i]      = seed_vector[i].hit2_z;


        }    



    }
    
    // Perform standard residual calculation
    if(square == false)
    {    
        int FLOPSPERCALC = 10;
        std::cout<<"Arrays: Residual"<<std::endl;
        tstart = dtime();  

        #pragma ivdep
        {
            //#pragma omp parallel for private(i)   
            #pragma omp parallel for simd private(i)
            for(int i = 0; i < seed_entries; i++)
            {
                float slope = (hit2_radius_array[i] - hit0_radius_array[i]) / (hit2_z_array[i] - hit0_z_array[i]);
                float inter =  hit2_radius_array[i] - hit2_z_array[i] * slope;

                residual_array[i] = fabs (hit1_radius_array[i] - slope * hit1_z_array[i] - inter) /
                sqrtf (1.0f + slope * slope);  
            }
        } 
        tstop = dtime();
        ttime = tstop - tstart;
        gflops_res = (double)( 1.0e-9*FLOPSPERCALC*seed_entries);

        if( ttime > 0.0)
        {
            std::cout<<" GFLOPS = " << gflops_res << ", Time(s) = " << ttime << ", GFLOPS/s = " << 
                gflops_res/ttime << std::endl << std::endl;
        }    
    }

    // Perform squared residual calculation
    else
    {
        int FLOPSPERCALC = 10;
        std::cout<<"Arrays: Residual Squared"<<std::endl;
        tstart = dtime();  

        #pragma ivdep
        {
            //#pragma omp parallel for private(i)
            #pragma omp parallel for simd private(i)
            for(int i = 0; i < seed_entries; i++)
            {
                float slope = (hit2_radius_array[i] - hit0_radius_array[i]) / (hit2_z_array[i] - hit0_z_array[i]);
                float inter =  hit2_radius_array[i] - hit2_z_array[i] * slope;
                
                residual_array[i] = hit1_radius_array[i] - slope * hit1_z_array[i] - inter;
                residual_array[i] = residual_array[i] * residual_array[i] / (1.0f + slope * slope);                
            }
        } 
        tstop = dtime();
        ttime = tstop - tstart;

        if( ttime > 0.0)
        {
            std::cout<<" GFLOPS = " << gflops_sqr << ", Time(s) = " << ttime << ", GFLOPS/s = " << 
                gflops_sqr/ttime << std::endl << std::endl;
        }    

    }        

    // Free the aligned memory blocks 
    _mm_free(hit0_z_array);
    _mm_free(hit1_z_array);
    _mm_free(hit2_z_array);
    _mm_free(hit0_radius_array);
    _mm_free(hit1_radius_array);
    _mm_free(hit2_radius_array);
    _mm_free(slope_array);
    _mm_free(intercept_array);
    _mm_free(residual_array);
}    
*/

/////////////////////////////////////
//          Main Method            //
/////////////////////////////////////
int main()
{   
    // Read in the vector   
    std::cout<<"Reading file.."<<std::endl << std::endl;   
    std::vector<seed_t> seed_vector = ReadFile();
    seed_vector  = ExtendVector(seed_vector, 100); 

    // Allocate memory for and create the Array of Structs (smallArrays_t)
    smallArrays_t * AOSOA;
    AOSOA  = MakeAOSOA(seed_vector, 16);



    return 0;
}   

