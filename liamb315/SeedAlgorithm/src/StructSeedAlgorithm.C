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
    double tseconds = 0.0;
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

// Function that returns the Residual of a Particular Seed
inline float CalcResidual(const seed_t & seed)
{
    float resid = 0;

    // Retrieve seed hit cooridinates in r-z plane
    float radius0    = seed.hit0_radius;
    float radius1    = seed.hit1_radius;
    float radius2    = seed.hit2_radius;

    float seed_hit0z = seed.hit0_z;
    float seed_hit1z = seed.hit1_z;
    float seed_hit2z = seed.hit2_z;

    // Approximate best-fit line parameters
    float slope     = (radius2 - radius0)/(seed_hit2z - seed_hit0z);
    float intercept = radius2 - seed_hit2z * slope;

    // Residual distance of third seed_hit coordinate from line
    resid = fabs(radius1 - slope*seed_hit1z - intercept)
          / sqrtf(slope*slope + 1);
    
    return resid;
}

// Function that returns the Residual of a Particular Seed
inline float CalcSquareResidual(const seed_t & seed)
{
    float sqresid = 0;

    // Retrieve seed hit cooridinates in r-z plane
    float radius0    = seed.hit0_radius;
    float radius1    = seed.hit1_radius;
    float radius2    = seed.hit2_radius;

    float seed_hit0z = seed.hit0_z;
    float seed_hit1z = seed.hit1_z;
    float seed_hit2z = seed.hit2_z;

    // Approximate best-fit line parameters
    float slope     = (radius2 - radius0)/(seed_hit2z - seed_hit0z);
    float intercept = radius2 - seed_hit2z * slope;

    // Square residual distance of third seed_hit coordinate from line
    sqresid = (radius1 - slope*seed_hit1z - intercept)
            * (radius1 - slope*seed_hit1z - intercept)
            / (slope*slope + 1);
        
    return sqresid;
}


/////////////////////////////////////
//          Main Method            //
/////////////////////////////////////
int main()
{
    double tstart0, tstop0, ttime0;
    double tstart0_imp, tstop0_imp, ttime0_imp;

    std::cout<<"Reading file.."<<std::endl << std::endl;   
    
    // Read in the vector   
    std::vector<seed_t> seed_vector = ReadFile();
    seed_vector  = ExtendVector(seed_vector, 100); 
    const int seed_entries = seed_vector.size();
  
    ////////////////////////////
    // Structs:  Residual 
    ////////////////////////////
    float * residuals;
    residuals = (float*) _mm_malloc(sizeof(float)*seed_entries,64);

    std::cout<<"Structs: Residual" << std::endl;
    int i = 0, num_threads, chunk;

    #pragma omp parallel
    #pragma omp master
    num_threads = omp_get_num_threads();
    chunk = seed_entries/num_threads;
    
    tstart0 = dtime();   
    #pragma ivdep
    {     
        //#pragma omp parallel for private(i)
        #pragma omp parallel for simd private(i)
        for(i = 0; i < seed_entries; i++)
        {
            residuals[i] = CalcResidual(seed_vector[i]);
        }
    }

    tstop0 = dtime();
    ttime0 = tstop0 - tstart0;
    std::cout<<" Elapsed time: " << ttime0 <<"s"<< std::endl<<std::endl;


    ////////////////////////////
    // Structs:  Square Residual 
    ////////////////////////////
    float * square_residuals;
    square_residuals = (float*) _mm_malloc(sizeof(float)*seed_entries,64);

    std::cout<<"Structs: Residual Squared" << std::endl;
    //int i, num_threads, chunk;

    tstart0_imp = dtime();   
    #pragma ivdep
    {     
        //#pragma omp parallel for private(i)
        #pragma omp parallel for simd private(i)
        for(i = 0; i < seed_entries; i++)
        {
            square_residuals[i] = CalcSquareResidual(seed_vector[i]);
        }
    }

    tstop0_imp = dtime();
    ttime0_imp = tstop0_imp - tstart0_imp;
    std::cout<<" Elapsed time: " << ttime0_imp <<"s"<< std::endl<<std::endl;

    
    ////////////////////////////
    // Error Checking 
    ////////////////////////////
    int tot_errs = 0;
    for(int i = 0; i < seed_entries; i++)
    {   
        float epsilon = 0.00001;
        
        if(fabs(residuals[i] - sqrt(square_residuals[i])) > epsilon)
        {    
            std::cout<<residuals[i] - sqrt(square_residuals[i]) << std::endl;
            tot_errs++;
        }
    }
    std::cout<<"Total errors for square method: " << tot_errs << std::endl;

    _mm_free(residuals);
    _mm_free(square_residuals);

    return 0;
}   

