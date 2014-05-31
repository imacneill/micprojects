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

// Struct to contain the data of the Root file
struct seed_t
{
    // Members of the struct
    float hit0_radius, hit1_radius, hit2_radius;
    float hit0_z, hit1_z, hit2_z;
};

// Reads the result of the .txt and outputs a vector of the seed structs
std::vector<seed_t> ReadFile()
{
    using namespace std;

    // Define the ifstream
    std::ifstream datafile;
    datafile.open("data/seedDataFile.txt", std::ifstream::in);

    // Vector of seed structs to store
    std::vector<seed_t>  seed_vector;

    std::string line;
    
    while(std::getline(datafile, line))
    {
        std::stringstream line_stream(line);
        seed_t seed;
        std::string field;

        // Read in a CSV format
        std::getline(line_stream, field, ','); seed.hit0_radius = atof(field.c_str());
        std::getline(line_stream, field, ','); seed.hit1_radius = atof(field.c_str());
        std::getline(line_stream, field, ','); seed.hit2_radius = atof(field.c_str());
        std::getline(line_stream, field, ','); seed.hit0_z = atof(field.c_str());
        std::getline(line_stream, field, ','); seed.hit1_z = atof(field.c_str());
        std::getline(line_stream, field, ','); seed.hit2_z = atof(field.c_str());
        
        // Add each seed struct into the vector
        seed_vector.push_back(seed);
    }
    
    datafile.close();
 
    return seed_vector;
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
          / sqrt(slope*slope + 1);
    
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

// Calculate the Residuals Using Arrays
void ComputeArrayResidual(const std::vector<seed_t> & seed_vector, bool square = false)
{
    if(seed_vector.empty() == true)
    {
        std::cout<< "Seed vector empty" << std::endl;
        return;    
    }        
    
    double tstart = 0, tstop = 0, ttime = 0;
    int i = 0, num_threads = 0, chunk = 0;

    int seed_entries = seed_vector.size();

    // Define array pointers
    float * hit0_z_array;
    float * hit1_z_array;
    float * hit2_z_array;
    float * hit0_radius_array;
    float * hit1_radius_array;
    float * hit2_radius_array;
    float * slope_array;
    float * intercept_array;
    float * residual_array;

    // Allocate memory 
    hit0_z_array      = (float*) _mm_malloc(sizeof(float)*seed_entries, 64); 
    hit1_z_array      = (float*) _mm_malloc(sizeof(float)*seed_entries, 64); 
    hit2_z_array      = (float*) _mm_malloc(sizeof(float)*seed_entries, 64); 
    hit0_radius_array = (float*) _mm_malloc(sizeof(float)*seed_entries, 64); 
    hit1_radius_array = (float*) _mm_malloc(sizeof(float)*seed_entries, 64); 
    hit2_radius_array = (float*) _mm_malloc(sizeof(float)*seed_entries, 64); 
    slope_array       = (float*) _mm_malloc(sizeof(float)*seed_entries, 64); 
    intercept_array   = (float*) _mm_malloc(sizeof(float)*seed_entries, 64); 
    residual_array    = (float*) _mm_malloc(sizeof(float)*seed_entries, 64); 

    #pragma omp parallel
    #pragma omp master
    num_threads = omp_get_num_threads();
    chunk = seed_entries/num_threads;

    // Pragma tells the compiler to ignore non-obvious dependencies
    //#pragma ivdep
    //#pragma omp parallel for
    {
        for(int i = 0; i < seed_entries; i++)
        {    
            hit0_radius_array[i] = seed_vector[i].hit0_radius; 
            hit1_radius_array[i] = seed_vector[i].hit1_radius; 
            hit2_radius_array[i] = seed_vector[i].hit2_radius; 
            hit0_z_array[i]      = seed_vector[i].hit0_z; 
            hit1_z_array[i]      = seed_vector[i].hit1_z; 
            hit2_z_array[i]      = seed_vector[i].hit2_z; 
        }
    }
    
    // Perform standard residual calculation
    if(square == false)
    {    
        std::cout<<"Arrays: Residual"<<std::endl;
        tstart = dtime();  

        //#pragma omp parallel shared(chunk) private(i)
        //{
            
            #pragma ivdep
            {
            //#pragma omp for schedule(dynamic, chunk) nowait  // **LOOP WAS NOT VECTORIZED:  NOT INNER LOOP**
            #pragma omp parallel for private(i)   // **LOOP WAS NOT VECTORIZED:  NOT INNER LOOP**
            for(int i = 0; i < seed_entries; i++)
            {
                residual_array[i]  = fabs(hit1_radius_array[i] - (hit2_radius_array[i] - hit0_radius_array[i])
                        / (hit2_z_array[i] - hit0_z_array[i])*hit1_z_array[i] - (hit2_radius_array[i] - hit2_z_array[i]*slope_array[i]))
                    / sqrt((  hit2_radius_array[i] - hit0_radius_array[i])/(hit2_z_array[i] - hit0_z_array[i]   ) 
                            * (hit2_radius_array[i] - hit0_radius_array[i])/(hit2_z_array[i] - hit0_z_array[i]) + 1); 
          //  }
            }

        } 
        tstop = dtime();
        ttime = tstop - tstart;
        std::cout<<" Elapsed time: " << ttime <<"s"<< std::endl<<std::endl;
    }

    // Perform squared residual calculation
    else
    {
        std::cout<<"Arrays: Residual Squared"<<std::endl;
        tstart = dtime();  

        // 1.  #pragma omp parallel shared(chunk) private(i)
        #pragma ivdep
        {
            // In the parallel region, query how many threads are available and divide the arrays
            //2.  num_threads = omp_get_num_threads();
            //3.  chunk = seed_entries/num_threads;

            //#pragma omp for schedule(dynamic, chunk) nowait  // **LOOP WAS NOT VECTORIZED:  NOT INNER LOOP**
            //4.  #pragma omp for schedule(dynamic, chunk)   // **LOOP WAS NOT VECTORIZED:  NOT INNER LOOP**
            #pragma omp parallel for private(i)
            for(int i = 0; i < seed_entries; i++)
            {
                slope_array[i]     = (hit2_radius_array[i] - hit0_radius_array[i])/(hit2_z_array[i] - hit0_z_array[i]); 
                intercept_array[i] = hit2_radius_array[i] - hit2_z_array[i]*slope_array[i];
                residual_array[i]  = (hit1_radius_array[i] - slope_array[i]*hit1_z_array[i] - intercept_array[i]);
                residual_array[i]  = residual_array[i] * residual_array[i]            
                                   / (slope_array[i] * slope_array[i] + 1); 
            }
        } 
        tstop = dtime();
        ttime = tstop - tstart;
        std::cout<<" Elapsed time: " << ttime <<"s"<< std::endl<<std::endl;
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


/////////////////////////////////////
//          Main Method            //
/////////////////////////////////////
int main()
{
    
    ////////////////////////////
    //   Set Up
    ////////////////////////////
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
        #pragma omp parallel for private(i)        
        for(i = 0; i < seed_entries; i++)
        {
            residuals[i] = CalcResidual(seed_vector[i]);
        }
    }
  
    /* 2.
    #pragma ivdep 
    {
    #pragma omp parallel shared(chunk) private(i) 
    //#pragma omp parallel shared(seed_vector, chunk) private(i) 
        { 
            // In the parallel region, query how many threads are available and divide the arrays
            num_threads = omp_get_num_threads();
            chunk = seed_entries/num_threads;
   
            #pragma omp for schedule(dynamic, chunk) nowait // **LOOP WAS NOT VECTORIZED:  NOT INNER LOOP**
            for(i = 0; i < seed_entries; i++)
            {
                residuals[i] = CalcResidual(seed_vector[i]);
            }
        }// end of parallel section
    }
    */

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
        #pragma omp parallel for private(i)        
        for(i = 0; i < seed_entries; i++)
        {
            square_residuals[i] = CalcSquareResidual(seed_vector[i]);
        }
    }


    /*  3.  
    #pragma ivdep 
    {
    #pragma omp parallel shared(chunk) private(i) 
    //#pragma omp parallel shared(seed_vector, chunk) private(i) 
        { 
            // In the parallel region, query how many threads are available and divide the arrays
            num_threads = omp_get_num_threads();
            chunk = seed_entries/num_threads;
    
            #pragma omp for schedule(dynamic, chunk) nowait  // **LOOP WAS NOT VECTORIZED:  NOT INNER LOOP**
            for(i = 0; i < seed_vector.size(); i++)
            {
                square_residuals[i] = CalcSquareResidual(seed_vector[i]);
            }
        }// end of parallel section
    }
    */
    tstop0_imp = dtime();
    ttime0_imp = tstop0_imp - tstart0_imp;
    std::cout<<" Elapsed time: " << ttime0_imp <<"s"<< std::endl<<std::endl;


    ////////////////////////////
    // Arrays:  Residual
    ////////////////////////////
    ComputeArrayResidual(seed_vector, false);
   
    
    ////////////////////////////
    // Arrays:  Square Residual
    ////////////////////////////
    ComputeArrayResidual(seed_vector, true);

    
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

