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

// User
#include "Seed.h"
#include "ReadSeedFile.h"


/////////////////////////////////////
//             Methods             //
/////////////////////////////////////

// Reads the result of the .txt and outputs a vector of the seed structs
std::vector<seed_t> ReadFile()
{
    using namespace std;

    // Define the ifstream
    std::ifstream datafile;
    datafile.open("../data/seedDataFile.txt", std::ifstream::in);
    
    // Check the validitiy of the ifstream
    if(!datafile.good())
    {
        std::cout<<"Error opening the seed data file" << std::endl;
    }    
        
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