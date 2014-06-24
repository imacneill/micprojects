#ifndef OBJECTBENCHMARKS_H
#define OBJECTBENCHMARKS_H

#include "Seeds.h"

 
void fillSeeds(const unsigned int NENTRIES, float *x0, float *y0, float *x1, float *y1, float *x2, float *y2, Seed2D *seed2D);
void copySeedToSeedsa(const unsigned int NENTRIES, float *x0, float *x1, float *x2, const unsigned int NENTRIES_sa, const unsigned int NENTRIES_sasize, Seed2Dsa *seed2Dsa); 



#endif
