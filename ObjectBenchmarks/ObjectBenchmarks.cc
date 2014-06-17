#include <iostream>
#include <cmath>
#include <omp.h>

#include "Seeds.h"
#include "BenchmarkTools.h"
#include "ObjectBenchmarks.h"



int main(){
  
  //// make test arrays
  //// one test with an array for each coordinate
  //// one test with an array that holds a struct that holds each coordinate
  //// one test with an array that holds a struct that holds pointers to each coordinate
  //// one test with an array that holds a struct that holds pointers to each coordinate that have been alligned with 64 bits

  const unsigned int NENTRIES = 10000000;
  const unsigned int ALLIGN = 64;
  float *x0 = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *y0 = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *x1 = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *y1 = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *x2 = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *y2 = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  
  Seed2D *seed2D = (Seed2D*) _mm_malloc(sizeof(Seed2D)*NENTRIES, ALLIGN);
  // Seed2Dp *seed2Dp = (Seed2Dp*) _mm_malloc(sizeof(Seed2Dp)*NENTRIES, ALLIGN);
  // Seed2D64a *seed2D64a = (Seed2D64a*) _mm_malloc(sizeof(Seed2D64a)*NENTRIES, ALLIGN);
  
  //  fillSeeds(NENTRIES, x0, y0, x1, y1, x2, y2, seed2D, seed2Dp, seed2D64a);
  fillSeeds(NENTRIES, x0, y0, x1, y1, x2, y2, seed2D);
  
  //// make arrays to hold the residuals that are calculated
  float *residd = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *resid2Dd = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  //
  float *resid = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *resid2D = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  // float *resid2Dp = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  // float *resid2D64a = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  //explicitly using inline functions for these
  float *residIn = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *resid2DIn = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  // float *resid2DpIn = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  // float *resid2D64aIn = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  
  
  //// make the loops
  /// no function
  stopwatch wall_bmad("Wall Array Direct", &get_wall_time);
  wall_bmad.timestart();
  //  #pragma simd
  #pragma omp parallel for simd
  for(unsigned int i = 0; i < NENTRIES; ++i){
	float slope = (y2[i] - y0[i])/(x2[i] - x0[i]);
	float intercept = y0[i] - slope * x0[i];
	residd[i] = std::abs( slope*x0[i] - y0[i] + intercept ) / sqrt( slope*slope + intercept*intercept );
  }
  wall_bmad.timestop();
  if( wall_bmad.done() ){ std::cout<<wall_bmad<<std::endl; }
  
  //struct
  stopwatch wall_bmsd("Wall Struct Direct", &get_wall_time);
  wall_bmsd.timestart();
  //  #pragma simd
#pragma omp parallel for simd
  for(unsigned int i = 0; i < NENTRIES; ++i){
	float slope = (seed2D[i].y2_ - seed2D[i].y0_)/(seed2D[i].x2_ - seed2D[i].x0_);
	float intercept = seed2D[i].y0_ - slope * seed2D[i].x0_;
	 resid2Dd[i] = std::abs( slope*seed2D[i].x0_ - seed2D[i].y0_ + intercept ) / sqrt( slope*slope + intercept*intercept );
  }
  wall_bmsd.timestop();
  if( wall_bmsd.done() ){ std::cout<<wall_bmsd<<std::endl; }


  ///no explicit inline
  //array
  stopwatch wall_bma("Wall Array", &get_wall_time);
  wall_bma.timestart();
  //  #pragma simd
  #pragma omp parallel for simd
  for(unsigned int i = 0; i < NENTRIES; ++i){
	resid[i] = calcResid(x0[i], y0[i], x1[i], y1[i], x2[i], y2[i]);
  }
  wall_bma.timestop();
  if( wall_bma.done() ){ std::cout<<wall_bma<<std::endl; }
  
  //struct
  stopwatch wall_bms("Wall Struct", &get_wall_time);
  wall_bms.timestart();
  //  #pragma simd
#pragma omp parallel for simd
  for(unsigned int i = 0; i < NENTRIES; ++i){
	resid2D[i] = calcResid(seed2D[i]);
  }
  wall_bms.timestop();
  if( wall_bms.done() ){ std::cout<<wall_bms<<std::endl; }
  
  // //struct pointer
  // stopwatch wall_bmp("Wall Struct Pointer", &get_wall_time);
  // wall_bmp.timestart();
  // for(unsigned int i = 0; i < NENTRIES; ++i){
  // 	resid2Dp[i] = calcResid(seed2Dp[i]);
  // }
  // wall_bmp.timestop();
  // if( wall_bmp.done() ){ std::cout<<wall_bmp<<std::endl; }
  
  // //struct pointer memory alligned
  // stopwatch wall_bm64a("Wall Struct Pointer Alligned", &get_wall_time);
  // wall_bm64a.timestart();
  // for(unsigned int i = 0; i < NENTRIES; ++i){
  // 	resid2D64a[i] = calcResid(seed2D64a[i]);
  // }
  // wall_bm64a.timestop();
  // if( wall_bm64a.done() ){ std::cout<<wall_bm64a<<std::endl; }

  /// explicit inline
  //array
  stopwatch wall_bmai("Wall Array Inline", &get_wall_time);
  wall_bmai.timestart();
  //  #pragma simd
#pragma omp parallel for simd
  for(unsigned int i = 0; i < NENTRIES; ++i){
	residIn[i] = calcResid(x0[i], y0[i], x1[i], y1[i], x2[i], y2[i]);
  }
  wall_bmai.timestop();
  if( wall_bmai.done() ){ std::cout<<wall_bmai<<std::endl; }
  
  //struct
  stopwatch wall_bmsi("Wall Struct Inline", &get_wall_time);
  wall_bmsi.timestart();
  //  #pragma simd
#pragma omp parallel for simd
  for(unsigned int i = 0; i < NENTRIES; ++i){
	resid2DIn[i] = calcResid(seed2D[i]);
  }
  wall_bmsi.timestop();
  if( wall_bmsi.done() ){ std::cout<<wall_bmsi<<std::endl; }
  
  // //struct pointer
  // stopwatch wall_bmpi("Wall Struct Pointer Inline", &get_wall_time);
  // wall_bmpi.timestart();
  // for(unsigned int i = 0; i < NENTRIES; ++i){
  // 	resid2DpIn[i] = calcResid(seed2Dp[i]);
  // }
  // wall_bmpi.timestop();
  // if( wall_bmpi.done() ){ std::cout<<wall_bmpi<<std::endl; }
  
  // //struct pointer memory alligned
  // stopwatch wall_bm64ai("Wall Struct Pointer Alligned Inline", &get_wall_time);
  // wall_bm64ai.timestart();
  // for(unsigned int i = 0; i < NENTRIES; ++i){
  // 	resid2D64aIn[i] = calcResid(seed2D64a[i]);
  // }
  // wall_bm64ai.timestop();
  // if( wall_bm64ai.done() ){ std::cout<<wall_bm64ai<<std::endl; }
  
  
  
  
  std::cout<<std::endl<<std::endl<<std::endl<<"Garbage: "
		   <<sumArray(residd, NENTRIES)
		   <<" "<<sumArray(resid2Dd, NENTRIES)
		   <<sumArray(resid, NENTRIES)
		   <<" "<<sumArray(resid2D, NENTRIES)//<<" "<<sumSeedArray(seed2Dp, NENTRIES)<<" "<<sumSeedArray(seed2D64a, NENTRIES)
		   <<" "<<sumArray(residIn, NENTRIES)
		   <<" "<<sumArray(resid2DIn, NENTRIES)
		   <<std::endl;
  
  
  
  _mm_free(x0);
  _mm_free(y0);
  _mm_free(x1);
  _mm_free(y1);
  _mm_free(x2);
  _mm_free(y2);
  
  _mm_free(seed2D);
  // _mm_free(seed2Dp);
  // _mm_free(seed2D64a);


  _mm_free(residd);
  _mm_free(resid2Dd);
  //
  _mm_free(resid);
  _mm_free(resid2D);
  //  _mm_free(resid2Dp);
  //  _mm_free(resid2D64a);
  //explicitly using inline functions for these
  _mm_free(residIn);
  _mm_free(resid2DIn);
  //  _mm_free(resid2DpIn);
  //  _mm_free(resid2D64aIn);
  
  return 0;
}




void fillSeeds(const unsigned int NENTRIES, float *x0, float *y0, float *x1, float *y1, float *x2, float *y2, 
			   Seed2D *seed2D, Seed2Dp *seed2Dp, Seed2D64a *seed2D64a){

  #pragma omp parallel for simd
  for(unsigned int i=0; i<NENTRIES; i++){
	generateSeeds((seed2D+i), (seed2Dp+i), (seed2D64a+i));
	x0[i]=seed2D[i].x0_;
	y0[i]=seed2D[i].y0_;
	x1[i]=seed2D[i].x1_;
	y1[i]=seed2D[i].y1_;
	x2[i]=seed2D[i].x2_;
	y2[i]=seed2D[i].y2_;
  }

}

void fillSeeds(const unsigned int NENTRIES, float *x0, float *y0, float *x1, float *y1, float *x2, float *y2, 
			   Seed2D *seed2D){

  #pragma omp parallel for simd
  for(unsigned int i=0; i<NENTRIES; i++){
	generateSeeds((seed2D+i));
	x0[i]=seed2D[i].x0_;
	y0[i]=seed2D[i].y0_;
	x1[i]=seed2D[i].x1_;
	y1[i]=seed2D[i].y1_;
	x2[i]=seed2D[i].x2_;
	y2[i]=seed2D[i].y2_;
  }
}
