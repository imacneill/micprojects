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
  //// one test with an array that holds a struct that holds pointers to each coordinate that have been alligned with 64 bits

  const  int NENTRIES = 10000000;
  const  int NENTRIES_sa = 156250;
  const  int NENTRIES_sasize = 16;
  const  int ALLIGN = 64;
  float *x0 = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *y0 = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *x1 = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *y1 = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *x2 = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *y2 = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  
  Seed2D *seed2D = (Seed2D*) _mm_malloc(sizeof(Seed2D)*NENTRIES, ALLIGN);
  Seed2Dsa *seed2Dsa = (Seed2Dsa*) _mm_malloc(sizeof(Seed2Dsa)*NENTRIES_sa, ALLIGN);
  Seed2Dhsa<NENTRIES_sasize> *seed2Dhsa = (Seed2Dhsa<NENTRIES_sasize>*) _mm_malloc(sizeof(Seed2Dhsa<NENTRIES_sasize>)*NENTRIES_sa, ALLIGN);

  
  fillSeeds(NENTRIES, x0, y0, x1, y1, x2, y2, seed2D);
  copySeedToSeedsa(NENTRIES, x0, x1, x2, NENTRIES_sa, NENTRIES_sasize, seed2Dsa);

  for(unsigned int i=0; i<NENTRIES_sa; ++i){
	if(NENTRIES > (i * NENTRIES_sasize)){
	  seed2Dhsa[i].setX(x0, x1, x2, i * NENTRIES_sasize, NENTRIES_sasize);
	}
  }


  

  
  //// make arrays to hold the residuals that are calculated
  float *residd = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *resid2Dd = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *resid2Dsad = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *resid2Dhsad = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);

  float *resid = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *resid2D = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *resid2Dsa = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *resid2Dhsa = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);

  //explicitly using inline functions for these
  float *residIn = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *resid2DIn = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *resid2DsaIn = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *resid2DhsaIn = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  
  



  //// make the loops
  /// no function
  stopwatch wall_bmad("Wall Array Direct", &get_wall_time);
  wall_bmad.timestart();
#pragma simd
  for(int i = 0; i < NENTRIES; ++i){
	float slope = (y2[i] - y0[i])/(x2[i] - x0[i]);
	float intercept = y0[i] - slope * x0[i];
	residd[i] = std::abs( slope*x0[i] - y0[i] + intercept ) / sqrt( slope*slope + intercept*intercept ); 
  }
  wall_bmad.timestop();
  if( wall_bmad.done() ){ std::cout<<wall_bmad<<std::endl; }
  
  //struct
  stopwatch wall_bmsd("Wall Struct Direct", &get_wall_time);
  wall_bmsd.timestart();
#pragma simd
  for(int i = 0; i < NENTRIES; ++i){
	float slope = (seed2D[i].y2_ - seed2D[i].y0_)/(seed2D[i].x2_ - seed2D[i].x0_);
	float intercept = seed2D[i].y0_ - slope * seed2D[i].x0_;
	 resid2Dd[i] = std::abs( slope*seed2D[i].x0_ - seed2D[i].y0_ + intercept ) / sqrt( slope*slope + intercept*intercept );
  }
  wall_bmsd.timestop();
  if( wall_bmsd.done() ){ std::cout<<wall_bmsd<<std::endl; }


  //struct short array
  stopwatch wall_bmsad("Wall Struct Short Array Direct", &get_wall_time);
  wall_bmsad.timestart();
  for(int i = 0; i < NENTRIES_sa; ++i){
#pragma simd
	for(int j = 0; j < NENTRIES_sasize; ++j){
	  float slope = (seed2Dsa[i].y2_[j] - seed2Dsa[i].y0_[j])/(seed2Dsa[i].x2_[j] - seed2Dsa[i].x0_[j]);
	  float intercept = seed2Dsa[i].y0_[j] - slope * seed2Dsa[i].x0_[j];
	  resid2Dsad[i*NENTRIES_sasize+j] = std::abs( slope*seed2Dsa[i].x0_[j] - seed2Dsa[i].y0_[j] + intercept ) / sqrt( slope*slope + intercept*intercept );
	}
  }
  wall_bmsad.timestop();
  if( wall_bmsad.done() ){ std::cout<<wall_bmsad<<std::endl; }


  //struct short array hard coded
  stopwatch wall_bmhsad("Wall Struct Hand Short Array Direct", &get_wall_time);
  wall_bmhsad.timestart();
  for(int i = 0; i < NENTRIES_sa; ++i){
#pragma simd
	for(int j = 0; j < NENTRIES_sasize; ++j){
	  float slope = (seed2Dhsa[i].y2_[j] - seed2Dhsa[i].y0_[j])/(seed2Dhsa[i].x2_[j] - seed2Dhsa[i].x0_[j]);
	  float intercept = seed2Dhsa[i].y0_[j] - slope * seed2Dhsa[i].x0_[j];
	  resid2Dhsad[i*NENTRIES_sasize+j] = std::abs( slope*seed2Dhsa[i].x0_[j] - seed2Dhsa[i].y0_[j] + intercept ) / sqrt( slope*slope + intercept*intercept );
	}
  }
  wall_bmhsad.timestop();
  if( wall_bmhsad.done() ){ std::cout<<wall_bmhsad<<std::endl; }


  std::cout<<std::endl;



  ///no explicit inline
  //array
  stopwatch wall_bma("Wall Array", &get_wall_time);
  wall_bma.timestart();
#pragma simd
  //  #pragma omp parallel for simd
  for(int i = 0; i < NENTRIES; ++i){
	resid[i] = calcResid(x0[i], y0[i], x1[i], y1[i], x2[i], y2[i]);
  }
  wall_bma.timestop();
  if( wall_bma.done() ){ std::cout<<wall_bma<<std::endl; }
  
  //struct
  stopwatch wall_bms("Wall Struct", &get_wall_time);
  wall_bms.timestart();
#pragma simd
  //#pragma omp parallel for simd
  for(int i = 0; i < NENTRIES; ++i){
	resid2D[i] = calcResid(seed2D[i]);
  }
  wall_bms.timestop();
  if( wall_bms.done() ){ std::cout<<wall_bms<<std::endl; }

  //struct short array
  stopwatch wall_bmsa("Wall Struct Short Array", &get_wall_time);
  wall_bmsa.timestart();
  //#pragma simd
  //#pragma omp parallel for simd
  for(int i = 0; i < NENTRIES_sa; ++i){
	//#pragma simd
	//	for(int j = 0; j < NENTRIES_sasize; ++j){
	//  resid2Dsad[i*NENTRIES_sasize+j] = calcResid(seed2Dsa[i], j);
	calcResid(seed2Dsa[i], resid2Dsa, i*NENTRIES_sasize, (i+1)*NENTRIES_sasize);

	//}
  }
  wall_bmsa.timestop();
  if( wall_bmsa.done() ){ std::cout<<wall_bmsa<<std::endl; }

  //struct short array hard coded
  stopwatch wall_bmhsa("Wall Struct Short Array Hand", &get_wall_time);
  wall_bmhsa.timestart();
  //#pragma simd
  //#pragma omp parallel for simd
  for(int i = 0; i < NENTRIES_sa; ++i){
	//#pragma simd
	//	for(int j = 0; j < NENTRIES_sasize; ++j){
	//  resid2Dsad[i*NENTRIES_sasize+j] = calcResid(seed2Dsa[i], j);
	calcResid(seed2Dhsa[i], resid2Dhsa, i*NENTRIES_sasize, (i+1)*NENTRIES_sasize);

	//}
  }
  wall_bmhsa.timestop();
  if( wall_bmhsa.done() ){ std::cout<<wall_bmhsa<<std::endl; }
  


  std::cout<<std::endl;



  /// explicit inline
  //array
  stopwatch wall_bmai("Wall Array Inline", &get_wall_time);
  wall_bmai.timestart();
#pragma simd
  //#pragma omp parallel for simd
  for(int i = 0; i < NENTRIES; ++i){
	residIn[i] = calcResidIn(x0[i], y0[i], x1[i], y1[i], x2[i], y2[i]);
  }
  wall_bmai.timestop();
  if( wall_bmai.done() ){ std::cout<<wall_bmai<<std::endl; }
  
  //struct
  stopwatch wall_bmsi("Wall Struct Inline", &get_wall_time);
  wall_bmsi.timestart();
#pragma simd
  //#pragma omp parallel for simd
  for(int i = 0; i < NENTRIES; ++i){
	resid2DIn[i] = calcResidIn(seed2D[i]);
  }
  wall_bmsi.timestop();
  if( wall_bmsi.done() ){ std::cout<<wall_bmsi<<std::endl; }

  //struct short array
  stopwatch wall_bmsai("Wall Struct Short Array Inline", &get_wall_time);
  wall_bmsai.timestart();
  for(int i = 0; i < NENTRIES_sa; ++i){
	//#pragma simd
	//	for(int j = 0; j < NENTRIES_sasize; ++j){
	//	  resid2DsaIn[i*NENTRIES_sasize+j] = calcResidIn(seed2Dsa[i], j);
	calcResidIn(seed2Dsa[i], resid2DsaIn, i*NENTRIES_sasize, (i+1)*NENTRIES_sasize);
	  //	}
  }
  wall_bmsai.timestop();
  if( wall_bmsai.done() ){ std::cout<<wall_bmsai<<std::endl; }
  
  
  
  
  
  std::cout<<std::endl<<std::endl<<std::endl<<"Garbage: "
		   <<sumArray(residd, NENTRIES)
		   <<" "<<sumArray(resid2Dd, NENTRIES)
		   <<" "<<sumArray(resid2Dsad, NENTRIES)
		   <<" "<<sumArray(resid2Dhsad, NENTRIES)
		   <<" "<<sumArray(resid, NENTRIES)
		   <<" "<<sumArray(resid2D, NENTRIES)
		   <<" "<<sumArray(resid2Dsa, NENTRIES)
		   <<" "<<sumArray(resid2Dhsa, NENTRIES)
		   <<" "<<sumArray(residIn, NENTRIES)
		   <<" "<<sumArray(resid2DIn, NENTRIES)
		   <<" "<<sumArray(resid2DsaIn, NENTRIES)
		   <<" "<<sumArray(resid2DhsaIn, NENTRIES)
		   <<std::endl;
  
  
  
  _mm_free(x0);
  _mm_free(y0);
  _mm_free(x1);
  _mm_free(y1);
  _mm_free(x2);
  _mm_free(y2);
  _mm_free(seed2D);
  _mm_free(seed2Dsa);
  _mm_free(seed2Dhsa);


  _mm_free(residd);
  _mm_free(resid2Dd);
  _mm_free(resid2Dsad);
  _mm_free(resid2Dhsad);

  _mm_free(resid);
  _mm_free(resid2D);
  _mm_free(resid2Dsa);
  _mm_free(resid2Dhsa);

  _mm_free(residIn);
  _mm_free(resid2DIn);
  _mm_free(resid2DsaIn);
  _mm_free(resid2DhsaIn);


  
  return 0;
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

void copySeedToSeedsa(const unsigned int NENTRIES, float *x0, float *x1, float *x2, const unsigned int NENTRIES_sa, const unsigned int NENTRIES_sasize, Seed2Dsa *seed2Dsa){
  for(unsigned int i=0; i<NENTRIES_sa; ++i){
	if(NENTRIES > (i * NENTRIES_sasize)){
	  seed2Dsa[i].setX(x0, x1, x2, i * NENTRIES_sasize, NENTRIES_sasize);
	}
  }
}


