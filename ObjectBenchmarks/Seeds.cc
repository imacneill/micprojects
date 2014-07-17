#include <random>
#include <iostream>
#include <cmath>
#include "Seeds.h"

std::default_random_engine generator(0xbeef0133);
std::normal_distribution<float> distribution(0.,1.);


float calcResid(const Seed2D& seed){
  float slope = (seed.y2_ - seed.y0_)/(seed.x2_ - seed.x0_);
  float intercept = seed.y0_ - slope * seed.x0_;
  return std::abs( slope*seed.x0_ - seed.y0_ + intercept ) / sqrt( slope*slope + intercept*intercept );
}

void calcResid(const Seed2Dsa& seed, float *resid, const int begin, const int end){
  #pragma simd
  for(int i = 0; i<(end-begin); ++i){
	float slope = ((seed.y2_[i]) - (seed.y0_[i]))/((seed.x2_[i]) - (seed.x0_[i]));
	float intercept = (seed.y0_[i]) - slope * (seed.x0_[i]);
	resid[i+begin] = std::abs( slope*(seed.x0_[i]) - (seed.y0_[i]) + intercept ) / sqrt( slope*slope + intercept*intercept );
  }
}

float calcResid(const float x0, const float y0, const float x1, const float y1, const float x2, const float y2){
  float slope = (y2 - y0)/(x2 - x0);
  float intercept = y0 - slope * x0;
  return std::abs( slope*x0 - y0 + intercept ) / sqrt( slope*slope + intercept*intercept );
}




void calcResidVectorized(Seed2D *seed, float *resid, const int begin, const int end){
#pragma simd
  for(int i = 0; i<(end-begin); ++i){
	float slope = ((seed[i].y2_) - (seed[i].y0_))/((seed[i].x2_) - (seed[i].x0_));
	float intercept = (seed[i].y0_) - slope * (seed[i].x0_);
	// resid[i+begin] = std::abs( slope*(seed[i].x0_) - (seed[i].y0_) + intercept ) / sqrt( slope*slope + intercept*intercept );
	resid[i+begin] = ( slope*(seed[i].x0_) - (seed[i].y0_) + intercept ) / ( slope*slope + intercept*intercept );
  }
}
void calcResidVectorized(const Seed2Dsa& seed, float *resid, const int begin, const int end){
#pragma simd
  for(int i = 0; i<(end-begin); ++i){
	float slope = ((seed.y2_[i]) - (seed.y0_[i]))/((seed.x2_[i]) - (seed.x0_[i]));
	float intercept = (seed.y0_[i]) - slope * (seed.x0_[i]);
	resid[i+begin] = std::abs( slope*(seed.x0_[i]) - (seed.y0_[i]) + intercept ) / sqrt( slope*slope + intercept*intercept );
  }
}
void calcResidVectorized(float *x0, float *y0, float *x1, float *y1, float *x2, float *y2, float *resid,  int begin,  int end){
  #pragma simd
  for(int i = 0; i<(end-begin); ++i){
	float slope = ((y2[i]) - (y0[i]))/((x2[i]) - (x0[i])); // 3 operations
	float intercept = (y0[i]) - slope * (x0[i]); // 2 operations
	//	resid[i+begin] = std::abs( slope*(x0[i]) - (y0[i]) + intercept ) / sqrt( slope*slope + intercept*intercept );
	resid[i+begin] = (slope*(x0[i]) - (y0[i]) + intercept ) / ( slope*slope + intercept*intercept ); // 7 operations
	// 12 operations total
  }
}

void calcResidVectorizedDouble(float *x0, float *y0, float *x1, float *y1, float *x2, float *y2, float *resid, const int inner, const int outer){
  for(int i=0; i<outer; ++i){
	#pragma simd
	for(int j=0; j<inner; ++j){
	  float slope = ((y2[i*inner+j]) - (y0[i*inner+j]))/((x2[i*inner+j]) - (x0[i*inner+j])); // 3 operations
	  float intercept = (y0[i*inner+j]) - slope * (x0[i*inner+j]); // 2 operations
	  //	resid[i+begin] = std::abs( slope*(x0[i*inner+j]) - (y0[i*inner+j]) + intercept ) / sqrt( slope*slope + intercept*intercept );
	  resid[i*inner+j] = (slope*(x0[i*inner+j]) - (y0[i*inner+j]) + intercept ) / ( slope*slope + intercept*intercept ); // 7 operations
	  // 12 operations total
	}
  }
}






void generateSeeds(Seed2D* seed2D){
  float tmp0 = distribution(generator);
  float tmp1 = distribution(generator);
  float tmp2 = distribution(generator);
  seed2D->setX(tmp0, tmp1, tmp2);
}





float sumSeedArray(Seed2D* seed2D, const int length){
  float tmp = 0;
  for(int i=0; i<length; ++i){
	tmp = tmp + seed2D[i].x0_ - seed2D[i].x1_ + seed2D[i].x2_ - seed2D[i].y0_ + seed2D[i].y1_ - seed2D[i].y2_;
  }
  return tmp;
}

float sumSeedArray(Seed2Dsa* seed2Dsa, const int length, const int arraysize){
  float tmp = 0;
  for(int i=0; i<length; ++i){
	for(int j=0; j<arraysize; ++j){
	  tmp = tmp + (seed2Dsa[i].x0_[j]) - (seed2Dsa[i].x1_[j]) + (seed2Dsa[i].x2_[j]) - (seed2Dsa[i].y0_[j]) + (seed2Dsa[i].y1_[j]) - (seed2Dsa[i].y2_[j]);
	}  
  }
  return tmp;
}


void addArrays(float *ina, float *inb, float *out, const int start, const int stop){
  #pragma simd
  for(int i = start; i < stop-start; ++i){
	out[start+i]=ina[start+i]*ina[start+i] + inb[start+i]*inb[start+i]-ina[start+i]/inb[start+i];
  }
}

void addArrays(float *ina, float *inb, float *inc, float *out, const int start, const int stop){
  #pragma simd
  for(int i = start; i < stop-start; ++i){
	out[start+i]=ina[start+i] + inb[start+i] + inc[start+1]*inc[start+i] - inc[start+i]/inb[start+i];
  }
}

void addArrays(float *ina, float *inb, float *inc, float *ind, float *out, const int start, const int stop){
  #pragma simd
  for(int i = start; i < stop-start; ++i){
	out[start+i]=ina[start+i] + inb[start+i] + inc[start+i] + ind[start+i]/inc[start+i] - ind[start+i];
  }
}

void addArrays(float *ina, float *inb, float *inc, float *ind, float *ine, float *out, const int start, const int stop){
  #pragma simd
  for(int i = start; i < stop-start; ++i){
	out[start+i]=ina[start+i] + inb[start+i] + inc[start+i] + ind[start+i]/ine[start+i] - ind[start+i];
  }
}

void addArrays(float *ina, float *inb, float *inc, float *ind, float *ine, float *inf, float *out, const int start, const int stop){
  #pragma simd
  for(int i = start; i < stop-start; ++i){
	out[start+i]=ina[start+i] + inb[start+i] + inc[start+i] + ind[start+i]/ine[start+i]- inf[start+i];
  }
}
