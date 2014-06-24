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
