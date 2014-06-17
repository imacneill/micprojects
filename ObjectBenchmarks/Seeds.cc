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

float calcResid(const Seed2Dp& seed){
  float slope = ((*seed.y2_) - (*seed.y0_))/((*seed.x2_) - (*seed.x0_));
  float intercept = (*seed.y0_) - slope * (*seed.x0_);
  return std::abs( slope*(*seed.x0_) - (*seed.y0_) + intercept ) / sqrt( slope*slope + intercept*intercept );
}

float calcResid(const Seed2D64a& seed){
  float slope = ((*seed.y2_) - (*seed.y0_))/((*seed.x2_) - (*seed.x0_));
  float intercept = (*seed.y0_) - slope * (*seed.x0_);
  return std::abs( slope*(*seed.x0_) - (*seed.y0_) + intercept ) / sqrt( slope*slope + intercept*intercept );
}

float calcResid(const float x0, const float y0, const float x1, const float y1, const float x2, const float y2){
  float slope = (y2 - y0)/(x2 - x0);
  float intercept = y0 - slope * x0;
  return std::abs( slope*x0 - y0 + intercept ) / sqrt( slope*slope + intercept*intercept );
}


float calcResidIn(const Seed2D& seed){
  float slope = (seed.y2_ - seed.y0_)/(seed.x2_ - seed.x0_);
  float intercept = seed.y0_ - slope * seed.x0_;
  return std::abs( slope*seed.x0_ - seed.y0_ + intercept ) / sqrt( slope*slope + intercept*intercept );
}

float calcResidIn(const Seed2Dp& seed){
  float slope = ((*seed.y2_) - (*seed.y0_))/((*seed.x2_) - (*seed.x0_));
  float intercept = (*seed.y0_) - slope * (*seed.x0_);
  return std::abs( slope*(*seed.x0_) - (*seed.y0_) + intercept ) / sqrt( slope*slope + intercept*intercept );
}

float calcResidIn(const Seed2D64a& seed){
  float slope = ((*seed.y2_) - (*seed.y0_))/((*seed.x2_) - (*seed.x0_));
  float intercept = (*seed.y0_) - slope * (*seed.x0_);
  return std::abs( slope*(*seed.x0_) - (*seed.y0_) + intercept ) / sqrt( slope*slope + intercept*intercept );
}

float calcResidIn(const float x0, const float y0, const float x1, const float y1, const float x2, const float y2){
  float slope = (y2 - y0)/(x2 - x0);
  float intercept = y0 - slope * x0;
  return std::abs( slope*x0 - y0 + intercept ) / sqrt( slope*slope + intercept*intercept );
}



void generateSeeds(Seed2D* seed2D, Seed2Dp* seed2Dp, Seed2D64a* seed2D64a){
  float tmp0 = distribution(generator);
  float tmp1 = distribution(generator);
  float tmp2 = distribution(generator);
  seed2D->setX(tmp0, tmp1, tmp2);
  seed2Dp->setX(tmp0, tmp1, tmp2);
  seed2D64a->setX(tmp0, tmp1, tmp2);

  // seed2D = generateSeed2D(tmp0, tmp1, tmp2);
  // seed2Dp = generateSeed2Dp(tmp0, tmp1, tmp2);
  // seed2D64a = generateSeed2D64a(tmp0, tmp1, tmp2);
}

void generateSeeds(Seed2D* seed2D){
  float tmp0 = distribution(generator);
  float tmp1 = distribution(generator);
  float tmp2 = distribution(generator);
  seed2D->setX(tmp0, tmp1, tmp2);
}

// Seed2D* generateSeed2D(const float x0, const float x1, const float x2){
//   return new Seed2D(x0,x1,x2); 
// }
// Seed2Dp* generateSeed2Dp(const float x0, const float x1, const float x2){
//   return new Seed2Dp(x0,x1,x2); 
// }
// Seed2D64a* generateSeed2D64a(const float x0, const float x1, const float x2){
//   return new Seed2D64a(x0,x1,x2); 
// }



float sumSeedArray(Seed2D* seed2D, const unsigned int length){
  float tmp = 0;
  for(unsigned int i=0; i<length; ++i){
	tmp = tmp + seed2D[i].x0_ - seed2D[i].x1_ + seed2D[i].x2_ - seed2D[i].y0_ + seed2D[i].y1_ - seed2D[i].y2_;
  }
  return tmp;
}

float sumSeedArray(Seed2Dp* seed2Dp, const unsigned int length){
  float tmp = 0;
  for(unsigned int i=0; i<length; ++i){
	tmp = tmp + (*seed2Dp[i].x0_) - (*seed2Dp[i].x1_) + (*seed2Dp[i].x2_) - (*seed2Dp[i].y0_) + (*seed2Dp[i].y1_) - (*seed2Dp[i].y2_);
  }
  return tmp;
}

float sumSeedArray(Seed2D64a* seed2D64a, const unsigned int length){
  float tmp = 0;
  for(unsigned int i=0; i<length; ++i){
	tmp = tmp + (*seed2D64a[i].x0_) - (*seed2D64a[i].x1_) + (*seed2D64a[i].x2_) - (*seed2D64a[i].y0_) + (*seed2D64a[i].y1_) - (*seed2D64a[i].y2_);
  }
  return tmp;
}
