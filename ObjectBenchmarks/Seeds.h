#ifndef SEEDS_H
#define SEEDS_H

#include <iostream>
#include<cstdlib>
#include<cstdio>

class Seed2D{
public:
  Seed2D() : y0_(0.), y1_(1.), y2_(2.), x0_(0), x1_(0), x2_(0) {}
  float x0_;
  float y0_;
  float x1_;
  float y1_;
  float x2_;
  float y2_;

  void setX(float x0, float x1, float x2){ x0_=x0; x1_=x1; x2_=x2; }
};

class Seed2Dp{
public:
  Seed2Dp(){
	y0_ = new float; (*y0_) = 0;
	y1_ = new float; (*y1_) = 1;
	y2_ = new float; (*y2_) = 2;
	x0_ = new float; (*x0_) = 0;
	x1_ = new float; (*x1_) = 0;
	x2_ = new float; (*x2_) = 0;
  }
  ~Seed2Dp(){
	if(x0_){delete x0_;}
	if(x1_){delete x1_;}
	if(x2_){delete x2_;}
	if(y0_){delete y0_;}
	if(x1_){delete y1_;}
	if(x2_){delete y2_;}
  }
  float *x0_;
  float *y0_;
  float *x1_;
  float *y1_;
  float *x2_;
  float *y2_;

  void setX(float x0, float x1, float x2){ (*x0_)=x0; (*x1_)=x1; (*x2_)=x2; }
};

class Seed2D64a{
public:
  Seed2D64a(){
	x0_=(float*) _mm_malloc(sizeof(float)*1, 64); *x0_=0;
	x1_=(float*) _mm_malloc(sizeof(float)*1, 64); *x1_=0;
	x2_=(float*) _mm_malloc(sizeof(float)*1, 64); *x2_=0;
	y0_=(float*) _mm_malloc(sizeof(float)*1, 64); *y0_=0;
	y1_=(float*) _mm_malloc(sizeof(float)*1, 64); *y1_=1;
	y2_=(float*) _mm_malloc(sizeof(float)*1, 64); *y2_=2;
  }
  ~Seed2D64a(){
	if(x0_){_mm_free(x0_);}
	if(x1_){_mm_free(x1_);}
	if(x2_){_mm_free(x2_);}
	if(y0_){_mm_free(y0_);}
	if(x1_){_mm_free(y1_);}
	if(x2_){_mm_free(y2_);}
  }
  float *x0_;
  float *y0_;
  float *x1_;
  float *y1_;
  float *x2_;
  float *y2_;

  void setX(float x0, float x1, float x2){ (*x0_)=x0; (*x1_)=x1; (*x2_)=x2; }
};

float calcResid(const Seed2D& seed);
float calcResid(const Seed2Dp& seed);
float calcResid(const Seed2D64a& seed);
float calcResid(const float x0, const float y0, const float x1, const float y1, const float x2, const float y2);

inline float calcResidIn(const Seed2D& seed);
inline float calcResidIn(const Seed2Dp& seed);
inline float calcResidIn(const Seed2D64a& seed);
inline float calcResidIn(const float x0, const float y0, const float x1, const float y1, const float x2, const float y2);

void generateSeeds(Seed2D* seed2D, Seed2Dp* seed2Dp, Seed2D64a* seed2D64a);
void generateSeeds(Seed2D* seed2D);
// Seed2D* generateSeed2D(const float x0, const float x1, const float x2);
// Seed2Dp* generateSeed2Dp(const float x0, const float x1, const float x2);
// Seed2D64a* generateSeed2D64a(const float x0, const float x1, const float x2);

float sumSeedArray(Seed2D* seed2D, const unsigned int length);
float sumSeedArray(Seed2Dp* seed2Dp, const unsigned int length);
float sumSeedArray(Seed2D64a* seed2D64a, const unsigned int length);


#endif
