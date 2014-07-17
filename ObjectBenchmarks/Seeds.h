#ifndef SEEDS_H
#define SEEDS_H

#include <iostream>
#include<cstdlib>
#include<cstdio>
#include <cmath>

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

//
//// Hard coded short arrays
// template it with # of entrie
template<int i>
class Seed2Dhsa{ //seed short array
public:
  Seed2Dhsa(){
	entries_=i;
  }
 
  ~Seed2Dhsa(){}

  int entries_;

  float x0_[i];
  float y0_[i];
  float x1_[i];
  float y1_[i];
  float x2_[i];
  float y2_[i];

  void setX(float *x0, float *x1, float *x2, int begin, int entries){ 
	entries_ = entries;
	#pragma simd
	for(int j=0; j < entries_; ++j){
	  x0_[j]=x0[begin+j]; 
	  x1_[j]=x1[begin+j];  
	  x2_[j]=x2[begin+j];
	  y0_[j]=0;
	  y1_[j]=1.;
	  y2_[j]=2.;
	}	
  }

  void print(){
	for(int j=0; j<entries_; ++j){
	  std::cout<<x0_[j]<<" ";
	}
	std::cout<<std::endl;
	for(int j=0; j<entries_; ++j){
	  std::cout<<x1_[j]<<" ";
	}
	std::cout<<std::endl;
	for(int j=0; j<entries_; ++j){
	  std::cout<<x2_[j]<<" ";
	}
	std::cout<<std::endl;
  }
};


class Seed2Dsa{ //seed short array
public:
  Seed2Dsa(){}
 
  ~Seed2Dsa(){
	if(x0_){_mm_free(x0_);}
	if(x1_){_mm_free(x1_);}
	if(x2_){_mm_free(x2_);}
	if(y0_){_mm_free(y0_);}
	if(x1_){_mm_free(y1_);}
	if(x2_){_mm_free(y2_);}
  }

  int entries_;

  float *x0_;
  float *y0_;
  float *x1_;
  float *y1_;
  float *x2_;
  float *y2_;

  void setX(float *x0, float *x1, float *x2, int begin, int entries){ 
	entries_ = entries;
	x0_=(float*) _mm_malloc(sizeof(float)*entries_, 64); 
	x1_=(float*) _mm_malloc(sizeof(float)*entries_, 64); 
	x2_=(float*) _mm_malloc(sizeof(float)*entries_, 64); 
	y0_=(float*) _mm_malloc(sizeof(float)*entries_, 64); 
	y1_=(float*) _mm_malloc(sizeof(float)*entries_, 64); 
	y2_=(float*) _mm_malloc(sizeof(float)*entries_, 64); 

	#pragma simd
	for(int i=0; i < entries_; ++i){
	  x0_[i]=x0[begin+i]; 
	  x1_[i]=x1[begin+i];  
	  x2_[i]=x2[begin+i];
	  y0_[i]=0;
	  y1_[i]=1.;
	  y2_[i]=2.;
	}	
  }

  void print(){
	for(int i=0; i<entries_; ++i){
	  std::cout<<x0_[i]<<" ";
	}
	std::cout<<std::endl;
	for(int i=0; i<entries_; ++i){
	  std::cout<<x1_[i]<<" ";
	}
	std::cout<<std::endl;
	for(int i=0; i<entries_; ++i){
	  std::cout<<x2_[i]<<" ";
	}
	std::cout<<std::endl;
  }
};


// class Seed2D64a{
// public:
//   Seed2D64a(){
// 	x0_=(float*) _mm_malloc(sizeof(float)*1, 64); *x0_=0;
// 	x1_=(float*) _mm_malloc(sizeof(float)*1, 64); *x1_=0;
// 	x2_=(float*) _mm_malloc(sizeof(float)*1, 64); *x2_=0;
// 	y0_=(float*) _mm_malloc(sizeof(float)*1, 64); *y0_=0;
// 	y1_=(float*) _mm_malloc(sizeof(float)*1, 64); *y1_=1;
// 	y2_=(float*) _mm_malloc(sizeof(float)*1, 64); *y2_=2;
//   }
//   ~Seed2D64a(){
// 	if(x0_){_mm_free(x0_);}
// 	if(x1_){_mm_free(x1_);}
// 	if(x2_){_mm_free(x2_);}
// 	if(y0_){_mm_free(y0_);}
// 	if(x1_){_mm_free(y1_);}
// 	if(x2_){_mm_free(y2_);}
//   }
//   float *x0_;
//   float *y0_;
//   float *x1_;
//   float *y1_;
//   float *x2_;
//   float *y2_;

//   void setX(float x0, float x1, float x2){ (*x0_)=x0; (*x1_)=x1; (*x2_)=x2; }
// };

// class Seed2Dp{
// public:
//   Seed2Dp(){
// 	y0_ = new float; (*y0_) = 0;
// 	y1_ = new float; (*y1_) = 1;
// 	y2_ = new float; (*y2_) = 2;
// 	x0_ = new float; (*x0_) = 0;
// 	x1_ = new float; (*x1_) = 0;
// 	x2_ = new float; (*x2_) = 0;
//   }
//   ~Seed2Dp(){
// 	if(x0_){delete x0_;}
// 	if(x1_){delete x1_;}
// 	if(x2_){delete x2_;}
// 	if(y0_){delete y0_;}
// 	if(x1_){delete y1_;}
// 	if(x2_){delete y2_;}
//   }
//   float *x0_;
//   float *y0_;
//   float *x1_;
//   float *y1_;
//   float *x2_;
//   float *y2_;

//   void setX(float x0, float x1, float x2){ (*x0_)=x0; (*x1_)=x1; (*x2_)=x2; }
// };

float calcResid(const Seed2D& seed);
//float calcResid(const Seed2Dp& seed);
//float calcResid(const Seed2D64a& seed);
void calcResid(const Seed2Dsa& seed, float *resid, const int begin, const int end);
//float calcResid(const Seed2Dsa& seed, const int i);
float calcResid(const float x0, const float y0, const float x1, const float y1, const float x2, const float y2);
template <int j> void calcResid(const Seed2Dhsa<j>& seed, float *resid, const int begin, const int end){
  #pragma simd
  for(int i = 0; i<(end-begin); ++i){
	float slope = ((seed.y2_[i]) - (seed.y0_[i]))/((seed.x2_[i]) - (seed.x0_[i]));
	float intercept = (seed.y0_[i]) - slope * (seed.x0_[i]);
	resid[i+begin] = std::abs( slope*(seed.x0_[i]) - (seed.y0_[i]) + intercept ) / sqrt( slope*slope + intercept*intercept );
  }
}

inline float calcResidIn(const Seed2D& seed){
  float slope = (seed.y2_ - seed.y0_)/(seed.x2_ - seed.x0_);
  float intercept = seed.y0_ - slope * seed.x0_;
  return std::abs( slope*seed.x0_ - seed.y0_ + intercept ) / sqrt( slope*slope + intercept*intercept );
}
//inline float calcResidIn(const Seed2Dp& seed);
//inline float calcResidIn(const Seed2D64a& seed);
// inline float calcResidIn(const Seed2Dsa& seed, const int i){
//   float slope = ((seed.y2_[i]) - (seed.y0_[i]))/((seed.x2_[i]) - (seed.x0_[i]));
//   float intercept = (seed.y0_[i]) - slope * (seed.x0_[i]);
//   return std::abs( slope*(seed.x0_[i]) - (seed.y0_[i]) + intercept ) / sqrt( slope*slope + intercept*intercept );
// }

inline void calcResidIn(const Seed2Dsa& seed, float *resid, const int begin, const int end){
  #pragma simd
  for(int i = 0; i<(end-begin); ++i){
	float slope = ((seed.y2_[i]) - (seed.y0_[i]))/((seed.x2_[i]) - (seed.x0_[i]));
	float intercept = (seed.y0_[i]) - slope * (seed.x0_[i]);
	resid[i+begin] = std::abs( slope*(seed.x0_[i]) - (seed.y0_[i]) + intercept ) / sqrt( slope*slope + intercept*intercept );
  }
}
// inline template<int j> float calcResidIn(const Seed2Dhsa<j>& seed, const int i){
//   float slope = ((seed.y2_[i]) - (seed.y0_[i]))/((seed.x2_[i]) - (seed.x0_[i]));
//   float intercept = (seed.y0_[i]) - slope * (seed.x0_[i]);
//   return std::abs( slope*(seed.x0_[i]) - (seed.y0_[i]) + intercept ) / sqrt( slope*slope + intercept*intercept );
// }
inline float calcResidIn(const float x0, const float y0, const float x1, const float y1, const float x2, const float y2){
  float slope = (y2 - y0)/(x2 - x0);
  float intercept = y0 - slope * x0;
  return std::abs( slope*x0 - y0 + intercept ) / sqrt( slope*slope + intercept*intercept );
}





void calcResidVectorized(Seed2D *seed, float *resid, const int begin, const int end);
void calcResidVectorized(const Seed2Dsa& seed, float *resid, const int begin, const int end);
void calcResidVectorized(float *x0, float *y0, float *x1, float *y1, float *x2, float *y2, float *resid, const int begin, const int end);
template <int j> void calcResidVectorized(Seed2Dhsa<j> *seed, float *resid, const int begin, const int end){
  int outer = (end-begin)/j; // will always work if end-begin is a power of 2
  for(int i = 0; i<outer; ++i){
#pragma simd
	for(int k = 0; k<j; ++k){ 
	  float slope = ((seed[i].y2_[k]) - (seed[i].y0_[k]))/((seed[i].x2_[k]) - (seed[i].x0_[k]));
	  float intercept = (seed[i].y0_[k]) - slope * (seed[i].x0_[k]);
	  //	  resid[i+begin] = std::abs( slope*(seed[i].x0_[k]) - (seed[i].y0_[k]) + intercept ) / sqrt( slope*slope + intercept*intercept );
	  resid[i+begin] = ( slope*(seed[i].x0_[k]) - (seed[i].y0_[k]) + intercept ) / ( slope*slope + intercept*intercept );
	}
  }
}


void calcResidVectorizedDouble(float *x0, float *y0, float *x1, float *y1, float *x2, float *y2, float *resid, const int inner, const int outer);
template <int j> void calcResidVectorizedDouble(Seed2Dhsa<j> *seed, float *resid, const int inner, const int outer){
  for(int i = 0; i<outer; ++i){
#pragma simd
	for(int k = 0; k<inner; ++k){ 
	  float slope = ((seed[i].y2_[k]) - (seed[i].y0_[k]))/((seed[i].x2_[k]) - (seed[i].x0_[k]));
	  float intercept = (seed[i].y0_[k]) - slope * (seed[i].x0_[k]);
	  //	  resid[i+begin] = std::abs( slope*(seed[i].x0_[k]) - (seed[i].y0_[k]) + intercept ) / sqrt( slope*slope + intercept*intercept );
	  resid[i*inner+k] = ( slope*(seed[i].x0_[k]) - (seed[i].y0_[k]) + intercept ) / ( slope*slope + intercept*intercept );
	}
  }
}






//void generateSeeds(Seed2D* seed2D, Seed2Dp* seed2Dp, Seed2D64a* seed2D64a);
void generateSeeds(Seed2D* seed2D);
// Seed2D* generateSeed2D(const float x0, const float x1, const float x2);
// Seed2Dp* generateSeed2Dp(const float x0, const float x1, const float x2);
// Seed2D64a* generateSeed2D64a(const float x0, const float x1, const float x2);

float sumSeedArray(Seed2D* seed2D, const int length);
//float sumSeedArray(Seed2Dp* seed2Dp, const int length);
//float sumSeedArray(Seed2D64a* seed2D64a, const int length);
float sumSeedArray(Seed2Dsa* seed2Dsa, const int length, const int arraysize);




void addArrays(float *ina, float *inb, float *out, const int start, const int stop);
void addArrays(float *ina, float *inb, float *inc, float *out, const int start, const int stop);
void addArrays(float *ina, float *inb, float *inc, float *ind, float *out, const int start, const int stop);
void addArrays(float *ina, float *inb, float *inc, float *ind, float *ine, float *out, const int start, const int stop);
void addArrays(float *ina, float *inb, float *inc, float *ind, float *ine, float *inf, float *out, const int start, const int stop);



#endif
