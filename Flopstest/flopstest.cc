#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<string.h>
#include<sys/time.h>


//#include "BenchmarkTools.h"

double dtime(){
  double tseconds = 0.0;
  struct timeval mytime;
  gettimeofday(&mytime,(struct timezone*)0);
  tseconds = (double)(mytime.tv_sec+mytime.tv_usec*1.0e-6);
  return (tseconds);
}





int main(){

  const int SIZE = 1024*1024;
  const int INNER = 16384;
  const int OUTER = 1000000;

  const double ncores = 1;
  
  const double thread_factor = 2; // divide by this if unpaired thread, need 2 threads for full utilization
  #ifdef __ISMIC__
  const double vec_size = 512;
  const double frequency = 1.238094; // GHz
  #endif
  #ifdef __ISHOST__
  const double vec_size = 256;
  const double frequency = 2.000003; // GHz
  #endif
  const double lanes = vec_size / 32.;

  const double theory_gflops = ncores * frequency * lanes / thread_factor;


  // float a[1048576] __attribute__((align(64)));
  // float b[1048576] __attribute__((align(64)));
  // float c[1048576] __attribute__((align(64)));
  float *a = (float*)_mm_malloc(sizeof(float)*SIZE,64);
  float *b = (float*)_mm_malloc(sizeof(float)*SIZE,64);
  float *c = (float*)_mm_malloc(sizeof(float)*SIZE,64);

							 
  for(int i = 0; i<SIZE; ++i){
	a[i] = 0.1*i;
	b[i] = 0.2*i;
	c[i]=0;
  }
  
  //  stopwatch wall_bm1("Wall a+b", &get_wall_time);
  double start = dtime();
  int j, k;  
  //  wall_bm1.timestart();
  for(j=0; j<OUTER; j++){
	for(k=0; k<INNER; k++){
	  c[k] += a[k] + b[k];
	}
  }
  double stop = dtime();
  //wall_bm1.timestop();
  //  if( wall_bm1.done() ){
	//std::cout<<wall_bm1<<std::endl;
	//std::cout<<wall_bm1.elapsed()<<" vs "<<stop-start<<std::endl;
	double exp_gflops = 1.0e-9 * (double) OUTER * (double) INNER;
	std::cout<<exp_gflops<<" Gflops in "<<stop-start<<" seconds."<<" Theory max "<<theory_gflops<<" flops/s. "<<exp_gflops/(stop-start)/theory_gflops * 100.<<"%"<<std::endl;
	//  }


  float tmp = 0;
  float sign = -1.;
  for(int i=0; i<INNER; i++){
	sign*=-1;
  	tmp+=(c[i]*sign+(float)i*sign);
  }
   std::cout<<tmp<<std::endl;

  _mm_free(a);
  _mm_free(b);
  _mm_free(c);


  return 0;
}
