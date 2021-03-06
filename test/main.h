#ifndef MAIN_H
#define MAIN_H

#include<cstdlib>
#include<cstdio>
#include<omp.h>

#include <sys/time.h>
#include <time.h>

enum FillType{ZERO,LINEAR};

template<typename T>
void fillArray(T*& array, unsigned int length, FillType filltype= LINEAR){
#pragma omp parallel for simd 
  for( unsigned int i=0; i<length; ++i ){
	if( filltype == LINEAR ){
	  array[i] = (T) i;
	}else if(filltype == ZERO){
	  array[i] = 0;
	}
  }
}

template<typename T>
void newArray(T*& array, unsigned int length, int allignBytes = -1, FillType filltype = LINEAR){
  if(allignBytes < 0){
	array = (T*) malloc(sizeof(T)*length);
  }else{
	array = (T*) _mm_malloc(sizeof(T)*length, allignBytes);
  }
  fillArray(array, length, filltype);
}

template<typename T>
T sumArray(T* array, unsigned int length){
  T tmp = 0;
  for(unsigned int i=0; i< length; ++i){
	tmp+=array[i];
  }
  return tmp;
}


double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,0)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}


#endif
