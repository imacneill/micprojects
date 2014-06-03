#ifndef BASICBENCHMARKS_H
#define BASICBENCHMARKS_H

#include<cstdlib>
#include<cstdio>
#include<omp.h>

enum FillType{ZERO,LINEAR};

// don't be a jackass and use this function for something other than an array with a number as a datatype
template<typename T>
void fillArray(T*& array, unsigned int length, FillType filltype= LINEAR){
  //#pragma ivdep
#pragma omp parallel for simd //private(i)
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



// // don't be a jackass and use this function for something other than an array with a number as a datatype
// //template<typename T>
// void fillArray(float*& array, unsigned int length, FillType filltype= LINEAR){
//   for( unsigned int i=0; i<length; ++i ){
// 	if( filltype == LINEAR ){
// 	  array[i] = (float) i;
// 	}else if(filltype == ZERO){
// 	  array[i] = 0;
// 	}
//   }
// }

// //template<typename T>
// void newArray(float*& array, unsigned int length, int allignBytes = -1, FillType filltype = LINEAR){
//   if(allignBytes < 0){
// 	array = (float*) malloc(sizeof(float)*length);
//   }else{
// 	array = (float*) _mm_malloc(sizeof(float)*length, allignBytes);
//   }
//   fillArray(array, length, filltype);
// }


#endif
