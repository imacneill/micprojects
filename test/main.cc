#include <iostream>

#include "main.h"


int main(){

  float *a, *b, *c;
  const unsigned int NENTRIES = 10000000;
  const int ALIGN = 64;
  newArray(a, NENTRIES, ALIGN, LINEAR);
  newArray(b, NENTRIES, ALIGN, LINEAR);
  newArray(c, NENTRIES, ALIGN, ZERO);
  
  //#pragma simd          //tells the compiler that the vectors are independent of each other, similar to #pragma ivdep
  //#pragma omp parallel for       //tells the compiler to run this in parallel on available threads
  //#pragma omp parallel for simd      //tells the compiler to run this in parallel on available threads and that the vectors are independent

  double starttime = get_wall_time();
  for(unsigned int i = 0; i<NENTRIES; ++i){
	c[i]=a[i]+b[i];
  }
  double endtime = get_wall_time();
  std::cout<<"Total time: "<<endtime-starttime<<" sec"<<std::endl<<sumArray(c, NENTRIES)<<std::endl;

  return 0;
}
