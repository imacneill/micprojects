#include <iostream>
#include <cmath>
#include <sys/time.h>
#include <time.h>
#include <omp.h>

#include "BasicBenchmarks.h"
#include "BenchmarkTools.h"


int main(){

  const unsigned int NENTRIES = 10000000;
  const int ALIGN = 64;
  // tests a+b, a-b, a*b, a/b, sqrt(a), a*a, a+b+c, a+b*c, (a+b)*c, a*scalar, a+b+c+d, a+b+c+d+e, a+b+c+d+e+f, a+b:c+d:e+f:ab+cd+ef

  float scalar = 5.;

  float *a, *b, *c, *d, *e, *f, 
	*tests, *test0, *test1, *test2, *test3, *test4, *test5, *test6, *test7, 
	*test8, *test9, *test10, *test11, *test12, *test13, *test14;

  newArray(a, NENTRIES, ALIGN, LINEAR);
  newArray(b, NENTRIES, ALIGN, LINEAR);
  newArray(c, NENTRIES, ALIGN, LINEAR);
  newArray(d, NENTRIES, ALIGN, LINEAR);
  newArray(e, NENTRIES, ALIGN, LINEAR);
  newArray(f, NENTRIES, ALIGN, LINEAR);

  newArray(tests, NENTRIES, ALIGN, ZERO);
  newArray(test0, NENTRIES, ALIGN, ZERO);
  newArray(test1, NENTRIES, ALIGN, ZERO);
  newArray(test2, NENTRIES, ALIGN, ZERO);
  newArray(test3, NENTRIES, ALIGN, ZERO);
  newArray(test4, NENTRIES, ALIGN, ZERO);
  newArray(test5, NENTRIES, ALIGN, ZERO);
  newArray(test6, NENTRIES, ALIGN, ZERO);
  newArray(test7, NENTRIES, ALIGN, ZERO);
  newArray(test8, NENTRIES, ALIGN, ZERO);
  newArray(test9, NENTRIES, ALIGN, ZERO);
  newArray(test10, NENTRIES, ALIGN, ZERO);
  newArray(test11, NENTRIES, ALIGN, ZERO);
  newArray(test12, NENTRIES, ALIGN, ZERO);
  newArray(test13, NENTRIES, ALIGN, ZERO);
  newArray(test14, NENTRIES, ALIGN, ZERO);


// #pragma omp parallel
// #pragma omp master




  //=scalar
  stopwatch wall_bms("Wall =scalar", &get_wall_time);
  stopwatch cpu_bms("CPU =scalar", &get_cpu_time);
  wall_bms.timestart();
  cpu_bms.timestart();
#ifdef _USEsimd_
#pragma simd
#endif
#ifdef _USEparallel_
#pragma omp parallel for
#endif
#ifdef _USEparallelsimd_
#pragma omp parallel for simd
#endif 
  for(unsigned int i = 0; i < NENTRIES; ++i){
	tests[i] = 10;
  }
  wall_bms.timestop();
  cpu_bms.timestop();
  if( wall_bms.done() && cpu_bms.done() ){
	std::cout<<wall_bms<<std::endl;
	//	std::cout<<cpu_bm0<<std::endl;
  }

  //=a
  stopwatch wall_bm0("Wall =a", &get_wall_time);
  stopwatch cpu_bm0("CPU =a", &get_cpu_time);
  wall_bm0.timestart();
  cpu_bm0.timestart();
#ifdef _USEsimd_
#pragma simd
#endif
#ifdef _USEparallel_
#pragma omp parallel for
#endif
#ifdef _USEparallelsimd_
#pragma omp parallel for simd
#endif 
  for(unsigned int i = 0; i < NENTRIES; ++i){
	test0[i] = a[i];
  }
  wall_bm0.timestop();
  cpu_bm0.timestop();
  if( wall_bm0.done() && cpu_bm0.done() ){
	std::cout<<wall_bm0<<std::endl;
	//	std::cout<<cpu_bm0<<std::endl;
  }


  //a+b
  stopwatch wall_bm1("Wall a+b", &get_wall_time);
  stopwatch cpu_bm1("CPU a+b", &get_cpu_time);
  wall_bm1.timestart();
  cpu_bm1.timestart();

#ifdef _USEsimd_
#pragma simd
#endif
#ifdef _USEparallel_
#pragma omp parallel for
#endif
#ifdef _USEparallelsimd_
#pragma omp parallel for simd
#endif 
  for(unsigned int i = 0; i < NENTRIES; ++i){
	test1[i] = a[i] + b[i];
  }
  wall_bm1.timestop();
  cpu_bm1.timestop();
  if( wall_bm1.done() && cpu_bm1.done() ){
	std::cout<<wall_bm1<<std::endl;
	//	std::cout<<cpu_bm1<<std::endl;
  }


  //a-b
  stopwatch wall_bm2("Wall a-b", &get_wall_time);
  stopwatch cpu_bm2("CPU a-b", &get_cpu_time);
  wall_bm2.timestart();
  cpu_bm2.timestart();
#ifdef _USEsimd_
#pragma simd
#endif
#ifdef _USEparallel_
#pragma omp parallel for
#endif
#ifdef _USEparallelsimd_
#pragma omp parallel for simd
#endif 
  for(unsigned int i = 0; i < NENTRIES; ++i){
	test2[i] = a[i] - b[i];
  }
  wall_bm1.timestop();
  cpu_bm1.timestop();
  if( wall_bm2.done() && cpu_bm2.done() ){
	std::cout<<wall_bm2<<std::endl;
	//	std::cout<<cpu_bm2<<std::endl;
  }

  //a*b
  stopwatch wall_bm3("Wall a*b", &get_wall_time);
  stopwatch cpu_bm3("CPU a*b", &get_cpu_time);
  wall_bm3.timestart();
  cpu_bm3.timestart();
#ifdef _USEsimd_
#pragma simd
#endif
#ifdef _USEparallel_
#pragma omp parallel for
#endif
#ifdef _USEparallelsimd_
#pragma omp parallel for simd
#endif 
  for(unsigned int i = 0; i < NENTRIES; ++i){
	test3[i] = a[i] * b[i];
  }
  wall_bm3.timestop();
  cpu_bm3.timestop();
  if( wall_bm3.done() && cpu_bm3.done() ){
	std::cout<<wall_bm3<<std::endl;
	//	std::cout<<cpu_bm3<<std::endl;
  }

  //a/b
  stopwatch wall_bm4("Wall a/b", &get_wall_time);
  stopwatch cpu_bm4("CPU a/b", &get_cpu_time);
  wall_bm4.timestart();
  cpu_bm4.timestart();
#ifdef _USEsimd_
#pragma simd
#endif
#ifdef _USEparallel_
#pragma omp parallel for
#endif
#ifdef _USEparallelsimd_
#pragma omp parallel for simd
#endif 
  for(unsigned int i = 0; i < NENTRIES; ++i){
	test4[i] = a[i] / b[i];
  }
  wall_bm4.timestop();
  cpu_bm4.timestop();
  if( wall_bm4.done() && cpu_bm4.done() ){
	std::cout<<wall_bm4<<std::endl;
	//	std::cout<<cpu_bm4<<std::endl;
  }

  //sqrt(a)
  stopwatch wall_bm5("Wall sqrt(a)", &get_wall_time);
  stopwatch cpu_bm5("CPU sqrt(a)", &get_cpu_time);
  wall_bm5.timestart();
  cpu_bm5.timestart();
#ifdef _USEsimd_
#pragma simd
#endif
#ifdef _USEparallel_
#pragma omp parallel for
#endif
#ifdef _USEparallelsimd_
#pragma omp parallel for simd
#endif 
  for(unsigned int i = 0; i < NENTRIES; ++i){
	test5[i] = sqrt(a[i]);
  }
  wall_bm5.timestop();
  cpu_bm5.timestop();
  if( wall_bm5.done() && cpu_bm5.done() ){
	std::cout<<wall_bm5<<std::endl;
	//	std::cout<<cpu_bm5<<std::endl;
  }

  //a*a
  stopwatch wall_bm6("Wall a*a", &get_wall_time);
  stopwatch cpu_bm6("CPU a*a", &get_cpu_time);
  wall_bm6.timestart();
  cpu_bm6.timestart();
#ifdef _USEsimd_
#pragma simd
#endif
#ifdef _USEparallel_
#pragma omp parallel for
#endif
#ifdef _USEparallelsimd_
#pragma omp parallel for simd
#endif 
  for(unsigned int i = 0; i < NENTRIES; ++i){
	test6[i] = a[i]*a[i];
  }
  wall_bm6.timestop();
  cpu_bm6.timestop();
  if( wall_bm6.done() && cpu_bm6.done() ){
	std::cout<<wall_bm6<<std::endl;
	//std::cout<<cpu_bm6<<std::endl;
  }

  //a+b+c
  stopwatch wall_bm7("Wall a+b+c", &get_wall_time);
  stopwatch cpu_bm7("CPU a+b+c", &get_cpu_time);
  wall_bm7.timestart();
  cpu_bm7.timestart();
#ifdef _USEsimd_
#pragma simd
#endif
#ifdef _USEparallel_
#pragma omp parallel for
#endif
#ifdef _USEparallelsimd_
#pragma omp parallel for simd
#endif 
  for(unsigned int i = 0; i < NENTRIES; ++i){
	test7[i] = a[i] + b[i] + c[i];
  }
  wall_bm7.timestop();
  cpu_bm7.timestop();
  if( wall_bm7.done() && cpu_bm7.done() ){
	std::cout<<wall_bm7<<std::endl;
	//	std::cout<<cpu_bm7<<std::endl;
  }

//   //a+b*c
//   stopwatch wall_bm8("Wall 8", &get_wall_time);
//   stopwatch cpu_bm8("CPU 8", &get_cpu_time);
//   wall_bm8.timestart();
//   cpu_bm8.timestart();
// #pragma omp parallel for simd //private(i)
//   for(unsigned int i = 0; i < NENTRIES; ++i){
// 	test8[i] = a[i] + b[i] * c[i];
//   }
//   wall_bm8.timestop();
//   cpu_bm8.timestop();
//   if( wall_bm8.done() && cpu_bm8.done() ){
// 	std::cout<<wall_bm8<<std::endl;
// 	std::cout<<cpu_bm8<<std::endl;
//   }

//   //(a+b)*c
//   stopwatch wall_bm9("Wall 9", &get_wall_time);
//   stopwatch cpu_bm9("CPU 9", &get_cpu_time);
//   wall_bm9.timestart();
//   cpu_bm9.timestart();
// #pragma omp parallel for simd //private(i)
//   for(unsigned int i = 0; i < NENTRIES; ++i){
// 	test9[i] = (a[i] + b[i]) * c[i];
//   }
//   wall_bm9.timestop();
//   cpu_bm9.timestop();
//   if( wall_bm9.done() && cpu_bm9.done() ){
// 	std::cout<<wall_bm9<<std::endl;
// 	std::cout<<cpu_bm9<<std::endl;
//   }

//   //a*scalar
//   stopwatch wall_bm10("Wall 10", &get_wall_time);
//   stopwatch cpu_bm10("CPU 10", &get_cpu_time);
//   wall_bm10.timestart();
//   cpu_bm10.timestart();
// #pragma omp parallel for simd //private(i)
//   for(unsigned int i = 0; i < NENTRIES; ++i){
// 	test10[i] = a[i] + b[i] * c[i];
//   }
//   wall_bm10.timestop();
//   cpu_bm10.timestop();
//   if( wall_bm10.done() && cpu_bm10.done() ){
// 	std::cout<<wall_bm10<<std::endl;
// 	std::cout<<cpu_bm10<<std::endl;
//   }
  
  //a+b+c+d
  stopwatch wall_bm11("Wall a+b+c+d", &get_wall_time);
  stopwatch cpu_bm11("CPU a+b+c+d", &get_cpu_time);
  wall_bm11.timestart();
  cpu_bm11.timestart();
#ifdef _USEsimd_
#pragma simd
#endif
#ifdef _USEparallel_
#pragma omp parallel for
#endif
#ifdef _USEparallelsimd_
#pragma omp parallel for simd
#endif 
  for(unsigned int i = 0; i < NENTRIES; ++i){
	test11[i] = a[i] + b[i] + c[i] + d[i];
  }
  wall_bm11.timestop();
  cpu_bm11.timestop();
  if( wall_bm11.done() && cpu_bm11.done() ){
	std::cout<<wall_bm11<<std::endl;
	//	std::cout<<cpu_bm11<<std::endl;
  }

  //a+b+c+d+e
  stopwatch wall_bm12("Wall a+b+c+d+e", &get_wall_time);
  stopwatch cpu_bm12("CPU a+b+c+d+e", &get_cpu_time);
  wall_bm12.timestart();
  cpu_bm12.timestart();
#ifdef _USEsimd_
#pragma simd
#endif
#ifdef _USEparallel_
#pragma omp parallel for
#endif
#ifdef _USEparallelsimd_
#pragma omp parallel for simd
#endif 
  for(unsigned int i = 0; i < NENTRIES; ++i){
	test12[i] = a[i] + b[i] + c[i] + d[i] + e[i];
  }
  wall_bm12.timestop();
  cpu_bm12.timestop();
  if( wall_bm12.done() && cpu_bm12.done() ){
	std::cout<<wall_bm12<<std::endl;
	//	std::cout<<cpu_bm12<<std::endl;
  }

  //a+b+c+d+e+f
  stopwatch wall_bm13("Wall a+b+c+d+e+f", &get_wall_time);
  stopwatch cpu_bm13("CPU a+b+c+d+e+f", &get_cpu_time);
  wall_bm13.timestart();
  cpu_bm13.timestart();
  //  #pragma ivdep
#ifdef _USEsimd_
#pragma simd
#endif
#ifdef _USEparallel_
#pragma omp parallel for
#endif
#ifdef _USEparallelsimd_
#pragma omp parallel for simd
#endif 
  for(unsigned int i = 0; i < NENTRIES; ++i){
	test13[i] = a[i] + b[i] + c[i] + d[i] + e[i] + f[i];
  }
  wall_bm13.timestop();
  cpu_bm13.timestop();
  if( wall_bm13.done() && cpu_bm13.done() ){
	std::cout<<wall_bm13<<std::endl;
	//	std::cout<<cpu_bm13<<std::endl;
  }


  //a+b:c+d:e+f
  stopwatch wall_bm14("Wall a+b:c+d:e+f", &get_wall_time);
  stopwatch cpu_bm14("CPU a+b:c+d:e+f", &get_cpu_time);
  wall_bm14.timestart();
  cpu_bm14.timestart();
  //  #pragma ivdep
#ifdef _USEsimd_
#pragma simd
#endif
#ifdef _USEparallel_
#pragma omp parallel for
#endif
#ifdef _USEparallelsimd_
#pragma omp parallel for simd
#endif 
  for(unsigned int i = 0; i < NENTRIES; ++i){
	float tmp1 = a[i] + b[i];
	float tmp2 = c[i] + d[i];
	float tmp3 = e[i] + f[i];
	test14[i] = tmp1+tmp2+tmp3;
  }
  wall_bm14.timestop();
  cpu_bm14.timestop();
  if( wall_bm14.done() && cpu_bm13.done() ){
	std::cout<<wall_bm14<<std::endl;
	//	std::cout<<cpu_bm13<<std::endl;
  }

  std::cout<<sumArray(test1, NENTRIES)<<sumArray(test2, NENTRIES)<<sumArray(test3, NENTRIES)<<sumArray(test4, NENTRIES)
			 <<sumArray(test5, NENTRIES)<<sumArray(test6, NENTRIES)<<sumArray(test7, NENTRIES)<<sumArray(test8, NENTRIES)
			 <<sumArray(test9, NENTRIES)<<sumArray(test10, NENTRIES)<<sumArray(test11, NENTRIES)<<sumArray(test12, NENTRIES)
			 <<sumArray(test13, NENTRIES)<<sumArray(test14, NENTRIES)<<std::endl;
 

  _mm_free(a);
  _mm_free(b);
  _mm_free(c);
  _mm_free(d);
  _mm_free(e);
  _mm_free(f);

  _mm_free(test1);
  _mm_free(test2);
  _mm_free(test3);
  _mm_free(test4);
  _mm_free(test5);
  _mm_free(test6);
  _mm_free(test7);
  _mm_free(test8);
  _mm_free(test9);
  _mm_free(test10);
  _mm_free(test11);
  _mm_free(test12);
  _mm_free(test13);
  _mm_free(test14);

  return 0;
}
