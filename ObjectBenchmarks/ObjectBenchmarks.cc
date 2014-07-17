#include <iostream>
#include <fstream>
#include <cmath>
#include <omp.h>
#include <string>
#include <sstream>

#include "Seeds.h"
#include "BenchmarkTools.h"
#include "ObjectBenchmarks.h"



int main(){
  
  //// make test arrays
  //// one test with an array for each coordinate
  //// one test with an array that holds a struct that holds each coordinate
  //// one test with an array that holds a struct that holds pointers to each coordinate that have been alligned with 64 bits

  const  int NENTRIES = 10000000;
  // useful pairs, 156250:16, 78125:128
  const  int NENTRIES_sa = 156250;
  const  int NENTRIES_sasize = 16;
  const  int ALLIGN = 64;
  float *x0 = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *y0 = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *x1 = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *y1 = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *x2 = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *y2 = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  
  Seed2D *seed2D = (Seed2D*) _mm_malloc(sizeof(Seed2D)*NENTRIES, ALLIGN);
  Seed2Dsa *seed2Dsa = (Seed2Dsa*) _mm_malloc(sizeof(Seed2Dsa)*NENTRIES_sa, ALLIGN);
  Seed2Dhsa<NENTRIES_sasize> *seed2Dhsa = (Seed2Dhsa<NENTRIES_sasize>*) _mm_malloc(sizeof(Seed2Dhsa<NENTRIES_sasize>)*NENTRIES_sa, ALLIGN);

  //  std::cout<<sizeof(Seed2D)<<" struct size"<<std::endl;
  fillSeeds(NENTRIES, x0, y0, x1, y1, x2, y2, seed2D);
  copySeedToSeedsa(NENTRIES, x0, x1, x2, NENTRIES_sa, NENTRIES_sasize, seed2Dsa);

  for(unsigned int i=0; i<NENTRIES_sa; ++i){
	if(NENTRIES > (i * NENTRIES_sasize)){
	  seed2Dhsa[i].setX(x0, x1, x2, i * NENTRIES_sasize, NENTRIES_sasize);
	}
  }


  

  //// make arrays to hold the residuals that are calculated
  float *residd = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *residdn = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *resid2Dd = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *resid2Ddn = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *resid2Dsad = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *resid2Dhsad = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);

  float *resid = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *resid2D = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *resid2Dsa = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *resid2Dhsa = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);

  //explicitly using inline functions for these
  float *residIn = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *resid2DIn = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *resid2DsaIn = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *resid2DhsaIn = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  
  





  const int TOTAL = 2097152 * 1000;
  float *residScan2_4 = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *residScan3_4 = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *residScan4_4 = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);

  float *residScan = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *resid2DScan = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  float *resid2DhsaScan = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);


  // std::ofstream out_2arrays_5calc;
  // out_2arrays_5calc.open("output/out_2arrays_5calc.txt");

  // std::ofstream out_3arrays_5calc;
  // out_3arrays_5calc.open("output/out_3arrays_5calc.txt");

  // std::ofstream out_4arrays_5calc;
  // out_4arrays_5calc.open("output/out_4arrays_5calc.txt");

  // std::ofstream resid_6arrays_12calc;
  // resid_6arrays_12calc.open("output/resid_6arrays_12calc.txt");

  // std::ofstream resid_s6floats_12calc;
  // resid_s6floats_12calc.open("output/resid_s6floats_12calc.txt");

  // std::ofstream resid_shsa_12calc;
  // resid_shsa_12calc.open("output/resid_shsa_12calc.txt");



  // //  for(int INNER = 8; INNER < (262144+1); INNER*=2){
  // for(int INNER = 8; INNER < (2097152+1); INNER*=2){

  // 	const int OUTER = TOTAL/INNER;


  // 	std::cout<<"Scan with reuse the same cache "<<INNER<<std::endl;
  // 	stopwatch wall_bm_a_2_4("Wall Array Reuse Cache", &get_wall_time);
  // 	wall_bm_a_2_4.timestart();
  // 	for(int i=0; i<OUTER; ++i){
  // 	  addArrays(x0, x1, residScan2_4, 0, INNER);
  // 	}
  // 	wall_bm_a_2_4.timestop();
  // 	if( wall_bm_a_2_4.done() ){ 
  // 	  std::cout<<wall_bm_a_2_4<<" GFlops = "<<5.*(float)OUTER*(float)INNER/1e9/wall_bm_a_2_4.elapsed()<<std::endl; 
  // 	  out_2arrays_5calc<<INNER<<" "<<5.*(float)OUTER*(float)INNER/1e9/wall_bm_a_2_4.elapsed()<<std::endl;
  // 	}

  // 	std::cout<<"Scan with reuse the same cache "<<INNER<<std::endl;
  // 	stopwatch wall_bm_a_3_4("Wall Array Reuse Cache", &get_wall_time);
  // 	wall_bm_a_3_4.timestart();
  // 	for(int i=0; i<OUTER; ++i){
  // 	  addArrays(x0, x1, x2, residScan3_4, 0, INNER);
  // 	}
  // 	wall_bm_a_3_4.timestop();
  // 	if( wall_bm_a_3_4.done() ){ 
  // 	  std::cout<<wall_bm_a_3_4<<" GFlops = "<<5.*(float)OUTER*(float)INNER/1e9/wall_bm_a_3_4.elapsed()<<std::endl; 
  // 	  out_3arrays_5calc<<INNER<<" "<<5.*(float)OUTER*(float)INNER/1e9/wall_bm_a_3_4.elapsed()<<std::endl;
  // 	}

  // 	std::cout<<"Scan with reuse the same cache "<<INNER<<std::endl;
  // 	stopwatch wall_bm_a_4_4("Wall Array Reuse Cache", &get_wall_time);
  // 	wall_bm_a_4_4.timestart();
  // 	for(int i=0; i<OUTER; ++i){
  // 	  addArrays(x0, x1, x2, y1, residScan4_4, 0, INNER);
  // 	}
  // 	wall_bm_a_4_4.timestop();
  // 	if( wall_bm_a_4_4.done() ){ 
  // 	  std::cout<<wall_bm_a_4_4<<" GFlops = "<<5.*(float)OUTER*(float)INNER/1e9/wall_bm_a_4_4.elapsed()<<std::endl; 
  // 	  out_4arrays_5calc<<INNER<<" "<<5.*(float)OUTER*(float)INNER/1e9/wall_bm_a_4_4.elapsed()<<std::endl;
  // 	}


  // 	std::cout<<"Scan with reuse the same cache "<<INNER<<std::endl;
  // 	stopwatch wall_bm_a_rv("Wall Array Reuse Cache, resid vectorized", &get_wall_time);
  // 	wall_bm_a_rv.timestart();
  // 	for(int i=0; i<OUTER; ++i){
  // 	  calcResidVectorized(x0, y0, x1, y1, x2, y2, residScan, 0, INNER);
  // 	}
  // 	wall_bm_a_rv.timestop();
  // 	if( wall_bm_a_rv.done() ){ 
  // 	  std::cout<<wall_bm_a_rv<<" GFlops = "<<12.*(float)OUTER*(float)INNER/1e9/wall_bm_a_rv.elapsed()<<std::endl; 
  // 	  resid_6arrays_12calc<<INNER<<" "<<12.*(float)OUTER*(float)INNER/1e9/wall_bm_a_rv.elapsed()<<std::endl;
  // 	}

  // 	std::cout<<"Scan with reuse the same cache "<<INNER<<std::endl;
  // 	stopwatch wall_bm_s_rv("Wall struct Reuse Cache, resid vectorized", &get_wall_time);
  // 	wall_bm_s_rv.timestart();
  // 	for(int i=0; i<OUTER; ++i){
  // 	  calcResidVectorized(seed2D, resid2DScan, 0, INNER);
  // 	}
  // 	wall_bm_s_rv.timestop();
  // 	if( wall_bm_s_rv.done() ){ 
  // 	  std::cout<<wall_bm_s_rv<<" GFlops = "<<12.*(float)OUTER*(float)INNER/1e9/wall_bm_s_rv.elapsed()<<std::endl; 
  // 	  resid_s6floats_12calc<<INNER<<" "<<12.*(float)OUTER*(float)INNER/1e9/wall_bm_s_rv.elapsed()<<std::endl;
  // 	}

  // 	std::cout<<"Scan with reuse the same cache "<<INNER<<std::endl;
  // 	stopwatch wall_bm_shsa_rv("Wall struct hsa Reuse Cache, resid vectorized", &get_wall_time);
  // 	wall_bm_shsa_rv.timestart();
  // 	for(int i=0; i<OUTER; ++i){
  // 	  calcResidVectorized(seed2Dhsa, resid2DhsaScan, 0, INNER);
  // 	}
  // 	wall_bm_shsa_rv.timestop();
  // 	if( wall_bm_shsa_rv.done() ){ 
  // 	  std::cout<<wall_bm_shsa_rv<<" GFlops = "<<12.*(float)OUTER*(float)INNER/1e9/wall_bm_shsa_rv.elapsed()<<std::endl; 
  // 	  resid_shsa_12calc<<INNER<<" "<<12.*(float)OUTER*(float)INNER/1e9/wall_bm_shsa_rv.elapsed()<<std::endl;
  // 	}


  // 	// std::cout<<"Scan with reuse the same cache "<<INNER<<std::endl;
  // 	// stopwatch wall_bm_a_c("Wall Array Reuse Cache", &get_wall_time);
  // 	// wall_bm_a_c.timestart();
  // 	// for(int i=0; i<TOTAL; ++i){
  // 	//   calcResidVectorized(x0, y0, x1, y1, x2, y2, residScan, 0, INNER);
  // 	// }
  // 	// wall_bm_a_c.timestop();
  // 	// if( wall_bm_a_c.done() ){ std::cout<<wall_bm_a_c<<" GFlops = "<<14*(float)TOTAL*(float)INNER/1e9/wall_bm_a_c.elapsed()<<std::endl; }

  // }

  // out_2arrays_5calc.close();
  // out_3arrays_5calc.close();
  // out_4arrays_5calc.close();
  // resid_6arrays_12calc.close();
  // resid_s6floats_12calc.close();
  // resid_shsa_12calc.close();



  std::ofstream resid_dl_6arrays_12calc;
  //  resid_dl_6arrays_12calc.open("output/resid_dl_6arrays_12calc.txt");

  std::ofstream resid_dl_shsa_12calc;
  //resid_dl_shsa_12calc.open("output/resid_dl_shsa_12calc.txt");

  for(int ENTRIES = 8; ENTRIES < (2097152+1); ENTRIES*=8){
   	std::cout<<ENTRIES<<std::endl;
	//  const int ENTRIES = 2097152;


	std::stringstream ss;
	ss<<'e'<<ENTRIES;;
	std::string details;
	ss>>details;
	std::string filea("output/resid_dl_6arrays_");
	filea+=details;
	filea+=std::string(".txt");
	std::cout<<filea<<std::endl;
	resid_dl_6arrays_12calc.open(filea.c_str());

	std::string fileshsa("output/resid_dl_shsa_");
	fileshsa+=details;
	fileshsa+=std::string(".txt");
	std::cout<<fileshsa<<std::endl;
	resid_dl_shsa_12calc.open(fileshsa.c_str());


  for(int VECTOR = 8; VECTOR < (ENTRIES+1); VECTOR*=2){
	const int INNER = ENTRIES / VECTOR;
	const int OUTER = 1000;



	
	std::cout<<"Scan with reuse the same cache , double loop "<<VECTOR<<std::endl;
	stopwatch wall_bm_a_d("Wall Array Reuse Cache, resid vectorized, double loop", &get_wall_time);
	wall_bm_a_d.timestart();
	for(int i=0; i<OUTER; ++i){
	  //	  calcResidVectorized(x0, y0, x1, y1, x2, y2, residScan, 0, INNER);
	  calcResidVectorizedDouble(x0, y0, x1, y1, x2, y2, residScan, VECTOR, INNER);
	}
	wall_bm_a_d.timestop();
	if( wall_bm_a_d.done() ){ 
	  std::cout<<wall_bm_a_d<<" GFlops = "<<12.*(float)OUTER*(float)VECTOR*(float)INNER/1e9/wall_bm_a_d.elapsed()<<std::endl; 
	  resid_dl_6arrays_12calc<<INNER<<" "<<12.*(float)OUTER*(float)VECTOR*(float)INNER/1e9/wall_bm_a_d.elapsed()<<std::endl;
	}

	std::cout<<"Scan with reuse the same cache , double loop "<<VECTOR<<std::endl;
	stopwatch wall_bm_sasa_d("Wall Array Reuse Cache, resid vectorized, double loop", &get_wall_time);
	wall_bm_sasa_d.timestart();
	for(int i=0; i<OUTER; ++i){
	  //	  calcResidVectorized(x0, y0, x1, y1, x2, y2, residScan, 0, INNER);
	  calcResidVectorizedDouble(x0, y0, x1, y1, x2, y2, residScan, VECTOR, INNER);
	}
	wall_bm_sasa_d.timestop();
	if( wall_bm_sasa_d.done() ){ 
	  std::cout<<wall_bm_sasa_d<<" GFlops = "<<12.*(float)OUTER*(float)VECTOR*(float)INNER/1e9/wall_bm_sasa_d.elapsed()<<std::endl; 
	  resid_dl_shsa_12calc<<INNER<<" "<<12.*(float)OUTER*(float)VECTOR*(float)INNER/1e9/wall_bm_sasa_d.elapsed()<<std::endl;
	}
	
  }

  resid_dl_6arrays_12calc.close();
  resid_dl_shsa_12calc.close();
  }
  //  resid_dl_6arrays_12calc.close();
  //  resid_dl_shsa_12calc.close();

  _mm_free(residScan);
  _mm_free(resid2DScan);
  _mm_free(resid2DhsaScan);
  _mm_free(residScan2_4);
  _mm_free(residScan3_4);
  _mm_free(residScan4_4);



  // //loop over the same data for max cache efficiency
  // const int REPETITIONS = 78125*10;
  // const int DATALENGTH = 128;
  // float *residRepeat = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  // float *resid2DRepeat = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  // float *resid2DsaRepeat = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);
  // float *resid2DhsaRepeat = (float*) _mm_malloc(sizeof(float)*NENTRIES, ALLIGN);



  // std::cout<<"reuse the same cache"<<std::endl;
  // stopwatch wall_bm_a_c("Wall Array Reuse Cache", &get_wall_time);
  // wall_bm_a_c.timestart();
  // for(int i=0; i<REPETITIONS; ++i){
  // 	calcResidVectorized(x0, y0, x1, y1, x2, y2, residRepeat, 0, 128);
  // }
  // wall_bm_a_c.timestop();
  // if( wall_bm_a_c.done() ){ std::cout<<wall_bm_a_c<<std::endl; }

  // stopwatch wall_bm_s_c("Wall Struct Reuse Cache", &get_wall_time);
  // wall_bm_s_c.timestart();
  // for(int i=0; i<REPETITIONS; ++i){
  // 	calcResidVectorized(seed2D, resid2DRepeat, 0, 128);
  // }
  // wall_bm_s_c.timestop();
  // if( wall_bm_s_c.done() ){ std::cout<<wall_bm_s_c<<std::endl; }

  // stopwatch wall_bm_sp_c("Wall Struct Short Array Reuse Cache", &get_wall_time);
  // wall_bm_sp_c.timestart();
  // for(int i=0; i<REPETITIONS; ++i){
  // 	calcResidVectorized(seed2Dsa[0], resid2DRepeat, 0, 128);
  // }
  // wall_bm_sp_c.timestop();
  // if( wall_bm_sp_c.done() ){ std::cout<<wall_bm_sp_c<<std::endl; }

  // stopwatch wall_bm_sh_c("Wall Struct Static Short Array Reuse Cache", &get_wall_time);
  // wall_bm_sh_c.timestart();
  // for(int i=0; i<REPETITIONS; ++i){
  // 	calcResidVectorized(seed2Dhsa[0], resid2DhsaRepeat, 0, 128);
  // }
  // wall_bm_sh_c.timestop();
  // if( wall_bm_sh_c.done() ){ std::cout<<wall_bm_sh_c<<std::endl; }

  // std::cout<<std::endl<<std::endl;








//   //// make the loops
//   /// no function
//   stopwatch wall_bmad("Wall Array Direct", &get_wall_time);
//   wall_bmad.timestart();
// #pragma simd
//   for(int i = 0; i < NENTRIES; ++i){
// 	float slope = (y2[i] - y0[i])/(x2[i] - x0[i]);
// 	float intercept = y0[i] - slope * x0[i];
// 	residd[i] = std::abs( slope*x0[i] - y0[i] + intercept ) / sqrt( slope*slope + intercept*intercept ); 
//   }
//   wall_bmad.timestop();
//   if( wall_bmad.done() ){ std::cout<<wall_bmad<<std::endl; }
//   //nested
//   stopwatch wall_bmadn("Wall Array Direct Nested", &get_wall_time);
//   wall_bmadn.timestart();
//   for(int i = 0; i < NENTRIES_sa; ++i){
// #pragma simd
// 	for(int j = 0; j < NENTRIES_sasize; ++j){
// 	  float slope = (y2[i*NENTRIES_sasize+j] - y0[i*NENTRIES_sasize+j])/(x2[i*NENTRIES_sasize+j] - x0[i*NENTRIES_sasize+j]);
// 	  float intercept = y0[j] - slope * x0[j];
// 	  residdn[i*NENTRIES_sasize+j] = std::abs( slope*x0[i*NENTRIES_sasize+j] - y0[i*NENTRIES_sasize+j] + intercept ) / sqrt( slope*slope + intercept*intercept ); 
// 	}
//   }
//   wall_bmadn.timestop();
//   if( wall_bmadn.done() ){ std::cout<<wall_bmadn<<std::endl; }
  
//   //struct
//   stopwatch wall_bmsd("Wall Struct Direct", &get_wall_time);
//   wall_bmsd.timestart();
// #pragma simd
//   for(int i = 0; i < NENTRIES; ++i){
// 	float slope = (seed2D[i].y2_ - seed2D[i].y0_)/(seed2D[i].x2_ - seed2D[i].x0_);
// 	float intercept = seed2D[i].y0_ - slope * seed2D[i].x0_;
// 	 resid2Dd[i] = std::abs( slope*seed2D[i].x0_ - seed2D[i].y0_ + intercept ) / sqrt( slope*slope + intercept*intercept );
//   }
//   wall_bmsd.timestop();
//   if( wall_bmsd.done() ){ std::cout<<wall_bmsd<<std::endl; }
//   //nested
//   stopwatch wall_bmsdn("Wall Struct Direct Nested", &get_wall_time);
//   wall_bmsdn.timestart();
//   for(int i = 0; i < NENTRIES_sa; ++i){
// #pragma simd
// 	for(int j = 0; j < NENTRIES_sasize; ++j){
// 	  float slope = (seed2D[i*NENTRIES_sasize+j].y2_ - seed2D[i*NENTRIES_sasize+j].y0_)/(seed2D[i*NENTRIES_sasize+j].x2_ - seed2D[i*NENTRIES_sasize+j].x0_);
// 	  float intercept = seed2D[i*NENTRIES_sasize+j].y0_ - slope * seed2D[i*NENTRIES_sasize+j].x0_;
// 	  resid2Ddn[i*NENTRIES_sasize+j] = std::abs( slope*seed2D[i*NENTRIES_sasize+j].x0_ - seed2D[i*NENTRIES_sasize+j].y0_ + intercept ) / sqrt( slope*slope + intercept*intercept );
// 	}
//   }
//   wall_bmsdn.timestop();
//   if( wall_bmsdn.done() ){ std::cout<<wall_bmsdn<<std::endl; }


//   //struct short array
//   stopwatch wall_bmsad("Wall Struct Short Array Direct", &get_wall_time);
//   wall_bmsad.timestart();
//   for(int i = 0; i < NENTRIES_sa; ++i){
// #pragma simd
// 	for(int j = 0; j < NENTRIES_sasize; ++j){
// 	  float slope = (seed2Dsa[i].y2_[j] - seed2Dsa[i].y0_[j])/(seed2Dsa[i].x2_[j] - seed2Dsa[i].x0_[j]);
// 	  float intercept = seed2Dsa[i].y0_[j] - slope * seed2Dsa[i].x0_[j];
// 	  resid2Dsad[i*NENTRIES_sasize+j] = std::abs( slope*seed2Dsa[i].x0_[j] - seed2Dsa[i].y0_[j] + intercept ) / sqrt( slope*slope + intercept*intercept );
// 	}
//   }
//   wall_bmsad.timestop();
//   if( wall_bmsad.done() ){ std::cout<<wall_bmsad<<std::endl; }


//   //struct short array hard coded
//   stopwatch wall_bmhsad("Wall Struct Hand Short Array Direct", &get_wall_time);
//   wall_bmhsad.timestart();
//   for(int i = 0; i < NENTRIES_sa; ++i){
// #pragma simd
// 	for(int j = 0; j < NENTRIES_sasize; ++j){
// 	  float slope = (seed2Dhsa[i].y2_[j] - seed2Dhsa[i].y0_[j])/(seed2Dhsa[i].x2_[j] - seed2Dhsa[i].x0_[j]);
// 	  float intercept = seed2Dhsa[i].y0_[j] - slope * seed2Dhsa[i].x0_[j];
// 	  resid2Dhsad[i*NENTRIES_sasize+j] = std::abs( slope*seed2Dhsa[i].x0_[j] - seed2Dhsa[i].y0_[j] + intercept ) / sqrt( slope*slope + intercept*intercept );
// 	}
//   }
//   wall_bmhsad.timestop();
//   if( wall_bmhsad.done() ){ std::cout<<wall_bmhsad<<std::endl; }


//   std::cout<<std::endl;



//   ///no explicit inline
//   //array
//   stopwatch wall_bma("Wall Array", &get_wall_time);
//   wall_bma.timestart();
// #pragma simd
//   //  #pragma omp parallel for simd
//   for(int i = 0; i < NENTRIES; ++i){
// 	resid[i] = calcResid(x0[i], y0[i], x1[i], y1[i], x2[i], y2[i]);
//   }
//   wall_bma.timestop();
//   if( wall_bma.done() ){ std::cout<<wall_bma<<std::endl; }
  
//   //struct
//   stopwatch wall_bms("Wall Struct", &get_wall_time);
//   wall_bms.timestart();
// #pragma simd
//   //#pragma omp parallel for simd
//   for(int i = 0; i < NENTRIES; ++i){
// 	resid2D[i] = calcResid(seed2D[i]);
//   }
//   wall_bms.timestop();
//   if( wall_bms.done() ){ std::cout<<wall_bms<<std::endl; }

//   //struct short array
//   stopwatch wall_bmsa("Wall Struct Short Array", &get_wall_time);
//   wall_bmsa.timestart();
//   //#pragma simd
//   //#pragma omp parallel for simd
//   for(int i = 0; i < NENTRIES_sa; ++i){
// 	//#pragma simd
// 	//	for(int j = 0; j < NENTRIES_sasize; ++j){
// 	//  resid2Dsad[i*NENTRIES_sasize+j] = calcResid(seed2Dsa[i], j);
// 	calcResid(seed2Dsa[i], resid2Dsa, i*NENTRIES_sasize, (i+1)*NENTRIES_sasize);

// 	//}
//   }
//   wall_bmsa.timestop();
//   if( wall_bmsa.done() ){ std::cout<<wall_bmsa<<std::endl; }

//   //struct short array hard coded
//   stopwatch wall_bmhsa("Wall Struct Short Array Hand", &get_wall_time);
//   wall_bmhsa.timestart();
//   //#pragma simd
//   //#pragma omp parallel for simd
//   for(int i = 0; i < NENTRIES_sa; ++i){
// 	//#pragma simd
// 	//	for(int j = 0; j < NENTRIES_sasize; ++j){
// 	//  resid2Dsad[i*NENTRIES_sasize+j] = calcResid(seed2Dsa[i], j);
// 	calcResid(seed2Dhsa[i], resid2Dhsa, i*NENTRIES_sasize, (i+1)*NENTRIES_sasize);

// 	//}
//   }
//   wall_bmhsa.timestop();
//   if( wall_bmhsa.done() ){ std::cout<<wall_bmhsa<<std::endl; }
  


//   std::cout<<std::endl;



//   /// explicit inline
//   //array
//   stopwatch wall_bmai("Wall Array Inline", &get_wall_time);
//   wall_bmai.timestart();
// #pragma simd
//   //#pragma omp parallel for simd
//   for(int i = 0; i < NENTRIES; ++i){
// 	residIn[i] = calcResidIn(x0[i], y0[i], x1[i], y1[i], x2[i], y2[i]);
//   }
//   wall_bmai.timestop();
//   if( wall_bmai.done() ){ std::cout<<wall_bmai<<std::endl; }
  
//   //struct
//   stopwatch wall_bmsi("Wall Struct Inline", &get_wall_time);
//   wall_bmsi.timestart();
// #pragma simd
//   //#pragma omp parallel for simd
//   for(int i = 0; i < NENTRIES; ++i){
// 	resid2DIn[i] = calcResidIn(seed2D[i]);
//   }
//   wall_bmsi.timestop();
//   if( wall_bmsi.done() ){ std::cout<<wall_bmsi<<std::endl; }

//   //struct short array
//   stopwatch wall_bmsai("Wall Struct Short Array Inline", &get_wall_time);
//   wall_bmsai.timestart();
//   for(int i = 0; i < NENTRIES_sa; ++i){
// 	//#pragma simd
// 	//	for(int j = 0; j < NENTRIES_sasize; ++j){
// 	//	  resid2DsaIn[i*NENTRIES_sasize+j] = calcResidIn(seed2Dsa[i], j);
// 	calcResidIn(seed2Dsa[i], resid2DsaIn, i*NENTRIES_sasize, (i+1)*NENTRIES_sasize);
// 	  //	}
//   }
//   wall_bmsai.timestop();
//   if( wall_bmsai.done() ){ std::cout<<wall_bmsai<<std::endl; }
  
  
  int pass = 0, fail = 0;

  float tmpsum =0;
  for(int i=0; i<NENTRIES; i++){
	if(
	   abs( residd[i]) < 1. &&
	   abs( residdn[i]) < 1. &&
	   abs( resid2Dd[i]) < 1. &&
	   abs( resid2Dsad[i]) < 1. &&
	   abs( resid2Dhsad[i]) < 1. &&
	   abs( resid[i]) < 1. &&
	   abs( resid2D[i]) < 1. &&
	   abs( resid2Dsa[i]) < 1. &&
	   abs( resid2Dhsa[i]) < 1. &&
	   abs( residIn[i]) < 1. &&
	   abs( resid2DIn[i]) < 1. &&
	   abs( resid2DsaIn[i]) < 1. &&
	   abs( resid2DhsaIn[i]) < 1.
	   ){
	  pass++;
	}else{
	  fail++;
	}
	
	tmpsum=tmpsum
	  + residd[i]
	  - residdn[i]
	  + resid2Dd[i]
	  - resid2Dsad[i]
	  + resid2Dhsad[i]
	  - resid[i]
	  + resid2D[i]
	  - resid2Dsa[i]
	  + resid2Dhsa[i]
	  - residIn[i]
	  + resid2DIn[i]
	  - resid2DsaIn[i]
	  + resid2DhsaIn[i];
  }
  
  std::cout<<std::endl<<std::endl<<std::endl<<"Garbage: "
		   <<sumArray(residd, NENTRIES)
		   <<" "<<sumArray(residdn, NENTRIES)
		   <<" "<<sumArray(resid2Dd, NENTRIES)
		   <<" "<<sumArray(resid2Dsad, NENTRIES)
		   <<" "<<sumArray(resid2Dhsad, NENTRIES)
		   <<" "<<sumArray(resid, NENTRIES)
		   <<" "<<sumArray(resid2D, NENTRIES)
		   <<" "<<sumArray(resid2Dsa, NENTRIES)
		   <<" "<<sumArray(resid2Dhsa, NENTRIES)
		   <<" "<<sumArray(residIn, NENTRIES)
		   <<" "<<sumArray(resid2DIn, NENTRIES)
		   <<" "<<sumArray(resid2DsaIn, NENTRIES)
		   <<" "<<sumArray(resid2DhsaIn, NENTRIES)
		   <<std::endl;
  
  
  
  _mm_free(x0);
  _mm_free(y0);
  _mm_free(x1);
  _mm_free(y1);
  _mm_free(x2);
  _mm_free(y2);
  _mm_free(seed2D);
  _mm_free(seed2Dsa);
  _mm_free(seed2Dhsa);


  _mm_free(residd);
  _mm_free(resid2Dd);
  _mm_free(resid2Dsad);
  _mm_free(resid2Dhsad);

  _mm_free(resid);
  _mm_free(resid2D);
  _mm_free(resid2Dsa);
  _mm_free(resid2Dhsa);

  _mm_free(residIn);
  _mm_free(resid2DIn);
  _mm_free(resid2DsaIn);
  _mm_free(resid2DhsaIn);


  
  return 0;
}






void fillSeeds(const unsigned int NENTRIES, float *x0, float *y0, float *x1, float *y1, float *x2, float *y2, 
			   Seed2D *seed2D){

  #pragma omp parallel for simd
  for(unsigned int i=0; i<NENTRIES; i++){
	generateSeeds((seed2D+i));
	x0[i]=seed2D[i].x0_;
	y0[i]=seed2D[i].y0_;
	x1[i]=seed2D[i].x1_;
	y1[i]=seed2D[i].y1_;
	x2[i]=seed2D[i].x2_;
	y2[i]=seed2D[i].y2_;
  }
}

void copySeedToSeedsa(const unsigned int NENTRIES, float *x0, float *x1, float *x2, const unsigned int NENTRIES_sa, const unsigned int NENTRIES_sasize, Seed2Dsa *seed2Dsa){
  for(unsigned int i=0; i<NENTRIES_sa; ++i){
	if(NENTRIES > (i * NENTRIES_sasize)){
	  seed2Dsa[i].setX(x0, x1, x2, i * NENTRIES_sasize, NENTRIES_sasize);
	}
  }
}


