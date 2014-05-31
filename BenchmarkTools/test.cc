#include <iostream>
#include "BenchmarkTools.h"

int main(){

  stopwatch wall_bm(&get_wall_time);
  stopwatch cpu_bm(&get_cpu_time);

  wall_bm.timestart();
  cpu_bm.timestart();
  
  double a=0;
  for(unsigned int i = 0; i < 100000; i++){
	for(unsigned int j = 0; j < 100000; j++){
	  a+=(double)i;
	  a-=(double)j;
	}
  }

  wall_bm.timestop();
  cpu_bm.timestop();

  if( wall_bm.done() && cpu_bm.done() ){
	std::cout<<"CPU:  "<<cpu_bm.elapsed()<<std::endl;
	std::cout<<"Wall: "<<wall_bm.elapsed()<<std::endl;
  }
  std::cout<<"a: "<<a<<std::endl; //must do something with a or it will compile away the loop
  //  a++;

  return 0;
}
