#include <iostream>
#include "BenchmarkTools.h"

int main(){

  stopwatch wall_bm("Wall", &get_wall_time);
  stopwatch cpu_bm("CPU", &get_cpu_time);

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
	std::cout<<wall_bm<<std::endl;
	std::cout<<cpu_bm<<std::endl;
  }
  std::cout<<"a: "<<a<<std::endl; //must do something with a or it will compile away the loop
  //  a++;

  wall_bm.timestart();
  cpu_bm.timestart();
  
  double b=0;
  for(unsigned int i = 0; i < 100000; i++){
	for(unsigned int j = 0; j < 100000; j++){
	  b+=(double)i;
	  b-=(double)j;
	}
  }

  wall_bm.timestop();
  cpu_bm.timestop();

  if( wall_bm.done() && cpu_bm.done() ){
	std::cout<<wall_bm<<std::endl;
	std::cout<<cpu_bm<<std::endl;
  }
  std::cout<<"a+b: "<<b<<std::endl; //must do something with a or it will compile away the loop


  return 0;
}
