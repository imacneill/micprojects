#include <sstream>
#include <iostream>

#include <sys/time.h>
#include <time.h>

#include "BenchmarkTools.h"

double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,0)){
        //  Handle error
        return 0;
    }
    return (double)(time.tv_sec + time.tv_usec * .000001);
}
double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}

int calcFlops(Processor type){
  float frequency = 0; //GHz
  float ncores = 0;
  float nlanes = 0;

  if(type == HOST){
	frequency = 2.000003; //GHz
	ncores = 6;
	nlanes = 8; //256byte vector unit? with 32 byte floats
  }
  else if(type == MIC){
	frequency = 1.238094; //GHz
	ncores = 61;
	nlanes = 16; //512byte vector unit? with 32 byte floats
  }else{
	std::cout<<"Unrecognized processor type in countFlops()"<<std::endl;
  }

  return frequency * ncores * nlanes *2;
}

void addArrays(float *ina, float *inb, float *out, const int start, const int stop);

///// stopwatch class functions /////
stopwatch::stopwatch(const std::string& name, double (*watch)()):start_(0),stop_(0),elapsed_(0),started_(false),stopped_(false),name_(name){watch_=watch;}
double stopwatch::start() const {return start_;}
double stopwatch::stop() const {return stop_;}
double stopwatch::elapsed() const {return elapsed_;}
std::string stopwatch::elapsedstr() const {
  std::stringstream ss;
  ss<<elapsed_;
  std::string tmp;
  ss>>tmp;
  return tmp;
}
bool stopwatch::started() const {return started_;}
bool stopwatch::stopped() const {return stopped_;}
bool stopwatch::done() const {return (started_ && stopped_);}
void stopwatch::timestart(){
  start_ = (*watch_)();
  started_ = true;
  stopped_ = false;
}
void stopwatch::timestop(){
  stop_ = (*watch_)(); 
  if(started_ && !stopped_){
	elapsed_ += (stop_ - start_);
	stopped_ = true;
  }
}
void stopwatch::timereset(){
  start_=0;
  stop_=0;
  elapsed_=0;
  started_=false;
  stopped_=false;
}

std::string stopwatch::name() const {return name_;}

std::string stopwatch::report() const{
  std::string tmp=std::string("Benchmark ")+name_;
  if(done()){
	tmp+=std::string(": ")+elapsedstr()+std::string(" secs");
  }else if(started() && !stopped()){
	tmp+=std::string(" started but not stopped. ")+elapsedstr()+std::string(" secs so far.");
  }else{
	tmp+=std::string(" is not being used properly.");
  }
  return tmp;
}

std::ostream& operator<<(std::ostream& os, const stopwatch& obj){
  os<<obj.report();
  return os;
}
