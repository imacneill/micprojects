#include <sstream>

#include <sys/time.h>
#include <time.h>

#include "BenchmarkTools.h"

double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,0)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}

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
