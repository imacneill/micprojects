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
stopwatch::stopwatch(double (*watch)()):start_(0),stop_(0),elapsed_(0),started_(false),stopped_(false){watch_=watch;}
double stopwatch::start(){return start_;}
double stopwatch::stop(){return stop_;}
double stopwatch::elapsed(){return elapsed_;}
bool stopwatch::started(){return started_;}
bool stopwatch::stopped(){return stopped_;}
bool stopwatch::done(){return (started_ && stopped_);}
void stopwatch::timestart(){
  start_ = (*watch_)();
  started_ = true;
}
void stopwatch::timestop(){
  stop_ = (*watch_)(); 
  stopped_ = true;
  if(started_){elapsed_ = stop_ - start_;}
}

