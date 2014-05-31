#ifndef BENCHMARKTOOLS_H
#define BENCHMARKTOOLS_H

double get_wall_time();
double get_cpu_time();



class stopwatch{
public:
  stopwatch(double (*watch)());
  double start();
  double stop();
  double elapsed();
  bool started();
  bool stopped();
  bool done();
  void timestart();
  void timestop();

private:
  double (*watch_)();
  double start_;
  double stop_;
  double elapsed_;
  bool started_;
  bool stopped_;
};









#endif
