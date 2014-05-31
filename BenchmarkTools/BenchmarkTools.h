#ifndef BENCHMARKTOOLS_H
#define BENCHMARKTOOLS_H

#include <string>
#include <ostream>

double get_wall_time();
double get_cpu_time();



class stopwatch{
public:
  stopwatch(const std::string& name, double (*watch)());
  double start() const;
  double stop() const;
  double elapsed() const;
  std::string elapsedstr() const;
  bool started() const;
  bool stopped() const;
  bool done() const;
  void timestart();
  void timestop();
  std::string name() const;
  std::string report() const;

private:
  double (*watch_)();
  const std::string name_;
  double start_;
  double stop_;
  double elapsed_;
  bool started_;
  bool stopped_;
};



std::ostream& operator<<(std::ostream& os, const stopwatch& obj);



#endif
