#ifndef SPHSOL_OUTPUT_H_
#define SPHSOL_OUTPUT_H_

#include"global.h"
#include"particle.h"


//--------------------------------------------------------------------------------------------------
// Output setup
//--------------------------------------------------------------------------------------------------
class Output
{
public:
  Output() {};
  ~Output() {};

  void Initialize(double total_time, int total_step, const std::vector<int>& output_flag);
  void AddOutputTime(const std::vector<double>& add_time);

  void PrintTime();

  void PrintParticles(const Particle* particles, int particles_num);
  void PrintFragments(const Particle* particles, int particles_num, std::vector<int>* fragments);
  void PrintNeighbors(std::vector<Neighbor>* neighbors, int particles_num);

  void PrintLogger(const char* format, ...);
  void PrintMyLogo();

  std::vector<double> output_time_;  // Output time list
    
  int output_step_;  // Find next output time in the list

  std::vector<int> output_flag_;

};


//--------------------------------------------------------------------------------------------------
// Get current date [year-month-day]
//--------------------------------------------------------------------------------------------------
inline Vec3i GetCurrentDate()
{
  time_t t = time(0);
  tm* now = localtime(&t);
  Vec3i current_date(now->tm_year + 1900, now->tm_mon + 1, now->tm_mday);
  return(current_date);
}
//--------------------------------------------------------------------------------------------------
// Get current time [hr:min:sec]
//--------------------------------------------------------------------------------------------------
inline Vec3i GetCurrentTime()
{
  time_t t = time(0);
  tm* now = localtime(&t);
  Vec3i current_time(now->tm_hour, now->tm_min, now->tm_sec);
  return(current_time);
}


//--------------------------------------------------------------------------------------------------
// Get time interval from start to now
//--------------------------------------------------------------------------------------------------
inline double GetTimeInterval(timeval* start, bool reset_start = false)
{
  timeval now ;
  gettimeofday(&now, NULL);
  
  double time_interval = (now.tv_sec - start->tv_sec) + (now.tv_usec - start->tv_usec) / 1000000.0;

  // Reset the start time to now
  if (reset_start)
  {
    gettimeofday(start, NULL);
  }

  return time_interval;
}


#endif
