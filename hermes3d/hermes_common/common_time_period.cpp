#include <math.h> 
#include <time.h> 
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <fstream>

#include "common_time_period.h"

using namespace std;

TimePeriod::TimePeriod(const char *name) : period_name(name == NULL ? "unnamed" : name) {
  //initialization
#ifdef WIN32 //Windows
  LARGE_INTEGER freq;
  if (QueryPerformanceFrequency(&freq))
    frequency = (double)freq.QuadPart;
  else
    frequency = -1;
#endif //Linux
  tick_reset();
}

TimePeriod::SysTime TimePeriod::get_time() const {
#ifdef WIN32 //Windows
  if (frequency > 0) {
    LARGE_INTEGER ticks;
    QueryPerformanceCounter(&ticks);
    return ticks.QuadPart;
  }
  else {
    return clock();
  }
#elif defined(__APPLE__) //Mac
  // FIXME: implement time measurement on Mac
  timespec tm;
  return tm;
#else //Linux
  timespec tm;
  clock_gettime(CLOCK_REALTIME, &tm);
  return tm;
#endif
}

double TimePeriod::period_in_seconds(const SysTime& begin, const SysTime& end) const {
#ifdef WIN32 //Windows
  uint64_t period = end - begin;
  if (frequency > 0)
    return period / frequency;
  else
    return period / (double)CLOCKS_PER_SEC;
#else //Linux
  int sec_corr = 0;
  long period_nsec = end.tv_nsec - begin.tv_nsec;
  if (period_nsec < 0) {
    sec_corr += -1;
    period_nsec += 1000000000UL;
  }
  long period_sec = (long)(end.tv_sec - begin.tv_sec) + sec_corr;
  return period_sec + (period_nsec / 1E9);
#endif
}

const TimePeriod& TimePeriod::tick(TimerPeriodTickType type) {
  SysTime cur_time = get_time();
  if (type == HERMES_ACCUMULATE) {
    double secs = period_in_seconds(last_time, cur_time);
    accum += secs;
    last_period = secs;
  }
  else {
    last_period = -1;
  }
  last_time = cur_time;
  return *this;
}

const TimePeriod& TimePeriod::tick_reset() {
  tick(HERMES_SKIP);
  reset();
  return *this;
}

const TimePeriod& TimePeriod::reset() {
  accum = 0;
  last_period = -1;
  return *this;
}

string TimePeriod::to_string(double secs) const {
  if (secs < 0)
    return "NO TIME";
  else {
  int hours = (int) secs / (3600);
  int mins = (int) fmod(secs, 3600) / 60;
  secs = fmod(secs, 60);

  stringstream str;
  if (hours > 0)
    str << hours << "h ";
  if (hours > 0 || mins > 0)
    str << mins << "m ";
  str << secs << "s";

  return str.str();
  }
}

ostream& operator<<(ostream& stream, const TimePeriod& period) {
  stream << period.accumulated_str();
  return stream;
}
