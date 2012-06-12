// This file is part of HermesCommon
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes2D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
/*! \file time_period.cpp
    \brief File containing the class TimePeriod and relevant definitions for measuring time.
*/
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <fstream>

#include "time_period.h"

namespace Hermes
{
  TimePeriod::TimePeriod(const char *name) : period_name(name == NULL ? "unnamed" : name)
  {
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

  TimePeriod::SysTime TimePeriod::get_time() const
  {
#ifdef WIN32 //Windows
    if (frequency > 0)
    {
      LARGE_INTEGER ticks;
      QueryPerformanceCounter(&ticks);
      return ticks.QuadPart;
    }
    else
    {
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

  double TimePeriod::period_in_seconds(const SysTime& begin, const SysTime& end) const
  {
#ifdef WIN32 //Windows
    uint64_t period = end - begin;
    if (frequency > 0)
      return period / frequency;
    else
      return period / (double)CLOCKS_PER_SEC;
#else //Linux
    int sec_corr = 0;
    long period_nsec = end.tv_nsec - begin.tv_nsec;
    if (period_nsec < 0)
    {
      sec_corr += -1;
      period_nsec += 1000000000UL;
    }
    long period_sec = (long)(end.tv_sec - begin.tv_sec) + sec_corr;
    return period_sec + (period_nsec / 1E9);
#endif
  }

  const TimePeriod& TimePeriod::tick(TimerPeriodTickType type)
  {
    SysTime cur_time = get_time();
    if (type == HERMES_ACCUMULATE)
    {
      double secs = period_in_seconds(last_time, cur_time);
      accum += secs;
      last_period = secs;
    }
    else
      last_period = 0.0;

    last_time = cur_time;
    return *this;
  }

  const std::string& TimePeriod::name() const 
  { 
    return period_name; 
  }

  double TimePeriod::accumulated() const 
  { 
    return accum; 
  }

  std::string TimePeriod::accumulated_str() const { 
    return to_string(accum); 
  }

  double TimePeriod::last() const 
  { 
    return last_period; 
  }

  std::string TimePeriod::last_str() const 
  {
    return to_string(last_period); 
  }


  const TimePeriod& TimePeriod::tick_reset()
  {
    tick(HERMES_SKIP);
    reset();
    return *this;
  }

  const TimePeriod& TimePeriod::reset()
  {
    accum = 0;
    last_time = get_time();
    last_period = 0.0;
    return *this;
  }

  std::string TimePeriod::to_string(double secs) const
  {
    if (secs < 0)
      return "NO TIME";
    else
    {
      int hours = (int) secs / (3600);
      int mins = (int) fmod(secs, 3600) / 60;
      secs = fmod(secs, 60);

      std::stringstream str;
      if (hours > 0)
        str << hours << "h ";
      if (hours > 0 || mins > 0)
        str << mins << "m ";
      str << secs << "s";

      return str.str();
    }
  }

  std::ostream& operator<<(std::ostream& stream, const TimePeriod& period)
  {
    stream << period.accumulated_str();
    return stream;
  }
}
