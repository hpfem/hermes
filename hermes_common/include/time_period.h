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
/*! \file time_period.h
    \brief File containing the class TimePeriod and relevant definitions for measuring time.
*/
#ifndef __HERMES_time_period_H
#define __HERMES_time_period_H

// For uint64_t type on windows.
#include "compat.h"

#ifdef _MSC_VER
#include <inttypes.h>
#endif

namespace Hermes
{
  /// Tick type. Used by the class Hermes::TimePeriod.
  enum TimerPeriodTickType
  {
    HERMES_ACCUMULATE, ///< Accumulate a period between ticks.
    HERMES_SKIP ///< Skip period between ticks, i.e., do not accumulate it.
  };

  /// A named time period measurement with accumulation.
  /** An instance of the timer should not be used across threads. The class is not thread-safe.
  *  \todo Measure time that CPU spent on the task instead of a global time. */
  class HERMES_API TimePeriod
  {
  public:
    TimePeriod(const char *name = NULL); ///< Constructs internal structures and starts measuring.

    const TimePeriod& reset(); ///< Resets accumulated time.
    const TimePeriod& tick_reset(); ///< Starts a new period and resets accumulated time.
    const TimePeriod& tick(TimerPeriodTickType type = HERMES_ACCUMULATE); ///< Starts/ends a new period.

    /// Returns a name of the time period if any.
    const std::string& name() const;

    /// Returns accumulated time (in seconds).
    double accumulated() const;

    /// Returns accumulated time in human readable form.
    std::string accumulated_str() const;

    /// Returns last measured period (in seconds).
    /** \return Returns the length of the last measured time period. -1 if period was skipped. */
    double last() const;

    /// Returns last measured period in human readable form.
    std::string last_str() const;

  private:
#ifdef WIN32 //Windows
    typedef uint64_t SysTime;
    double frequency; ///< Frequency of the performance timer. If zero, no hi-res timer is supported. (Win32 only)
#else //Linux
    typedef timespec SysTime;
#endif
    const std::string period_name; ///< Name of the timer (can be empty)
    double last_period; ///< Time of the last measured period.
    SysTime last_time; ///< Time when the timer was started/resumed (in platform-dependent units).
    double accum; ///< Time accumulator (in seconds).

    SysTime get_time() const; ///< Returns current time (in platform-dependent units).
    double period_in_seconds(const SysTime& begin, const SysTime& end) const; ///< Calculates distance between times (in platform specific units) and returns it in seconds.
    std::string to_string(const double time) const; ///< Converts time from seconds to human readable form.

    friend std::ofstream& operator<<(std::ofstream& stream, const Hermes::TimePeriod& period);
  };
  extern std::ostream& operator<<(std::ostream& stream, const Hermes::TimePeriod& period);
}
#endif
