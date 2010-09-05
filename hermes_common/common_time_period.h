// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef __HERMES_COMMON_TIME_PERIOD_H
#define __HERMES_COMMON_TIME_PERIOD_H

// For uint64_t type on windows.
#ifdef _MSC_VER
#include <inttypes.h>
#endif

/// Tick type. Used by the class TimePeriod.
enum TimerPeriodTickType {
  HERMES_ACCUMULATE, ///< Accumulate a period between ticks.
  HERMES_SKIP ///< Skip period between ticks, i.e., do not accumulate it.
};

/// A named time period measurement with accumulation.
/** An instance of the timer should not be used across threads. The class is not thread-safe.
 *  \todo Measure time that CPU spent on the task instead of a global time. */
class TimePeriod {
public:
  TimePeriod(const char *name = NULL); ///< Constructs internal structures and starts measuring.

  const TimePeriod& reset(); ///< Resets accumulated time.
  const TimePeriod& tick_reset(); ///< Starts a new period and resets accumulated time.
  const TimePeriod& tick(TimerPeriodTickType type = HERMES_ACCUMULATE); ///< Starts/ends a new period.

  /// Returns a name of the time period if any.
  const std::string& name() const { return period_name; }

  /// Returns accumulated time (in seconds).
  double accumulated() const { return accum; };

  /// Returns accumulated time in human readable form.
  std::string accumulated_str() const { return to_string(accum); };

  /// Returns last measured period (in seconds).
  /** \return Returns the length of the last measured time period. -1 if period was skipped. */
  double last() const { return last_period; };

  /// Returns last measured period in human readable form.
  std::string last_str() const { return to_string(last_period); };

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

  friend std::ofstream& operator<<(std::ofstream& stream, const TimePeriod& period);
};
extern std::ostream& operator<<(std::ostream& stream, const TimePeriod& period);

#endif
