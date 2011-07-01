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
/*! \file hermes_logging.h
\brief Functions and support for logging of events.
*/
/** \addtogroup g_logging Event Logging
*  \{
*  \brief Functions and support for logging of events.
*
*
*  By default, all logs are written into a file 'hermes.log' created in the current directory.
*  A logs created by the application can be directed to a separate file. The file is specified through
*  a compiler directive ::HERMES_REPORT_FILE, i.e.,\code
#define HERMES_REPORT_FILE "application.log"
*  \endcode
*  In a case, the application is a test (i.e., a directive ::HERMES_TEST is defined),
*  the output is directed to a file 'test.log'.
*  The output to a file can be suppressed specifying ::HERMES_REPORT_NO_FILE.
*
*  \section s_example Supported Directives
*  The following list contains directives that controls even logging.
*  - ::HERMES_REPORT_WARN: Define to allow warning.
*  - ::HERMES_REPORT_WARN_INTRO: Define to allow warning about integration issues.
*  - ::HERMES_REPORT_INFO: Define to allow info.
*  - ::HERMES_REPORT_VERBOSE: Define to allow verbose.
*  - ::HERMES_REPORT_TIME: Define to allow time measurement reports.
*  - ::HERMES_REPORT_TRACE: Define to allow execution tracing.
*  - ::HERMES_REPORT_FILE "my_file.log": Define to direct output of a file \c my_file_log.
*  - ::HERMES_REPORT_NO_FILE: Define to avoid any output file. It always overrides ::HERMES_REPORT_FILE.
*  - ::HERMES_REPORT_REPORT_ALL: Define to allow all events to be reported except integration warnings (::HERMES_REPORT_WARN_INTRO). It overrides all settings.
*  - ::HERMES_REPORT_RUNTIME_CONTROL: Define to allow controling of event logging through boolean variables.
*    Notice this will enforce evaluation of all parameters of logging macros even though a logging of a given event
*    is not enabled.
*  - ::HERMES_NO_LOGO: Define to disable logo message. This directive has to be defined at the compilation time of Hermes library.
*
*  \section s_usage Usage Guidelines
*  - Do \b not put any computation (e.g., \c it++) into parameters if the result of the computation is
*    used outside the macro. If a given event is not logged, no code may be generated and therefore
*    your computation will never be executed.
*  - Do \b not use new line characters (i.e. \c \\n or \c \\r) inside the message. Used a space or an exclamation
*    mark instead, see below.
*  - Use an exclamation mark at the beginning of the message to emphasize the message.
*  - Use a space at the beginning of the message to generate a sub-info.
*  - Use the event logging wisely since every line is flushed to a file if file output is enabled. This allows
*    to obtain all logged vents in a case a SIGSEGV happens on a remote machine, e.g., if an application is
*    executed using a batch system on a cluster.
*  - The format of the message is similar to the function printf().
*
*  \section s_example Examples of Use
*  The following code \code
#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include <hermes.h>
#include <solver_umfpack.h>
*  \endcode will enable logging of events warning, info and verbose in the application and it
*  copies the output to a file \c application.log.
*
*  The following code \code
info("Result is %d", 32);
info(" Probability of error is %g", 0.1);
trace("Computation is done.");
info("!Done");
*  \endcode will generate \verbatim
I Result is 32
Probability of error is 0.1
R Computation is done.

I Done. \endverbatim on screen if all events are enabled.
*/

#ifndef __HERMES_COMMON_LOGGING_H_
#define __HERMES_COMMON_LOGGING_H_

#include "compat.h"
#include <pthread.h>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <map>
#include <cstdio>
#include <stdarg.h>
#include <sstream>

namespace Hermes
{
  namespace Logging
  {
    /// Writes a fancy formatted text to a console. \internal \ingroup g_logging
    /** \param[in] code An event code, e.g., ::HERMES_EC_ERROR.
    *  \param[in] emphasize True if the message should be emphasized.
    *  \param[in] text A message. A C-style string.
    *  \return True if the message was written. False if it failed due to some reasone. */
    HERMES_API bool write_console(const char code, const bool emphasize, const char* text);

    /// Info about a log record. Used for output log function. \internal
    class HERMES_API HermesLogEventInfo 
    {
    public:
      const char code;          ///< An event code character. For defails see event characters, e.g., ::HERMES_EC_ERROR
      const char* log_file;     ///< Log file name.
      const char* src_function; ///< A name of a function/method at which the event was generated.
      const char* src_file;     ///< A source file at which the event was generated.
      const int src_line;       ///< A line in the source file at which the event was generated.
      HermesLogEventInfo(const char code, const char* log_file, const char* src_function, const char* src_file, const int src_line);

    };

    /// Exits the application if the condition is true. \internal
    /** Used by macros error() and error_if().
    *  \param[in] cond True if the function should exit.
    *  \param[in] code Exit code returned by the application throught exit(). */
    HERMES_API void hermes_exit_if(bool cond, int code = -1);

    /// Logs an event if the condition is true. \internal
    /** Used by all even logging macros. Since this function returns a copy of the parameter cond,
    *  it can be used to call a function hermes2d_exit_if() or a function(). Thanks to that, the macro
    *  behaves as a function rather than a block of code. Also, this allows a debugger to a particular
    *  code.
    *  \param[in] cond True if the event should be logged.
    *  \param[in] info Info about the event.
    *  \param[in] msg A message or prinf-like formatting string.
    *  \return A value of the parameter cond. */
    HERMES_API bool hermes_log_message_if(bool cond, const HermesLogEventInfo& info, const char* msg, ...);

    /* file operations */
    void HERMES_API hermes_fwrite(const void* ptr, size_t size, size_t nitems, FILE* stream);
    void HERMES_API hermes_fread(void* ptr, size_t size, size_t nitems, FILE* stream);

    /// Logging output monitor. \internal \ingroup g_logging
    /** This class protects a logging function __hermes_log_message_if() in multithreded environment. */
    class LoggerMonitor 
    {
      pthread_mutexattr_t mutex_attr; ///< Mutext attributes.
      pthread_mutex_t mutex; ///< Mutex that protects monitor.

    public:
      /// Constructor. Creates a mutex.
      LoggerMonitor() 
      {
        pthread_mutexattr_init(&mutex_attr);
        pthread_mutexattr_settype(&mutex_attr, PTHREAD_MUTEX_RECURSIVE);
        pthread_mutex_init(&mutex, &mutex_attr);
      };
      /// Destructor. Deletes a mutex.
      ~LoggerMonitor() 
      {
        pthread_mutex_destroy(&mutex);
        pthread_mutexattr_destroy(&mutex_attr);
      };

      /// Enters protected section.
      void enter() { pthread_mutex_lock(&mutex); };

      /// Leaves protected section.
      void leave() { pthread_mutex_unlock(&mutex); };
    };
  }
}

// Preprocessor directives follow.
using namespace Hermes::Logging;

/* event codes */
#define HERMES_EC_ERROR 'E' ///< An event code: errors. \internal
#define HERMES_EC_ASSERT 'X' ///< An event code: asserts. \internal
#define HERMES_EC_WARNING 'W' ///< An event code: warnings. \internal
#define HERMES_EC_INFO 'I' ///< An event code: info about results. \internal
#define HERMES_EC_VERBOSE 'V' ///< An event code: more details about details. \internal
#define HERMES_EC_TRACE 'R' ///< An event code: execution tracing. \internal
#define HERMES_EC_TIME 'T' ///< An event code: time measurements. \internal
#define HERMES_EC_DEBUG 'D' ///< An event code: general debugging messages. \internal

/// A size of a delimiter in a log file. \internal \ingroup g_logging
#define HERMES_LOG_FILE_DELIM_SIZE 80

/* function name */
/** \def __CURRENT_FUNCTION
*  \brief A platform-dependent string defining a current function. \internal */
#ifdef _WIN32 //Win32
# ifdef __MINGW32__ //MinGW
#   define __CURRENT_FUNCTION __func__
# else //MSVC and other compilers
#   define __CURRENT_FUNCTION __FUNCTION__
# endif
#else //Linux and Mac
# define __CURRENT_FUNCTION __PRETTY_FUNCTION__
#endif

/* log file */
#undef HERMES_LOG_FILE
#ifdef HERMES_REPORT_NO_FILE
#  define HERMES_LOG_FILE NULL
#else
# ifdef HERMES_REPORT_FILE
#  define HERMES_LOG_FILE HERMES_REPORT_FILE
# else
#  ifndef HERMES_TEST
#    define HERMES_LOG_FILE "hermes.log" // default filename for a library
#  else
#    define HERMES_LOG_FILE "test.log" // default filename for a library test
#  endif
# endif
#endif

// Builds info about an event. \internal
#define HERMES_BUILD_LOG_INFO(__event) HermesLogEventInfo(__event, HERMES_LOG_FILE, __CURRENT_FUNCTION, __FILE__, __LINE__)

/* error and assert macros */
#define error(...) hermes_exit_if(hermes_log_message_if(true, HERMES_BUILD_LOG_INFO(HERMES_EC_ERROR), __VA_ARGS__))
#define error_if(__cond, ...) hermes_exit_if(hermes_log_message_if(__cond, HERMES_BUILD_LOG_INFO(HERMES_EC_ERROR), __VA_ARGS__))
#ifndef NDEBUG
# define assert_msg(__cond, ...) assert(!hermes_log_message_if(!(__cond), HERMES_BUILD_LOG_INFO(HERMES_EC_ASSERT), __VA_ARGS__))
#else
# define assert_msg(__cond, ...)
#endif

/* reporting macros */
#ifdef HERMES_REPORT_ALL
# undef HERMES_REPORT_WARNING
# define HERMES_REPORT_WARNING
# undef HERMES_REPORT_NO_INTR_WARNING
# undef HERMES_REPORT_INFO
# define HERMES_REPORT_INFO
# undef HERMES_REPORT_VERBOSE
# define HERMES_REPORT_VERBOSE
# undef HERMES_REPORT_TRACE
# define HERMES_REPORT_TRACE
# undef HERMES_REPORT_TIME
# define HERMES_REPORT_TIME
#endif

#if defined(_DEBUG) || !defined(NDEBUG)
# define __HERMES_REP_DEBG true
#else
# define __HERMES_REP_DEBG false
#endif

#if defined(HERMES_REPORT_WARNING)
# define warn(...) hermes_log_message_if(true, HERMES_BUILD_LOG_INFO(HERMES_EC_WARNING), __VA_ARGS__)
# define warn_if(__cond, ...) hermes_log_message_if((__cond), HERMES_BUILD_LOG_INFO(HERMES_EC_WARNING), __VA_ARGS__)
#else
# define warn(...)
# define warn_if(__cond, ...)
#endif
#if defined(HERMES_REPORT_INTR_WARNING)
# define warn_intr(...) hermes_log_message_if(true, HERMES_BUILD_LOG_INFO(HERMES_EC_WARNING), __VA_ARGS__)
#else
# define warn_intr(...)
#endif
#if defined(HERMES_REPORT_INFO)
# define info(...) hermes_log_message_if(true, HERMES_BUILD_LOG_INFO(HERMES_EC_INFO), __VA_ARGS__)
# define info_if(__cond, ...) hermes_log_message_if((__cond), HERMES_BUILD_LOG_INFO(HERMES_EC_INFO), __VA_ARGS__)
#else
# define info(...)
# define info_if(__cond, ...)
#endif
#if defined(HERMES_REPORT_VERBOSE)
# define verbose(...) hermes_log_message_if(true, HERMES_BUILD_LOG_INFO(HERMES_EC_VERBOSE), __VA_ARGS__)
#else
# define verbose(...)
#endif
#if defined(HERMES_REPORT_TRACE)
# define trace(...) hermes_log_message_if(true, HERMES_BUILD_LOG_INFO(HERMES_EC_TRACE), __VA_ARGS__)
#else
# define trace(...)
#endif
#if defined(HERMES_REPORT_TIME)
# define report_time(...) hermes_log_message_if(true, HERMES_BUILD_LOG_INFO(HERMES_EC_TIME), __VA_ARGS__)
#else
# define report_time(...)
#endif
#if defined(_DEBUG) || !defined(NDEBUG)
# define debug_log(...) hermes_log_message_if(true, HERMES_BUILD_LOG_INFO(HERMES_EC_DEBUG), __VA_ARGS__)
#else
# define debug_log(...)
#endif
#endif