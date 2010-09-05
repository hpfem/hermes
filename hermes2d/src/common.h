// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __H2D_COMMON_H
#define __H2D_COMMON_H

#include "config.h"

// common headers
#include <stdexcept>
#include <cstdarg>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h> //allows to use offsetof
#include <string.h>
#include <assert.h>
#include <pthread.h>
#include <math.h>
#include <time.h>

#include <float.h>

#include <cmath>

// Matrix solvers
enum MatrixSolverType 
{
   SOLVER_UMFPACK, 
   SOLVER_PETSC, 
   SOLVER_MUMPS
};

// STL stuff
#include <algorithm>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <queue>
#include <sstream>
#include <fstream>

// platform compatibility stuff
#include "compat.h"

// others
#include <Judy.h>
#include "auto_local_array.h"
#include "common_time_period.h"
#include "../../hermes_common/tuple.h"

// Enabling second derivatives in weak forms. Turned off by default. Second
// derivatives are employed, among others, by stabilization methods for
// transport equations. For usage see the example linear-convection-diffusion.
#define H2D_SECOND_DERIVATIVES_ENABLED

enum // node types
{
  H2D_TYPE_VERTEX = 0,
  H2D_TYPE_EDGE = 1
};

#define H2D_NUM_MODES 2 ///< A number of modes, see enum ElementMode.

enum ElementMode { // element modes
  H2D_MODE_TRIANGLE = 0,
  H2D_MODE_QUAD = 1
};

// default projection norm is H1 norm
// FIXME: this global variable should be declared here but 
// doing so leads to compilation problems. That's why it 
// is temporarily in linsystem.cpp.
//const int H2D_DEFAULT_PROJ_NORM = 1;

const int H2D_ANY = -1234;

// how many bits the order number takes
const int H2D_ORDER_BITS = 5;
const int H2D_ORDER_MASK = (1 << H2D_ORDER_BITS) - 1;
const int H2D_DEFAULT_PROJ_TYPE = 1;

// macros for combining quad horizontal and vertical orders
#define H2D_MAKE_QUAD_ORDER(h_order, v_order) (((v_order) << H2D_ORDER_BITS) + (h_order))
#define H2D_GET_H_ORDER(order) ((order) & H2D_ORDER_MASK)
#define H2D_GET_V_ORDER(order) ((order) >> H2D_ORDER_BITS)
extern H2D_API const std::string get_quad_order_str(const int quad_order); ///< Returns string representation of the quad order: used for debugging purposses.

extern H2D_API int make_edge_order(int edge, int encoded_order, int mode); ///< Returns the correct axial order for given edge.

#include "scalar.h"

typedef int int2[2];
typedef int int3[3];
typedef int int4[4];
typedef int int5[5];

typedef double double2[2];
typedef double double3[3];
typedef double double4[4];
typedef double double2x2[2][2];
typedef double double3x2[3][2];

typedef scalar scalar2[2];
typedef scalar scalar3[3];


inline int sqr(int x) { return x*x; }
inline double sqr(double x) { return x*x; }
#ifdef H2D_COMPLEX
inline double sqr(cplx x) { return std::norm(x); }
#endif

inline double magn(double x) { return fabs(x); }
#ifdef H2D_COMPLEX
inline double magn(cplx x) { return std::abs(x); }
#endif

inline double conj(double a) {  return a; }
#ifdef H2D_COMPLEX
inline cplx conj(cplx a) { return std::conj(a); }
#endif

#define H2D_IS_INT(x) ((int) (x) == (x))

/* logging functions */
/** \addtogroup g_logging Event Logging
 *  \{
 *  \brief Functions and support for logging of events.
 *
 *  Hermes2D controls even logging through:
 *  - compiler directives (e.g. ::H2D_REPORT_INFO). Directives has to be included prior including Hermes2D header files.
 *  - boolean variables (e.g. ::__h2d_report_info). These variables can be set anytime but their direct use is discouraged
 *    because they are integeded to be used by Python wreappers. Initial settings of these variables is given by the compiler
 *    directives.
 *
 *  By default, all logs are written into a file 'hermes2d.log' created in the current directory.
 *  A logs created by the application can be directed to a separate file. The file is specified through
 *  a compiler directive ::H2D_REPORT_FILE, i.e.,\code
    #define H2D_REPORT_FILE "application.log"
 *  \endcode
 *  In a case, the application is a test (i.e., a directive ::H2D_TEST is defined),
 *  the output is directed to a file 'test.log'.
 *  The output to a file can be suppressed specifying ::H2D_REPORT_NO_FILE.
 *
 *  \section s_example Supported Directives
 *  The following list contains directives that controls even logging.
 *  - ::H2D_REPORT_WARN: Define to allow warning.
 *  - ::H2D_REPORT_WARN_INTRO: Define to allow warning about integration issues.
 *  - ::H2D_REPORT_INFO: Define to allow info.
 *  - ::H2D_REPORT_VERBOSE: Define to allow verbose.
 *  - ::H2D_REPORT_TIME: Define to allow time measurement reports.
 *  - ::H2D_REPORT_TRACE: Define to allow execution tracing.
 *  - ::H2D_REPORT_FILE "my_file.log": Define to direct output of a file \c my_file_log.
 *  - ::H2D_REPORT_NO_FILE: Define to avoid any output file. It always overrides ::H2D_REPORT_FILE.
 *  - ::H2D_REPORT_REPORT_ALL: Define to allow all events to be reported except integration warnings (::H2D_REPORT_WARN_INTRO). It overrides all settings.
 *  - ::H2D_REPORT_RUNTIME_CONTROL: Define to allow controling of event logging through boolean variables.
 *    Notice this will enforce evaluation of all parameters of logging macros even though a logging of a given event
 *    is not enabled.
 *  - ::H2D_NO_LOGO: Define to disable logo message. This directive has to be defined at the compilation time of Hermes2D library.
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
    #define H2D_REPORT_WARN
    #define H2D_REPORT_INFO
    #define H2D_REPORT_VERBOSE
    #define H2D_REPORT_FILE "application.log"
    #include <hermes2d.h>
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

/* basic logging functions */
/// Info about a log record. Used for output log function. \internal
struct H2D_API Hermes2DLogEventInfo {
  const char code; ///< An event code character. For defails see event characters, e.g., ::H2D_EC_ERROR
  const char* log_file; ///< Log file name.
  const char* src_function; ///< A name of a function/method at which the event was generated.
  const char* src_file; ///< A source file at which the event was generated.
  const int src_line; ///< A line in the source file at which the event was generated.
  Hermes2DLogEventInfo(const char code, const char* log_file, const char* src_function, const char* src_file, const int src_line)
    : code(code), log_file(log_file), src_function(src_function), src_file(src_file), src_line(src_line) {};
};
/// Builds info about an event. \internal
#define H2D_BUILD_LOG_INFO(__event) Hermes2DLogEventInfo(__event, H2D_LOG_FILE, __CURRENT_FUNCTION, __FILE__, __LINE__)

/// Exits the application if the condition is true. \internal
/** Used by macros error() and error_if().
 *  \param[in] cond True if the function should exit.
 *  \param[in] code Exit code returned by the application throught exit(). */
extern H2D_API void hermes2d_exit_if(bool cond, int code = -1);

/// Logs an event if the condition is true. \internal
/** Used by all even logging macros. Since this function returns a copy of the parameter cond,
 *  it can be used to call a function hermes2d_exit_if() or a function(). Thanks to that, the macro
 *  behaves as a function rather than a block of code. Also, this allows a debugger to a particular
 *  code.
 *  \param[in] cond True if the event should be logged.
 *  \param[in] info Info about the event.
 *  \param[in] msg A message or prinf-like formatting string.
 *  \return A value of the parameter cond. */
extern H2D_API bool hermes2d_log_message_if(bool cond, const Hermes2DLogEventInfo& info, const char* msg, ...);

/* log file */
#undef H2D_LOG_FILE
#ifdef H2D_REPORT_NO_FILE
#  define H2D_LOG_FILE NULL
#else
# ifdef H2D_REPORT_FILE
#  define H2D_LOG_FILE H2D_REPORT_FILE
# else
#  ifndef H2D_TEST
#    define H2D_LOG_FILE "hermes2d.log" // default filename for a library
#  else
#    define H2D_LOG_FILE "test.log" // default filename for a library test
#  endif
# endif
#endif

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

/* event codes */
#define H2D_EC_ERROR 'E' ///< An event code: errors. \internal
#define H2D_EC_ASSERT 'X' ///< An event code: asserts. \internal
#define H2D_EC_WARNING 'W' ///< An event code: warnings. \internal
#define H2D_EC_INFO 'I' ///< An event code: info about results. \internal
#define H2D_EC_VERBOSE 'V' ///< An event code: more details about details. \internal
#define H2D_EC_TRACE 'R' ///< An event code: execution tracing. \internal
#define H2D_EC_TIME 'T' ///< An event code: time measurements. \internal
#define H2D_EC_DEBUG 'D' ///< An event code: general debugging messages. \internal

/* error and assert macros */
#define error(...) hermes2d_exit_if(hermes2d_log_message_if(true, H2D_BUILD_LOG_INFO(H2D_EC_ERROR), __VA_ARGS__))
#define error_if(__cond, ...) hermes2d_exit_if(hermes2d_log_message_if(__cond, H2D_BUILD_LOG_INFO(H2D_EC_ERROR), __VA_ARGS__))
#ifndef NDEBUG
# define assert_msg(__cond, ...) assert(!hermes2d_log_message_if(!(__cond), H2D_BUILD_LOG_INFO(H2D_EC_ASSERT), __VA_ARGS__))
#else
# define assert_msg(__cond, ...)
#endif

/* reporting macros */
#ifdef H2D_REPORT_ALL
# undef H2D_REPORT_WARNING
# define H2D_REPORT_WARNING
# undef HERMED2D_REPORT_NO_INTR_WARNING
# undef H2D_REPORT_INFO
# define H2D_REPORT_INFO
# undef H2D_REPORT_VERBOSE
# define H2D_REPORT_VERBOSE
# undef H2D_REPORT_TRACE
# define H2D_REPORT_TRACE
# undef H2D_REPORT_TIME
# define H2D_REPORT_TIME
#endif
/** \def H2D_RCTR(__var)
 *  \brief Defines a condition that can control whether logging of a given event is enabled. \internal
 *  An argument \a __var spefies a variable which can control a logging of a given event during
 *  runtime if runtime control is enabled through a preprocessor directive ::H2D_REPORT_RUNTIME_CONTROL. */
#ifdef H2D_REPORT_RUNTIME_CONTROL
# define H2D_RCTR(__var) __var /* reports will be controled also by runtime report control variables */
extern H2D_API bool __h2d_report_warn;
extern H2D_API bool __h2d_report_warn_intr;
extern H2D_API bool __h2d_report_info;
extern H2D_API bool __h2d_report_verbose;
extern H2D_API bool __h2d_report_trace;
extern H2D_API bool __h2d_report_time;
extern H2D_API bool __h2d_report_debug;
#else
# define H2D_RCTR(__var) true /* reports will be controled strictly by preprocessor directives */
#endif

#if defined(H2D_REPORT_WARNING) || defined(H2D_REPORT_RUNTIME_CONTROL)
# define warn(...) hermes2d_log_message_if(true && H2D_RCTR(__h2d_report_warn), H2D_BUILD_LOG_INFO(H2D_EC_WARNING), __VA_ARGS__)
# define warn_if(__cond, ...) hermes2d_log_message_if((__cond) && H2D_RCTR(__h2d_report_warn), H2D_BUILD_LOG_INFO(H2D_EC_WARNING), __VA_ARGS__)
#else
# define warn(...)
# define warn_if(__cond, ...)
#endif
#if defined(H2D_REPORT_INTR_WARNING) || defined(H2D_REPORT_RUNTIME_CONTROL)
# define warn_intr(...) hermes2d_log_message_if(H2D_RCTR(__h2d_report_warn_intr), H2D_BUILD_LOG_INFO(H2D_EC_WARNING), __VA_ARGS__)
#else
# define warn_intr(...)
#endif
#if defined(H2D_REPORT_INFO) || defined(H2D_REPORT_RUNTIME_CONTROL)
# define info(...) hermes2d_log_message_if(true  && H2D_RCTR(__h2d_report_info), H2D_BUILD_LOG_INFO(H2D_EC_INFO), __VA_ARGS__)
# define info_if(__cond, ...) hermes2d_log_message_if((__cond) && H2D_RCTR(__h2d_report_warn), H2D_BUILD_LOG_INFO(H2D_EC_INFO), __VA_ARGS__)
#else
# define info(...)
# define info_if(__cond, ...)
#endif
#if defined(H2D_REPORT_VERBOSE) || defined(H2D_REPORT_RUNTIME_CONTROL)
# define verbose(...) hermes2d_log_message_if(true && H2D_RCTR(__h2d_report_verbose), H2D_BUILD_LOG_INFO(H2D_EC_VERBOSE), __VA_ARGS__)
#else
# define verbose(...)
#endif
#if defined(H2D_REPORT_TRACE) || defined(H2D_REPORT_RUNTIME_CONTROL)
# define trace(...) hermes2d_log_message_if(true && H2D_RCTR(__h2d_report_trace), H2D_BUILD_LOG_INFO(H2D_EC_TRACE), __VA_ARGS__)
#else
# define trace(...)
#endif
#if defined(H2D_REPORT_TIME) || defined(H2D_REPORT_RUNTIME_CONTROL)
# define report_time(...) hermes2d_log_message_if(true && H2D_RCTR(__h2d_report_time), H2D_BUILD_LOG_INFO(H2D_EC_TIME), __VA_ARGS__)
#else
# define report_time(...)
#endif
#if defined(_DEBUG) || !defined(NDEBUG) || defined(H2D_REPORT_RUNTIME_CONTROL)
# define debug_log(...) hermes2d_log_message_if(true && H2D_RCTR(__h2d_report_debug), H2D_BUILD_LOG_INFO(H2D_EC_DEBUG), __VA_ARGS__)
#else
# define debug_log(...)
#endif

/** \def error(...)
 *  \brief Logs an error and quits the application. For details see \ref s_usage. */
/** \def error_if(__cond, ...)
 *  \brief If \a __cond is true, it logs an error and quits the application.
 *
 *  For usage guidelines, see \ref s_usage. */

/** \def assert_msg(__cond, ...)
 *  \brief If \a __cond is false, it logs a message and invokes assert().
 *
 *  Similar to the function assert() if ::NDEBUG is defined, no code may be generated for the macro.
 *  For usage guidelines, see \ref s_usage. */

/** \def warn(...)
 *  \brief Logs a warning.
 *
 *  If ::H2D_REPORT_WARN is not defined, no code may be generated for the macro.
 *  For usage guidelines, see \ref s_usage. */

/** \def warn_if(__cond, ...)
 *  \brief If \a __cond is true, it logs a warning.
 *
 *  If ::H2D_REPORT_WARN is not defined, no code may be generated for the macro.
 *  For usage guidelines, see \ref s_usage. */

/** \def warn_intr(...)
 *  \brief Logs an warning about integration. This is used to report integration issues which
 *  may occur frequently.
 *
 *  If ::H2D_REPORT_INTR_WARNING is not defined, no code may be generated for the macro.
 *  For usage guidelines, see \ref s_usage. */

/** \def info(...)
 *  \brief Logs info about a result of an operation.
 *
 *  If ::H2D_REPORT_INFO is not defined, no code may be generated for the macro.
 *  For usage guidelines, see \ref s_usage. */

/** \def info_if(__cond, ...)
 *  \brief If \a __cond is true, it logs info about a result of na operation.
 *
 *  If ::H2D_REPORT_WARN is not defined, no code may be generated for the macro.
 *  For usage guidelines, see \ref s_usage. */

/** \def verbose(...)
 *  \brief Logs detailed info about a result of an operation. It should be used as a second level
 *  of the macro info().
 *
 *  If ::H2D_REPORT_VERBOSE is not defined, no code may be generated for the macro.
 *  For usage guidelines, see \ref s_usage. */

/** \def report_time(...)
 *  \brief Logs information about measured time.
 *
 *  If ::H2D_REPORT_TIME is not defined, no code may be generated for the macro.
 *  For usage guidelines, see \ref s_usage. */

/** \def trace(...)
 *  \brief Logs information about executed code portions.
 *
 *  If ::H2D_REPORT_TRACE is not defined, no code may be generated for the macro.
 *  For usage guidelines, see \ref s_usage. */

/** \def debug_log(...)
 *  \brief Logs a general debugging information.
 *
 *  This macro should be used for debugging outputs and it is suggested to remove
 *  almost all of them an debugged issue is solved.
 *  If ::NDEBUG is defined, no code may be generated for the macro.
 *  For usage guidelines, see \ref s_usage. */

/** \} */

/* file operations */
void __hermes2d_fwrite(const void* ptr, size_t size, size_t nitems, FILE* stream, const Hermes2DLogEventInfo& err_info);
void __hermes2d_fread(void* ptr, size_t size, size_t nitems, FILE* stream, const Hermes2DLogEventInfo& err_info);

#define hermes2d_fwrite(ptr, size, nitems, stream) \
      __hermes2d_fwrite((ptr), (size), (nitems), (stream), H2D_BUILD_LOG_INFO(H2D_EC_ERROR))

#define hermes2d_fread(ptr, size, nitems, stream) \
      __hermes2d_fread((ptr), (size), (nitems), (stream), H2D_BUILD_LOG_INFO(H2D_EC_ERROR))

/* python support */
/// Throws an exception std::runtime_error. Used by Python wrappers.
/** \param[in] text A text (a cause) of the exception. */
extern H2D_API void throw_exception(char *text);

/* Uncomment this line to disable internal mesh compatibility 
   tests in Traverse:begin(). */ 
//#define H2D_DISABLE_MULTIMESH_TESTS

#endif

