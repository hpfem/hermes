#ifndef __COMMON_LOGGING_H_
#define __COMMON_LOGGING_H_

#include "compat.h"
#include <pthread.h>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <map>
#include <cstdio>
#include <stdarg.h>
#include <sstream>

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
extern HERMES_API void hermes_exit_if(bool cond, int code = -1);

/// Logs an event if the condition is true. \internal
/** Used by all even logging macros. Since this function returns a copy of the parameter cond,
 *  it can be used to call a function hermes2d_exit_if() or a function(). Thanks to that, the macro
 *  behaves as a function rather than a block of code. Also, this allows a debugger to a particular
 *  code.
 *  \param[in] cond True if the event should be logged.
 *  \param[in] info Info about the event.
 *  \param[in] msg A message or prinf-like formatting string.
 *  \return A value of the parameter cond. */
extern HERMES_API bool hermes_log_message_if(bool cond, const HermesLogEventInfo& info, const char* msg, ...);

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

/* file operations */
void __hermes_fwrite(const void* ptr, size_t size, size_t nitems, FILE* stream, const HermesLogEventInfo& err_info);
void __hermes_fread(void* ptr, size_t size, size_t nitems, FILE* stream, const HermesLogEventInfo& err_info);

#endif
