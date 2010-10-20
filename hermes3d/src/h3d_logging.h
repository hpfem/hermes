#ifndef __H3D_LOGGING_H_
#define __H3D_LOGGING_H_

#include "../../hermes_common/logging.h"

/* log file */
#undef H3D_LOG_FILE
#ifdef H3D_REPORT_NO_FILE
#  define H3D_LOG_FILE NULL
#else
# ifdef H3D_REPORT_FILE
#  define H3D_LOG_FILE H3D_REPORT_FILE
# else
#  ifndef H3D_TEST
#    define H3D_LOG_FILE "hermes3d.log" // default filename for a library
#  else
#    define H3D_LOG_FILE "test.log" // default filename for a library test
#  endif
# endif
#endif

/// Builds info about an event. \internal
#define H3D_BUILD_LOG_INFO(__event) HermesLogEventInfo(__event, H3D_LOG_FILE, __CURRENT_FUNCTION, __FILE__, __LINE__)

/* error and assert macros */
#define error(...) hermes_exit_if(hermes_log_message_if(true, H3D_BUILD_LOG_INFO(HERMES_EC_ERROR), __VA_ARGS__))
#define error_if(__cond, ...) hermes_exit_if(hermes_log_message_if(__cond, H3D_BUILD_LOG_INFO(HERMES_EC_ERROR), __VA_ARGS__))
#ifndef NDEBUG
# define assert_msg(__cond, ...) assert(!hermes_log_message_if(!(__cond), H3D_BUILD_LOG_INFO(HERMES_EC_ASSERT), __VA_ARGS__))
#else
# define assert_msg(__cond, ...)
#endif

/* reporting macros */
#define H3D_REPORT_TIME	 // report time by default

#ifdef H3D_REPORT_ALL
# undef H3D_REPORT_WARNING
# define H3D_REPORT_WARNING
# undef H3D_REPORT_NO_INTR_WARNING
# undef H3D_REPORT_INFO
# define H3D_REPORT_INFO
# undef H3D_REPORT_VERBOSE
# define H3D_REPORT_VERBOSE
# undef H3D_REPORT_TRACE
# define H3D_REPORT_TRACE
# undef H3D_REPORT_TIME
# define H3D_REPORT_TIME
#endif
/** \def H3D_RCTR(__var)
 *  \brief Defines a condition that can control whether logging of a given event is enabled. \internal
 *  An argument \a __var spefies a variable which can control a logging of a given event during
 *  runtime if runtime control is enabled through a preprocessor directive ::H3D_REPORT_RUNTIME_CONTROL. */
#ifdef H3D_REPORT_RUNTIME_CONTROL
# define H3D_RCTR(__var) __var /* reports will be controled also by runtime report control variables */
extern HERMES_API bool __H3D_report_warn;
extern HERMES_API bool __H3D_report_warn_intr;
extern HERMES_API bool __H3D_report_info;
extern HERMES_API bool __H3D_report_verbose;
extern HERMES_API bool __H3D_report_trace;
extern HERMES_API bool __H3D_report_time;
extern HERMES_API bool __H3D_report_debug;
#else
# define H3D_RCTR(__var) true /* reports will be controled strictly by preprocessor directives */
#endif

#if defined(H3D_REPORT_WARNING) || defined(H3D_REPORT_RUNTIME_CONTROL)
# define warn(...) hermes_log_message_if(true && H3D_RCTR(__H3D_report_warn), H3D_BUILD_LOG_INFO(HERMES_EC_WARNING), __VA_ARGS__)
# define warn_if(__cond, ...) hermes_log_message_if((__cond) && H3D_RCTR(__H3D_report_warn), H3D_BUILD_LOG_INFO(HERMES_EC_WARNING), __VA_ARGS__)
#else
# define warn(...)
# define warn_if(__cond, ...)
#endif
#if defined(H3D_REPORT_INTR_WARNING) || defined(H3D_REPORT_RUNTIME_CONTROL)
# define warn_intr(...) hermes_log_message_if(H3D_RCTR(__H3D_report_warn_intr), H3D_BUILD_LOG_INFO(HERMES_EC_WARNING), __VA_ARGS__)
#else
# define warn_intr(...)
#endif
#if defined(H3D_REPORT_INFO) || defined(H3D_REPORT_RUNTIME_CONTROL)
# define info(...) hermes_log_message_if(true  && H3D_RCTR(__H3D_report_info), H3D_BUILD_LOG_INFO(HERMES_EC_INFO), __VA_ARGS__)
# define info_if(__cond, ...) hermes_log_message_if((__cond) && H3D_RCTR(__H3D_report_warn), H3D_BUILD_LOG_INFO(HERMES_EC_INFO), __VA_ARGS__)
#else
# define info(...)
# define info_if(__cond, ...)
#endif
#if defined(H3D_REPORT_VERBOSE) || defined(H3D_REPORT_RUNTIME_CONTROL)
# define verbose(...) hermes_log_message_if(true && H3D_RCTR(__H3D_report_verbose), H3D_BUILD_LOG_INFO(HERMES_EC_VERBOSE), __VA_ARGS__)
#else
# define verbose(...)
#endif
#if defined(H3D_REPORT_TRACE) || defined(H3D_REPORT_RUNTIME_CONTROL)
# define trace(...) hermes_log_message_if(true && H3D_RCTR(__H3D_report_trace), H3D_BUILD_LOG_INFO(HERMES_EC_TRACE), __VA_ARGS__)
#else
# define trace(...)
#endif
#if defined(H3D_REPORT_TIME) || defined(H3D_REPORT_RUNTIME_CONTROL)
# define report_time(...) hermes_log_message_if(true && H3D_RCTR(__H3D_report_time), H3D_BUILD_LOG_INFO(HERMES_EC_TIME), __VA_ARGS__)
#else
# define report_time(...)
#endif
#if defined(_DEBUG) || !defined(NDEBUG) || defined(H3D_REPORT_RUNTIME_CONTROL)
# define debug_log(...) hermes_log_message_if(true && H3D_RCTR(__H3D_report_debug), H3D_BUILD_LOG_INFO(HERMES_EC_DEBUG), __VA_ARGS__)
#else
# define debug_log(...)
#endif

#endif
