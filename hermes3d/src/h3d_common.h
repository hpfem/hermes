// This file is part of Hermes3D.
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes3D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 2 of the License,
// or (at your option) any later version.
//
// Hermes3D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes3D. If not, see <http://www.gnu.org/licenses/>.

#ifndef __H3D_COMMON_H_
#define __H3D_COMMON_H_

#include "../../hermes_common/common.h"

// H3D-specific error codes.
#define H3D_ERR_FACE_INDEX_OUT_OF_RANGE         "Face index out of range."
#define H3D_ERR_EDGE_INDEX_OUT_OF_RANGE         "Edge index out of range."
#define H3D_ERR_TETRA_NOT_COMPILED              "hermes3d was not built with tetra elements."
#define H3D_ERR_HEX_NOT_COMPILED                "hermes3d was not built with hex elements."
#define H3D_ERR_PRISM_NOT_COMPILED              "hermes3d was not built with prism elements."

// 1D element modes
enum ElementMode1D {
	MODE_LINE = 0
};

// 2D element modes
enum ElementMode2D {
	MODE_TRIANGLE = 0,
	MODE_QUAD = 1
};

// 3D element modes
enum ElementMode3D {
	MODE_TETRAHEDRON = 0,
	MODE_HEXAHEDRON = 1,
	MODE_PRISM = 2
};

enum ESpaceType {
	H1 = 1,
	Hcurl = 2,
	Hdiv = 3,
	L2 = 4
};

// points
struct HERMES_API Point1D {
	double x;			// coordinates of a point
};

struct HERMES_API Point2D {
	double x, y;		// coordinates of a point
};

struct HERMES_API Point3D {
	double x, y, z;		// coordinates of a point
};

inline double dot_product(const Point3D &a, const Point3D &b) { return a.x * b.x + a.y * b.y + a.z * b.z;}

inline Point3D cross_product(Point3D a, Point3D b) {
	Point3D r = {
		a.y * b.z - a.z * b.y,
		a.z * b.x - a.x * b.z,
		a.x * b.y - a.y * b.x
	};
	return r;
}

inline Point3D lin_comb(Point3D a, double coef_a, Point3D b, double coef_b) {
	Point3D r = {
		coef_a * a.x + coef_b * b.x,
		coef_a * a.y + coef_b * b.y,
		coef_a * a.z + coef_b * b.z,
	};
	return r;
}

inline double norm(const Point3D &pt) { return sqrt(dot_product(pt, pt)); }

inline Point3D normalize(const Point3D &pt) {
	double n = norm(pt);
	Point3D res = { pt.x / n, pt.y / n, pt.z / n };
	return res;
}


struct Vector3D {
	scalar x, y, z;		// coordinates of a point

	Vector3D() {
		x = y = z = 0;
	}

	Vector3D(double x, double y, double z) {
		this->x = x;
		this->y = y;
		this->z = z;
	}

	void set(double x, double y, double z) {
		this->x = x;
		this->y = y;
		this->z = z;
	}

	scalar dot_product(Vector3D vec2) {return x * vec2.x + y * vec2.y + z * vec2.z;};
	scalar dot_product(Point3D vec2) {return x * vec2.x + y * vec2.y + z * vec2.z;};
	double norm() { return sqrt(REAL(dot_product(*this)));};
	void cross_product(Vector3D a, Vector3D b) {
		x = a.y * b.z - a.z * b.y;
		y = a.z * b.x - a.x * b.z;
		z = a.x * b.y - a.y * b.x;
	}
	void cross_product(Point3D a, Vector3D b) {
		x = a.y * b.z - a.z * b.y;
		y = a.z * b.x - a.x * b.z;
		z = a.x * b.y - a.y * b.x;
	}
	void cross_product(Vector3D a, Point3D b) {
		x = a.y * b.z - a.z * b.y;
		y = a.z * b.x - a.x * b.z;
		z = a.x * b.y - a.y * b.x;
	}

	void normalize(){
		double n = norm();
		x /= n;
		y /= n;
		z /= n;
	}

	void subtract(Vector3D b){
		x -= b.x;
		y -= b.y;
		z -= b.z;
	}
};

// maximal polynomial order of elements
#define H3D_MAX_ELEMENT_ORDER							10

#define H3D_EC_TIME 'T' ///< An event code: time measurements. \internal
#define H3D_REPORT_TIME


/* basic logging functions */
/// Info about a log record. Used for output log function. \internal
struct HERMES_API Hermes3DLogEventInfo {
  const char code; ///< An event code character. For defails see event characters, e.g., ::H3D_EC_ERROR
  const char* log_file; ///< Log file name.
  const char* src_function; ///< A name of a function/method at which the event was generated.
  const char* src_file; ///< A source file at which the event was generated.
  const int src_line; ///< A line in the source file at which the event was generated.
  Hermes3DLogEventInfo(const char code, const char* log_file, const char* src_function, const char* src_file, const int src_line)
    : code(code), log_file(log_file), src_function(src_function), src_file(src_file), src_line(src_line) {};
};

/// Builds info about an event. \internal
#define H3D_BUILD_LOG_INFO(__event) Hermes3DLogEventInfo(__event, H3D_LOG_FILE, __CURRENT_FUNCTION, __FILE__, __LINE__)


/// Exits the application if the condition is true. \internal
/** Used by macros error() and error_if().
 *  \param[in] cond True if the function should exit.
 *  \param[in] code Exit code returned by the application throught exit(). */
void HERMES_API hermes3d_exit_if(bool cond, int code = -1);


/// Logs an event if the condition is true. \internal
/** Used by all even logging macros. Since this function returns a copy of the parameter cond,
 *  it can be used to call a function hermes2d_exit_if() or a function(). Thanks to that, the macro
 *  behaves as a function rather than a block of code. Also, this allows a debugger to a particular
 *  code.
 *  \param[in] cond True if the event should be logged.
 *  \param[in] info Info about the event.
 *  \param[in] msg A message or prinf-like formatting string.
 *  \return A value of the parameter cond. */
bool HERMES_API hermes3d_log_message_if(bool cond, const Hermes3DLogEventInfo& info, const char* msg, ...);


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
#define H3D_EC_ERROR 'E' ///< An event code: errors. \internal
#define H3D_EC_ASSERT 'X' ///< An event code: asserts. \internal
#define H3D_EC_WARNING 'W' ///< An event code: warnings. \internal
#define H3D_EC_INFO 'I' ///< An event code: info about results. \internal
#define H3D_EC_VERBOSE 'V' ///< An event code: more details about details. \internal
#define H3D_EC_TRACE 'R' ///< An event code: execution tracing. \internal
#define H3D_EC_TIME 'T' ///< An event code: time measurements. \internal
#define H3D_EC_DEBUG 'D' ///< An event code: general debugging messages. \internal

/* error and assert macros */
#define error(...) hermes3d_exit_if(hermes3d_log_message_if(true, H3D_BUILD_LOG_INFO(H3D_EC_ERROR), __VA_ARGS__))
#define error_if(__cond, ...) hermes3d_exit_if(hermes3d_log_message_if(__cond, H3D_BUILD_LOG_INFO(H3D_EC_ERROR), __VA_ARGS__))
#ifndef NDEBUG
# define assert_msg(__cond, ...) assert(!hermes3d_log_message_if(!(__cond), H3D_BUILD_LOG_INFO(H3D_EC_ASSERT), __VA_ARGS__))
#else
# define assert_msg(__cond, ...)
#endif

/* reporting macros */
#ifdef H3D_REPORT_ALL
# undef H3D_REPORT_WARNING
# define H3D_REPORT_WARNING
# undef HERMED2D_REPORT_NO_INTR_WARNING
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
# define warn(...) hermes3d_log_message_if(true && H3D_RCTR(__H3D_report_warn), H3D_BUILD_LOG_INFO(H3D_EC_WARNING), __VA_ARGS__)
# define warn_if(__cond, ...) hermes3d_log_message_if((__cond) && H3D_RCTR(__H3D_report_warn), H3D_BUILD_LOG_INFO(H3D_EC_WARNING), __VA_ARGS__)
#else
# define warn(...)
# define warn_if(__cond, ...)
#endif
#if defined(H3D_REPORT_INTR_WARNING) || defined(H3D_REPORT_RUNTIME_CONTROL)
# define warn_intr(...) hermes3d_log_message_if(H3D_RCTR(__H3D_report_warn_intr), H3D_BUILD_LOG_INFO(H3D_EC_WARNING), __VA_ARGS__)
#else
# define warn_intr(...)
#endif
#if defined(H3D_REPORT_INFO) || defined(H3D_REPORT_RUNTIME_CONTROL)
# define info(...) hermes3d_log_message_if(true  && H3D_RCTR(__H3D_report_info), H3D_BUILD_LOG_INFO(H3D_EC_INFO), __VA_ARGS__)
# define info_if(__cond, ...) hermes3d_log_message_if((__cond) && H3D_RCTR(__H3D_report_warn), H3D_BUILD_LOG_INFO(H3D_EC_INFO), __VA_ARGS__)
#else
# define info(...)
# define info_if(__cond, ...)
#endif
#if defined(H3D_REPORT_VERBOSE) || defined(H3D_REPORT_RUNTIME_CONTROL)
# define verbose(...) hermes3d_log_message_if(true && H3D_RCTR(__H3D_report_verbose), H3D_BUILD_LOG_INFO(H3D_EC_VERBOSE), __VA_ARGS__)
#else
# define verbose(...)
#endif
#if defined(H3D_REPORT_TRACE) || defined(H3D_REPORT_RUNTIME_CONTROL)
# define trace(...) hermes3d_log_message_if(true && H3D_RCTR(__H3D_report_trace), H3D_BUILD_LOG_INFO(H3D_EC_TRACE), __VA_ARGS__)
#else
# define trace(...)
#endif
#if defined(H3D_REPORT_TIME) || defined(H3D_REPORT_RUNTIME_CONTROL)
# define report_time(...) hermes3d_log_message_if(true && H3D_RCTR(__H3D_report_time), H3D_BUILD_LOG_INFO(H3D_EC_TIME), __VA_ARGS__)
#else
# define report_time(...)
#endif
#if defined(_DEBUG) || !defined(NDEBUG) || defined(H3D_REPORT_RUNTIME_CONTROL)
# define debug_log(...) hermes3d_log_message_if(true && H3D_RCTR(__H3D_report_debug), H3D_BUILD_LOG_INFO(H3D_EC_DEBUG), __VA_ARGS__)
#else
# define debug_log(...)
#endif



#endif
