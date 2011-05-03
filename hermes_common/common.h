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

#ifndef __HERMES_COMMON_COMMON_H
#define __HERMES_COMMON_COMMON_H

// Include
//
// common headers

#ifdef _POSIX_C_SOURCE
# undef _POSIX_C_SOURCE	// typeinfo defines it
#endif
#ifdef _XOPEN_SOURCE
# undef _XOPEN_SOURCE	// typeinfo defines it
#endif
#include <typeinfo>

#include <stdexcept>
#include <cstdarg>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h> //allows to use offsetof
#include <assert.h>
#include <pthread.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <errno.h>
#include <cmath>
//
// STL stuff
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <queue>
#include <sstream>
#include <fstream>
#include <cstring>
//
// commonly used functions from hermes_common
#include "hermes_logging.h"       // logging
#include "common_time_period.h"   // timing utilities
#include "compat.h"               // platform compatibility stuff
#include "callstack.h"            // error tracing
#include "error.h"
//
#include "vector.h"
#include "tables.h"

#define HERMES  "Hermes"

// Decide which version of Hermes is being compiled and import
// the build options from the corresponding config.h file.
#ifndef CONFIG_H_INCLUDED
  #include "config.h"
  #define CONFIG_H_INCLUDED
#endif

// Error codes
#define HERMES_ERR_NOT_IMPLEMENTED                 "Not yet implemened."
#define HERMES_ERR_UNKNOWN_MODE                    "Unknown mode (mode = %d)."
#define HERMES_ERR_UNKNOWN_REFINEMENT_TYPE         "Unknown refinement type (refinement = %d)."

// Matrix solvers (maybe we could move it to solver/solver.h)
enum MatrixSolverType 
{
   SOLVER_UMFPACK = 0, 
   SOLVER_PETSC, 
   SOLVER_MUMPS,
   SOLVER_SUPERLU,
   SOLVER_AMESOS,
   SOLVER_AZTECOO
};

// Should be in the same order as MatrixSolverTypes above, so that the
// names may be accessed by the same enumeration variable.
const std::string MatrixSolverNames[6] = {
  "UMFPACK",
  "PETSc",
  "MUMPS",
  "SuperLU",
  "Trilinos/Amesos",
  "Trilinos/AztecOO"
};

#define UMFPACK_NOT_COMPILED  HERMES " was not built with UMFPACK support."
#define PETSC_NOT_COMPILED    HERMES " was not built with PETSC support."
#define MUMPS_NOT_COMPILED    HERMES " was not built with MUMPS support."
#define SUPERLU_NOT_COMPILED  HERMES " was not built with SUPERLU support."
#define NOX_NOT_COMPILED      HERMES " was not built with NOX support."
#define AMESOS_NOT_COMPILED   HERMES " was not built with AMESOS support."
#define AZTECOO_NOT_COMPILED  HERMES " was not built with AZTECOO support."
#define EPETRA_NOT_COMPILED   HERMES " was not built with EPETRA support."
#define IFPACK_NOT_COMPILED   HERMES " was not built with IFPACK support."
#define ML_NOT_COMPILED       HERMES " was not built with ML support."

// Spaces.
enum ESpaceType {
  HERMES_H1_SPACE = 0,
  HERMES_HCURL_SPACE = 1,
  HERMES_HDIV_SPACE = 2,
  HERMES_L2_SPACE = 3,
  HERMES_INVALID_SPACE = -9999
};


// Solutions.
enum ESolutionType {
  HERMES_UNDEF = -1,
  HERMES_SLN = 0,
  HERMES_EXACT = 1,
  HERMES_CONST = 2
};


// Projection norms.
enum ProjNormType
{
  HERMES_L2_NORM, 
  HERMES_H1_NORM, 
  HERMES_H1_SEMINORM, 
  HERMES_HCURL_NORM, 
  HERMES_HDIV_NORM,
  // Used for passing to projecting functions.
  HERMES_UNSET_NORM
};

// Splines.
// (In weak forms, NULL spline is translated into a constant spline with value 1.0.)
#define HERMES_DEFAULT_SPLINE      NULL

// Functions.
// (In weak forms, NULL function is translated into a constant function with value 1.0.)
#define HERMES_DEFAULT_FUNCTION    NULL

#ifdef HERMES_COMMON_COMPLEX

  #include <complex>

  typedef std::complex<double> cplx;
  typedef cplx complex2[2];
  typedef cplx scalar;

  #define CONJ(a)       (std::conj(a))
  #define REAL(a)       (std::real(a))
  #define IMAG(a)       (std::imag(a))
  #define ABS(a)        (std::abs(a))
  #define SCALAR_FMT      "(%lf, %lf)"
  #define SCALAR(a)     std::real(a), std::imag(a)
  
  inline double sqr(cplx x)   { return std::norm(x); }
  inline double magn(cplx x)  { return std::abs(x); }
  inline cplx conj(cplx a)    { return std::conj(a); }
  
  #ifdef WITH_BLAS         // always true for Hermes3D
  // BLAS-related functions

    #ifdef __cplusplus
    extern "C" {
    #endif

    extern int zscal_(int *, cplx *, cplx *, int *);
    extern int zaxpy_(int *, cplx *, cplx *, int *, cplx *, int *);
    extern int zcopy_(int *, cplx *, int *, cplx *, int *);

    #ifdef __cplusplus
    }
    #endif

    /// x <- alpha * x
    inline void blas_scal(int n, cplx alpha, cplx *x, int incx) { zscal_(&n, &alpha, x, &incx); }
    /// y <- alpha * x + y
    inline void blas_axpy(int n, cplx alpha, cplx *x, int incx, cplx *y, int incy) { zaxpy_(&n, &alpha, x, &incx, y, &incy); }
    /// y <- x
    inline void blas_copy(int n, cplx *x, int incx, cplx *y, int incy) { zcopy_(&n, x, &incx, y, &incx); }
  #endif

#else

  typedef double scalar;
  
  #define CONJ(a)       (a)
  #define REAL(a)       (a)
  #define IMAG(a)       (0)
  #define ABS(a)        (fabs(a))
  #define SCALAR_FMT      "%lf"
  #define SCALAR(a)     (a)

#endif

enum // node types
{
  HERMES_TYPE_VERTEX = 0,
  HERMES_TYPE_EDGE = 1
};

// 1D element modes
enum ElementMode1D {
	HERMES_MODE_LINE = 0
};

// 2D element modes
enum ElementMode2D {
	HERMES_MODE_TRIANGLE = 0,
	HERMES_MODE_QUAD = 1
};

// 3D element modes
enum ElementMode3D {
	HERMES_MODE_TET = 0,
	HERMES_MODE_HEX = 1,
	HERMES_MODE_PRISM = 2
};

// Points and vectors.
struct HERMES_API Point1D {
	double x;		// coordinates of a point
};

struct HERMES_API Point2D {
	double x, y;		// coordinates of a point
};

struct HERMES_API Point3D {
	double x, y, z;		// coordinates of a point
};

struct HERMES_API SplineCoeff {
  double a, b, c, d;		// four coefficients of a cubic spline.
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

// Default HERMES projection norm is the H1 norm.
const ProjNormType HERMES_DEFAULT_PROJ_NORM = HERMES_H1_NORM;

inline int sqr(int x) { return x*x; }
inline double sqr(double x) { return x*x; }
inline double magn(double x) { return fabs(x); }
inline double conj(double a) {  return a; }

#ifdef WITH_BLAS         // always true for Hermes3D
// BLAS-related functions
  #ifdef __cplusplus
  extern "C" {
    #endif
    
    extern int dscal_(int *, double *, double *, int *);
    extern int daxpy_(int *, double *, double *, int *, double *, int *);
    extern int dcopy_(int *,           double *, int *, double *, int *);
    
    #ifdef __cplusplus
  }
  #endif

  /// x <- alpha * x
  inline void blas_scal(int n, double alpha, double *x, int incx) { dscal_(&n, &alpha, x, &incx); }
  /// y <- alpha * x + y
  inline void blas_axpy(int n, double alpha, double *x, int incx, double *y, int incy) { daxpy_(&n, &alpha, x, &incx, y, &incy); }
  /// y <- x
  inline void blas_copy(int n, double *x, int incx, double *y, int incy) { dcopy_(&n, x, &incx, y, &incx); }
#endif

typedef int int2[2];
typedef int int3[3];
typedef int int4[4];
typedef int int5[5];

typedef double double2[2];
typedef double double3[3];
typedef double double4[4];
typedef double double2x2[2][2];
typedef double double3x2[3][2];
typedef double double3x3[3][3];

struct scalar2 
{ 
  scalar val[2]; 

 public:
  scalar2(scalar v1, scalar v2) 
  { 
    val[0] = v1; val[1] = v2; 
  }

  scalar& operator[] (int idx) 
  { 
    assert(idx >= 0 && idx < 2);
    return val[idx];
  }
};

struct scalar3
{ 
  scalar val[3]; 

 public:
  scalar3(scalar v1, scalar v2, scalar v3) 
  { 
    val[0] = v1; val[1] = v2, val[2] = v3; 
  }

  scalar& operator[] (int idx) 
  { 
    assert(idx >= 0 && idx < 3);
    return val[idx];
  }
};

typedef unsigned long long int uint64;


// Other Hermes macros and constants.
#define HERMES_IS_INT(x) ((int) (x) == (x))
#define countof(a) (sizeof(a)/sizeof(a[0]))

// Pi.
#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

const std::string HERMES_ANY = "-1234";
// For internal use.
const int HERMES_ANY_INT = -1234;
/// This defines the edge types used by discontinuous Galerkin weak forms.
const std::string H2D_DG_BOUNDARY_EDGE = "-12345";  ///< This is to be used by weak forms on the boundary. 
                                    ///< It complements H2D_ANY in that it ensures the forms are evaluated also on non-natural
                                    ///< boundaries (essential conditions may be enforced weakly in some DG methods).
const std::string H2D_DG_INNER_EDGE = "-1234567";    ///< This is to be used by weak forms specifying numerical flux through interior edges.
                                    ///< Forms with this identifier will receive DiscontinuousFunc representations of shape
                                    ///< and ext. functions, which they may query for values on either side of given interface.

// For internal use.
const int H2D_DG_INNER_EDGE_INT = -1234567;
const int H2D_DG_BOUNDARY_EDGE_INT = -12345;

const int HERMES_DIRICHLET_DOF = -1; // Dirichlet lift is a special DOF with number -1.

// For internal use (within Geom<Ord>).
const int HERMES_DUMMY_ELEM_MARKER = -9999;
const int HERMES_DUMMY_EDGE_MARKER = -8888;

/// This class makes command line arguments available to any other method in Hermes.
class HERMES_API CommandLineArgs
{
  public:  
    CommandLineArgs() {};

    int m_argc;
    char** m_argv;
  
    void set(int argc_in, char** argv_in) { 
      m_argc = argc_in; 
      m_argv = argv_in; 
    }
    bool check() { 
      return (m_argc > 0);
    }
    void missing_error() {
      error("Command line arguments have not been set."); 
    }
    int& get_argc() { 
      return m_argc; 
    }
    char**& get_argv() { 
      return m_argv; 
    }
};

/* python support */
/// Throws an exception std::runtime_error. Used by Python wrappers.
/** \param[in] text A text (a cause) of the exception. */
extern HERMES_API void throw_exception(char *text);


/// Int types handling.
#ifdef JU_WIN
typedef __int8           int8_t;
typedef __int16          int16_t;
typedef __int32          int32_t;
typedef __int64          int64_t;

typedef unsigned __int8  uint8_t;
typedef unsigned __int16 uint16_t;
typedef unsigned __int32 uint32_t;
typedef unsigned __int64 uint64_t;

#else
#include <inttypes.h>
#endif

#endif

