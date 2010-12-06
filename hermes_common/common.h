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

#ifndef __HERMES_COMMON_H
#define __HERMES_COMMON_H

// Include
//
// common headers
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
//.
#include "tuple.h"

#define HERMES  "Hermes"

// Decide which version of Hermes is being compiled and import
// the build options from the corresponding config.h file.
#if defined(H1D_REAL) || defined(H1D_COMPLEX)
  #ifndef CONFIG_H_INCLUDED
    #include "../hermes1d/src/config.h"
    #define CONFIG_H_INCLUDED
  #endif
#elif defined(H2D_REAL) || defined(H2D_COMPLEX)
  #ifndef CONFIG_H_INCLUDED
    #include "../hermes2d/src/config.h"
    #define CONFIG_H_INCLUDED
  #endif
#elif defined(H3D_REAL) || defined(H3D_COMPLEX)
  #ifndef CONFIG_H_INCLUDED
    #include "../hermes3d/src/config.h"
    #define CONFIG_H_INCLUDED
  #endif
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
   SOLVER_PARDISO,
   SOLVER_SUPERLU,
   SOLVER_AMESOS,
   SOLVER_AZTECOO
};

// Should be in the same order as MatrixSolverTypes above, so that the
// names may be accessed by the same enumeration variable.
const std::string MatrixSolverNames[7] = {
  "UMFPACK",
  "PETSc",
  "MUMPS",
  "Pardiso",
  "SuperLU",
  "Trilinos/Amesos",
  "Trilinos/AztecOO"
};

#define UMFPACK_NOT_COMPILED  HERMES " was not built with UMFPACK support."
#define PETSC_NOT_COMPILED    HERMES " was not built with PETSC support."
#define MUMPS_NOT_COMPILED    HERMES " was not built with MUMPS support."
#define PARDISO_NOT_COMPILED  HERMES " was not built with PARDISO support."
#define SUPERLU_NOT_COMPILED  HERMES " was not built with SUPERLU support."
#define NOX_NOT_COMPILED      HERMES " was not built with NOX support."
#define AMESOS_NOT_COMPILED   HERMES " was not built with AMESOS support."
#define AZTECOO_NOT_COMPILED  HERMES " was not built with AZTECOO support."
#define EPETRA_NOT_COMPILED   HERMES " was not built with EPETRA support."
#define IFPACK_NOT_COMPILED   HERMES " was not built with IFPACK support."
#define ML_NOT_COMPILED       HERMES " was not built with ML support."


// Projection norms:
enum ProjNormType
{
  HERMES_L2_NORM, 
  HERMES_H1_NORM, 
  HERMES_H1_SEMINORM, 
  HERMES_HCURL_NORM, 
  HERMES_HDIV_NORM
};

// Default HERMES projection norm is the H1 norm.
const ProjNormType HERMES_DEFAULT_PROJ_NORM = HERMES_H1_NORM;

inline int sqr(int x) { return x*x; }
inline double sqr(double x) { return x*x; }
inline double magn(double x) { return fabs(x); }
inline double conj(double a) {  return a; }

#if defined(H2D_COMPLEX) || defined(H3D_COMPLEX)

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

typedef scalar scalar2[2];
typedef scalar scalar3[3];

typedef unsigned long long int uint64;


// Other Hermes macros and constants.
#define HERMES_IS_INT(x) ((int) (x) == (x))
#define countof(a) (sizeof(a)/sizeof(a[0]))

// Pi.
#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

const int HERMES_ANY = -1234;
/// This defines the edge types used by discontinuous Galerkin weak forms.
enum DG_EdgeType
{
	H2D_DG_BOUNDARY_EDGE = -12345,  ///< This is to be used by weak forms on the boundary. 
                                    ///< It complements H2D_ANY in that it ensures the forms are evaluated also on non-natural
                                    ///< boundaries (essential conditions may be enforced weakly in some DG methods).
	H2D_DG_INNER_EDGE = -1234567    ///< This is to be used by weak forms specifying numerical flux through interior edges.
                                    ///< Forms with this identifier will receive DiscontinuousFunc representations of shape
                                    ///< and ext. functions, which they may query for values on either side of given interface.
};

const int HERMES_DIRICHLET_DOF = -1; // Dirichlet lift is a special DOF with number -1.

/// This class makes command line arguments available to any other method in Hermes.
class HERMES_API CommandLineArgs
{
  static int m_argc;
  static char** m_argv;
  
  CommandLineArgs() {};
  
  public:  
    static void set(int argc_in, char** argv_in) { 
      m_argc = argc_in; 
      m_argv = argv_in; 
    }
    static bool check() { 
      return (m_argc > 0);
    }
    static void missing_error() {
      error("Command line arguments have not been set."); 
    }
    static int& get_argc() { 
      return m_argc; 
    }
    static char**& get_argv() { 
      return m_argv; 
    }
};

#endif

