// This file is part of HermesCommon.
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

/*! \file common.h
    \brief File containing common definitions, and basic global enums etc. for HermesCommon.
*/
#ifndef __HERMES_COMMON_COMMON_H
#define __HERMES_COMMON_COMMON_H

#include <complex>

#include <stdexcept>
#include <cstdarg>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <errno.h>
#include <cmath>

#include <pthread.h>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <map>
#include <cstdio>
#include <stdarg.h>
#include <sstream>

#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <queue>
#include <sstream>
#include <fstream>
#include <cstring>
#include <iostream>
#include <omp.h>
#include <signal.h>

#include "config.h"

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

template<typename Scalar>
class Scalar2
{
  Scalar val[2];

public:
  Scalar2(Scalar v1, Scalar v2)
  {
    val[0] = v1; val[1] = v2;
  }

  Scalar& operator[] (int idx)
  {
    assert(idx >= 0 && idx < 2);
    return val[idx];
  }
};

template<typename Scalar>
class Scalar3
{
  Scalar val[3];

public:
  Scalar3(Scalar v1, Scalar v2, Scalar v3)
  {
    val[0] = v1; val[1] = v2, val[2] = v3;
  }

  Scalar& operator[] (int idx)
  {
    assert(idx >= 0 && idx < 3);
    return val[idx];
  }
};

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

// Pi.
#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

namespace Hermes
{
  inline double sqr(int x) { return x*x; }
  inline double sqr(double x) { return x*x; }
  inline double sqrt(double x) { return std::sqrt(x); }
  inline double sqr(std::complex<double> x)   { return std::norm(x); }
  inline double magn(double x) { return fabs(x); }
  inline double conj(double a) {  return a; }
  inline double cos(double x) { return std::cos(x); }
  inline double sin(double x) { return std::sin(x); }
  inline double atan(double x) { return std::atan(x); }
  inline double atan2(double x, double y) { return std::atan2(x, y); }
  inline double abs(double x) { return std::abs(x); }
  inline double pow(double x, double y) { return std::pow(x, y); }
  inline double log(double x) { return std::log(x); }

  // Represents "any" part of the boundary when deciding where (on which elements) to assemble the form at hand.
  const std::string HERMES_ANY = "-1234";
  // For internal use.
  const int HERMES_ANY_INT = -1234;
  /// This is to be used by weak forms specifying numerical flux through interior edges.
  /// Forms with this identifier will receive DiscontinuousFunc representations of shape
  /// and ext. functions, which they may query for values on either side of given interface.
  const std::string H2D_DG_INNER_EDGE = "-1234567";    
  // For internal use.
  const int H2D_DG_INNER_EDGE_INT = -1234567;

  namespace Helpers
  {
    inline void fprint_num(FILE*f, double x)
    {
      fprintf(f, "%lf", x);
    }

    inline void fprint_num(FILE*f, std::complex<double> x)
    {
      fprintf(f, "(%lf, %lf)", x.real(), x.imag());
    }
  }

  namespace BLAS
  {
#ifdef WITH_BLAS
    // BLAS-related functions

    // Complex part.
#ifdef __cplusplus
    extern "C" 
    {
#endif
      extern int zscal_(int *, std::complex<double> *, std::complex<double> *, int *);
      extern int zaxpy_(int *, std::complex<double> *, std::complex<double> *, int *, std::complex<double> *, int *);
      extern int zcopy_(int *, std::complex<double> *, int *, std::complex<double> *, int *);

      extern int dscal_(int *, double *, double *, int *);
      extern int daxpy_(int *, double *, double *, int *, double *, int *);
      extern int dcopy_(int *,           double *, int *, double *, int *);

#ifdef __cplusplus
    }
#endif

    /// x <- alpha * x
    inline void blas_scal(int n, std::complex<double> alpha, std::complex<double> *x, int incx) { zscal_(&n, &alpha, x, &incx); }
    /// y <- alpha * x + y
    inline void blas_axpy(int n, std::complex<double> alpha, std::complex<double> *x, int incx, std::complex<double> *y, int incy) { zaxpy_(&n, &alpha, x, &incx, y, &incy); }
    /// y <- x
    inline void blas_copy(int n, std::complex<double> *x, int incx, std::complex<double> *y, int incy) { zcopy_(&n, x, &incx, y, &incx); }

    /// x <- alpha * x
    inline void blas_scal(int n, double alpha, double *x, int incx) { dscal_(&n, &alpha, x, &incx); }
    /// y <- alpha * x + y
    inline void blas_axpy(int n, double alpha, double *x, int incx, double *y, int incy) { daxpy_(&n, &alpha, x, &incx, y, &incy); }
    /// y <- x
    inline void blas_copy(int n, double *x, int incx, double *y, int incy) { dcopy_(&n, x, &incx, y, &incx); }

#endif
  }
}
#endif

/**

\mainpage

This manual documents the source code of hermes_common. It is intended for the developers of
the library. If you are only interested in using hermes_common together with hermesNd, please refer to the User's Manual.

*/
