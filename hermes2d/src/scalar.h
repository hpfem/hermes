// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef __HERMES_COMMON_SCALAR_H
#define __HERMES_COMMON_SCALAR_H

#include <complex>
typedef std::complex<double> cplx;
typedef cplx complex2[2];

#ifdef H2D_COMPLEX
  typedef cplx scalar;
#else
  typedef double scalar;
#endif

#endif
