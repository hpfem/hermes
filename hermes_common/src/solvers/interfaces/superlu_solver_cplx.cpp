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
/*! \file superlu_solver.cpp
\brief SuperLU solver interface.
*/
#include "config.h"
#ifdef WITH_SUPERLU
#include "superlu_solver.h"
#include "callstack.h"
#include <slu_zdefs.h>

namespace Hermes
{
  namespace Solvers
  {
    template <>
    void SuperLUSolver<std::complex<double> >::create_csc_matrix (SuperMatrix *A, int m, int n, int nnz,
      SuperLuType<std::complex<double> >::Scalar *nzval,
      int *rowind, int *colptr, Stype_t stype, Dtype_t dtype, Mtype_t mtype)
    {
      zCreate_CompCol_Matrix (A, m, n, nnz, (doublecomplex*) nzval, rowind, colptr, stype, dtype, mtype);
    }

    template<>
    void SuperLUSolver<std::complex<double> >::create_dense_matrix (SuperMatrix *X, int m, int n, SuperLuType<std::complex<double> >::Scalar *x,
      int ldx, Stype_t stype, Dtype_t dtype, Mtype_t mtype)
    {
      zCreate_Dense_Matrix (X, m, n, (doublecomplex *) x, ldx, stype, dtype, mtype);
    }

    template <>
    void  SuperLUSolver<std::complex<double> >::solver_driver (superlu_options_t *options, SuperMatrix *A, int *perm_c, int *perm_r, int *etree, char *equed, double *R,
      double *C, SuperMatrix *L, SuperMatrix *U, void *work, int lwork, SuperMatrix *B, SuperMatrix *X,
      double *recip_pivot_growth, double *rcond, double *ferr, double *berr, slu_memusage_t *mem_usage, SuperLUStat_t *stat,
      int *info)
    {
      zgssvx(options, A, perm_c, perm_r, etree, equed, R, C, L, U, work, lwork, B, X, recip_pivot_growth, rcond, ferr, berr, (mem_usage_t*) mem_usage, stat, info);
    }
  }
}
#endif