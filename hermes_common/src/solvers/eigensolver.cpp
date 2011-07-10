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

#include "eigensolver.h"
#include "umfpack_solver.h"

#ifdef WITH_PYTHON

using Teuchos::RCP;
using Teuchos::Ptr;
using Teuchos::rcp;
using Teuchos::ptr;
using Teuchos::null;
using Teuchos::rcp_dynamic_cast;

namespace Hermes 
{
  PyMODINIT_FUNC initeigen(void); /* proto*/

  template<typename Scalar>
  EigenSolver<Scalar>::EigenSolver(const RCP<Matrix<Scalar> > &A, const RCP<Matrix<Scalar> > &B) 
  {
    this->A = A;
    this->B = B;
    this->n_eigs=0;

    initeigen();
  }

  template<typename Scalar>
  static void wrap_CSC(const Ptr<Python> p, const std::string name,
    const RCP<CSCMatrix<double> > A)
  {
    p->push_numpy_int_inplace("_IA", A->get_Ai(), A->get_nnz());
    p->push_numpy_int_inplace("_JA", A->get_Ap(), A->get_size()+1);
    p->push_numpy_double_inplace("_A", A->get_Ax(), A->get_nnz());
    p->push_int("n", A->get_size());
    p->exec("from scipy.sparse import csc_matrix\n");
    p->exec(name + " = csc_matrix((_A, _IA, _JA), shape=(n, n))");
  }

  template<typename Scalar>
  static void wrap_CSC(const Ptr<Python> p, const std::string name,
    const RCP<CSCMatrix<std::complex<double> > > A)
  {
    p->push_numpy_int_inplace("_IA", A->get_Ai(), A->get_nnz());
    p->push_numpy_int_inplace("_JA", A->get_Ap(), A->get_size()+1);
    throw std::runtime_error("Eigenproblem with complex numbers is not supported.");
    p->push_int("n", A->get_size());
    p->exec("from scipy.sparse import csc_matrix\n");
    p->exec(name + " = csc_matrix((_A, _IA, _JA), shape=(n, n))");
  }

  template<typename Scalar>
  void EigenSolver<Scalar>::solve(int n_eigs, double target_value, double tol,
    int max_iter) 
  {
    // Support CSCMatrix<Scalar> only for now:
    RCP<CSCMatrix<Scalar> > A = rcp_dynamic_cast<CSCMatrix<Scalar> >(this->A, true);
    RCP<CSCMatrix<Scalar> > B = rcp_dynamic_cast<CSCMatrix<Scalar> >(this->B, true);
    wrap_CSC<Scalar>(ptr(&p), "A", A);
    wrap_CSC<Scalar>(ptr(&p), "B", B);
    this->p.exec("from eigen import solve_eig_pysparse");
    this->p.push_double("target_value", target_value);
    this->p.push_int("n_eigs", n_eigs);
    this->p.push_double("jdtol", tol);
    this->p.push_int("max_iter", max_iter);

    printf("Solving the system A * x = lambda * B * x\n");
    this->p.exec("eigs = solve_eig_pysparse(A, B, target_value=target_value, n_eigs=n_eigs, jdtol=jdtol, max_iter=max_iter)");
    this->p.exec("n_eigs = len(eigs)");
    this->n_eigs = this->p.pull_int("n_eigs");
  }

  template<typename Scalar>
  double EigenSolver<Scalar>::get_eigenvalue(int i)
  {
    if (i >= 0 && i < this->n_eigs) 
    {
      this->p.push_int("i", i);
      this->p.exec("E = eigs[i][0]");
      return this->p.pull_double("E");
    } else
      throw std::runtime_error("'i' must obey 0 <= i < n_eigs");
  }

  template<typename Scalar>
  void EigenSolver<Scalar>::get_eigenvector(int i, double **vec, int *n)
  {
    if (i >= 0 && i < this->n_eigs) 
    {
      this->p.push_int("i", i);
      this->p.exec("vec = eigs[i][1]");
      this->p.pull_numpy_double_inplace("vec", vec, n);
    } else
      throw std::runtime_error("'i' must obey 0 <= i < n_eigs");
  }

  template class HERMES_API EigenSolver<double>;
  template class HERMES_API EigenSolver<std::complex<double> >;
}

#endif
