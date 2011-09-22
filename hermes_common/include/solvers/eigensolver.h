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
/*! \file eigensolver.h
    \brief class for solving Eigenproblems. Currently using Python.
*/
#ifndef __HERMES_EIGENSOLVER_H
#define __HERMES_EIGENSOLVER_H

#include "config.h"
#ifdef WITH_PYTHON

#include "matrix.h"

#include "python_API/python_api.h"

// RCP
#ifndef WITH_TRILINOS
#include "Teuchos_RCP.hpp"
#else
#include "Teuchos_RCP.hpp"
#endif

namespace Hermes
{

using Teuchos::RCP;
using Teuchos::Ptr;
using Teuchos::rcp;
using Teuchos::null;


template <typename Scalar>
class HERMES_API EigenSolver
{
public:
    EigenSolver(const RCP<Algebra::Matrix<Scalar> > &A, const RCP<Algebra::Matrix<Scalar> > &B);


    // Solves for 'n_eigs' eigenvectors, around the 'target_value'. Use
    // 'get_eigenvalue' and 'get_eigenvector' to retrieve the
    // eigenvalues/eigenvectors:
    void solve(int n_eigs = 4, double target_value = -1, double tol = 1e-6,
            int max_iter = 150);

    // Returns the number of calculated eigenvalues
    int get_n_eigs()
    {
        return this->n_eigs;
    }
    // Returns the i-th eigenvalue
    double get_eigenvalue(int i);
    // Returns the i-th eigenvector. A pointer will be returned into an
    // internal array, as well as the size of the vector. You don't own the
    // memory and it will be deallocated once the EigenSolver() class is
    // deleted. You need to make a copy of it if you want to store it
    // permanently.
    void get_eigenvector(int i, double **vec, int *n);

    void print_eigenvalues()
    {
        printf("Eigenvalues:\n");
        for (int i = 0; i < this->get_n_eigs(); i++)
            printf("%3d: %f\n", i, this->get_eigenvalue(i));
    }

private:
    RCP<Algebra::Matrix<Scalar> > A, B;
    int n_eigs;
    Python p;
};

}

#endif

#endif
