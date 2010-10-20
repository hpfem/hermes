// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "matrix.h"
#include "solvers.h"

#include <coord_double.h>
#include <compcol_double.h>
#include <mvvd.h>
#include <ilupre_double.h>
#include <bicg.h>
#include <cg.h>
#include <cgs.h>
#include <bicgstab.h>
#include <cheby.h>
#include <gmres.h>
#include <ir.h>
#include <qmr.h>

bool CommonSolverSparseLib::_solve(Matrix *mat, double *res)
{
    printf("SparseLib++ solver\n");

    CSCMatrix *Acsc = NULL;

    if (CooMatrix *mcoo = dynamic_cast<CooMatrix*>(mat))
        Acsc = new CSCMatrix(mcoo);
    else if (DenseMatrix *mden = dynamic_cast<DenseMatrix *>(mat))
        Acsc = new CSCMatrix(mden);
    else if (CSCMatrix *mcsc = dynamic_cast<CSCMatrix*>(mat))
        Acsc = mcsc;
    else if (CSRMatrix *mcsr = dynamic_cast<CSRMatrix*>(mat))
        Acsc = new CSCMatrix(mcsr);
    else
        _error("Matrix type not supported.");

    int nnz = Acsc->get_nnz();
    int size = Acsc->get_size();

    CompCol_Mat_double Acc = CompCol_Mat_double(size, size, nnz,
                                                Acsc->get_Ax(), Acsc->get_Ai(), Acsc->get_Ap());

    // rhs
    VECTOR_double rhs(res, size);

    // preconditioner
    CompCol_ILUPreconditioner_double ILU(Acc);
    VECTOR_double xv = ILU.solve(rhs);

    // method
    int result = -1;
    switch (method)
    {
    case HERMES_CommonSolverSparseLibSolver_ConjugateGradientSquared:
        result = CGS(Acc, xv, rhs, ILU, maxiter, tolerance);
        break;
    case CommonSolverSparseLibSolver_RichardsonIterativeRefinement:
        result = IR(Acc, xv, rhs, ILU, maxiter, tolerance);
        break;
    default:
        _error("SparseLib++ error. Method is not defined.");
    }

    if (result == 0)
        printf("SparseLib++ solver: maxiter: %i, tol: %e\n", maxiter, tolerance);
    else
        _error("SparseLib++ error.");

    double *x;
    x = (double*) malloc(size * sizeof(double));

    for (int i = 0 ; i < xv.size() ; i++)
        x[i] = xv(i);

    memcpy(res, x, size*sizeof(double));
    delete[] x;

    if (!dynamic_cast<CSCMatrix*>(mat))
        delete Acsc;
}

bool CommonSolverSparseLib::_solve(Matrix *mat, cplx *res)
{
    _error("CommonSolverSparseLib::solve(Matrix *mat, cplx *res) not implemented.");
}
