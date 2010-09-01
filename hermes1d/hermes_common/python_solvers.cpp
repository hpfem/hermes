// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "matrix.h"
#include "solvers.h"

#ifdef COMMON_WITH_SCIPY
#include "python_api.h"

bool CommonSolverNumPy::_solve(Matrix *mat, double *res)
{
  //printf("NumPy solver\n");

    CSRMatrix M(mat);
    Python *p = new Python();
    p->push("m", c2py_CSRMatrix(&M));
    p->push("rhs", c2numpy_double_inplace(res, mat->get_size()));
    p->exec("A = m.to_scipy_csr().todense()");
    p->exec("from numpy.linalg import solve");
    p->exec("x = solve(A, rhs)");
    double *x;
    int n;
    numpy2c_double_inplace(p->pull("x"), &x, &n);
    memcpy(res, x, n*sizeof(double));
    delete p;
		return true;
}

bool CommonSolverNumPy::_solve(Matrix *mat, cplx *res)
{
  //printf("NumPy solver - cplx\n");

    CSRMatrix M(mat);
    Python *p = new Python();
    p->push("m", c2py_CSRMatrix(&M));
    p->push("rhs", c2numpy_double_complex_inplace(res, mat->get_size()));
    p->exec("A = m.to_scipy_csr().todense()");
    p->exec("from numpy.linalg import solve");
    p->exec("x = solve(A, rhs)");
    cplx *x;
    int n;
    numpy2c_double_complex_inplace(p->pull("x"), &x, &n);
    memcpy(res, x, n*sizeof(cplx));
    delete p;
		return true;
}

bool CommonSolverSciPyUmfpack::_solve(Matrix *mat, double *res)
{
  //printf("SciPy UMFPACK solver\n");

    CSCMatrix M(mat);
    Python *p = new Python();
    p->push("m", c2py_CSCMatrix(&M));
    p->push("rhs", c2numpy_double_inplace(res, mat->get_size()));
    p->exec("A = m.to_scipy_csc()");
    p->exec("from scipy.sparse.linalg import spsolve");
    // Turn off warnings in spsolve (only there)
    p->exec("import warnings");
    p->exec("with warnings.catch_warnings():\n"
            "    warnings.simplefilter('ignore')\n"
            "    x = spsolve(A, rhs)");
    double *x;
    int n;
    numpy2c_double_inplace(p->pull("x"), &x, &n);
    memcpy(res, x, n*sizeof(double));
    delete p;
    return true;
}

bool CommonSolverSciPyUmfpack::_solve(Matrix *mat, cplx *res)
{
  //printf("SciPy UMFPACK solver - cplx\n");

    CSCMatrix M(mat);
    Python *p = new Python();
    p->push("m", c2py_CSCMatrix(&M));
    p->push("rhs", c2numpy_double_complex_inplace(res, mat->get_size()));
    p->exec("A = m.to_scipy_csc()");
    p->exec("from scipy.sparse.linalg import spsolve");
    // Turn off warnings in spsolve (only there)
    p->exec("import warnings");
    p->exec("with warnings.catch_warnings():\n"
            "    warnings.simplefilter('ignore')\n"
            "    x = spsolve(A, rhs)");
    cplx *x;
    int n;
    numpy2c_double_complex_inplace(p->pull("x"), &x, &n);
    memcpy(res, x, n*sizeof(cplx));
    delete p;
    return true;
}

bool CommonSolverSciPyCG::_solve(Matrix *mat, double *res)
{
  //printf("SciPy CG solver\n");

    CSRMatrix M(mat);
    Python *p = new Python();
    p->push("m", c2py_CSRMatrix(&M));
    p->push("rhs", c2numpy_double_inplace(res, mat->get_size()));
    p->exec("A = m.to_scipy_csr()");
    p->exec("from scipy.sparse.linalg import cg");
    p->exec("x, res = cg(A, rhs)");
    double *x;
    int n;
    numpy2c_double_inplace(p->pull("x"), &x, &n);
    memcpy(res, x, n*sizeof(double));
    delete p;
		return true;
}

bool CommonSolverSciPyCG::_solve(Matrix *mat, cplx *res)
{
    _error("CommonSolverSciPyCG::solve(Matrix *mat, cplx *res) not implemented.");
}

bool CommonSolverSciPyGMRES::_solve(Matrix *mat, double *res)
{
  //printf("SciPy GMRES solver\n");

    CSRMatrix M(mat);
    Python *p = new Python();
    p->push("m", c2py_CSRMatrix(&M));
    p->push("rhs", c2numpy_double_inplace(res, mat->get_size()));
    p->exec("A = m.to_scipy_csr()");
    p->exec("from scipy.sparse.linalg import gmres");
    p->exec("x, res = gmres(A, rhs)");
    double *x;
    int n;
    numpy2c_double_inplace(p->pull("x"), &x, &n);
    memcpy(res, x, n*sizeof(double));
    delete p;
		return true;
}

bool CommonSolverSciPyGMRES::_solve(Matrix *mat, cplx *res)
{
    _error("CommonSolverSciPyGMRES::solve(Matrix *mat, cplx *res) not implemented.");
}

#else

bool CommonSolverNumPy::_solve(Matrix *mat, double *res)
{
    _error("CommonSolverNumPy::solve(Matrix *mat, double *res) not implemented.");
}

bool CommonSolverNumPy::_solve(Matrix *mat, cplx *res)
{
    _error("CommonSolverNumPy::solve(Matrix *mat, cplx *res) not implemented.");
}

bool CommonSolverSciPyUmfpack::_solve(Matrix *mat, double *res)
{
    _error("CommonSolverSciPyUmfpack::solve(Matrix *mat, double *res) not implemented.");
}

bool CommonSolverSciPyUmfpack::_solve(Matrix *mat, cplx *res)
{
    _error("CommonSolverSciPyUmfpack::solve(Matrix *mat, cplx *res) not implemented.");
}

bool CommonSolverSciPyCG::_solve(Matrix *mat, double *res)
{
    _error("CommonSolverSciPyCG::solve(Matrix *mat, double *res) not implemented.");
}

bool CommonSolverSciPyCG::_solve(Matrix *mat, cplx *res)
{
    _error("CommonSolverSciPyCG::solve(Matrix *mat, cplx *res) not implemented.");
}

bool CommonSolverSciPyGMRES::_solve(Matrix *mat, double *res)
{
    _error("CommonSolverSciPyGMRES::solve(Matrix *mat, double *res) not implemented.");
}

bool CommonSolverSciPyGMRES::_solve(Matrix *mat, cplx *res)
{
    _error("CommonSolverSciPyGMRES::solve(Matrix *mat, cplx *res) not implemented.");
}


#endif
