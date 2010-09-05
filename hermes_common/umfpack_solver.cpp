// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "matrix.h"
#include "solvers.h"

#ifdef COMMON_WITH_UMFPACK
#include <umfpack.h>

double control_array[UMFPACK_CONTROL];
double info_array[UMFPACK_INFO];

void print_status(int status)
{
    switch (status)
    {
    case UMFPACK_OK:
        break;
    case UMFPACK_WARNING_singular_matrix:
        _error("UMFPACK: singular stiffness matrix!");
        break;

    case UMFPACK_ERROR_out_of_memory:           _error("UMFPACK: out of memory!");
    case UMFPACK_ERROR_argument_missing:        _error("UMFPACK: argument missing");
    case UMFPACK_ERROR_invalid_Symbolic_object: _error("UMFPACK: invalid Symbolic object");
    case UMFPACK_ERROR_invalid_Numeric_object:  _error("UMFPACK: invalid Numeric object");
    case UMFPACK_ERROR_different_pattern:       _error("UMFPACK: different pattern");
    case UMFPACK_ERROR_invalid_system:          _error("UMFPACK: invalid system");
    case UMFPACK_ERROR_n_nonpositive:           _error("UMFPACK: n nonpositive");
    case UMFPACK_ERROR_invalid_matrix:          _error("UMFPACK: invalid matrix");
    case UMFPACK_ERROR_internal_error:          _error("UMFPACK: internal error");
    default:                                    _error("UMFPACK: unknown error");
    }
}

bool CommonSolverUmfpack::_solve(Matrix *mat, double *res)
{
    printf("UMFPACK solver\n");

    CSCMatrix *Acsc = NULL;

    if (CooMatrix *mcoo = dynamic_cast<CooMatrix*>(mat))
        Acsc = new CSCMatrix(mcoo);
    else if (CSCMatrix *mcsc = dynamic_cast<CSCMatrix*>(mat))
        Acsc = mcsc;
    else if (CSRMatrix *mcsr = dynamic_cast<CSRMatrix*>(mat))
        Acsc = new CSCMatrix(mcsr);
    else
        _error("Matrix type not supported.");

    int nnz = Acsc->get_nnz();
    int size = Acsc->get_size();

    // solve
    umfpack_di_defaults(control_array);

    /* symbolic analysis */
    void *symbolic, *numeric;
    int status_symbolic = umfpack_di_symbolic(size, size,
                                              Acsc->get_Ap(), Acsc->get_Ai(), NULL, &symbolic,
                                              control_array, info_array);
    print_status(status_symbolic);

    /* LU factorization */
    int status_numeric = umfpack_di_numeric(Acsc->get_Ap(), Acsc->get_Ai(), Acsc->get_Ax(), symbolic, &numeric,
                                            control_array, info_array);
    print_status(status_numeric);

    umfpack_di_free_symbolic(&symbolic);

    double *x;
    x = (double*) malloc(size * sizeof(double));

    /* solve system */
    int status_solve = umfpack_di_solve(UMFPACK_A,
                                        Acsc->get_Ap(), Acsc->get_Ai(), Acsc->get_Ax(), x, res, numeric,
                                        control_array, info_array);

    print_status(status_solve);

    umfpack_di_free_numeric(&numeric);

    if (symbolic) umfpack_di_free_symbolic(&symbolic);
    if (numeric) umfpack_di_free_numeric(&numeric);

    memcpy(res, x, size*sizeof(double));
    delete[] x;

    if (!dynamic_cast<CSCMatrix*>(mat))
        delete Acsc;
}

bool CommonSolverUmfpack::_solve(Matrix *mat, cplx *res)
{
    printf("UMFPACK solver - cplx\n");

    CSCMatrix *Acsc = NULL;

    if (CooMatrix *mcoo = dynamic_cast<CooMatrix*>(mat))
        Acsc = new CSCMatrix(mcoo);
    else if (CSCMatrix *mcsc = dynamic_cast<CSCMatrix*>(mat))
        Acsc = mcsc;
    else if (CSRMatrix *mcsr = dynamic_cast<CSRMatrix*>(mat))
        Acsc = new CSCMatrix(mcsr);
    else
        _error("Matrix type not supported.");

    int nnz = Acsc->get_nnz();
    int size = Acsc->get_size();

    // complex components
    double *Axr = new double[nnz];
    double *Axi = new double[nnz];
    cplx *Ax = Acsc->get_Ax_cplx();
    for (int i = 0; i < nnz; i++)
    {
        Axr[i] = Ax[i].real();
        Axi[i] = Ax[i].imag();
    }

    umfpack_zi_defaults(control_array);

    /* symbolic analysis */
    void *symbolic, *numeric;
    int status_symbolic = umfpack_zi_symbolic(size, size,
                                              Acsc->get_Ap(), Acsc->get_Ai(), NULL, NULL, &symbolic,
                                              control_array, info_array);
    print_status(status_symbolic);

    /* LU factorization */
    int status_numeric = umfpack_zi_numeric(Acsc->get_Ap(), Acsc->get_Ai(), Axr, Axi, symbolic, &numeric,
                                            control_array, info_array);
    print_status(status_numeric);

    umfpack_zi_free_symbolic(&symbolic);

    double *xr, *xi;
    xr = (double*) malloc(size * sizeof(double));
    xi = (double*) malloc(size * sizeof(double));

    double *resr = new double[size];
    double *resi = new double[size];

    for (int i = 0; i < size; i++)
    {
        resr[i] = res[i].real();
        resi[i] = res[i].imag();
    }

    /* solve system */
    int status_solve = umfpack_zi_solve(UMFPACK_A,
                                        Acsc->get_Ap(), Acsc->get_Ai(), Axr, Axi, xr, xi, resr, resi, numeric,
                                        control_array, info_array);

    print_status(status_solve);

    umfpack_zi_free_numeric(&numeric);

    delete[] resr;
    delete[] resi;
    delete[] Axr;
    delete[] Axi;

    if (symbolic) umfpack_di_free_symbolic(&symbolic);
    if (numeric) umfpack_di_free_numeric(&numeric);

    for (int i = 0; i < Acsc->get_size(); i++)
        res[i] = cplx(xr[i], xi[i]);

    delete xr;
    delete xi;

    if (!dynamic_cast<CSCMatrix*>(mat))
        delete Acsc;
}

#else

bool CommonSolverUmfpack::_solve(Matrix *mat, double *res)
{
    _error("CommonSolverUmfpack::solve(Matrix *mat, double *res) not implemented.");
}

bool CommonSolverUmfpack::_solve(Matrix *mat, cplx *res)
{
    _error("CommonSolverUmfpack::solve(Matrix *mat, cplx *res) not implemented.");
}

#endif
