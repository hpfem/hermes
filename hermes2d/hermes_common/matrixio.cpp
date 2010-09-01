// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/
// Harwell-Boeing format

#include "matrixio.h"
#include "iohb.h"

CSCMatrix *read_hb_csc(const char *filename)
{
    int M, N, nnz, nrhs;
    readHB_info(filename, &M, &N, &nnz, &nrhs);

    // square matrix
    int size = N;

    // allocate matrix
    int *Ap = new int[size + 1];
    int *Ai = new int[nnz];
    double *Ax = new double[nnz];

    readHB_mat_double(filename, Ap, Ai, Ax);
    CSCMatrix *Acsc = new CSCMatrix(size, nnz, Ap, Ai, Ax);

    return Acsc;
}

CSRMatrix *read_hb_csr(const char *filename)
{
    CSCMatrix *Acsc = read_hb_csc(filename);
    CSRMatrix *Acsr = new CSRMatrix(Acsc);
    delete Acsc;

    return Acsr;
}

CooMatrix *read_hb_coo(const char *filename)
{
    CSCMatrix *Acsc = read_hb_csc(filename);
    Acsc->print();
    CooMatrix *Acoo = new CooMatrix(Acsc);
    delete Acsc;

    return Acoo;
}

double *read_rhs(const char *filename, int j)
{
    int M, N, nnz, nrhs;
    readHB_info(filename, &M, &N, &nnz, &nrhs);

    double *rhs = new double[nnz];

    if (j >= 0 && j < nrhs)
        readHB_rhs_double(filename, rhs, j);
    else
        printf("Error: HB rhs #%d in file '%s' not found.\n", j, filename);

    return rhs;
}

void write_hb_csc(const char *filename, CSCMatrix *A, double *rhs)
{
    const char *title = "Hermes Common";
    const char *key = "Matrix";

    writeHB_mat_double(filename, A->get_size(), A->get_size(), A->get_nnz(),
        A->get_Ap(), A->get_Ai(), A->get_Ax(), 1, rhs, title, key);
}

void write_hb_csr(const char *filename, CSRMatrix *A, double *rhs)
{
    CSCMatrix *Acsc = new CSCMatrix(A);

    write_hb_csc(filename, Acsc, rhs);
}

void write_hb_coo(const char *filename, CooMatrix *A, double *rhs)
{
    CSCMatrix *Acsc = new CSCMatrix(A);

    write_hb_csc(filename, Acsc, rhs);
}
