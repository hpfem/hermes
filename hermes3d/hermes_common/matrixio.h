// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "matrix.h"

// Harwell-Boeing format
CSCMatrix *read_hb_csc(const char *filename);
CSRMatrix *read_hb_csr(const char *filename);
CooMatrix *read_hb_coo(const char *filename);
void write_hb_csc(const char *filename, CSCMatrix *A, double *rhs);
void write_hb_csr(const char *filename, CSRMatrix *A, double *rhs);
void write_hb_coo(const char *filename, CooMatrix *A, double *rhs);

double *read_rhs(const char *filename, int j = 0);
