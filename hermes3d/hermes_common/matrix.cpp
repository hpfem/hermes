// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "matrix.h"

// print vector - int
void print_vector(const char *label, int *value, int size) {
    printf("%s [", label);
    for (int i = 0; i < size; i++) {
        if (i < size-1)
            printf("%i, ", value[i]);
        else
            printf("%i", value[i]);
    }
    printf("]\n");
}

// print vector - double
void print_vector(const char *label, double *value, int size) {
    printf("%s [", label);
    for (int i = 0; i < size; i++) {
        if (i < size-1)
            printf("%f, ", value[i]);
        else
            printf("%f", value[i]);
    }
    printf("]\n");
}

// print vector - cplx
void print_vector(const char *label, cplx *value, int size) {
    printf("%s [", label);
    for (int i = 0; i < size; i++) {
        if (i < size-1)
            printf("(%f, %f), ", value[i].real(), value[i].imag());
        else
            printf("(%f, %f)", value[i].real(), value[i].imag());
    }
    printf("]\n");
}

/// Transposes an m by n matrix. If m != n, the array matrix in fact has to be
/// a square matrix of the size max(m, n) in order for the transpose to fit inside it.
template<typename T>
void transpose(T** matrix, int m, int n)
{
    int min = std::min(m, n);
    for (int i = 0; i < min; i++)
        for (int j = i+1; j < min; j++)
            std::swap(matrix[i][j], matrix[j][i]);

    if (m < n)
        for (int i = 0; i < m; i++)
            for (int j = m; j < n; j++)
                matrix[j][i] = matrix[i][j];
    else if (n < m)
        for (int i = n; i < m; i++)
            for (int j = 0; j < n; j++)
                matrix[j][i] = matrix[i][j];
}

// *********************************************************************************************************************

DenseMatrix::DenseMatrix(Matrix *m)
{
    this->size = m->get_size();
    this->complex = m->is_complex();
    init();

    if (dynamic_cast<CooMatrix *>(m))
        this->add_from_coo((CooMatrix *)m);
    else
        _error("Matrix type not supported.");
}

DenseMatrix::DenseMatrix(CooMatrix *m)
{
    this->size = m->get_size();
    this->complex = m->is_complex();
    init();

    this->add_from_coo(m);
}

DenseMatrix::DenseMatrix(int size, bool is_complex)
{
    this->complex = is_complex;
    this->size = size;

    init();
}

DenseMatrix::~DenseMatrix()
{
    free_data();
}

void DenseMatrix::free_data()
{
    if (this->A != NULL) { 
      delete[] A; 
      this->A = NULL; 
    };
    if (this->A_cplx != NULL) { 
      delete[] this->A_cplx; 
      this->A_cplx = NULL; 
    };

    this->size = 0;
}

void DenseMatrix::init()
{
    this->A = NULL;
    this->A_cplx = NULL;

    if (complex)
        this->A_cplx = _new_matrix<cplx>(this->size, this->size);
    else
        this->A = _new_matrix<double>(this->size, this->size);

    if (complex)
    {
        for (int i = 0; i<size; i++)
            for (int j = 0; j<size; j++)
                this->A_cplx[i][j] = 0;
    }
    else
    {
        for (int i = 0; i<size; i++)
            for (int j = 0; j<size; j++)
                this->A[i][j] = 0;
    }
}

void DenseMatrix::add_from_coo(CooMatrix *m)
{
    int nnz = m->get_nnz();

    // get data
    int *row = new int[nnz];
    int *col = new int[nnz];

    if (complex)
    {
        cplx *data = new cplx[nnz];
        m->get_row_col_data(row, col, data);

        for (int i = 0; i < nnz; i++)
            A_cplx[row[i]][col[i]] = data[i];

        if (data != NULL) delete[] data;
    }
    else
    {
        double *data = new double[nnz];
        m->get_row_col_data(row, col, data);

        for (int i = 0; i < nnz; i++)
            A[row[i]][col[i]] = data[i];

        if (data != NULL) delete[] data;
    }

    // free data
    if (row != NULL) delete[] row;
    if (col != NULL) delete[] col;
}

int DenseMatrix::get_nnz()
{
    int nnz = 0;
    for (int i = 0; i < this->size; i++)
    {
        for (int j = 0; j < this->size; j++)
        {
            if (complex)
                if (fabs(A[i][j]) > 1e-12)
                    nnz++;
            else
                if (std::abs(A_cplx[i][j]) > 1e-12)
                    nnz++;
        }
    }
    return nnz;
}

void DenseMatrix::set_zero()
{
    if (this->complex)
        for (int i = 0; i<size; i++)
            for (int j = 0; j<size; j++)
                this->A_cplx[i][j] = 0;
    else
        for (int i = 0; i<size; i++)
            for (int j = 0; j<size; j++)
                this->A[i][j] = 0;
}

void DenseMatrix::print()
{
    for (int i = 0; i < this->size; i++)
    {
        for (int j = 0; j < this->size; j++)
        {
            if (complex)
                printf("(%i, %i): (%f, %f)\n",
                       i, j, A_cplx[i][j].real(), A_cplx[i][j].imag());
            else
                printf("(%i, %i): %f\n",
                       i, j, A[i][j]);
        }
    }
}

// *********************************************************************************************************************

CooMatrix::CooMatrix(bool complex) : Matrix()
{
    init();

    this->complex = complex;
    this->size = 0;
}

CooMatrix::CooMatrix(int size, bool complex) : Matrix()
{
    init();

    this->complex = complex;
    this->size = size;
}

CooMatrix::CooMatrix(Matrix *m)
{
    init();

    if (dynamic_cast<CSCMatrix*>(m))
        this->add_from_csc((CSCMatrix*)m);
    else if (dynamic_cast<CSRMatrix*>(m))
        this->add_from_csr((CSRMatrix*)m);
    else
        _error("Matrix type not supported.");
}

CooMatrix::CooMatrix(CSRMatrix *m)
{
    init();
    this->add_from_csr(m);
}

CooMatrix::CooMatrix(CSCMatrix *m)
{
    init();
    this->add_from_csc(m);
}

CooMatrix::~CooMatrix()
{
    this->free_data();
}
void CooMatrix::set_zero()
{
    this->free_data();
}

void CooMatrix::free_data()
{
  A_cplx.clear();
  A.clear();
  this->size = 0;
}

void CooMatrix::add_from_csr(CSRMatrix *m)
{
    free_data();

    this->complex = m->is_complex();

    int *Ap = m->get_Ap();
    int *Ai = m->get_Ai();
    double *Ax = m->get_Ax();
    cplx *Ax_cplx = m->get_Ax_cplx();

    int count = 0;
    // loop through rows...
    for (int i = 1; i <= m->get_size(); i++)
    {
        for (int j = count; j < Ap[i]; j++)
        {
            if (is_complex())
                add(i-1, Ai[count], Ax_cplx[count]);
            else
                add(i-1, Ai[count], Ax[count]);
            count++;
        }
    }
}

void CooMatrix::add_from_csc(CSCMatrix *m)
{
    free_data();

    this->complex = m->is_complex();

    int *Ap = m->get_Ap();
    int *Ai = m->get_Ai();
    double *Ax = m->get_Ax();
    cplx *Ax_cplx = m->get_Ax_cplx();

    int count = 0;
    // loop through columns...
    for (int i = 1; i <= m->get_size(); i++)
    {
        for (int j = count; j < Ap[i]; j++)
        {
            if (is_complex())
                add(Ai[count], i-1, Ax_cplx[count]);
            else
                add(Ai[count], i-1, Ax[count]);
            count++;
        }
    }
}

void CooMatrix::add(int m, int n, double v)
{
    if (this->complex)
        _error("can't use add(int, int, double) for complex matrix");

    // adjusting size if necessary
    if (m+1 > this->size) this->size = m+1;
    if (n+1 > this->size) this->size = n+1;

    // add new
    A[m][n] += v;
}

void CooMatrix::add(int m, int n, cplx v)
{
    if (!(this->complex))
        _error("can't use add(int, int, cplx) for real matrix");

    // adjusting size if necessary
    if (m+1 > this->size) this->size = m+1;
    if (n+1 > this->size) this->size = n+1;

    A_cplx[m][n] += v;
}

void CooMatrix::copy_into(Matrix *m)
{
    m->free_data();

    int index = 0;
    if (this->complex)
    {
        for(std::map<size_t, std::map<size_t, cplx> >::const_iterator it_row = A_cplx.begin(); it_row != A_cplx.end(); ++it_row)
        {
            for(std::map<size_t, cplx>::const_iterator it_col = it_row->second.begin(); it_col != it_row->second.end(); ++it_col)
            {
                m->add(it_row->first,
                       it_col->first,
                       cplx(it_col->second.real(), it_col->second.imag()));

                index++;
            }
        }
    }
    else
    {
        for(std::map<size_t, std::map<size_t, double> >::const_iterator it_row = A.begin(); it_row != A.end(); ++it_row)
        {
            for(std::map<size_t, double>::const_iterator it_col = it_row->second.begin(); it_col != it_row->second.end(); ++it_col)
            {
                m->add(it_row->first,
                       it_col->first,
                       it_col->second);

                index++;
            }
        }
    }
}

void CooMatrix::get_row_col_data(int *row, int *col, double *data)
{
    int index = 0;
    for(std::map<size_t, std::map<size_t, double> >::const_iterator it_row = A.begin(); it_row != A.end(); ++it_row)
    {
        for(std::map<size_t, double>::const_iterator it_col = it_row->second.begin(); it_col != it_row->second.end(); ++it_col)
        {
            row[index] = it_row->first;
            col[index] = it_col->first;
            data[index] = it_col->second;

            index++;
        }
    }
}

void CooMatrix::get_row_col_data(int *row, int *col, cplx *data)
{
    int index = 0;
    for(std::map<size_t, std::map<size_t, cplx> >::const_iterator it_row = A_cplx.begin(); it_row != A_cplx.end(); ++it_row)
    {
        for(std::map<size_t, cplx>::const_iterator it_col = it_row->second.begin(); it_col != it_row->second.end(); ++it_col)
        {
            row[index] = it_row->first;
            col[index] = it_col->first;
            data[index] = it_col->second;

            index++;
        }
    }
}

void CooMatrix::get_row_col_data(int *row, int *col, double *data_real, double *data_imag)
{
    int index = 0;
    for(std::map<size_t, std::map<size_t, cplx> >::const_iterator it_row = A_cplx.begin(); it_row != A_cplx.end(); ++it_row)
    {
        for(std::map<size_t, cplx>::const_iterator it_col = it_row->second.begin(); it_col != it_row->second.end(); ++it_col)
        {
            row[index] = it_row->first;
            col[index] = it_col->first;
            data_real[index] = it_col->second.real();
            data_imag[index] = it_col->second.imag();

            index++;
        }
    }
}

int CooMatrix::get_nnz()
{
    int nnz = 0;
    if (complex)
        for(std::map<size_t, std::map<size_t, cplx> >::const_iterator it_row = A_cplx.begin(); it_row != A_cplx.end(); ++it_row)
            nnz += it_row->second.size();
    else
        for(std::map<size_t, std::map<size_t, double> >::const_iterator it_row = A.begin(); it_row != A.end(); ++it_row)
            nnz += it_row->second.size();

    return nnz;
}

void CooMatrix::times_vector(double* vec, double* result, int rank)
{
    for (int i=0; i < rank; i++) result[i] = 0;

    for(std::map<size_t, std::map<size_t, double> >::const_iterator it_row = A.begin(); it_row != A.end(); ++it_row)
    {
        for(std::map<size_t, double>::const_iterator it_col = it_row->second.begin(); it_col != it_row->second.end(); ++it_col)
        {
            result[it_row->first] += it_col->second * vec[it_col->first];
        }
    }
}

void CooMatrix::print()
{
    printf("\nCoo Matrix:\n");

    if (is_complex())
    {
        for(std::map<size_t, std::map<size_t, cplx> >::const_iterator it_row = A_cplx.begin(); it_row != A_cplx.end(); ++it_row)
        {
            for(std::map<size_t, cplx>::const_iterator it_col = it_row->second.begin(); it_col != it_row->second.end(); ++it_col)
            {
                printf("(%u, %u): (%g, %g)\n",
                       it_row->first,
                       it_col->first,
                       ((cplx) it_col->second).real(),
                       ((cplx) it_col->second).imag());
            }
        }
    }
    else
    {
        for(std::map<size_t, std::map<size_t, double> >::const_iterator it_row = A.begin(); it_row != A.end(); ++it_row)
        {
            for(std::map<size_t, double>::const_iterator it_col = it_row->second.begin(); it_col != it_row->second.end(); ++it_col)
            {
                printf("(%u, %u): %g\n",
                       it_row->first,
                       it_col->first,
                       it_col->second);
            }
        }
    }
}

// *********************************************************************************************************************
CSRMatrix::CSRMatrix(int size) : Matrix()
{
    init();
    this->size = size;
}

CSRMatrix::CSRMatrix(CooMatrix *m) : Matrix()
{
    init();
    this->add_from_coo(m);
}

CSRMatrix::CSRMatrix(CSCMatrix *m) : Matrix()
{
    init();
    this->add_from_csc(m);
}

CSRMatrix::CSRMatrix(DenseMatrix *m) : Matrix()
{
    init();
    this->add_from_dense(m);
}

CSRMatrix::CSRMatrix(Matrix *m) : Matrix()
{
    init();

    if (dynamic_cast<CooMatrix*>(m))
        this->add_from_coo((CooMatrix*)m);
    else if (dynamic_cast<CSCMatrix*>(m))
        this->add_from_csc((CSCMatrix*)m);
    else if (dynamic_cast<DenseMatrix*>(m))
        this->add_from_dense((DenseMatrix*)m);
    else
        _error("Matrix type not supported.");
}

CSRMatrix::~CSRMatrix()
{
    free_data();
}

void CSRMatrix::init()
{
    this->complex = false;
    this->size = 0;
    this->nnz = 0;

    this->Ax = NULL;
    this->Ax_cplx = NULL;
    this->Ap = NULL;
    this->Ai = NULL;
}

void CSRMatrix::free_data()
{
    if (this->Ap != NULL) { delete[] this->Ap; this->Ap = NULL; }
    if (this->Ai != NULL) { delete[] this->Ai; this->Ai = NULL; }
    if (this->Ax != NULL) { delete[] this->Ax; this->Ax = NULL; }
    if (this->Ax_cplx != NULL) { delete[] this->Ax_cplx; this->Ax_cplx = NULL; }

    this->size = 0;
    this->nnz = 0;
}

void CSRMatrix::add_from_dense(DenseMatrix *m)
{
    this->size = m->get_size();

    // count nnz
    this->nnz = 0;
    for(int i = 0; i < this->size; i++)
    {
        for(int j = 0; j < this->size; j++)
        {
            double v = m->get(i, j);
            if (fabs(v) > 1e-12)
                this->nnz++;
        }
    }

    // allocate arrays
    this->Ax = new double[this->nnz];
    this->Ap = new int[this->size+1];
    this->Ai = new int[this->nnz];

    int count = 0;
    this->Ap[0] = 0;
    for(int i = 0; i < this->size; i++)
    {
        for(int j = 0; j < this->size; j++)
        {
            double v = m->get(i, j);
            if (fabs(v) > 1e-12)
            {
                this->Ax[count] = v;
                this->Ai[count] = j;
                count++;
            }
        }
        this->Ap[i+1] = count;
    }
}

void CSRMatrix::add_from_coo(CooMatrix *m)
{
    free_data();

    this->size = m->get_size();
    this->nnz = m->get_nnz();
    this->complex = m->is_complex();

    // allocate data
    this->Ap = new int[this->size + 1];
    this->Ai = new int[this->nnz];
    if (is_complex())
        this->Ax_cplx = new cplx[this->nnz];
    else
        this->Ax = new double[this->nnz];

    // get data
    int *row = new int[this->nnz];
    int *col = new int[this->nnz];
    if (is_complex())
    {
        cplx *data = new cplx[this->nnz];
        m->get_row_col_data(row, col, data);
        coo_to_csr(this->size, this->nnz, row, col, data, Ap, Ai, Ax_cplx);
        if (data != NULL) delete[] data;
    }
    else
    {
        double *data = new double[this->nnz];
        m->get_row_col_data(row, col, data);
        coo_to_csr(this->size, this->nnz, row, col, data, Ap, Ai, Ax);
        if (data != NULL) delete[] data;
    }

    // free data
    if (row != NULL) delete[] row;
    if (col != NULL) delete[] col;
}

void CSRMatrix::add_from_csc(CSCMatrix *m)
{
    free_data();

    this->size = m->get_size();
    this->nnz = m->get_nnz();
    this->complex = m->is_complex();

    // allocate data
    this->Ap = new int[this->size + 1];
    this->Ai = new int[this->nnz];
    if (is_complex())
        this->Ax_cplx = new cplx[this->nnz];
    else
        this->Ax = new double[this->nnz];

    if (is_complex())
    {
        csc_to_csr(this->size, this->nnz, m->get_Ap(), m->get_Ai(), m->get_Ax_cplx(), Ap, Ai, Ax_cplx);
    }
    else
    {
        csc_to_csr(this->size, this->nnz, m->get_Ap(), m->get_Ai(), m->get_Ax(), Ap, Ai, Ax);
    }
}

void CSRMatrix::print()
{
    printf("\nCSR Matrix:\n");
    printf("size: %i\n", this->size);
    printf("nzz: %i\n", this->nnz);

    print_vector("row_ptr", this->Ap, this->size+1);
    print_vector("col_ind", this->Ai, this->nnz);
    if (is_complex())
        print_vector("data", this->Ax_cplx, this->nnz);
    else
        print_vector("data", this->Ax, this->nnz);
}

// *********************************************************************************************************************

CSCMatrix::CSCMatrix(int size) : Matrix()
{
    init();
    this->size = size;
}

CSCMatrix::CSCMatrix(CooMatrix *m) : Matrix()
{
    init();
    this->add_from_coo(m);
}

CSCMatrix::CSCMatrix(DenseMatrix *m) : Matrix()
{
    init();
    this->add_from_dense(m);
}

CSCMatrix::CSCMatrix(CSRMatrix *m) : Matrix()
{
    init();
    this->add_from_csr(m);
}

CSCMatrix::CSCMatrix(Matrix *m) : Matrix()
{
    init();

    if (dynamic_cast<CooMatrix*>(m))
        this->add_from_coo((CooMatrix *) m);
    else if (dynamic_cast<DenseMatrix *>(m))
        this->add_from_dense((DenseMatrix *) m);
    else if (dynamic_cast<CSRMatrix *>(m))
        this->add_from_csr((CSRMatrix *) m);
    else
        _error("Matrix type not supported.");
}

CSCMatrix::CSCMatrix(int size, int nnz, int *Ap, int *Ai, double *Ax)
{
    this->size = size;
    this->nnz = nnz;
    this->complex = false;

    this->Ap = Ap;
    this->Ai = Ai;
    this->Ax = Ax;
}

CSCMatrix::CSCMatrix(int size, int nnz, int *Ap, int *Ai, cplx *Ax_cplx)
{
    this->size = size;
    this->nnz = nnz;
    this->complex = true;

    this->Ap = Ap;
    this->Ai = Ai;
    this->Ax_cplx = Ax_cplx;
}


CSCMatrix::~CSCMatrix()
{
    free_data();
}

void CSCMatrix::init()
{
    this->complex = false;
    this->size = 0;
    this->nnz = 0;

    this->Ax = NULL;
    this->Ax_cplx = NULL;
    this->Ap = NULL;
    this->Ai = NULL;
}

void CSCMatrix::free_data()
{
    if (this->Ap != NULL) { delete[] this->Ap; this->Ap = NULL; }
    if (this->Ai != NULL) { delete[] this->Ai; this->Ai = NULL; }
    if (this->Ax != NULL) { delete[] this->Ax; this->Ax = NULL; }
    if (this->Ax_cplx != NULL) { delete[] this->Ax_cplx; this->Ax_cplx = NULL; }

    size = 0;
    nnz = 0;
}

void CSCMatrix::add_from_dense(DenseMatrix *m)
{
    free_data();

    this->size = m->get_size();
    this->nnz = m->get_nnz();
    this->complex = m->is_complex();

    // allocate data
    this->Ap = new int[this->size + 1];
    this->Ai = new int[this->nnz];
    if (is_complex())
        this->Ax_cplx = new cplx[this->nnz];
    else
        this->Ax = new double[this->nnz];

    // get data
    int *row = new int[this->nnz];
    int *col = new int[this->nnz];
    if (is_complex())
    {
        cplx *data = new cplx[this->nnz];
        dense_to_coo(size, nnz, m->get_A_cplx(), row, col, data);
        coo_to_csc(this->size, this->nnz, row, col, data, Ap, Ai, Ax_cplx);
        if (data != NULL) delete[] data;
    }
    else
    {
        double *data = new double[this->nnz];
        dense_to_coo(size, nnz, m->get_A(), row, col, data);
        print_vector("row", row, this->nnz);
        print_vector("col", col, this->nnz);
        print_vector("Ax", data, this->nnz);
        coo_to_csc(this->size, this->nnz, row, col, data, Ap, Ai, Ax);
        if (data != NULL) delete[] data;
    }

    // free data
    if (row != NULL) delete[] row;
    if (col != NULL) delete[] col;
}

void CSCMatrix::add_from_coo(CooMatrix *m)
{
    free_data();

    this->size = m->get_size();
    this->nnz = m->get_nnz();
    this->complex = m->is_complex();

    // allocate data
    this->Ap = new int[this->size + 1];
    this->Ai = new int[this->nnz];
    if (is_complex())
        this->Ax_cplx = new cplx[this->nnz];
    else
        this->Ax = new double[this->nnz];

    // get data
    int *row = new int[this->nnz];
    int *col = new int[this->nnz];
    if (is_complex())
    {
        cplx *data = new cplx[this->nnz];
        m->get_row_col_data(row, col, data);
        coo_to_csc(this->size, this->nnz, row, col, data, Ap, Ai, Ax_cplx);
        if (data != NULL) delete[] data;
    }
    else
    {
        double *data = new double[this->nnz];
        m->get_row_col_data(row, col, data);
        coo_to_csc(this->size, this->nnz, row, col, data, Ap, Ai, Ax);
        if (data != NULL) delete[] data;
    }

    // free data
    if (row != NULL) delete[] row;
    if (col != NULL) delete[] col;
}

void CSCMatrix::add_from_csr(CSRMatrix *m)
{
    free_data();

    this->size = m->get_size();
    this->nnz = m->get_nnz();
    this->complex = m->is_complex();

    // allocate data
    this->Ap = new int[this->size + 1];
    this->Ai = new int[this->nnz];
    if (is_complex())
        this->Ax_cplx = new cplx[this->nnz];
    else
        this->Ax = new double[this->nnz];

    if (is_complex())
    {
        csr_to_csc(this->size, this->nnz, m->get_Ap(), m->get_Ai(), m->get_Ax_cplx(), Ap, Ai, Ax_cplx);
    }
    else
    {
        csr_to_csc(this->size, this->nnz, m->get_Ap(), m->get_Ai(), m->get_Ax(), Ap, Ai, Ax);
    }
}

void CSCMatrix::print()
{
    printf("\nCSC Matrix:\n");
    printf("size: %i\n", this->size);
    printf("nzz: %i\n", this->nnz);

    print_vector("col_ptr", this->Ap, this->size+1);
    print_vector("row_ind", this->Ai, this->nnz);
    if (is_complex())
        print_vector("data", this->Ax_cplx, this->nnz);
    else
        print_vector("data", this->Ax, this->nnz);
}

// ******************************************************************************************************************************

template<typename T>
void dense_to_coo(int size, int nnz, T **Ad, int *row, int *col, T *A)
{
    int count = 0;
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            if (std::abs(Ad[i][j]) > 1e-12)
            {
                row[count] = i;
                col[count] = j;
                A[count] = Ad[i][j];

                count++;
            }
        }
    }
}

template<typename T>
void coo_to_csr(int size, int nnz, int *row, int *col, T *A, int *Ap, int *Ai, T *Ax)
{
    std::fill(Ap, Ap + size, 0);

    for (int n = 0; n < nnz; n++)
    {
        Ap[row[n]]++;
    }

    // cumsum the nnz per row to get this->row_ptr[]
    for(int i = 0, cumsum = 0; i < size; i++)
    {
        int temp = Ap[i];
        Ap[i] = cumsum;
        cumsum += temp;
    }
    Ap[size] = nnz;

    // write Aj, Ax into Bj, Bx
    for(int n = 0; n < nnz; n++)
    {
        int index  = row[n];
        int dest = Ap[index];

        Ai[dest] = col[n];
        Ax[dest] = A[n];

        Ap[index]++;
    }

    for(int i = 0, last = 0; i <= size; i++)
    {
        int temp = Ap[i];
        Ap[i]  = last;
        last   = temp;
    }
}

template<typename T>
void coo_to_csc(int size, int nnz, int *row, int *col, T *A, int *Ap, int *Ai, T *Ax)
{
    coo_to_csr(size, nnz, col, row, A, Ap, Ai, Ax);
}

template<typename T>
void csr_to_csc(int size, int nnz, int *Arp, int *Ari, T *Arx, int *Acp, int *Aci, T *Acx)
{
    //compute number of non-zero entries per column of A
    std::fill(Acp, Acp + size, 0);

    for (int n = 0; n < nnz; n++)
        Acp[Ari[n]]++;

    // cumsum the nnz per column to get Bp[]
    for(int col = 0, cumsum = 0; col < size; col++)
    {
        int temp  = Acp[col];
        Acp[col] = cumsum;
        cumsum += temp;
    }
    Acp[size] = nnz;

    for(int row = 0; row < size; row++)
    {
        for(int jj = Arp[row]; jj < Arp[row+1]; jj++)
        {
            int col  = Ari[jj];
            int dest = Acp[col];

            Aci[dest] = row;
            Acx[dest] = Arx[jj];

            Acp[col]++;
        }
    }

    for(int col = 0, last = 0; col <= size; col++)
    {
        int temp  = Acp[col];
        Acp[col] = last;
        last = temp;
    }
}

template<typename T>
void csc_to_csr(int size, int nnz, int *Acp, int *Aci, T *Acx, int *Arp, int *Ari, T *Arx)
{
    csr_to_csc(size, nnz, Acp, Aci, Acx, Arp, Ari, Arx);
}

template<typename T>
void csc_to_coo(int size, int nnz, int *Ap, int *Ai, T *Ax, int *row, int *col, T *A)
{
    int count = 0;
    // loop through columns...
    for (int i = 1; i <= size; i++)
    {
        for (int j = count; j < Ap[i]; j++)
        {
            A[count] = Ax[count];
            row[count] = Ai[count];
            col[count] = i-1;

            count++;
        }
    }
}

template<typename T>
void csr_to_coo(int size, int nnz, int *Ap, int *Ai, T *Ax, int *row, int *col, T *A)
{
    int count = 0;
    // loop through rows...
    for (int i = 1; i <= size; i++)
    {
        for (int j = count; j < Ap[i]; j++)
        {
            A[count] = Ax[count];
            col[count] = Ai[count];
            row[count] = i-1;

            count++;
        }
    }
}

// matrix vector multiplication
void mat_dot(Matrix *A, double *x, double *result, int n_dof)
{
    A->times_vector(x, result, n_dof);
}

// vector vector multiplication
double vec_dot(double *r, double *s, int n_dof)
{
    double result = 0;
    for (int i=0; i < n_dof; i++) result += r[i]*s[i];
    return result;
}

/// Solves the set of n linear equations A*x = b, where a is a positive-definite symmetric matrix.
/// a[n][n] and p[n] are input as the output of the routine choldc. Only the lower
/// subdiagonal portion of a is accessed. b[n] is input as the right-hand side vector. The
/// solution vector is returned in x[n]. a, n, and p are not modified and can be left in place
/// for successive calls with different right-hand sides b. b is not modified unless you identify b and
/// x in the calling sequence, which is allowed. The right-hand side b can be complex, in which case
/// the solution x is also complex.
void ludcmp(double** a, int n, int* indx, double* d)
{
    int i, imax = 0, j, k;
    double big, dum, sum, temp;
    double* vv = new double[n];

    *d = 1.0;
    for (i = 0; i < n; i++)
    {
        big=0.0;
        for (j = 0; j < n; j++)
            if ((temp = fabs(a[i][j])) > big)
                big = temp;
        if (big == 0.0) _error("Singular matrix!");
        vv[i] = 1.0 / big;
    }
    for (j = 0; j < n; j++)
    {
        for (i = 0; i < j; i++)
        {
            sum = a[i][j];
            for (k = 0; k < i; k++) sum -= a[i][k]*a[k][j];
            a[i][j] = sum;
        }
        big = 0.0;
        for (i = j; i < n; i++)
        {
            sum = a[i][j];
            for (k = 0; k < j; k++) sum -= a[i][k]*a[k][j];
            a[i][j] = sum;
            if ((dum = vv[i]*fabs(sum)) >= big)
            {
                big = dum;
                imax = i;
            }
        }
        if (j != imax)
        {
            for (k = 0; k < n; k++)
            {
                dum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = dum;
            }
            *d = -(*d);
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if (a[j][j] == 0.0) a[j][j] = 1e-20;
        if (j != n-1)
        {
            dum = 1.0 / (a[j][j]);
            for (i = j+1; i < n; i++) a[i][j] *= dum;
        }
    }
    if (vv != NULL) delete [] vv;
}

/// Solves the set of n linear equations AX = B. Here a[n][n] is input, not as the matrix
/// A but rather as its LU decomposition, determined by the routine ludcmp. indx[n] is input
/// as the permutation vector returned by ludcmp. b[n] is input as the right-hand side vector
/// B, and returns with the solution vector X. a, n, and indx are not modified by this routine
/// and can be left in place for successive calls with different right-hand sides b. This routine takes
/// into account the possibility that b will begin with many zero elements, so it is efficient for use
/// in matrix inversion.
void lubksb(double** a, int n, int* indx, double* b)
{
    int i, ip, j;
    double sum;

    for (i = 0; i < n; i++)
    {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        for (j = 0; j < i; j++) sum -= a[i][j]*b[j];
        b[i] = sum;
    }
    for (i = n-1; i >= 0; i--)
    {
        sum = b[i];
        for (j = i+1; j < n; j++) sum -= a[i][j]*b[j];
        b[i] = sum / a[i][i];
    }
}
