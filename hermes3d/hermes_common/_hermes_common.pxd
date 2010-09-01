# Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
# Distributed under the terms of the BSD license (see the LICENSE
# file for the exact terms).
# Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

from ref cimport PyObject, Py_INCREF, Py_DECREF

ctypedef double complex cplx


cdef extern from *:
    ctypedef char* char_p       "char*"
    ctypedef char* const_char_p "const char*"

    # This is just the C++ "delete" statement
    void delete(...)

cdef extern from "Python.h":
    ctypedef struct RichPyTypeObject "PyTypeObject":
        object tp_new(RichPyTypeObject*, object, PyObject*)


cdef extern from "math.h":
    double c_sqr "sqr"(double x)
    double c_sqrt "sqrt"(double x)
    double c_atan "atan"(double x)
    double c_pi "M_PI"


cdef extern from "stdlib.h":
    ctypedef unsigned long size_t
    void *malloc (size_t size)
    void free(void *mem)
    void *memcpy(void *dst, void *src, long n)

    void exit(int exit_code)

cdef extern from "arrayobject.h":
    cdef enum NPY_TYPES:
        NPY_BOOL
        NPY_BYTE
        NPY_UBYTE
        NPY_SHORT
        NPY_USHORT
        NPY_INT
        NPY_UINT
        NPY_LONG
        NPY_ULONG
        NPY_LONGLONG
        NPY_ULONGLONG
        NPY_FLOAT
        NPY_DOUBLE
        NPY_LONGDOUBLE
        NPY_CFLOAT
        NPY_CDOUBLE
        NPY_CLONGDOUBLE
        NPY_OBJECT
        NPY_STRING
        NPY_UNICODE
        NPY_VOID
        NPY_NTYPES
        NPY_NOTYPE

        NPY_COMPLEX32
        NPY_COMPLEX64
        NPY_COMPLEX128

    ctypedef int npy_intp

    ctypedef extern class numpy.ndarray [object PyArrayObject]:
        cdef char *data
        cdef int nd
        cdef npy_intp *dimensions
        cdef npy_intp *strides
        cdef int flags

    object PyArray_SimpleNewFromData(int nd, npy_intp* dims, int typenum,
            void* data)
    void import_array()


cdef extern from "stdcython.h":
    void throw_exception(char *msg)

cdef extern from "matrix.h":

    cdef struct c_Matrix "Matrix":
        int get_size()
        int is_complex()
        void add(int m, int n, double v)
        void add_cplx "add"(int m, int n, cplx v)

    cdef struct c_CooMatrix "CooMatrix":
        int get_nnz()
        void get_row_col_data(int *row, int *col, double *data)
        void get_row_col_data_cplx "get_row_col_data"(int *row, int *col,
                cplx *data)
    c_CooMatrix *new_CooMatrix "new CooMatrix" (int size, int is_complex)

    cdef struct c_CSCMatrix "CSCMatrix"

    cdef struct c_CSRMatrix "CSRMatrix":
        CSRMatrix(int size)
        int *get_Ap()
        int *get_Ai()
        double *get_Ax()
        cplx *get_Ax_cplx()
        int get_nnz()
    c_CSRMatrix *new_CSRMatrix_size "new CSRMatrix" (int size)
    c_CSRMatrix *new_CSRMatrix_coo_matrix "new CSRMatrix" (c_CooMatrix *m)
    c_CSRMatrix *new_CSRMatrix_csc_matrix "new CSRMatrix" (c_CSCMatrix *m)

    cdef struct c_CSCMatrix "CSCMatrix":
        CSCMatrix(int size)
        int *get_Ap()
        int *get_Ai()
        double *get_Ax()
        cplx *get_Ax_cplx()
        int get_nnz()
    c_CSCMatrix *new_CSCMatrix_size "new CSCMatrix" (int size)
    c_CSCMatrix *new_CSCMatrix_coo_matrix "new CSCMatrix" (c_CooMatrix *m)
    c_CSCMatrix *new_CSCMatrix_csr_matrix "new CSCMatrix" (c_CSRMatrix *m)


cdef api object c2numpy_int(int *A, int len)
cdef api object c2numpy_double(double *A, int len)
cdef api void numpy2c_double_inplace(object A_n, double **A_c, int *n)

cdef inline PY_NEW(T)

