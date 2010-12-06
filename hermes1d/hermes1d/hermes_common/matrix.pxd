from cpp cimport matrix

cdef class Matrix:
    cdef matrix.Matrix *thisptr

cdef class SparseMatrix(Matrix):

    cdef matrix.SparseMatrix *getptr(self)

cdef class Vector:
    cdef matrix.Vector *thisptr
