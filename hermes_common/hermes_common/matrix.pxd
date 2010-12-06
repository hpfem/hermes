cimport hermes_common.cpp.matrix

cdef class Matrix:
    cdef hermes_common.cpp.matrix.Matrix *thisptr

cdef class SparseMatrix(Matrix):

    cdef hermes_common.cpp.matrix.SparseMatrix *getptr(self)

cdef class Vector:
    cdef hermes_common.cpp.matrix.Vector *thisptr
