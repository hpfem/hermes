from cpp cimport matrix

cdef class Matrix:
    cdef matrix.Matrix *thisptr

cdef class SparseMatrix(Matrix):

    cdef matrix.SparseMatrix* as_SparseMatrix(self)

cdef class CSCMatrix(SparseMatrix):
    pass

cdef class Vector:
    cdef matrix.Vector *thisptr

cdef class AVector(Vector):
    pass
