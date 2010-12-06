from hermes1d.hermes_common.cpp.matrix cimport (Matrix as cMatrix,
        SparseMatrix as cSparseMatrix, UMFPackMatrix as cUMFPackMatrix,
        UMFPackVector as cUMFPackVector, Vector as cVector)

cdef class Matrix:
    cdef cMatrix *thisptr

    cpdef int get_size(self)

cdef class SparseMatrix(Matrix):

    cdef cSparseMatrix* as_SparseMatrix(self)

cdef class CSCMatrix(SparseMatrix):

    cdef cUMFPackMatrix* as_UMFPackMatrix(self)

cdef class Vector:
    cdef cVector *thisptr

cdef class AVector(Vector):

    cdef cUMFPackVector* as_UMFPackVector(self)
