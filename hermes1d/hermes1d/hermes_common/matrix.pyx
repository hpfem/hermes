cdef class Matrix:
    pass

cdef class SparseMatrix(Matrix):

    cdef matrix.SparseMatrix* as_SparseMatrix(self):
        return <matrix.SparseMatrix*> self.thisptr

cdef class UMFPackMatrix(SparseMatrix):

    def __cinit__(self):
        self.thisptr = new matrix.UMFPackMatrix()

    def __dealloc__(self):
        del self.thisptr

cdef class Vector:
    pass

cdef class UMFPackVector(Vector):

    def __cinit__(self):
        self.thisptr = new matrix.UMFPackVector()

    def __dealloc__(self):
        del self.thisptr
