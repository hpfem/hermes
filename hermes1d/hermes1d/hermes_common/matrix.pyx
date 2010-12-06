cdef class Matrix:
    pass

cdef class SparseMatrix(Matrix):

    cdef matrix.SparseMatrix* as_SparseMatrix(self):
        return <matrix.SparseMatrix*> self.thisptr

cdef class CSCMatrix(SparseMatrix):
    """
    Represents a CSC matrix.
    """

    def __cinit__(self):
        self.thisptr = new matrix.UMFPackMatrix()

    def __dealloc__(self):
        del self.thisptr

    def to_scipy_csc(self):
        return "1"

cdef class Vector:
    pass

cdef class AVector(Vector):
    """
    A Vector, represented by a C array internally.
    """

    def __cinit__(self):
        self.thisptr = new matrix.UMFPackVector()

    def __dealloc__(self):
        del self.thisptr
