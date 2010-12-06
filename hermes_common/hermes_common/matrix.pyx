cdef class Matrix:
    pass

cdef class SparseMatrix(Matrix):

    #def __cinit__(self):
    #    self.thisptr = new matrix.SparseMatrix()

    #def __dealloc__(self):
    #    del self.thisptr

    cdef matrix.SparseMatrix *getptr(self):
        return <matrix.SparseMatrix *>(self.thisptr)

cdef class Vector:

    #def __cinit__(self):
    #    self.thisptr = new matrix.Vector()

    #def __dealloc__(self):
    #    del self.thisptr
    pass
