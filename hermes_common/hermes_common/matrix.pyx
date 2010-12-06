cimport hermes_common.cpp.matrix

cdef class Matrix:

    def __cinit__(self):
        self.thisptr = new hermes_common.cpp.matrix.Matrix()

    def __dealloc__(self):
        del self.thisptr

cdef class SparseMatrix(Matrix):

    def __cinit__(self):
        self.thisptr = new hermes_common.cpp.matrix.SparseMatrix()

    def __dealloc__(self):
        del self.thisptr

    cdef hermes_common.cpp.matrix.SparseMatrix *getptr(self):
        return <hermes_common.cpp.matrix.SparseMatrix *>(self.thisptr)

cdef class Vector:

    def __cinit__(self):
        self.thisptr = new hermes_common.cpp.matrix.Vector()

    def __dealloc__(self):
        del self.thisptr
