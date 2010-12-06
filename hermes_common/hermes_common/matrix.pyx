cimport hermes_common.cpp.matrix

cdef class Matrix:

    def __cinit__(self):
        self.thisptr = new hermes_common.cpp.matrix.Matrix()

    def __dealloc__(self):
        del self.thisptr
