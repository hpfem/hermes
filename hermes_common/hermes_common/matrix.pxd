cimport hermes_common.cpp.matrix

cdef class Matrix:
    cdef hermes_common.cpp.matrix.Matrix *thisptr
