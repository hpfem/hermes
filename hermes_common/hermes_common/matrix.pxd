cimport hermes_common.cpp.matrix

cdef class Matrix:
    cdef hermes_common.cpp.matrix.Matrix *thisptr

cdef class Vector:
    cdef hermes_common.cpp.matrix.Vector *thisptr
