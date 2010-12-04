cimport hermes1d

cdef class Mesh:
    cdef hermes1d.Space *thisptr
    cdef object delptr
