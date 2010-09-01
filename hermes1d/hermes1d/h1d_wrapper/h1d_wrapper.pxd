cimport hermes1d

cdef class Mesh:
    cdef hermes1d.Mesh *thisptr
    cdef object delptr
