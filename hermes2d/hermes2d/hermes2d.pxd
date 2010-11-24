cimport hermes2d_defs

cdef class MeshFunction:
    cdef hermes2d_defs.MeshFunction *thisptr

cdef class Solution(MeshFunction):

    cdef hermes2d_defs.Solution *getptr(self)
