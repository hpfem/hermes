from numpy cimport ndarray
from libcpp.vector cimport vector
from libcpp.pair cimport pair

cimport schroedinger_defs
from hermes2d.hermes2d cimport Solution

cdef class ModuleSchroedinger:
    cdef schroedinger_defs.ModuleSchroedinger *thisptr

    def __init__(self):
        self.thisptr = new schroedinger_defs.ModuleSchroedinger()

    def __dealloc__(self):
        del self.thisptr
