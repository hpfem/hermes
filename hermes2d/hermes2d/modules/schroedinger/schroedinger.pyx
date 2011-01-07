from numpy cimport ndarray
from libcpp.vector cimport vector
from libcpp.pair cimport pair

cimport schroedinger_defs
from hermes2d.hermes2d cimport Solution

cdef class Potential:
    cdef schroedinger_defs.Potential *thisptr

    def get_value(self, double x, double y):
        return self.thisptr.get_value(x, y)

cdef class PotentialHarmonicOscillator(Potential):

    def __init__(self):
        self.thisptr = new schroedinger_defs.PotentialHarmonicOscillator()

    def __dealloc__(self):
        del self.thisptr

    def set_omega(self, double omega):
        cdef schroedinger_defs.PotentialHarmonicOscillator *s
        s = <schroedinger_defs.PotentialHarmonicOscillator *> self.thisptr
        s.set_omega(omega)

cdef class ModuleSchroedinger:
    cdef schroedinger_defs.ModuleSchroedinger *thisptr

    def __init__(self):
        self.thisptr = new schroedinger_defs.ModuleSchroedinger()

    def __dealloc__(self):
        del self.thisptr

    def set_potential(self, Potential potential):
        # TODO: how can we convert "potential" to "RCP<Potential>"?
        cdef schroedinger_defs.RCP[Potential] p
        self.thisptr.set_potential(p)
