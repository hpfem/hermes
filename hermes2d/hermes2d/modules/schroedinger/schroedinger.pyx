from numpy cimport ndarray
from libcpp.vector cimport vector
from libcpp.pair cimport pair

cimport schroedinger_defs
from hermes2d.hermes2d cimport Solution
from hermes2d.hermes_common.matrix cimport Matrix

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
    cdef schroedinger_defs.RCP[schroedinger_defs.ModuleSchroedinger] *thisptr

    def __init__(self):
        #self.thisptr = rcp(new schroedinger_defs.ModuleSchroedinger())
        pass

    def set_potential(self, Potential potential):
        pass
        #self.thisptr.set_potential(potential.thisptr)

    def assemble(self, Matrix A, Matrix B):
        pass
        #self.thisptr.assemble(A.thisptr.ptr(), B.this.ptr.ptr())
