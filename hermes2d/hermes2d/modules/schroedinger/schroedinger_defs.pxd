from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.pair cimport pair

from hermes2d.hermes2d_defs cimport Solution

cdef extern from "schroedinger.h":

    cdef cppclass Potential:
        double get_value(double x, double y)
        void get_values(int n, double *x, double *y, double *values)

    cdef cppclass PotentialHarmonicOscillator(Potential):
        PotentialHarmonicOscillator()
        void set_omega(double omega)

    cdef cppclass ModuleSchroedinger:
        ModuleSchroedinger()
        #~ModuleSchroedinger()
        #void set_potential(const RCP<Potential> &potential)
        #void assemble(const Ptr<SparseMatrix> &A, const Ptr<SparseMatrix> &B)
