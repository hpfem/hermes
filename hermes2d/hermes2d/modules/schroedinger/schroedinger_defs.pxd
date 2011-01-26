from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.pair cimport pair

from hermes2d.hermes2d_defs cimport Solution
from hermes2d.hermes_common.matrix cimport SparseMatrix

cdef extern from "hermes2d.h" namespace "Teuchos":

    cdef cppclass RCP[T]:
        pass

    cdef cppclass Ptr[T]:
        pass

cdef extern from "schroedinger.h" namespace "Schroedinger":

    cdef cppclass Potential:
        double get_value(double x, double y)
        void get_values(int n, double *x, double *y, double *values)

    cdef cppclass PotentialHarmonicOscillator(Potential):
        PotentialHarmonicOscillator()
        void set_omega(double omega)

    cdef cppclass ModuleSchroedinger:
        ModuleSchroedinger()
        void set_potential(RCP[Potential] &potential)
        void assemble(Ptr[SparseMatrix] &A, Ptr[SparseMatrix] &B)
