from hermes1d.h1d_wrapper.hermes1d cimport Mesh
from hermes_common._hermes_common cimport c_Matrix as Matrix

cdef extern from "forms.h":
    void assemble_dirac(Mesh *mesh, Matrix *A, Matrix *B, int _kappa,
        int _Z) except +
