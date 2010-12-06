from hermes1d.h1d_wrapper.hermes1d cimport Mesh
from hermes_common._hermes_common cimport c_Matrix

cdef extern from "forms.h":
    int c_eqn_type_R "eqn_type_R"
    int c_eqn_type_rR "eqn_type_rR"
    void c_assemble_schroedinger "assemble_schroedinger"(Mesh *mesh,
            c_Matrix *A, c_Matrix *B, int l, int Z, int equation_type) except +
