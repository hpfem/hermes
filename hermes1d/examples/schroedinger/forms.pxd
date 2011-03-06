from hermes1d.h1d_wrapper.hermes1d cimport Space
from hermes_common.cpp.matrix cimport SparseMatrix

cdef extern from "forms.h":
    int eqn_type_R
    int eqn_type_rR
    void assemble_schroedinger(Space *mesh, SparseMatrix *A, SparseMatrix *B,
            int l, int Z, int equation_type) except +
