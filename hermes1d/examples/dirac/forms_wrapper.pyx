from hermes1d.h1d_wrapper.h1d_wrapper cimport Mesh
from hermes_common._hermes_common cimport Matrix
cimport forms

def assemble_dirac(Mesh mesh, Matrix A, Matrix B, int kappa, int Z):
    forms.assemble_dirac(mesh.thisptr, A.thisptr, B.thisptr, kappa, Z)
