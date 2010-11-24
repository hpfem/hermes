from numpy cimport ndarray
from libcpp.vector cimport vector

cimport electrostatics_defs

cdef vector[int] array2vector_int(a):
    cdef vector[int] v
    for i in range(len(a)):
        v.push_back(a[i])
    return v

cdef vector[int] array2vector_double(a):
    cdef vector[int] v
    for i in range(len(a)):
        v.push_back(a[i])
    return v

cdef class Electrostatics:
    cdef electrostatics_defs.Electrostatics *thisptr

    def __init__(self):
        self.thisptr = new electrostatics_defs.Electrostatics()

    def __dealloc__(self):
        del self.thisptr

    def set_mesh_filename(self, filename):
        return self.thisptr.set_mesh_filename(filename)

    def set_initial_mesh_refinement(self, int init_ref_num):
        self.thisptr.set_initial_mesh_refinement(init_ref_num)

    def set_initial_poly_degree(self, int p):
        self.thisptr.set_initial_poly_degree(p)

    def set_material_markers(self, mat_markers):
        self.thisptr.set_material_markers(array2vector_int(mat_markers))
