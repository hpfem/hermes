from numpy cimport ndarray
from libcpp.vector cimport vector

cimport electrostatics_defs
from hermes2d.hermes2d cimport Solution

cdef vector[int] array2vector_int(a):
    cdef vector[int] v
    for i in range(len(a)):
        v.push_back(a[i])
    return v

cdef vector[double] array2vector_double(a):
    cdef vector[double] v
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

    def set_permittivity_array(self, p_array):
        self.thisptr.set_permittivity_array(array2vector_double(p_array))

    def set_charge_density_array(self, cd_array):
        self.thisptr.set_charge_density_array(array2vector_double(cd_array))

    def set_boundary_markers_value(self, bdy_markers_val):
        self.thisptr.set_boundary_markers_value(array2vector_int(bdy_markers_val))

    def set_boundary_values(self, bc_val):
        self.thisptr.set_boundary_values(array2vector_double(bc_val))

    def set_boundary_markers_derivative(self, bdy_markers_der):
        self.thisptr.set_boundary_markers_derivative(array2vector_int(bdy_markers_der))

    def set_boundary_derivatives(self, bc_der):
        self.thisptr.set_boundary_derivatives(array2vector_double(bc_der))

    def calculate(self):
        s = Solution()
        r = self.thisptr.calculate(s.getptr())
        return r, s
