from numpy cimport ndarray
from libcpp.vector cimport vector

cimport basic_defs
from hermes2d.hermes2d cimport double_pair
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

cdef vector[double_pair] array2vector_double_pair(a):
    cdef vector[double_pair] v
    for i in range(len(a)):
        v.push_back(a[i])
    return v

cdef class basic:
    cdef basic_defs.Basic *thisptr

    def __init__(self):
        self.thisptr = new basic_defs.Basic()

    def __dealloc__(self):
        del self.thisptr

    def set_mesh_str(self, mesh):
        self.thisptr.set_mesh_str(mesh)

    def set_initial_mesh_refinement(self, int init_ref_num):
        self.thisptr.set_initial_mesh_refinement(init_ref_num)

    def set_initial_poly_degree(self, int p):
        self.thisptr.set_initial_poly_degree(p)

    def set_material_markers(self, mat_markers):
        self.thisptr.set_material_markers(array2vector_int(mat_markers))

    def set_c1_array(self, c1_array):
        self.thisptr.set_c1_array(array2vector_double(c1_array))

    def set_c2_array(self, c2_array):
        self.thisptr.set_c2_array(array2vector_double(c2_array))

    def set_c3_array(self, c3_array):
        self.thisptr.set_c3_array(array2vector_double(c3_array))

    def set_c4_array(self, c4_array):
        self.thisptr.set_c4_array(array2vector_double(c4_array))

    def set_c5_array(self, c5_array):
        self.thisptr.set_c5_array(array2vector_double(c5_array))

    def set_dirichlet_markers(self, bdy_markers_dirichlet):
        self.thisptr.set_dirichlet_markers(array2vector_int(bdy_markers_dirichlet))

    def set_dirichlet_values(self, bdy_values_dirichlet):
        self.thisptr.set_dirichlet_values(array2vector_double(bdy_values_dirichlet))

    def set_neumann_markers(self, bdy_markers_neumann):
        self.thisptr.set_neumann_markers(array2vector_int(bdy_markers_neumann))

    def set_neumann_values(self, bdy_values_neumann):
        self.thisptr.set_neumann_values(array2vector_double(bdy_values_neumann))

    def set_newton_markers(self, bdy_markers_newton):
        self.thisptr.set_newton_markers(array2vector_int(bdy_markers_newton))

    def set_neumann_values(self, bdy_values_newton):
        self.thisptr.set_newton_values(array2vector_double_pair(bdy_values_newton))

    def calculate(self):
        s = Solution()
        r = self.thisptr.calculate(s.getptr())
        return r, s
