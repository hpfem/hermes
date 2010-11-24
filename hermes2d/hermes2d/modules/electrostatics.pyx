cimport electrostatics_defs

cdef class Electrostatics:
    cdef electrostatics_defs.Electrostatics *thisptr

    def __init__(self):
        self.thisptr = new electrostatics_defs.Electrostatics()

    def __dealloc__(self):
        del self.thisptr

    def set_initial_mesh_refinement(self, int init_ref_num):
        self.thisptr.set_initial_mesh_refinement(init_ref_num)

    def set_initial_poly_degree(self, int p):
        self.thisptr.set_initial_poly_degree(p)
