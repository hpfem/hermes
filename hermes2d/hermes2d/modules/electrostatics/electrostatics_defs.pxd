from libcpp cimport bool
from libcpp.vector cimport vector

from hermes2d.hermes2d_defs cimport Solution

cdef extern from "electrostatics.h":

    cdef cppclass Electrostatics:
        void set_mesh_str(char *mesh)
        void set_initial_mesh_refinement(int init_ref_num)
        void set_initial_poly_degree(int p)
        void set_material_markers(vector[int] &mat_markers)
        void set_permittivity_array(vector[double] &p_array)
        void set_charge_density_array(vector[double] &cd_array)
        void set_boundary_markers_value(vector[int] &bdy_markers_val)
        void set_boundary_values(vector[double] &bc_val)
        void set_boundary_markers_derivative(vector[int] &bdy_markers_der)
        void set_boundary_derivatives(vector[double] &bc_der)
        bool calculate(Solution* phi)
