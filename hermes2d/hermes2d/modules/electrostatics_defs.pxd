from libcpp cimport bool
from libcpp.vector cimport vector

cdef extern from "electrostatics.h":

    cdef cppclass Electrostatics:
        bool set_mesh_filename(char *filename)
        void set_initial_mesh_refinement(int init_ref_num)
        void set_initial_poly_degree(int p)
        void set_material_markers(vector[int] mat_markers)
        #void set_permittivity_array(const std::vector<double> p_array)
        #void set_charge_density_array(const std::vector<double> cd_array)
        #void set_boundary_markers_value(const std::vector<int> &bdy_markers_val)
        #void set_boundary_values(const std::vector<double> bc_val)
        #void set_boundary_markers_derivative(const std::vector<int> bdy_markers_der)
        #void set_boundary_derivatives(const std::vector<double> bc_der)
        #bool calculate(Solution* phi)

