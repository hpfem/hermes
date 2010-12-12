from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.pair cimport pair

from hermes2d.hermes2d_defs cimport Solution

# FIXME; Space needs wrapping
# from hermes2d.hermes2d_defs cimport Space

cdef extern from "basic.h":

    cdef cppclass ModuleBasic:
        void set_mesh_str(char *mesh)
        void set_initial_mesh_refinement(int init_ref_num)
        void set_matrix_solver(char *meshsolver_name)
        void set_initial_poly_degree(int p)
        void set_material_markers(vector[int] &mat_markers)
        void set_c1_array(vector[double] &c1_array)
        void set_c2_array(vector[double] &c2_array)
        void set_c3_array(vector[double] &c3_array)
        void set_c4_array(vector[double] &c4_array)
        void set_c5_array(vector[double] &c5_array)
        void set_dirichlet_markers(vector[int] &bdy_markers_dirichlet)
        void set_dirichlet_values(vector[int] &bdy_markers_dirichlet, vector[double] &bdy_values_dirichlet)
        void set_neumann_markers(vector[int] &bdy_markers_neumann)
        void set_neumann_values(vector[double] &bdy_values_neumann)
        void set_newton_markers(vector[int] &bdy_markers_newton)
        void set_newton_values(vector[pair[double, double]] &bdy_values_newton)
        double get_assembly_time()
        double get_solver_time()
        void get_solution(Solution* s)
        #void get_space(Space* s)
        bool calculate()
