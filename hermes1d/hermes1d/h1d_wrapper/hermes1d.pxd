# Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
# Distributed under the terms of the BSD license (see the LICENSE
# file for the exact terms).
# Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

cdef extern from "hermes1d.h":

    ctypedef double double4[4]
    ctypedef double double3[3]
    ctypedef int int3[3]
    ctypedef int int2[2]

    cdef cppclass Element:
        double x1, x2
        int p
        int *dof

    cdef cppclass Mesh:
        Mesh(double a, double b, int n_elem, int p_init, int n_eq, int
                n_sln, int print_banner)
        Mesh(int n_macro_elem, double *pts_array, int *p_array, int *m_array,
             int *div_array, int n_eq, int n_sln, int print_banner)
        void create(double A, double B, int n)
        int get_n_base_elems()
        int get_n_dofs()
        void set_poly_orders(int poly_order)
        int assign_dofs()
        Element *get_base_elems()
        void set_dirichlet_bc_left(int eq_n, double val)
        void set_dirichlet_bc_right(int eq_n, double val)
        void copy_vector_to_mesh(double *y, int sln)
        void copy_mesh_to_vector(double *y, int sln)
        void plot(char* filename)

    cdef cppclass Linearizer:
        Linearizer(Mesh *mesh)
        void plot_solution(char *out_filename,
                int plotting_elem_subdivision)
        void get_xy_mesh(int comp, int plotting_elem_subdivision,
                double **x, double **y, int *n)

    cdef cppclass Iterator:
        Iterator(Mesh *mesh)
        void reset()
        Element *first_active_element()
        Element *next_active_element()
        Element *last_active_element()
