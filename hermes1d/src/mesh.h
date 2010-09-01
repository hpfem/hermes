// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef _MESH_H_
#define _MESH_H_

#include "common.h"
#include "legendre.h"
#include "lobatto.h"

class Element {
public:
    Element();
    Element(double x_left, double x_right, int level, int deg, 
            int n_eq, int n_sln, int marker);
    void free_element() {
        if (this->sons[0] != NULL) delete this->sons[0];
        if (this->sons[1] != NULL) delete this->sons[1];
    }
    ~Element() {
        this->free_element();
    }
    void init(double x1, double x2, int p_init, 
	      int id, int active, int level, int n_eq, int n_sln, int marker);
    void copy_into(Element *e_trg);
    void copy_recursively_into(Element *e_trg);
    double get_x_phys(double x_ref); // gets physical coordinate of a reference poin
    double calc_elem_norm_squared(int norm);
    void get_coeffs_from_vector(double *y, int sln=0);
    void copy_coeffs_to_vector(double *y, int sln=0);
    void copy_dofs(int sln_src, int sln_trg);
    void get_solution_quad(int flag, int quad_order, 
                           double val_phys[MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
			   double der_phys[MAX_EQN_NUM][MAX_QUAD_PTS_NUM], int sln=0);
    void get_solution_plot(double x_phys[MAX_PLOT_PTS_NUM], int pts_num,
         double val_phys[MAX_EQN_NUM][MAX_PLOT_PTS_NUM], 
			   double der_phys[MAX_EQN_NUM][MAX_PLOT_PTS_NUM], int sln=0);
    void get_solution_point(double x_phys, 
			    double val[MAX_EQN_NUM], double der[MAX_EQN_NUM], int sln=0);
    int create_cand_list(int adapt_type, int p_ref_left, int p_ref_right, int3 *cand_list);
    void print_cand_list(int num_cand, int3 *cand_list);
    void refine(int3 cand);
    void refine(int type, int p_left, int p_right);
    unsigned is_active();
    unsigned active;   // flag used by assembling algorithm
    double x1, x2;     // endpoints
    int p;             // poly degrees
    int marker;        // can be used to distinguish between material parameters
    int n_eq;          // number of equations (= number of solution components)
    int n_sln;         // number of solution copies
    int dof[MAX_EQN_NUM][MAX_P + 1];   // connectivity array of length p+1 
                                       // for every solution component
    double coeffs[MAX_SLN_NUM][MAX_EQN_NUM][MAX_P + 1];   // solution coefficient array of length p+1 
                                                          // for every component and every solution 
    int id;
    unsigned level;    // refinement level (zero for initial mesh elements) 
    Element *sons[2];  // for refinement
};

typedef Element* ElemPtr2[2];

class Mesh;

void copy_mesh_to_vector(Mesh *mesh, double *y, int sln=0);
void copy_vector_to_mesh(double *y, Mesh *mesh, int sln=0);

class Mesh {
    public:
        Mesh();
        // Creates equidistant mesh with uniform polynomial degree of elements.
        // All elements will have the same (zero) marker.
        Mesh(double a, double b, int n_elem, int p_init=1, int n_eq=1, int
                n_sln=1, bool print_banner=true);
        // Creates a general mesh (used, e.g., in example "neutronics").
        // n_macro_elem... number of macro elements
        // pts_array[]...  array of macroelement grid points
        // p_array[]...    array of macroelement poly degrees
        // m_array[]...    array of macroelement material markers
        // div_array[]...  array of macroelement equidistant divisions
        Mesh(int n_macro_elem, double *pts_array, int *p_array, int *m_array, 
             int *div_array, int n_eq=1, int n_sln=1, bool print_banner=true);
        ~Mesh() {
            if (this->base_elems != NULL) {
                delete[] this->base_elems;
            }
        }
        void free_elements() {
            if (this->base_elems != NULL) {
                delete[] this->base_elems;
            }
        }
        int assign_dofs();
        Element *get_base_elems() {
            return this->base_elems;
        }
        int get_n_base_elem() {
            return this->n_base_elem;
        }
        void set_n_base_elem(int n_base_elem) {
            this->n_base_elem = n_base_elem;
        }
        int get_n_active_elem() {
            return this->n_active_elem;
        }
        void set_n_active_elem(int n) {
            this->n_active_elem = n;
        }
        int get_n_dof() {
            return this->n_dof;
        }
        void set_n_dof(int n) {
            this->n_dof = n;
        }
        int get_n_eq() {
            return this->n_eq;
        }
        int set_n_eq(int n_eq) {
            this->n_eq = n_eq;
        }
        int get_n_sln() {
            return this->n_sln;
        }
        void set_n_sln(int n_sln) {
            this->n_sln = n_sln;
        }
        double get_left_endpoint() {
            return this->left_endpoint; 
        }
        void set_left_endpoint(double a) {
            this->left_endpoint = a; 
        }
        double get_right_endpoint() {
            return this->right_endpoint; 
        }
        void set_right_endpoint(double b) {
            this->right_endpoint = b; 
        }
        Element* first_active_element();
        Element* last_active_element();
        void set_bc_left_dirichlet(int eqn, double val);
        void set_bc_right_dirichlet(int eqn, double val);
        void refine_single_elem(int id, int3 cand);
        void refine_elems(int elem_num, int *id_array, int3 *cand_array);
        void reference_refinement(int start_elem_id, int elem_num);
        Mesh *replicate(); 
        void plot(const char* filename); // plots the mesh and polynomial degrees of elements
        void plot_element_error_p(int norm, FILE *f, Element *p, Element *e_ref,  
                                  int subdivision = 20); // plots error wrt. reference solution
        void plot_element_error_hp(int norm, FILE *f, Element *p, 
                                   Element *e_ref_left, Element *e_ref_right, 
                                   int subdivision = 20); // plots error wrt. reference solution
                                                          // if ref. refinement was hp-refinement
        void plot_element_error_exact(int norm, FILE *f, Element *p, 
				   exact_sol_type exact_sol,
                                   int subdivision = 20); // plots error wrt. reference solution
                                                          // if ref. refinement was hp-refinement
        void plot_error_estimate(int norm, Mesh* mesh_ref, const char *filename, 
                        int subdivision = 500);  // plots error wrt. reference solution
        void plot_error_estimate(int norm, ElemPtr2* elem_ref_pairs, const char *filename,  
                        int subdivision = 500);  // plots error wrt. reference solution
        void plot_error_exact(int norm, exact_sol_type exact_sol, const char *filename,  
                        int subdivision = 500); // plots error wrt. exact solution
        int assign_elem_ids();
        int n_active_elem;
        void copy_vector_to_mesh(double *y, int sln=0) {
            ::copy_vector_to_mesh(y, this, sln);
        }
        void copy_mesh_to_vector(double *y, int sln=0) {
            ::copy_mesh_to_vector(this, y, sln);
        }

    private:
        double left_endpoint, right_endpoint;
        int n_eq;            // number of equations in the system
        int n_sln;           // number of solution copies
        int n_base_elem;     // number of elements in the base mesh
        int n_dof;           // number of DOF (in each solution copy)
        Element *base_elems; // base mesh

};

// Returns updated coarse and reference meshes, with the last 
// coarse and reference mesh solutions on them, respectively. 
// The coefficient vectors and numbers of degrees of freedom 
// on both meshes are also updated. 
void adapt(int norm, int adapt_type, double threshold, 
           double *err_squared_array,
           Mesh* &mesh, Mesh* &mesh_ref);

// Returns updated coarse mesh, with the last 
// coarse solution on it. 
// The coefficient vector and number of degrees of freedom 
// also is updated. 
void adapt(int norm, int adapt_type, double threshold, 
           double *err_array, 
           Mesh* &mesh, ElemPtr2 *ref_elem_pairs);

void adapt_plotting(Mesh *mesh, Mesh *mesh_ref,
                    int norm, int exact_sol_provided, 
                    exact_sol_type exact_sol); 

void adapt_plotting(Mesh *mesh, ElemPtr2* ref_elem_pairs,
                    int norm, int exact_sol_provided, 
                    exact_sol_type exact_sol); 
#endif
