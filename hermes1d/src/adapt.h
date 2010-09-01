// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef _ADAPT_H_
#define _ADAPT_H_

#include "common.h"
#include "legendre.h"
#include "lobatto.h"
#include "mesh.h"
#include "matrix.h"
#include "iterator.h"

// Calculates the square in L2 or H1 norm of the difference between 
// the coarse and fine mesh solution, for all solution components.
// Assumes that the element was not refined in space for the 
// reference solution. 
// FIXME: to be moved to the Element class
double calc_elem_est_error_squared_p(int norm, Element *e, Element *e_ref);

// Calculates the square in L2 or H1 norm of the difference between 
// the coarse and fine mesh solution, for all solution components.
// Assumes that the element was refined in space for the 
// reference solution.
// FIXME: to be moved to the Element class
double calc_elem_est_error_squared_hp(int norm, Element *e, 
                        Element *e_ref_left, Element *e_ref_right);

// Calculates L2 or H1 norm of the difference between the coarse
// and reference solutions in all active elements of 'mesh'. Total
// error is returned.
double calc_error_estimate(int norm, Mesh* mesh, Mesh* mesh_ref, 
			   double *err_array);

// Calculates L2 or H1 norm of the difference between the coarse
// and reference solutions in all active elements of 'mesh'. Total
// error is returned.
double calc_error_estimate(int norm, Mesh* mesh, 
                           ElemPtr2* ref_element_pairs);

// Can be used for both the coarse and reference solutions
double calc_solution_norm(int norm, Mesh* mesh);
double calc_solution_norm(int norm, Mesh* mesh, ElemPtr2* ref_element_pairs);

// Sort err_array[] and returning array of sorted element indices
void sort_element_errors(int n, double *err_array, int *id_array); 

// Assumes that reference solution is defined on two half-elements 'e_ref_left'
// and 'e_ref_right'. The reference solution is projected onto the space of 
// (discontinuous) polynomials of degree 'p_left' on 'e_ref_left'
// and degree 'p_right' on 'e_ref_right'. Returned is projection error and
// number of dof added by this candidate.
void check_cand_coarse_hp_fine_hp(Element *e, Element *e_ref_left, Element *e_ref_right, 
                                  int p_left, int p_right,
                                  double &err, int &dof);

// Assumes that reference solution is defined on one single element 'e_ref' = 'e'. 
// The reference solution is projected onto the space of (discontinuous) 
// polynomials of degree 'p_left' on the left half of 'e' and degree 
// 'p_right' on the right half of 'e'. Returned is projection error and
// number of dof added by this candidate.
void check_cand_coarse_hp_fine_p(Element *e, Element *e_ref,
                                 int p_left, int p_right,
                                 double &err, int &dof);

// Assumes that reference solution is defined on two half-elements 'e_ref_left'
// and 'e_ref_right'. The reference solution is projected onto the space of 
// polynomials of degree 'p' on 'e'. Returned is projection error and
// number of dof added by this candidate.
void check_cand_coarse_p_fine_hp(Element *e, Element *e_ref_left, Element *e_ref_right, 
                                 int p, double &err, int &dof);

// Assumes that reference solution is defined on one single element 
// 'e_ref' (reference refinement did not split the element in space). 
// The reference solution is projected onto the space of 
// polynomials of degree 'p' on 'e'. Returned is projection error and
// number of dof added by this candidate.
void check_cand_coarse_p_fine_p(Element *e, Element *e_ref, int p,
                                double &err, int &dof);

// Error wrt. exact solution (if provided) on element 'e' 
double calc_elem_exact_error_squared(int norm, exact_sol_type exact_sol,
                                     Element *e, int order);

// Error wrt. exact solution (if provided) on the entire interval (A, B) 
double calc_error_exact(int norm, Mesh *mesh, 
                        exact_sol_type exact_sol); 

// Calculates L2 or H1 norm of function exact_sol in interval (A, B)
double calc_solution_norm(int norm, exact_sol_type exact_sol, 
                           int n_eq, 
                           double A, double B, int subdivision, 
                           int order);

// Selects best hp-refinement from the given list (distinguishes whether 
// the reference refinement on that element was p- or hp-refinement). 
// Each refinement candidate is a triple of integers. First one means 
// p-refinement (0) or hp-refinement (1). Second and/or third number are 
// the new proposed polynomial degrees.
int select_hp_refinement(Element *e, Element *e_ref, Element * e_ref2,
                         int num_cand, int3 *cand_list, int ref_sol_type, 
                         int norm); 


// Refine coarse mesh elements whose id_array >= 0, and 
// adjust the reference mesh accordingly.  
// Returns updated coarse and reference meshes, with the last 
// coarse and reference mesh solutions on them, respectively. 
// The coefficient vectors and numbers of degrees of freedom 
// on both meshes are also updated. 
void adapt(int norm, int adapt_type, double threshold, 
           double *err_squared_array,
           Mesh* &mesh, Mesh* &mesh_ref, 
           double * &y_prev, double* &y_prev_ref, 
           int &n_dof, int &n_dof_ref);

void adapt_plotting(Mesh *mesh, Mesh *mesh_ref, 
              double *y_prev, double *y_prev_ref,
              int norm, int exact_sol_provided, 
              exact_sol_type exact_sol); 




#endif
