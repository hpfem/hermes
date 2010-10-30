// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef _ADAPT_H_
#define _ADAPT_H_

#include "../../hermes_common/common.h"
#include "../../hermes_common/matrix.h"
#include "iterator.h"

// Calculates the square in L2 or H1 norm of the difference between 
// the coarse and fine mesh solution, for all solution components.
// Assumes that the element was not refined in space for the 
// reference solution. 
// FIXME: to be moved to the Element class
double HERMES_API calc_elem_est_error_squared_p(int norm, Element *e, Element *e_ref,
        int sln=0);

// Calculates the square in L2 or H1 norm of the difference between 
// the coarse and fine mesh solution, for all solution components.
// Assumes that the element was refined in space for the 
// reference solution.
// FIXME: to be moved to the Element class
double HERMES_API calc_elem_est_error_squared_hp(int norm, Element *e, 
                        Element *e_ref_left, Element *e_ref_right,
                        int sln=0);

// Calculates L2 or H1 norm of the difference between the coarse
// and reference solutions in all active elements of 'space'. Total
// error is returned.
double HERMES_API calc_err_est(int norm, Space* space, Space* space_ref, 
			   double *err_array = NULL, int sln=0);

// Calculates L2 or H1 norm of the difference between the coarse
// and reference solutions in all active elements of 'space'. Total
// error is returned.
double HERMES_API calc_err_est(int norm, Space* space, 
                           ElemPtr2* ref_element_pairs);

// Can be used for both the coarse and reference solutions
double HERMES_API calc_solution_norm(int norm, Space* space);
double HERMES_API calc_solution_norm(int norm, Space* space, ElemPtr2* ref_element_pairs);

// Sort err_array[] and returning array of sorted element indices
void HERMES_API sort_element_errors(int n, double *err_array, int *id_array); 

// Assumes that reference solution is defined on two half-elements 'e_ref_left'
// and 'e_ref_right'. The reference solution is projected onto the space of 
// (discontinuous) polynomials of degree 'p_left' on 'e_ref_left'
// and degree 'p_right' on 'e_ref_right'. Returned is projection error and
// number of dof added by this candidate.
void HERMES_API check_cand_coarse_hp_fine_hp(Element *e, Element *e_ref_left, Element *e_ref_right, 
                                  int p_left, int p_right,
                                  double &err, int &dof);

// Assumes that reference solution is defined on one single element 'e_ref' = 'e'. 
// The reference solution is projected onto the space of (discontinuous) 
// polynomials of degree 'p_left' on the left half of 'e' and degree 
// 'p_right' on the right half of 'e'. Returned is projection error and
// number of dof added by this candidate.
void HERMES_API check_cand_coarse_hp_fine_p(Element *e, Element *e_ref,
                                 int p_left, int p_right,
                                 double &err, int &dof);

// Assumes that reference solution is defined on two half-elements 'e_ref_left'
// and 'e_ref_right'. The reference solution is projected onto the space of 
// polynomials of degree 'p' on 'e'. Returned is projection error and
// number of dof added by this candidate.
void HERMES_API check_cand_coarse_p_fine_hp(Element *e, Element *e_ref_left, Element *e_ref_right, 
                                 int p, double &err, int &dof);

// Assumes that reference solution is defined on one single element 
// 'e_ref' (reference refinement did not split the element in space). 
// The reference solution is projected onto the space of 
// polynomials of degree 'p' on 'e'. Returned is projection error and
// number of dof added by this candidate.
void HERMES_API check_cand_coarse_p_fine_p(Element *e, Element *e_ref, int p,
                                double &err, int &dof);

// Error wrt. exact solution (if provided) on element 'e' 
double HERMES_API calc_elem_exact_error_squared(int norm, exact_sol_type exact_sol,
                                     Element *e, int order);

// Error wrt. exact solution (if provided) on the entire interval (A, B) 
double HERMES_API calc_err_exact(int norm, Space *space, 
                        exact_sol_type exact_sol, int neq, double a, double b); 

// Calculates L2 or H1 norm of function exact_sol in interval (A, B)
double HERMES_API calc_solution_norm(int norm, exact_sol_type exact_sol, 
                           int n_eq, 
                           double A, double B, int subdivision, 
                           int order);

// Selects best hp-refinement from the given list (distinguishes whether 
// the reference refinement on that element was p- or hp-refinement). 
// Each refinement candidate is a triple of integers. First one means 
// p-refinement (0) or hp-refinement (1). Second and/or third number are 
// the new proposed polynomial degrees.
int HERMES_API select_hp_refinement(Element *e, Element *e_ref, Element * e_ref2,
                         int num_cand, int3 *cand_list, int ref_sol_type, 
                         int norm); 


#endif
