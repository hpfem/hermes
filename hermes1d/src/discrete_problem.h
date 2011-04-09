// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef _DISCRETE_H_
#define _DISCRETE_H_

#include <vector>

#include "space.h"
#include "quad_std.h"
#include "legendre.h"
#include "lobatto.h"
#include "weakform.h"
#include "../../hermes_common/matrix.h"
#include "../../hermes_common/common.h"
#include "../../hermes_common/solver/solver.h"
#include "../../hermes_common/solver/dpinterface.h"
#include "../../hermes_common/tables.h"
#include "iterator.h"

class HERMES_API DiscreteProblem : public DiscreteProblemInterface {
public:
  DiscreteProblem(WeakForm* wf, Space* space, bool is_linear = true);

  void process_surf_forms(SparseMatrix *mat, Vector *res, int bdy_index);
  
  void process_vol_forms(SparseMatrix *mat, Vector *res);

  void assemble(scalar* coeff_vec, SparseMatrix *mat, Vector *rhs = NULL, 
                bool force_diagonal_blocks = false,
                bool add_dir_lift = true, Table* block_weights = NULL);

  int get_num_dofs();
  
  bool is_matrix_free();
  
  void create_sparse_structure(SparseMatrix* matrix, Vector* rhs = NULL, 
                               bool force_diagonal_blocks = false, Table* block_weights = NULL);
                               
  void invalidate_matrix() { return; }
private:
  WeakForm* wf;
  Space* space;
  bool is_linear;
};

// return coefficients for all shape functions on the element m,
// for all solution components
void calculate_elem_coeffs(Element *e, double **coeffs, int n_eq);

void element_quadrature(double a, double b, 
                        int order, double *pts, double *weights, int *num);

void element_shapefn(double a, double b, 
		     int k, int order, double *val, double *der);

void element_shapefn_point(double x_ref, double a, double b, 
			   int k, double &val, double &der);

void HERMES_API jfnk_cg(DiscreteProblem *dp, Space *space,
             double matrix_solver_tol, int matrix_solver_maxiter,  
	     double jfnk_epsilon, double jfnk_tol, int jfnk_maxiter, 
             MatrixSolverType matrix_solver, bool verbose = true);

// Splits the indicated elements and 
// increases poly degree in sons by one.
// Solution is transfered to new elements.
HERMES_API Space* construct_refined_space(Space* space, int order_increase = 1);

HERMES_API double get_l2_norm(Vector* vec);

HERMES_API void set_zero(scalar *vec, int length);


#endif
