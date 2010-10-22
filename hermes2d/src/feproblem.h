// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.


#ifndef __H2D_FEPROBLEM_H
#define __H2D_FEPROBLEM_H

#include "../../hermes_common/matrix.h"
#include "adapt/adapt.h"
#include "graph.h"
#include "forms.h"
#include "weakform.h"
#include "views/view.h"
#include "views/scalar_view.h"
#include "views/vector_view.h"
#include "views/order_view.h"
#include "ref_selectors/selector.h"
#include <map>

class Space;
class PrecalcShapeset;
class WeakForm;
class Matrix;
class SparseMatrix;
class Vector;
class Solver;

/// Instantiated template. It is used to create a clean Windows DLL interface.
HERMES_API_USED_TEMPLATE(Tuple<ProjNormType>);
HERMES_API_USED_TEMPLATE(Tuple<Space*>);
HERMES_API_USED_TEMPLATE(Tuple<MeshFunction*>);
HERMES_API_USED_TEMPLATE(Tuple<Solution*>);
HERMES_API_USED_TEMPLATE(Tuple<PrecalcShapeset*>);


/// Discrete problem class
///
/// This class does assembling into external matrix / vactor structures.
///
class HERMES_API DiscreteProblem 
{
public:
  DiscreteProblem(WeakForm* wf, Tuple<Space *> spaces, bool is_linear = false);
  virtual ~DiscreteProblem();
  void free();

  // Get pointer to n-th space.
  Space* get_space(int n) {  return this->spaces[n];  }

  // This is different from H3D.
  PrecalcShapeset* get_pss(int n) {  return this->pss[n];  }

  // Precalculate matrix sparse structure.
  void create(SparseMatrix* mat, Vector* rhs = NULL, bool rhsonly = false);

  // General assembling procedure for nonlinear problems. coeff_vec is the 
  // previous Newton vector.
  void assemble(scalar* coeff_vec, SparseMatrix* mat, Vector* rhs, bool rhsonly = false);

  // Assembling for linear problems. Same as the previous functions, but 
  // does not need the coeff_vector.
  void assemble(SparseMatrix* mat, Vector* rhs, bool rhsonly = false);

  // Get the number of unknowns.
  int get_num_dofs();

  bool is_matrix_free() { return wf->is_matrix_free(); }

  void invalidate_matrix() { have_matrix = false; }

protected:
  WeakForm* wf;

  bool is_linear;

  int ndof;
  int *sp_seq;
  int wf_seq;
  Tuple<Space *> spaces;

  scalar** matrix_buffer;                /// buffer for holding square matrix (during assembling)
  int matrix_buffer_dim;                 /// dimension of the matrix held by 'matrix_buffer'
  inline scalar** get_matrix_buffer(int n);

  bool have_spaces;
  bool have_matrix;

  bool values_changed;
  bool struct_changed;
  bool is_up_to_date();

  PrecalcShapeset** pss;    // This is different from H3D.
  int num_user_pss;         // This is different from H3D.

  ExtData<Ord>* init_ext_fns_ord(std::vector<MeshFunction *> &ext);
  ExtData<Ord>* init_ext_fns_ord(std::vector<MeshFunction *> &ext, int edge);
  ExtData<scalar>* init_ext_fns(std::vector<MeshFunction *> &ext, RefMap *rm, const int order);
  Func<double>* get_fn(PrecalcShapeset *fu, RefMap *rm, const int order);

  // Caching transformed values for element
  std::map<PrecalcShapeset::Key, Func<double>*, PrecalcShapeset::Compare> cache_fn;
  Geom<double>* cache_e[g_max_quad + 1 + 4 * g_max_quad + 4];
  double* cache_jwt[g_max_quad + 1 + 4 * g_max_quad + 4];

  void init_cache();
  void delete_cache();

  scalar eval_form(WeakForm::MatrixFormVol *mfv, Tuple<Solution *> u_ext, 
         PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv);
  scalar eval_form(WeakForm::VectorFormVol *vfv, Tuple<Solution *> u_ext, 
         PrecalcShapeset *fv, RefMap *rv);
  scalar eval_form(WeakForm::MatrixFormSurf *mfv, Tuple<Solution *> u_ext, 
         PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv, SurfPos* surf_pos);
  scalar eval_form(WeakForm::VectorFormSurf *vfv, Tuple<Solution *> u_ext, 
         PrecalcShapeset *fv, RefMap *rv, SurfPos* surf_pos);

};

/// Calls the required (de)initialization routines of the selected matrix solver.
HERMES_API bool initialize_solution_environment(MatrixSolverType matrix_solver, int argc, char* argv[]);
HERMES_API bool finalize_solution_environment(MatrixSolverType matrix_solver);

/// Selects the appropriate linear solver.
HERMES_API Vector* create_vector(MatrixSolverType matrix_solver);
HERMES_API SparseMatrix* create_matrix(MatrixSolverType matrix_solver);
HERMES_API Solver* create_linear_solver(MatrixSolverType matrix_solver, Matrix* matrix, Vector* rhs);

// Create globally refined space.
HERMES_API Tuple<Space *>* construct_refined_spaces(Tuple<Space *> coarse, int order_increase = 1);
HERMES_API Space* construct_refined_space(Space* coarse, int order_increase = 1);

HERMES_API double get_l2_norm(Vector* vec); 

#endif



