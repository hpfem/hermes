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

#include "adapt.h"
#include "matrix.h"
#include "graph.h"
#include "forms.h"
#include "weakform.h"
#include "views/view.h"
#include "views/scalar_view.h"
#include "views/vector_view.h"
#include "views/order_view.h"
#include "ref_selectors/selector.h"
#include <map>

typedef enum {H2D_L2_NORM, H2D_H1_NORM, H2D_HCURL_NORM, H2D_HDIV_NORM} ProjNormType;

class Space;
class PrecalcShapeset;
class WeakForm;
class Matrix;
class SparseMatrix;
class Vector;
class Solver;

// Default H2D projection norm in H1 norm.
extern int H2D_DEFAULT_PROJ_NORM;

/// Instantiated template. It is used to create a clean Windows DLL interface.
H2D_API_USED_TEMPLATE(Tuple<int>);
H2D_API_USED_TEMPLATE(Tuple<Space*>);
H2D_API_USED_TEMPLATE(Tuple<MeshFunction*>);
H2D_API_USED_TEMPLATE(Tuple<Solution*>);
H2D_API_USED_TEMPLATE(Tuple<PrecalcShapeset*>);

/// For projection, the user may provide bi/linear forms for each solution component stored
/// in tuples of following types
typedef Tuple< std::pair<WeakForm::matrix_form_val_t, WeakForm::matrix_form_ord_t> > matrix_forms_tuple_t;
typedef Tuple< std::pair<WeakForm::vector_form_val_t, WeakForm::vector_form_ord_t> > vector_forms_tuple_t;


/// Finite Element problem class
///
/// This class does assembling into passed-in structures.
///
class H2D_API FeProblem {
public:
  FeProblem(WeakForm *wf, Tuple<Space *> spaces, bool is_linear = false);
  virtual ~FeProblem();
  void free();

  Space* get_space(int n) {  return this->spaces[n];  }
  PrecalcShapeset* get_pss(int n) {  return this->pss[n];  }

  void create(SparseMatrix* mat, Vector* rhs = NULL, bool rhsonly = false);

  // General assembling procedure for nonlinear problems. coeff_vec is the 
  // previous Newton vector.
  void assemble(scalar* coeff_vec, SparseMatrix* mat, Vector* rhs,
                bool rhsonly = false);
  // Assembling for linear problems. Same as the previous functions, but 
  // does not need the coeff_vector.
  void assemble(SparseMatrix* mat, Vector* rhs, bool rhsonly = false);

  int get_num_dofs();
  bool is_matrix_free() { return wf->is_matrix_free(); }
  void invalidate_matrix() { have_matrix = false; }

protected:
  WeakForm *wf;

  bool is_linear;

  int ndof;
  int *sp_seq;
  int wf_seq;
  Tuple<Space *> spaces;
  PrecalcShapeset** pss;

  int num_user_pss;
  bool values_changed;
  bool struct_changed;
  bool have_spaces;
  bool have_matrix;
  bool is_up_to_date();

  scalar** buffer;
  int mat_size;

  scalar** get_matrix_buffer(int n)
  {
    if (n <= mat_size) return buffer;
    if (buffer != NULL) delete [] buffer;
    return (buffer = new_matrix<scalar>(mat_size = n));
  }

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
         PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv, EdgePos* ep);
  scalar eval_form(WeakForm::VectorFormSurf *vfv, Tuple<Solution *> u_ext, 
         PrecalcShapeset *fv, RefMap *rv, EdgePos* ep);

};

H2D_API int get_num_dofs(Tuple<Space *> spaces);

// Underlying function for global orthogonal projection.
// Not intended for the user. NOTE: the weak form here must be 
// a special projection weak form, which is different from 
// the weak form of the PDE. If you supply a weak form of the 
// PDE, the PDE will just be solved. 
void project_internal(Tuple<Space *> spaces, WeakForm *proj_wf, scalar* target_vec);

H2D_API void project_global(Tuple<Space *> spaces, Tuple<int> proj_norms, Tuple<MeshFunction *> source_meshfns, 
                    scalar* target_vec);

H2D_API void project_global(Tuple<Space *> spaces, Tuple<int> proj_norms, Tuple<Solution*> sols_src, Tuple<Solution*> sols_dest);

H2D_API void project_global(Tuple<Space *> spaces, matrix_forms_tuple_t proj_biforms, 
                    vector_forms_tuple_t proj_liforms, Tuple<MeshFunction*> source_meshfns, 
                    scalar* target_vec);

H2D_API void project_global(Space *space, 
                    std::pair<WeakForm::matrix_form_val_t, WeakForm::matrix_form_ord_t> proj_biform,
                    std::pair<WeakForm::vector_form_val_t, WeakForm::vector_form_ord_t> proj_liform,
                    ExactFunction source_fn, scalar* target_vec);

H2D_API void project_global(Space *space, ExactFunction2 source_fn, scalar* target_vec);

/// Selects the appropriate linear solver.
H2D_API Vector* create_vector(MatrixSolverType matrix_solver);
H2D_API SparseMatrix* create_matrix(MatrixSolverType matrix_solver);
H2D_API Solver* create_solver(MatrixSolverType matrix_solver, Matrix* matrix, Vector* rhs);

/// Basic Newton's loop. Takes a coefficient vector, delivers a coefficient vector (in the 
/// same variable "init_coeff_vector").
H2D_API bool solve_newton(Tuple<Space *> spaces, WeakForm* wf, scalar* coeff_vec, 
                          MatrixSolverType matrix_solver, double newton_tol, 
                          int newton_max_iter, bool verbose);

// Solve a typical linear problem (without automatic adaptivity).
// Feel free to adjust this function for more advanced applications.
H2D_API bool solve_linear(Tuple<Space *> spaces, WeakForm* wf, MatrixSolverType matrix_solver, 
                          Tuple<Solution *> solutions, scalar*coeff_vec = NULL);

// Do initial checks for the adaptivity process.
H2D_API void lin_adapt_begin(Tuple<Space *> spaces, Tuple<RefinementSelectors::Selector *> selectors, Tuple<int> proj_norms, TimePeriod * cpu_time);

// Create globally refined space.
H2D_API Tuple<Space *>* construct_refined_spaces(Tuple<Space *> coarse, int order_increase = 1);
H2D_API Space* construct_refined_space(Space* coarse, int order_increase = 1);
#endif



