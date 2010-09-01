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

#ifndef __H2D_DISCRETE_PROBLEM_H
#define __H2D_DISCRETE_PROBLEM_H

#include "common.h"
#include "matrix.h"
#include "matrix_old.h"
#include "forms.h"
#include "weakform.h"
#include <map>
#include "views/scalar_view.h"
#include "views/vector_view.h"
#include "views/order_view.h"
#include "function.h"
#include "ref_selectors/selector.h"
#include "graph.h"
#include "adapt.h"
#include "norm.h"

typedef enum {H2D_L2_NORM, H2D_H1_NORM, H2D_HCURL_NORM, H2D_HDIV_NORM} ProjNormType;

class Space;
class PrecalcShapeset;
class WeakForm;
class CommonSolver;

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

///
///
///
///
///
class H2D_API DiscreteProblem
{
public:

  DiscreteProblem();
  DiscreteProblem(WeakForm* wf_);
  DiscreteProblem(WeakForm* wf_, Space* s_);
  DiscreteProblem(WeakForm* wf_, Tuple<Space*> spaces_);
  virtual ~DiscreteProblem();

  void init(WeakForm* wf, CommonSolver* solver);
  void init_spaces(Tuple<Space*> spaces);
  void init_space(Space* s);         // single equation case
  void set_spaces(Tuple<Space*> spaces);
  void set_pss(Tuple<PrecalcShapeset*> pss);
  void set_pss(PrecalcShapeset* p);  // single equation case
  void copy(DiscreteProblem* sys);
  Space* get_space(int n) {
      if (n < 0 || n >= this->wf->neq) error("Bad index of space.");
      return this->spaces[n];
  }
  Space* get_space() {
      if (this->wf->neq != 1) error("get_space() can be used in single PDE case only.");
      return this->spaces[0];
  }
  Mesh* get_mesh(int n) {
      if (n < 0 || n >= this->wf->neq) error("Bad index of mesh.");
      return this->spaces[n]->mesh;
  }
  PrecalcShapeset* get_pss(int n) {
      if (n < 0 || n >= this->wf->neq) error("Bad index of precalc shapeset.");
      return this->pss[n];
  }
  PrecalcShapeset* get_pss() {
      return this->pss[0];
  }

  /// Assembles the matrix A and vectors Vec, Dir and RHS, and exposes them to the user. 
  /// Everything must be allocated in advance when assemble() is called. This is the generic 
  /// functionality to be used for linear problems, nonlinear problems, and eigenproblems.
  /// Soon this will be extended to assemble an arbitrary number of matrix and vector
  /// weak forms. 
  virtual void assemble(Vector* init_vec, Matrix* mat_ext, Vector* dir_ext, Vector* rhs_ext, 
                        bool rhsonly = false, bool is_complex = false);

  /// Basic function that just solves the matrix problem. The right-hand
  /// side enters through "vec" and the result is stored in "vec" as well. 
  bool solve_matrix_problem(Matrix* mat, Vector* vec); 

  /// Solves the matrix problem with "mat" and "rhs", and adds the result 
  /// to the vector "vec".
  virtual bool solve(Matrix* mat, Vector* rhs, Vector* vec);

  /// Frees the stiffness matrix, coefficient vectors, and matrix solver data.
  virtual void free();

  /// Saves the stiffness matrix in various formats.

  int get_num_dofs();
  int get_num_dofs(int i) {
    if (this->spaces[i] == NULL) error("spaces[%d] is NULL in DiscreteProblem::get_num_dofs().", i);
    return this->spaces[i]->get_num_dofs();
  }
  int get_num_spaces() { return this->wf->neq; };
  int get_matrix_size(Matrix* mat_ext);
  //void get_matrix(int*& Ap, int*& Ai, scalar*& Ax, int& size);
  //void get_rhs(scalar*& RHS, int& size) { RHS = this->RHS; size=this->get_num_dofs(); }
  //void get_solution_vector(scalar*& sln_vector, int& sln_vector_len)
  //     { sln_vector = Vec; sln_vector_len = this->get_num_dofs(); }

  /// Returns a copy of the solution vector.
  void get_solution_vector(std::vector<scalar>& sln_vector_out);

  /// Assigning DOF = enumerating basis functions in the FE spaces.
  int assign_dofs();  // all spaces

  /// Needed for problems where BC depend on time.
  void update_essential_bc_values();

  /// Frees spaces. Called automatically on destruction.
  void free_spaces();

  Tuple<Space *> spaces;
  WeakForm* wf;
  bool have_spaces;

  /// Matrix solvers.
  CommonSolver* solver;
  CommonSolver* solver_default;

protected:

  PrecalcShapeset** pss;

  bool mat_sym; ///< true if symmetric - then only upper half stored

  int RHS_length;
  int Vec_length;
  int Dir_length;

  void insert_block(Matrix *A, scalar** mat, int* iidx, int* jidx,
          int ilen, int jlen);

  ExtData<Ord>* init_ext_fns_ord(std::vector<MeshFunction *> &ext);
  ExtData<Ord>* init_ext_fns_ord(std::vector<MeshFunction *> &ext, int edge);
  ExtData<scalar>* init_ext_fns(std::vector<MeshFunction *> &ext, RefMap *rm, const int order);
  Func<double>* get_fn(PrecalcShapeset *fu, RefMap *rm, const int order);

  // Key for caching transformed function values on elements
  struct Key
  {
    int index;
    int order;
    int sub_idx;
    int shapeset_type;

    Key(int index, int order, int sub_idx, int shapeset_type)
    {
      this->index = index;
      this->order = order;
      this->sub_idx = sub_idx;
      this->shapeset_type = shapeset_type;
    }
  };

  struct Compare
  {
    bool operator()(Key a, Key b) const
    {
      if (a.index < b.index) return true;
      else if (a.index > b.index) return false;
      else
      {
        if (a.order < b.order) return true;
        else if (a.order > b.order) return false;
        else
        {
          if (a.sub_idx < b.sub_idx) return true;
          else if (a.sub_idx > b.sub_idx) return false;
          else
          {
            if (a.shapeset_type < b.shapeset_type) return true;
            else return false;
          }
        }
      }
    }
  };

  // Caching transformed values for element
  std::map<Key, Func<double>*, Compare> cache_fn;
  Geom<double>* cache_e[g_max_quad + 1 + 4 * g_max_quad + 4];
  double* cache_jwt[g_max_quad + 1 + 4 * g_max_quad + 4];

  void init_cache();
  void delete_cache();

  // evaluation of forms, general case
  scalar eval_form(WeakForm::MatrixFormVol *bf, Tuple<Solution *> sln, PrecalcShapeset *fu, 
                   PrecalcShapeset *fv, RefMap *ru, RefMap *rv);
  scalar eval_form(WeakForm::VectorFormVol *lf, Tuple<Solution *> sln, PrecalcShapeset *fv, 
                   RefMap *rv);
  scalar eval_form(WeakForm::MatrixFormSurf *bf, Tuple<Solution *> sln, PrecalcShapeset *fu, 
                   PrecalcShapeset *fv, RefMap *ru, RefMap *rv, EdgePos* ep);
  scalar eval_form(WeakForm::VectorFormSurf *lf, Tuple<Solution *> sln, PrecalcShapeset *fv, 
                   RefMap *rv, EdgePos* ep);

  scalar** get_matrix_buffer(int n)
  {
    if (n <= mat_size) return buffer;
    if (buffer != NULL) delete [] buffer;
    return (buffer = new_matrix<scalar>(mat_size = n));
  }

  scalar** buffer;
  int mat_size;

  int* sp_seq;
  int wf_seq;
  int num_user_pss;
  bool values_changed;
  bool struct_changed;
};

H2D_API int get_num_dofs(Tuple<Space *> spaces);

H2D_API void init_matrix_solver(MatrixSolverType matrix_solver, int ndof, 
                        Matrix* &mat, Vector* &rhs, 
                        CommonSolver* &solver, bool is_complex = false);

// Underlying function for global orthogonal projection.
// Not intended for the user. NOTE: the weak form here must be 
// a special projection weak form, which is different from 
// the weak form of the PDE. If you supply a weak form of the 
// PDE, the PDE will just be solved. 
void project_internal(Tuple<Space *> spaces, WeakForm *proj_wf, 
                    Tuple<Solution*> target_slns = Tuple<Solution*>(), Vector* target_vec = NULL, bool is_complex = false);

H2D_API void project_global(Tuple<Space *> spaces, Tuple<int> proj_norms, Tuple<MeshFunction *> source_meshfns, 
                    Tuple<Solution*> target_slns = Tuple<Solution*>(), Vector* target_vec = NULL, bool is_complex = false);

H2D_API void project_global(Tuple<Space *> spaces, matrix_forms_tuple_t proj_biforms, 
                    vector_forms_tuple_t proj_liforms, Tuple<MeshFunction*> source_meshfns, 
                    Tuple<Solution*> target_slns = Tuple<Solution*>(),
                    Vector* target_vec = NULL, bool is_complex = false);

H2D_API void project_global(Space *space, 
                    std::pair<WeakForm::matrix_form_val_t, WeakForm::matrix_form_ord_t> proj_biform,
                    std::pair<WeakForm::vector_form_val_t, WeakForm::vector_form_ord_t> proj_liform,
                    ExactFunction source_fn, Solution* target_sln = NULL,
                    Vector* target_vec = NULL, bool is_complex = false);

H2D_API void project_global(Space *space, ExactFunction2 source_fn, Solution* target_sln = NULL, Vector* target_vec = NULL, 
                    bool is_complex = false);

/// Basic Newton's loop. Takes a coefficient vector, delivers a coefficient vector (in the 
/// same variable "init_coeff_vector").
H2D_API bool solve_newton(Tuple<Space *> spaces, WeakForm* wf, Vector* init_coeff_vec,
                  MatrixSolverType matrix_solver, double newton_tol = 1e-5, 
                  int newton_max_iter = 100, bool verbose = false, bool is_complex = false);

// Solve a typical nonlinear problem using the Newton's method and 
// automatic adaptivity. 
// Feel free to adjust this function for more advanced applications.
H2D_API bool solve_newton_adapt(Tuple<Space *> spaces, WeakForm* wf, Vector *coeff_vec, 
                        MatrixSolverType matrix_solver, Tuple<int>proj_norms, 
                        Tuple<Solution *> slns, Tuple<Solution *> ref_slns, 
                        Tuple<WinGeom *> sln_win_geom, Tuple<WinGeom *> mesh_win_geom, 
                        Tuple<RefinementSelectors::Selector *> selectors, AdaptivityParamType* apt,
                        double newton_tol_coarse, double newton_tol_fine, int newton_max_iter, 
                        bool verbose = false, Tuple<ExactSolution *> exact_slns = Tuple<ExactSolution *>(), 
                        bool is_complex = false);

#endif
