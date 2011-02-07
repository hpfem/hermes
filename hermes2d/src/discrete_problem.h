/// This file is part of Hermes2D.
///
/// Hermes2D is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 2 of the License, or
/// (at your option) any later version.
///
/// Hermes2D is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with Hermes2D. If not, see <http:///www.gnu.org/licenses/>.

#define HERMES_REPORT_INFO
#define HERMES_REPORT_WARN

#ifndef __H2D_FEPROBLEM_H
#define __H2D_FEPROBLEM_H

#include "../../hermes_common/matrix.h"
#include "../../hermes_common/solver/solver.h"
#include "adapt/adapt.h"
#include "graph.h"
#include "weakform/forms.h"
#include "weakform/weakform.h"
#include "views/view.h"
#include "views/scalar_view.h"
#include "views/vector_view.h"
#include "views/order_view.h"
#include "function/function.h"
#include "neighbor.h"
#include "tables.h"
#include "ref_selectors/selector.h"
#include <map>

class Space;
class PrecalcShapeset;
class WeakForm;
class Matrix;
class SparseMatrix;
class Vector;
class Solver;

/// Discrete problem class.
///
/// This class does assembling into external matrix / vactor structures.
///
class HERMES_API DiscreteProblem
{
public:
  DiscreteProblem(WeakForm* wf, Hermes::vector<Space *> spaces, bool is_linear = false);
  virtual ~DiscreteProblem();
  void free();

  /// Get pointer to n-th space.
  Space* get_space(int n) { return this->spaces[n]; }

  /// Get whether the DiscreteProblem is linear.
  bool get_is_linear() { return is_linear;};

  /// Get the weak forms.
  WeakForm* get_weak_formulation() { return this->wf;};

  /// Get all spaces as a Hermes::vector.
  Hermes::vector<Space *> get_spaces() {return this->spaces;}

  /// Get the number of spaces.
  int get_num_spaces() {return this->spaces.size();}

  /// This is different from H3D.
  PrecalcShapeset* get_pss(int n) {  return this->pss[n];  }

  /// Precalculate matrix sparse structure.
  /// If force_diagonal_block == true, then (zero) matrix
  /// antries are created in diagonal blocks even if corresponding matrix weak
  /// forms do not exist. This is useful if the matrix is later to be merged with
  /// a matrix that has nonzeros in these blocks. The Table serves for optional
  /// weighting of matrix blocks in systems.
  void create_sparse_structure(SparseMatrix* mat, Vector* rhs = NULL, bool rhsonly = false,
              bool force_diagonal_blocks = false, Table* block_weights = NULL);

  /// Check whether it is sane to assemble.
  /// Throws errors if not.
  void assemble_sanity_checks(Table* block_weights);

  /// Converts coeff_vec to u_ext.
  /// If the supplied coeff_vec is NULL, it supplies the appropriate number of NULL pointers.
  void convert_coeff_vec(scalar* coeff_vec, Hermes::vector<Solution *> & u_ext, bool add_dir_lift);

  /// Initializes psss.
  void initialize_psss(Hermes::vector<PrecalcShapeset *>& spss);

  /// Initializes refmaps.
  void initialize_refmaps(Hermes::vector<RefMap *>& refmap);

  /// General assembling procedure for nonlinear problems. coeff_vec is the
  /// previous Newton vector. If force_diagonal_block == true, then (zero) matrix
  /// antries are created in diagonal blocks even if corresponding matrix weak
  /// forms do not exist. This is useful if the matrix is later to be merged with
  /// a matrix that has nonzeros in these blocks. The Table serves for optional
  /// weighting of matrix blocks in systems. The parameter add_dir_lift decides 
  /// whether Dirichlet lift will be added while coeff_vec is converted into 
  /// Solutions.
  void assemble(scalar* coeff_vec, SparseMatrix* mat, Vector* rhs = NULL, bool rhsonly = false,
		bool force_diagonal_blocks = false, bool add_dir_lift = true, Table* block_weights = NULL);

  /// Assembling for linear problems. Same as the previous functions, but
  /// does not need the coeff_vector.
  void assemble(SparseMatrix* mat, Vector* rhs = NULL, bool rhsonly = false, 
                bool force_diagonal_blocks = false, Table* block_weights = NULL);

  /// Assemble one stage.
  void assemble_one_stage(WeakForm::Stage& stage, 
			  SparseMatrix* mat, Vector* rhs, bool rhsonly, bool force_diagonal_blocks,
                          Table* block_weights,
                          Hermes::vector<PrecalcShapeset *>& spss, Hermes::vector<RefMap *>& refmap, 
                          Hermes::vector<Solution *>& u_ext);

  /// Initialize a state, returns a non-NULL Element.
  Element* init_state(WeakForm::Stage& stage, Hermes::vector<PrecalcShapeset *>& spss, 
                      Hermes::vector<RefMap *>& refmap, Element** e, Hermes::vector<bool>& isempty, 
                      Hermes::vector<AsmList *>& al);

  /// Assemble one state.
  void assemble_one_state(WeakForm::Stage& stage, 
                          SparseMatrix* mat, Vector* rhs, bool rhsonly, bool force_diagonal_blocks, 
                          Table* block_weights, Hermes::vector<PrecalcShapeset *>& spss, 
                          Hermes::vector<RefMap *>& refmap, Hermes::vector<Solution *>& u_ext, Element** e, 
                          bool* bnd, SurfPos* surf_pos, Element* trav_base);

  /// Calculates the integration order for assemble_volume_matrix_forms().
  int calc_order_matrix_form_vol(WeakForm::MatrixFormVol *mfv, Hermes::vector<Solution *> u_ext,
                                 PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv);

  /// Assemble volume matrix forms.
  void assemble_volume_matrix_forms(WeakForm::Stage& stage, 
       SparseMatrix* mat, Vector* rhs, bool rhsonly, bool force_diagonal_blocks, Table* block_weights,
       Hermes::vector<PrecalcShapeset *>& spss, Hermes::vector<RefMap *>& refmap, Hermes::vector<Solution *>& u_ext, 
       Hermes::vector<bool>& isempty, int marker, Hermes::vector<AsmList *>& al);

  /// Calculates the integration order for assemble_volume_vector_forms().
  int calc_order_vector_form_vol(WeakForm::VectorFormVol *mfv, Hermes::vector<Solution *> u_ext,
                                 PrecalcShapeset *fv, RefMap *rv);

  /// Assemble volume vector forms.
  void assemble_volume_vector_forms(WeakForm::Stage& stage, 
       SparseMatrix* mat, Vector* rhs, bool rhsonly, bool force_diagonal_blocks, Table* block_weights,
       Hermes::vector<PrecalcShapeset *>& spss, Hermes::vector<RefMap *>& refmap, Hermes::vector<Solution *>& u_ext, 
       Hermes::vector<bool>& isempty, int marker, Hermes::vector<AsmList *>& al);

  /// Assemble surface and DG forms.
  void assemble_surface_integrals(WeakForm::Stage& stage, 
       SparseMatrix* mat, Vector* rhs, bool rhsonly, bool force_diagonal_blocks, Table* block_weights,
       Hermes::vector<PrecalcShapeset *>& spss, Hermes::vector<RefMap *>& refmap, Hermes::vector<Solution *>& u_ext, 
       Hermes::vector<bool>& isempty, int marker, Hermes::vector<AsmList *>& al, bool bnd, SurfPos& surf_pos, 
       Hermes::vector<bool>& nat, int isurf, Element** e, Element* trav_base, Element* rep_element);

  /// Calculates the integration order for assemble_surface_matrix_forms().
  int calc_order_matrix_form_surf(WeakForm::MatrixFormSurf *mfs, Hermes::vector<Solution *> u_ext,
                                  PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv, SurfPos* surf_pos);
  /// Assemble surface matrix forms.
  void assemble_surface_matrix_forms(WeakForm::Stage& stage, 
       SparseMatrix* mat, Vector* rhs, bool rhsonly, bool force_diagonal_blocks, Table* block_weights,
       Hermes::vector<PrecalcShapeset *>& spss, Hermes::vector<RefMap *>& refmap, Hermes::vector<Solution *>& u_ext, 
       Hermes::vector<bool>& isempty, int marker, Hermes::vector<AsmList *>& al, bool bnd, SurfPos& surf_pos, 
       Hermes::vector<bool>& nat, int isurf, Element** e, Element* trav_base);

  /// Calculates the integration order for assemble_surface_vector_forms().
  int calc_order_vector_form_surf(WeakForm::VectorFormSurf *vfs, Hermes::vector<Solution *> u_ext,
                                  PrecalcShapeset *fv, RefMap *rv, SurfPos* surf_pos);
  /// Assemble surface vector forms.
  void assemble_surface_vector_forms(WeakForm::Stage& stage, 
       SparseMatrix* mat, Vector* rhs, bool rhsonly, bool force_diagonal_blocks, Table* block_weights,
       Hermes::vector<PrecalcShapeset *>& spss, Hermes::vector<RefMap *>& refmap, Hermes::vector<Solution *>& u_ext, 
       Hermes::vector<bool>& isempty, int marker, Hermes::vector<AsmList *>& al, bool bnd, SurfPos& surf_pos, 
       Hermes::vector<bool>& nat, int isurf, Element** e, Element* trav_base);

  /// Assemble DG matrix forms.
  void assemble_DG_matrix_forms(WeakForm::Stage& stage, 
       SparseMatrix* mat, Vector* rhs, bool rhsonly, bool force_diagonal_blocks, Table* block_weights,
       Hermes::vector<PrecalcShapeset *>& spss, Hermes::vector<RefMap *>& refmap, Hermes::vector<Solution *>& u_ext, 
       Hermes::vector<bool>& isempty, int marker, Hermes::vector<AsmList *>& al, bool bnd, SurfPos& surf_pos, 
       Hermes::vector<bool>& nat, int isurf, Element** e, Element* trav_base, Element* rep_element);

  /// Assemble DG vector forms.
  void assemble_DG_vector_forms(WeakForm::Stage& stage, 
       SparseMatrix* mat, Vector* rhs, bool rhsonly, bool force_diagonal_blocks, Table* block_weights,
       Hermes::vector<PrecalcShapeset *>& spss, Hermes::vector<RefMap *>& refmap, Hermes::vector<Solution *>& u_ext, 
       Hermes::vector<bool>& isempty, int marker, Hermes::vector<AsmList *>& al, bool bnd, SurfPos& surf_pos, 
       Hermes::vector<bool>& nat, int isurf, Element** e, Element* trav_base, Element* rep_element);

  /// Get the number of unknowns.
  int get_num_dofs();

  bool is_matrix_free() { return wf->is_matrix_free(); }

  void invalidate_matrix() { have_matrix = false; }

  void set_fvm() {this->is_fvm = true;}

  /// Experimental caching of vector valued (vector) forms.
  struct SurfVectorFormsKey
  {
    WeakForm::vector_form_val_t vfs;
    int element_id, isurf, shape_fn;
#ifdef _MSC_VER
    UINT64 sub_idx;
    SurfVectorFormsKey(WeakForm::vector_form_val_t vfs, int element_id, int isurf, int shape_fn, UINT64 sub_idx)
      : vfs(vfs), element_id(element_id), isurf(isurf), shape_fn(shape_fn), sub_idx(sub_idx) {};
#else
    unsigned int sub_idx;
    SurfVectorFormsKey(WeakForm::vector_form_val_t vfs, int element_id, int isurf, int shape_fn, unsigned int sub_idx)
      : vfs(vfs), element_id(element_id), isurf(isurf), shape_fn(shape_fn), sub_idx(sub_idx) {};
#endif
  };

  struct SurfVectorFormsKeyCompare
  {
    bool operator()(SurfVectorFormsKey a, SurfVectorFormsKey b) const
    {
      if (a.vfs < b.vfs)
        return true;
      else if (a.vfs > b.vfs)
        return false;
      else
          if (a.element_id < b.element_id)
            return true;
          else if (a.element_id > b.element_id)
            return false;
          else
              if (a.isurf < b.isurf)
                return true;
              else if (a.isurf > b.isurf)
                return false;
              else
                  if (a.sub_idx < b.sub_idx)
                    return true;
                  else if (a.sub_idx > b.sub_idx)
                    return false;
                  else
                    return (a.shape_fn < b.shape_fn);
    }
  };

  /// Cache of SurfaceVectorForms values.
  static std::map<SurfVectorFormsKey, double*, SurfVectorFormsKeyCompare> surf_forms_cache;

  /// Static key to surf_forms_cache.
  static SurfVectorFormsKey surf_forms_key;

  struct VolVectorFormsKey
  {
    WeakForm::vector_form_val_t vfv;
    int element_id, shape_fn;
    VolVectorFormsKey(WeakForm::vector_form_val_t vfv, int element_id, int shape_fn)
      : vfv(vfv), element_id(element_id), shape_fn(shape_fn) {};
  };

  struct VolVectorFormsKeyCompare
  {
    bool operator()(VolVectorFormsKey a, VolVectorFormsKey b) const
    {
      if (a.vfv < b.vfv)
        return true;
      else if (a.vfv > b.vfv)
        return false;
      else
          if (a.element_id < b.element_id)
            return true;
          else if (a.element_id > b.element_id)
            return false;
          else
            return (a.shape_fn < b.shape_fn);
    }
  };

  /// Cache of SurfaceVectorForms values.
  static std::map<VolVectorFormsKey, double*, VolVectorFormsKeyCompare> vol_forms_cache;

  /// Static key to vol_forms_cache.
  static VolVectorFormsKey vol_forms_key;

  /// Method to empty the above caches.
  static void empty_form_caches();

  void use_vector_valued_forms() { vector_valued_forms = true; };

protected:
  WeakForm* wf;

  /// If the problem has only constant test functions, there is no need for order calculation,
  /// which saves time.
  bool is_fvm;

  /// Experimental caching of vector valued forms.
  bool vector_valued_forms;

  bool is_linear;

  int ndof;
  int *sp_seq;
  int wf_seq;
  Hermes::vector<Space *> spaces;

  scalar** matrix_buffer;                // buffer for holding square matrix (during assembling)
  int matrix_buffer_dim;                 // dimension of the matrix held by 'matrix_buffer'
  scalar** get_matrix_buffer(int n);

  bool have_spaces;
  bool have_matrix;

  bool values_changed;
  bool struct_changed;
  bool is_up_to_date();

  PrecalcShapeset** pss;    // This is different from H3D.
  int num_user_pss;         // This is different from H3D.

  ExtData<Ord>* init_ext_fns_ord(Hermes::vector<MeshFunction *> &ext);
  ExtData<Ord>* init_ext_fns_ord(Hermes::vector<MeshFunction *> &ext, int edge);
  ExtData<Ord>* init_ext_fns_ord(Hermes::vector<MeshFunction *> &ext, NeighborSearch* nbs);
  ExtData<scalar>* init_ext_fns(Hermes::vector<MeshFunction *> &ext, RefMap *rm, const int order);
  ExtData<scalar>* init_ext_fns(Hermes::vector<MeshFunction *> &ext, NeighborSearch* nbs);
  Func<double>* get_fn(PrecalcShapeset *fu, RefMap *rm, const int order);
  Func<Ord>* get_fn_ord(const int order);

  Geom<double>* cache_e[g_max_quad + 1 + 4 * g_max_quad + 4];
  double* cache_jwt[g_max_quad + 1 + 4 * g_max_quad + 4];

  void init_cache();
  void delete_cache();

  scalar eval_form(WeakForm::MatrixFormVol *mfv, Hermes::vector<Solution *> u_ext,
         PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv);
  scalar eval_form(WeakForm::VectorFormVol *vfv, Hermes::vector<Solution *> u_ext,
         PrecalcShapeset *fv, RefMap *rv);
  scalar eval_form(WeakForm::MatrixFormSurf *mfv, Hermes::vector<Solution *> u_ext,
         PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv, SurfPos* surf_pos);
  scalar eval_form(WeakForm::VectorFormSurf *vfv, Hermes::vector<Solution *> u_ext,
         PrecalcShapeset *fv, RefMap *rv, SurfPos* surf_pos);

  /// Evaluation of forms, discontinuous Galerkin case.
  scalar eval_dg_form(WeakForm::MatrixFormSurf* mfs, Hermes::vector<Solution *> sln,
                      NeighborSearch* nbs_u, NeighborSearch* nbs_v, ExtendedShapeFnPtr efu, ExtendedShapeFnPtr efv,
                      SurfPos* ep);
  scalar eval_dg_form(WeakForm::VectorFormSurf* vfs, Hermes::vector<Solution *> sln,
                      NeighborSearch* nbs_v, PrecalcShapeset* fv, RefMap* rv,
                      SurfPos* ep);
  /// Class handling various caches used in assembling.
  class AssemblingCaches {
  public:
    /// Basic constructor and destructor.
    AssemblingCaches();
    ~AssemblingCaches();

    /// Key for caching precalculated shapeset values on transformed elements with constant
    /// jacobians.
    struct KeyConst
    {
      int index;
      int order;
#ifdef _MSC_VER
      UINT64 sub_idx;
#else
      unsigned int sub_idx;
#endif
      int shapeset_type;
      double inv_ref_map[2][2];
#ifdef _MSC_VER
      KeyConst(int index, int order, UINT64 sub_idx, int shapeset_type, double2x2* inv_ref_map) {
        this->index = index;
        this->order = order;
        this->sub_idx = sub_idx;
        this->shapeset_type = shapeset_type;
        this->inv_ref_map[0][0] = (*inv_ref_map)[0][0];
        this->inv_ref_map[0][1] = (*inv_ref_map)[0][1];
        this->inv_ref_map[1][0] = (*inv_ref_map)[1][0];
        this->inv_ref_map[1][1] = (*inv_ref_map)[1][1];
      }
#else
      KeyConst(int index, int order, unsigned int sub_idx, int shapeset_type, double2x2* inv_ref_map) {
        this->index = index;
        this->order = order;
        this->sub_idx = sub_idx;
        this->shapeset_type = shapeset_type;
        this->inv_ref_map[0][0] = (*inv_ref_map)[0][0];
        this->inv_ref_map[0][1] = (*inv_ref_map)[0][1];
        this->inv_ref_map[1][0] = (*inv_ref_map)[1][0];
        this->inv_ref_map[1][1] = (*inv_ref_map)[1][1];
      }
#endif
    };

    /// Functor that compares two above keys (needed e.g. to create a std::map indexed by these keys);
    struct CompareConst {
      bool operator()(KeyConst a, KeyConst b) const {
        if(a.inv_ref_map[0][0] < b.inv_ref_map[0][0]) return true;
        else if(a.inv_ref_map[0][0] > b.inv_ref_map[0][0]) return false;
        else
          if(a.inv_ref_map[0][1] < b.inv_ref_map[0][1]) return true;
          else if(a.inv_ref_map[0][1] > b.inv_ref_map[0][1]) return false;
          else
            if(a.inv_ref_map[1][0] < b.inv_ref_map[1][0]) return true;
            else if(a.inv_ref_map[1][0] > b.inv_ref_map[1][0]) return false;
            else
              if(a.inv_ref_map[1][1] < b.inv_ref_map[1][1]) return true;
              else if(a.inv_ref_map[1][1] > b.inv_ref_map[1][1]) return false;
              else
                if (a.index < b.index) return true;
                else if (a.index > b.index) return false;
                else
                  if (a.order < b.order) return true;
                  else if (a.order > b.order) return false;
                  else
                    if (a.sub_idx < b.sub_idx) return true;
                    else if (a.sub_idx > b.sub_idx) return false;
                    else
                      if (a.shapeset_type < b.shapeset_type) return true;
                      else return false;
      };
    };

    /// PrecalcShapeset stored values for Elements with constant jacobian of the reference mapping.
    /// For triangles.
    std::map<KeyConst, Func<double>*, CompareConst> const_cache_fn_triangles;
    /// For quads
    std::map<KeyConst, Func<double>*, CompareConst> const_cache_fn_quads;
    
    /// The same setup for elements with non-constant jacobians.
    /// This cache is deleted with every change of the state in assembling.
    struct KeyNonConst {
      int index;
      int order;
#ifdef _MSC_VER
      UINT64 sub_idx;
#else
      unsigned int sub_idx;
#endif
      int shapeset_type;
#ifdef _MSC_VER
      KeyNonConst(int index, int order, UINT64 sub_idx, int shapeset_type) {
        this->index = index;
        this->order = order;
        this->sub_idx = sub_idx;
        this->shapeset_type = shapeset_type;
      }
#else
      KeyNonConst(int index, int order, unsigned int sub_idx, int shapeset_type) {
        this->index = index;
        this->order = order;
        this->sub_idx = sub_idx;
        this->shapeset_type = shapeset_type;
      }
#endif
    };

    /// Functor that compares two above keys (needed e.g. to create a std::map indexed by these keys);
    struct CompareNonConst {
      bool operator()(KeyNonConst a, KeyNonConst b) const {
        if (a.index < b.index) return true;
        else if (a.index > b.index) return false;
        else {
          if (a.order < b.order) return true;
          else if (a.order > b.order) return false;
          else {
            if (a.sub_idx < b.sub_idx) return true;
            else if (a.sub_idx > b.sub_idx) return false;
            else {
              if (a.shapeset_type < b.shapeset_type) return true;
              else return false;
            }
          }
        }
      }
    };
    
    /// PrecalcShapeset stored values for Elements with constant jacobian of the reference mapping.
    /// For triangles.
    std::map<KeyNonConst, Func<double>*, CompareNonConst> cache_fn_triangles;
    /// For quads
    std::map<KeyNonConst, Func<double>*, CompareNonConst> cache_fn_quads;

    LightArray<Func<Ord>*> cache_fn_ord;
  };
  AssemblingCaches assembling_caches;
};

/// Create globally refined space.
HERMES_API Hermes::vector<Space *>* construct_refined_spaces(Hermes::vector<Space *> coarse, int order_increase = 1);
HERMES_API Space* construct_refined_space(Space* coarse, int order_increase = 1);

HERMES_API double get_l2_norm(Vector* vec);

/// New interface, still in developement
/// HERMES_API bool solve_newton(scalar* coeff_vec, DiscreteProblem* dp, Solver* solver, SparseMatrix* matrix,
///		               Vector* rhs, double NEWTON_TOL, int NEWTON_MAX_ITER, bool verbose,
///                             unsigned int stop_condition = NEWTON_WATCH_RESIDUAL);

HERMES_API bool solve_newton(scalar* coeff_vec, DiscreteProblem* dp, Solver* solver, SparseMatrix* matrix,
			     Vector* rhs, double NEWTON_TOL, int NEWTON_MAX_ITER, bool verbose = false,
                             bool residual_as_function = false,
                             double damping_coeff = 1.0, double max_allowed_residual_norm = 1e6);

HERMES_API bool solve_picard(WeakForm* wf, Space* space, Solution* sln_prev_iter,
                             MatrixSolverType matrix_solver, double picard_tol,
			     int picard_max_iter, bool verbose);

#endif
