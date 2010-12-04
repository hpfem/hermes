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
#include "../../hermes_common/solver/solver.h"
#include "adapt/adapt.h"
#include "graph.h"
#include "forms.h"
#include "weakform.h"
#include "views/view.h"
#include "views/scalar_view.h"
#include "views/vector_view.h"
#include "views/order_view.h"
#include "function.h"
#include "neighbor.h"
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
HERMES_API_USED_TEMPLATE(Hermes::Tuple<ProjNormType>);
HERMES_API_USED_TEMPLATE(Hermes::Tuple<Space*>);
HERMES_API_USED_TEMPLATE(Hermes::Tuple<MeshFunction*>);
HERMES_API_USED_TEMPLATE(Hermes::Tuple<Solution*>);
HERMES_API_USED_TEMPLATE(Hermes::Tuple<PrecalcShapeset*>);


/// Discrete problem class
///
/// This class does assembling into external matrix / vactor structures.
///
class HERMES_API DiscreteProblem 
{
public:
  DiscreteProblem(WeakForm* wf, Hermes::Tuple<Space *> spaces, bool is_linear = false);
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
  void assemble(scalar* coeff_vec, SparseMatrix* mat, Vector* rhs = NULL, bool rhsonly = false);

  // Assembling for linear problems. Same as the previous functions, but 
  // does not need the coeff_vector.
  void assemble(SparseMatrix* mat, Vector* rhs = NULL, bool rhsonly = false);

  // Get the number of unknowns.
  int get_num_dofs();

  bool is_matrix_free() { return wf->is_matrix_free(); }

  void invalidate_matrix() { have_matrix = false; }

  void set_fvm() {this->is_fvm = true;}  

  // Experimental caching of vector valued (vector) forms.
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

  // If the problem has only constant test functions, there is no need for order calculation, 
  // which saves time.
  bool is_fvm;

  // Experimental caching of vector valued forms.
  bool vector_valued_forms;

  bool is_linear;

  int ndof;
  int *sp_seq;
  int wf_seq;
  Hermes::Tuple<Space *> spaces;

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
  ExtData<Ord>* init_ext_fns_ord(std::vector<MeshFunction *> &ext, NeighborSearch* nbs);
  ExtData<scalar>* init_ext_fns(std::vector<MeshFunction *> &ext, RefMap *rm, const int order);
  ExtData<scalar>* init_ext_fns(std::vector<MeshFunction *> &ext, NeighborSearch* nbs);
  Func<double>* get_fn(PrecalcShapeset *fu, RefMap *rm, const int order);

  // Caching transformed values for element
  std::map<PrecalcShapeset::Key, Func<double>*, PrecalcShapeset::Compare> cache_fn;
  Geom<double>* cache_e[g_max_quad + 1 + 4 * g_max_quad + 4];
  double* cache_jwt[g_max_quad + 1 + 4 * g_max_quad + 4];

  void init_cache();
  void delete_cache();

  scalar eval_form(WeakForm::MatrixFormVol *mfv, Hermes::Tuple<Solution *> u_ext, 
         PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv);
  scalar eval_form(WeakForm::VectorFormVol *vfv, Hermes::Tuple<Solution *> u_ext, 
         PrecalcShapeset *fv, RefMap *rv);
  scalar eval_form(WeakForm::MatrixFormSurf *mfv, Hermes::Tuple<Solution *> u_ext, 
         PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv, SurfPos* surf_pos);
  scalar eval_form(WeakForm::VectorFormSurf *vfv, Hermes::Tuple<Solution *> u_ext, 
         PrecalcShapeset *fv, RefMap *rv, SurfPos* surf_pos);

  // Evaluation of forms, discontinuous Galerkin case.
  scalar eval_dg_form(WeakForm::MatrixFormSurf* mfs, Hermes::Tuple<Solution *> sln, 
                      NeighborSearch* nbs_u, NeighborSearch* nbs_v, ExtendedShapeFnPtr efu, ExtendedShapeFnPtr efv,
                      SurfPos* ep);
  scalar eval_dg_form(WeakForm::VectorFormSurf* vfs, Hermes::Tuple<Solution *> sln,
                      NeighborSearch* nbs_v, PrecalcShapeset* fv, RefMap* rv,
                      SurfPos* ep);
};

// Create globally refined space.
HERMES_API Hermes::Tuple<Space *>* construct_refined_spaces(Hermes::Tuple<Space *> coarse, int order_increase = 1);
HERMES_API Space* construct_refined_space(Space* coarse, int order_increase = 1);

HERMES_API double get_l2_norm(Vector* vec); 

HERMES_API bool solve_newton(scalar* coeff_vec, DiscreteProblem* dp, Solver* solver, SparseMatrix* matrix,
			     Vector* rhs, double NEWTON_TOL, int NEWTON_MAX_ITER, bool verbose);

#endif
