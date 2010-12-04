// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes3D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes3D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes3D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#ifndef _DISCRETE_PROBLEM_H_
#define _DISCRETE_PROBLEM_H_

#include "h3d_common.h"
#include "weakform.h"
#include "tuple.h"
#include "../../hermes_common/array.h"
#include "../../hermes_common/solver/solver.h"

class Space;
class Matrix;
class SparseMatrix;
class Vector;
class SurfPos;

/// Discrete problem class
///
/// This class does assembling into external matrix / vactor structures.
///
class HERMES_API DiscreteProblem {
public:
        DiscreteProblem(WeakForm *wf, Hermes::Tuple<Space *> sp, bool is_linear = false);
	virtual ~DiscreteProblem();
	void free();

        // Get pointer to n-th space.
        Space* get_space(int n) {  return this->spaces[n];  }

        // Precalculate matrix sparse structure.
	void create(SparseMatrix* mat, Vector* rhs = NULL, bool rhsonly = false);

        // General assembling procedure for nonlinear problems. coeff_vec is the 
        // previous Newton vector.
	void assemble(SparseMatrix* mat, Vector* rhs = NULL, bool rhsonly = false);

        // Assembling for linear problems. Same as the previous functions, but 
        // does not need the coeff_vector.
	void assemble(scalar* coeff_vec, SparseMatrix* mat, Vector* rhs = NULL,
                      bool rhsonly = false);

        // Get the number of unknowns.
	int get_num_dofs();

	bool is_matrix_free() { return wf->is_matrix_free(); }
  
        void invalidate_matrix() { have_matrix = false; }

protected:
	WeakForm* wf;

        bool is_linear;

	int ndof;				/// number of DOF
	int* sp_seq;				/// sequence numbers of spaces
        int wf_seq;
        Hermes::Tuple<Space *> spaces;

	scalar** matrix_buffer;		/// buffer for holding square matrix (during assembling)
	int matrix_buffer_dim;		/// dimension of the matrix held by 'matrix_buffer'
	inline scalar** get_matrix_buffer(int n);

	bool have_spaces;
	bool have_matrix;

        bool values_changed;
        bool struct_changed;
	bool is_up_to_date();

	// pre-transforming and fn. caching
	struct fn_key_t {
		int index;
		int order;
		int sub_idx;
		int ss_id;			// shapeset id

		fn_key_t(int index, int order, int sub_idx, int ss_id = -1) {
			this->index = index;
			this->order = order;
			this->sub_idx = sub_idx;
			this->ss_id = ss_id;
		}
	};

	struct FnCache {
		Array<double *> jwt;			// jacobian x weight
		Array<Geom<double> > e;		// geometries
		Map<fn_key_t, sFunc*> fn;		// shape functions
		Map<fn_key_t, mFunc*> ext;		// external functions
		Map<fn_key_t, mFunc*> sln;		// sln from prev iter

		~FnCache();
		void free();
	} fn_cache;

	scalar eval_form(WeakForm::MatrixFormVol *mfv, Hermes::Tuple<Solution *> u_ext, ShapeFunction *fu,
	                 ShapeFunction *fv, RefMap *ru, RefMap *rv);
	scalar eval_form(WeakForm::VectorFormVol *vfv, Hermes::Tuple<Solution *> u_ext, ShapeFunction *fv, RefMap *rv);
	scalar eval_form(WeakForm::MatrixFormSurf *mfs, Hermes::Tuple<Solution *> u_ext, ShapeFunction *fu,
	                 ShapeFunction *fv, RefMap *ru, RefMap *rv, SurfPos *surf_pos);
	scalar eval_form(WeakForm::VectorFormSurf *vfs, Hermes::Tuple<Solution *> u_ext, ShapeFunction *fv, RefMap *rv,
	                 SurfPos *surf_pos);

	sFunc *get_fn(ShapeFunction *fu, int order, RefMap *rm, const int np, const QuadPt3D *pt);
	sFunc *get_fn(ShapeFunction *fu, int order, RefMap *rm, int iface, const int np,
	              const QuadPt3D *pt);
	mFunc *get_fn(Solution *fu, int order, RefMap *rm, const int np, const QuadPt3D *pt);

	void init_ext_fns(ExtData<Ord> &fake_ud, std::vector<MeshFunction *> &ext);
	void init_ext_fns(ExtData<scalar> &ud, std::vector<MeshFunction *> &ext, int order,
	                  RefMap *rm, const int np, const QuadPt3D *pt);
};

HERMES_API Hermes::Tuple<Space *> * construct_refined_spaces(Hermes::Tuple<Space *> coarse, int order_increase, int refinement);
HERMES_API Space* construct_refined_space(Space* coarse, int order_increase, int refinement);

#endif /* _DISCRETE_PROBLEM_H_ */
