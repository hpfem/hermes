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

#include <common/array.h>
#include "weakform.h"
#include "tuple.h"

class Space;
class Matrix;
class SparseMatrix;
class Vector;
class FacePos;

/// Finite Element problem class
///
/// This class does assembling into passed-in structures.
///
class DiscreteProblem {
public:
	DiscreteProblem(WeakForm *wf);
	virtual ~DiscreteProblem();
	void free();

	void set_spaces(Tuple<Space *> sp);
	void set_space(Space* sp);

	void create(SparseMatrix *mat);
	void assemble(Vector* rhs, Matrix* jac, Vector* x = NULL);
	void assemble(Vector* rhs, Matrix* jac, Tuple<Solution*> u_ext =  Tuple<Solution*> ());

	int get_num_dofs();
	bool is_matrix_free() { return wf->is_matrix_free(); }

protected:
	WeakForm *wf;

	int ndof;					/// number of DOFs
	int *sp_seq;				/// sequence numbers of spaces
	Space **spaces;
	Solution **slns;
	bool have_spaces;

	scalar **matrix_buffer;		/// buffer for holding square matrix (during assembling)
	int matrix_buffer_dim;		/// dimension of the matrix held by 'matrix_buffer'
	inline scalar **get_matrix_buffer(int n);
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
		Array<geom_t<double> > e;		// geometries
		Map<fn_key_t, sfn_t*> fn;		// shape functions
		Map<fn_key_t, mfn_t*> ext;		// external functions
		Map<fn_key_t, mfn_t*> sln;		// sln from prev iter

		~FnCache();
		void free();
	} fn_cache;

	scalar eval_form(WeakForm::MatrixFormVol *mfv, Tuple<Solution *> u_ext, ShapeFunction *fu,
	                 ShapeFunction *fv, RefMap *ru, RefMap *rv);
	scalar eval_form(WeakForm::VectorFormVol *vfv, Tuple<Solution *> u_ext, ShapeFunction *fv, RefMap *rv);
	scalar eval_form(WeakForm::MatrixFormSurf *mfs, Tuple<Solution *> u_ext, ShapeFunction *fu,
	                 ShapeFunction *fv, RefMap *ru, RefMap *rv, FacePos *fp);
	scalar eval_form(WeakForm::VectorFormSurf *vfs, Tuple<Solution *> u_ext, ShapeFunction *fv, RefMap *rv,
	                 FacePos *fp);

	sfn_t *get_fn(ShapeFunction *fu, int order, RefMap *rm, const int np, const QuadPt3D *pt);
	sfn_t *get_fn(ShapeFunction *fu, int order, RefMap *rm, int iface, const int np,
	              const QuadPt3D *pt);
	mfn_t *get_fn(Solution *fu, int order, RefMap *rm, const int np, const QuadPt3D *pt);

	void init_ext_fns(user_data_t<ord_t> &fake_ud, std::vector<MeshFunction *> &ext);
	void init_ext_fns(user_data_t<scalar> &ud, std::vector<MeshFunction *> &ext, int order,
	                  RefMap *rm, const int np, const QuadPt3D *pt);
};

#endif /* _LINPROBLEM_H_ */
