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

#include "common.h"
#include "linear_problem.h"
#include "matrix.h"
#include "traverse.h"
#include <common/error.h>
#include <common/callstack.h>
#include <common/timer.h>


LinearProblem::LinearProblem(WeakForm *wf)
{
	_F_
	this->wf = wf;

	spaces = new Space *[wf->neq];
	sp_seq = new int[wf->neq];
	memset(sp_seq, -1, sizeof(int) * wf->neq);

	matrix_buffer = NULL;
	matrix_buffer_dim = 0;
}

LinearProblem::~LinearProblem()
{
	_F_
	delete [] spaces;
	delete [] sp_seq;
}

void LinearProblem::free()
{
	memset(sp_seq, -1, sizeof(int) * wf->neq);
}

void LinearProblem::set_spaces(Tuple<Space *> sp)
{
	_F_
	int n = sp.size();
	if (n <= 0 || n > wf->neq) EXIT("Bad number of spaces.");
        if (n != this->wf->neq) 
          error("Number of spaces must match the number of equations in LinearProblem::set_spaces()");
	for (int i = 0; i < n; i++) spaces[i] = sp[i];
	memset(sp_seq, -1, sizeof(int) * wf->neq);
}

scalar **LinearProblem::get_matrix_buffer(int n)
{
	_F_
	if (n <= matrix_buffer_dim) return matrix_buffer;
	if (matrix_buffer != NULL) delete [] matrix_buffer;
	matrix_buffer_dim = n;
	return (matrix_buffer = new_matrix<scalar>(n, n));
}

bool LinearProblem::assemble(Matrix *matrix, Vector *rhs)
{
	_F_
	Timer tmr;
	tmr.start();

	create(matrix, rhs);
	if (ndofs == 0) return false;

	bool bnd[10];						// FIXME: magic number - maximal possible number of faces
	FacePos fp[10];
	bool nat[wf->neq], isempty[wf->neq];

	AsmList al[wf->neq];
	AsmList *am, *an;
	ShapeFunction base_fn[wf->neq];
	ShapeFunction test_fn[wf->neq];
	ShapeFunction *fu, *fv;
	RefMap refmap[wf->neq];
	for (int i = 0; i < wf->neq; i++) {
		base_fn[i].set_shapeset(spaces[i]->get_shapeset());
		test_fn[i].set_shapeset(spaces[i]->get_shapeset());
		refmap[i].set_mesh(spaces[i]->get_mesh());
	}

	matrix_buffer = NULL;
	get_matrix_buffer(10);

	// obtain a list of assembling stages
	std::vector<WeakForm::Stage> stages;
	wf->get_stages(spaces, stages, matrix == NULL);

	// Loop through all assembling stages -- the purpose of this is increased performance
	// in multi-mesh calculations, where, e.g., only the right hand side uses two meshes.
	// In such a case, the bilinear forms are assembled over one mesh, and only the rhs
	// traverses through the union mesh. On the other hand, if you don't use multi-mesh
	// at all, there will always be only one stage in which all forms are assembled as usual.
	Traverse trav;
	for (unsigned ss = 0; ss < stages.size(); ss++) {
		WeakForm::Stage *s = &stages[ss];
		for (unsigned i = 0; i < s->idx.size(); i++)
			s->fns[i] = &base_fn[s->idx[i]];
		trav.begin(s->meshes.size(), &(s->meshes.front()), &(s->fns.front()));

		// assemble one stage
		Element **e;
		while ((e = trav.get_next_state(bnd, fp)) != NULL) {
			// find a non-NULL e[i]
			Element *e0;
			for (unsigned i = 0; i < s->idx.size(); i++)
				if ((e0 = e[i]) != NULL) break;
			if (e0 == NULL) continue;

			// obtain assembly lists for the element at all spaces
			memset(isempty, 0, sizeof(bool) * wf->neq);
			for (unsigned i = 0; i < s->idx.size(); i++) {
				int j = s->idx[i];
				if (e[i] == NULL) { isempty[j] = true; continue; }

				// TODO: do not obtain again if the element was not changed
				spaces[j]->get_element_assembly_list(e[i], al + j);
				test_fn[j].set_active_element(e[i]);
				test_fn[j].set_transform(base_fn + j);

				refmap[j].set_active_element(e[i]);
				refmap[j].force_transform(base_fn[j].get_transform(), base_fn[j].get_ctm());
			}
			int marker = e0->marker;

			fn_cache.free();
			// assemble volume bilinear forms //////////////////////////////////////
			for (unsigned ww = 0; ww < s->mfvol.size(); ww++) {
				WeakForm::MatrixFormVol *mfv = s->mfvol[ww];
				if (isempty[mfv->i] || isempty[mfv->j]) continue;
				if (mfv->area != ANY && !wf->is_in_area(marker, mfv->area)) continue;
				int m = mfv->i; fv = test_fn + m; am = al + m;
				int n = mfv->j; fu = base_fn + n; an = al + n;
				bool tra = (m != n) && (mfv->sym != UNSYM);
				bool sym = (m == n) && (mfv->sym == SYM);

				// assemble the local stiffness matrix for the form mfv
				scalar **mat = get_matrix_buffer(std::max(am->cnt, an->cnt));

				for (int i = 0; i < am->cnt; i++) {
					int k = am->dof[i];
					if (!tra && k == H3D_DIRICHLET_DOF) continue;
					fv->set_active_shape(am->idx[i]);

					if (!sym) { // unsymmetric block
						for (int j = 0; j < an->cnt; j++) {
							fu->set_active_shape(an->idx[j]);
							scalar bi = eval_form(mfv, NULL, fu, fv, refmap + n, refmap + m)
								* an->coef[j] * am->coef[i];
							if (an->dof[j] == H3D_DIRICHLET_DOF) rhs->add(k, -bi);
							else mat[i][j] = bi;
						}
					}
					else { // symmetric block
						for (int j = 0; j < an->cnt; j++) {
							if (j < i && an->dof[j] >= 0) continue;
							fu->set_active_shape(an->idx[j]);
							scalar bi = eval_form(mfv, NULL, fu, fv, refmap + n, refmap + m)
								* an->coef[j] * am->coef[i];
							if (an->dof[j] == H3D_DIRICHLET_DOF) rhs->add(k, -bi);
							else mat[i][j] = mat[j][i] = bi;
						}
					}
				}

				// insert the local stiffness matrix into the global one
				matrix->add(am->cnt, an->cnt, mat, am->dof, an->dof);

				// insert also the off-diagonal (anti-)symmetric block, if required
				if (tra) {
					if (mfv->sym < 0) chsgn(mat, am->cnt, an->cnt);
					transpose(mat, am->cnt, an->cnt);
					matrix->add(an->cnt, am->cnt, mat, an->dof, am->dof);

					// we also need to take care of the RHS...
					for (int j = 0; j < am->cnt; j++)
						if (am->dof[j] == H3D_DIRICHLET_DOF)
							for (int i = 0; i < an->cnt; i++)
								if (an->dof[i] != H3D_DIRICHLET_DOF)
									rhs->add(an->dof[i], -mat[i][j]);
				}
			}

			// assemble volume linear forms ////////////////////////////////////////
			for (unsigned ww = 0; ww < s->vfvol.size(); ww++) {
				WeakForm::VectorFormVol *vfv = s->vfvol[ww];
				if (isempty[vfv->i]) continue;
				if (vfv->area != ANY && !wf->is_in_area(marker, vfv->area)) continue;
				int m = vfv->i;  fv = test_fn + m;  am = al + m;

				for (int i = 0; i < am->cnt; i++) {
					if (am->dof[i] == H3D_DIRICHLET_DOF) continue;
					fv->set_active_shape(am->idx[i]);
					rhs->add(am->dof[i], eval_form(vfv, NULL, fv, refmap + m) * am->coef[i]);
				}
			}

			// assemble surface integrals now: loop through boundary faces of the element
			for (int iface = 0; iface < e[0]->get_num_faces(); iface++) {
				fn_cache.free();
				if (!bnd[iface]/* || !e0->en[edge]->bnd*/) continue;
				int marker = fp[iface].marker;

				// obtain the list of shape functions which are nonzero on this edge
				for (unsigned i = 0; i < s->idx.size(); i++) {
					if (e[i] == NULL) continue;
					int j = s->idx[i];
					if ((nat[j] = (spaces[j]->bc_type_callback(marker) == BC_NATURAL)))
						spaces[j]->get_boundary_assembly_list(e[i], iface, al + j);
				}

				// assemble surface bilinear forms ///////////////////////////////////
				for (unsigned ww = 0; ww < s->mfsurf.size(); ww++) {
					WeakForm::MatrixFormSurf *mfs = s->mfsurf[ww];
					if (isempty[mfs->i] || isempty[mfs->j]) continue;
					if (mfs->area != ANY && !wf->is_in_area(marker, mfs->area)) continue;
					int m = mfs->i; fv = test_fn + m; am = al + m;
					int n = mfs->j; fu = base_fn + n; an = al + n;

					if (!nat[m] || !nat[n]) continue;
					fp[iface].base = trav.get_base();
					fp[iface].space_v = spaces[m];
					fp[iface].space_u = spaces[n];

					scalar **mat = get_matrix_buffer(std::max(am->cnt, an->cnt));
					for (int i = 0; i < am->cnt; i++) {
						int k = am->dof[i];
						if (k == H3D_DIRICHLET_DOF) continue;
						fv->set_active_shape(am->idx[i]);
						for (int j = 0; j < an->cnt; j++) {
							fu->set_active_shape(an->idx[j]);
							scalar bi = eval_form(mfs, NULL, fu, fv, refmap + n, refmap + m, fp + iface)
								* an->coef[j] * am->coef[i];
							if (an->dof[j] != H3D_DIRICHLET_DOF) mat[i][j] = bi;
							else rhs->add(k, -bi);
						}
					}
					matrix->add(am->cnt, an->cnt, mat, am->dof, an->dof);
				}

				// assemble surface linear forms /////////////////////////////////////
				for (unsigned ww = 0; ww < s->vfsurf.size(); ww++) {
					WeakForm::VectorFormSurf *vfs = s->vfsurf[ww];
					if (isempty[vfs->i]) continue;
					if (vfs->area != ANY && !wf->is_in_area(marker, vfs->area)) continue;
					int m = vfs->i; fv = test_fn + m; am = al + m;

					if (!nat[m]) continue;
					fp[iface].base = trav.get_base();
					fp[iface].space_v = spaces[m];

					for (int i = 0; i < am->cnt; i++) {
						if (am->dof[i] == H3D_DIRICHLET_DOF) continue;
						fv->set_active_shape(am->idx[i]);
						rhs->add(am->dof[i], eval_form(vfs, NULL, fv, refmap + m, fp + iface)
						         * am->coef[i]);
					}
				}
			}
		}
		trav.finish();
	}

	SparseMatrix *mat = dynamic_cast<SparseMatrix *>(matrix);
	if (mat != NULL) mat->finish();
	if (rhs != NULL) rhs->finish();

	delete [] matrix_buffer;
	matrix_buffer = NULL;
	matrix_buffer_dim = 0;

	tmr.stop();
	time = tmr.get_seconds();

	return true;
}

void LinearProblem::create(Matrix *matrix, Vector *rhs)
{
	_F_

	// check if we can reuse the matrix structure
	bool up_to_date = true;
	for (int i = 0; i < wf->neq; i++)
		if (spaces[i]->get_seq() != sp_seq[i]) {
			up_to_date = false;
			break;
		}

	if (up_to_date) {
		// Reusing matrix structure
		if (matrix != NULL) matrix->zero();
		if (rhs != NULL) rhs->zero();
		return;
	}

	ndofs = 0;
	for (int i = 0; i < wf->neq; i++)
		ndofs += spaces[i]->get_dof_count();
	if (!ndofs) return;


	if (matrix != NULL) {
		SparseMatrix *mat = dynamic_cast<SparseMatrix *>(matrix);
		if (mat == NULL) {
			warning("Calling LinearProblem::create() with a matrix that is not sparse matrix.");
			return;
		}

		mat->free();
		mat->prealloc(ndofs);

		AsmList al[wf->neq];
		Mesh *meshes[wf->neq];
		bool **blocks = wf->get_blocks();

		// init multi-mesh traversal
		for (int i = 0; i < wf->neq; i++)
			meshes[i] = spaces[i]->get_mesh();

		Traverse trav;
		trav.begin(wf->neq, meshes);

		// loop through all elements
		Element **e;
		while ((e = trav.get_next_state(NULL, NULL)) != NULL) {
			// obtain assembly lists for the element at all spaces
			for (int i = 0; i < wf->neq; i++)
				// TODO: do not get the assembly list again if the element was not changed
				if (e[i] != NULL)
					spaces[i]->get_element_assembly_list(e[i], al + i);

			// go through all equation-blocks of the local stiffness matrix
			for (int m = 0; m < wf->neq; m++)
				for (int n = 0; n < wf->neq; n++)
					if (blocks[m][n] && e[m] != NULL && e[n] != NULL) {
						AsmList *am = al + m;
						AsmList *an = al + n;

						// pretend assembling of the element stiffness matrix
						// register nonzero elements
						for (int i = 0; i < am->cnt; i++)
							if (am->dof[i] != H3D_DIRICHLET_DOF)
								for (int j = 0; j < an->cnt; j++)
									if (an->dof[j] != H3D_DIRICHLET_DOF)
										mat->pre_add_ij(am->dof[i], an->dof[j]);
					}
		}

		trav.finish();
		delete [] blocks;

		mat->alloc();
	}

	if (rhs != NULL) rhs->alloc(ndofs);

	// save space seq numbers, so we can detect their changes
	for (int i = 0; i < wf->neq; i++)
		sp_seq[i] = spaces[i]->get_seq();
}

void LinearProblem::init_ext_fns(user_data_t<scalar> &ud, std::vector<MeshFunction *> &ext, int order,
                              RefMap *rm, const int np, const QuadPt3D *pt)
{
	_F_
	ud.nf = ext.size();
	mfn_t *ext_fn = new mfn_t[ud.nf];
	for (int i = 0; i < ud.nf; i++) {
		fn_key_t key(ext[i]->seq, order, ext[i]->get_transform());
		mfn_t *efn = NULL;
		if (!fn_cache.ext.lookup(key, efn)) {
			efn = init_fn(ext[i], rm, np, pt);
			fn_cache.ext.set(key, efn);
		}
		assert(efn != NULL);
		ext_fn[i] = *efn;
	}
	ud.ext = ext_fn;
}

void LinearProblem::init_ext_fns(user_data_t<ord_t> &fake_ud, std::vector<MeshFunction *> &ext)
{
	_F_
	fake_ud.nf = ext.size();
	fn_t<ord_t> *fake_ext_fn = new fn_t<ord_t>[fake_ud.nf];
	for (int i = 0; i < fake_ud.nf; i++) {
		fake_ext_fn[i] = init_fn(ext[i]->get_fn_order());
	}
	fake_ud.ext = fake_ext_fn;
}

sfn_t *LinearProblem::get_fn(ShapeFunction *fu, int order, RefMap *rm, const int np,
                          const QuadPt3D *pt)
{
	_F_
	fn_key_t key(fu->get_active_shape(), order, fu->get_transform(), fu->get_shapeset()->id);
	sfn_t *u = NULL;
	if (!fn_cache.fn.lookup(key, u)) {
		u = init_fn(fu, rm, np, pt);
		fn_cache.fn.set(key, u);
	}
	return u;
}

sfn_t *LinearProblem::get_fn(ShapeFunction *fu, int order, RefMap *rm, int iface, const int np,
                          const QuadPt3D *pt)
{
	fn_key_t key(fu->get_active_shape(), order, fu->get_transform(), fu->get_shapeset()->id);
	sfn_t *u = NULL;
	if (!fn_cache.fn.lookup(key, u)) {
		u = init_fn(fu, rm, iface, np, pt);
		fn_cache.fn.set(key, u);
	}
	return u;
}

scalar LinearProblem::eval_form(WeakForm::MatrixFormVol *mfv, fn_t<scalar> *u_ext[], ShapeFunction *fu, ShapeFunction *fv,
                             RefMap *ru, RefMap *rv)
{
	_F_
	Element *elem = fu->get_active_element();

	// determine the integration order
	fn_t<ord_t> ou = init_fn(fu->get_fn_order());
	fn_t<ord_t> ov = init_fn(fv->get_fn_order());

	user_data_t<ord_t> fake_ud;
	init_ext_fns(fake_ud, mfv->ext);

	double fake_wt = 1.0;
	geom_t<ord_t> fake_e = init_geom(elem->marker);
	ord_t o = mfv->ord(1, &fake_wt, NULL, &ou, &ov, &fake_e, &fake_ud);
	order3_t order = ru->get_inv_ref_order();
	switch (order.type) {
		case MODE_TETRAHEDRON: order += order3_t(o.get_order()); break;
		case MODE_HEXAHEDRON: order += order3_t(o.get_order(), o.get_order(), o.get_order()); break;
	}
	order.limit();
	int ord_idx = order.get_idx();

	free_fn(&ou);
	free_fn(&ov);

	// eval the form
	Quad3D *quad = get_quadrature(elem->get_mode());
	int np = quad->get_num_points(order);
	QuadPt3D *pt = quad->get_points(order);

	double *jwt = NULL;
	geom_t<double> e;
	if (!fn_cache.e.exists(ord_idx)) {
		fn_cache.jwt[ord_idx] = ru->get_jacobian(np, pt);
		fn_cache.e[ord_idx] = init_geom(elem->marker, ru, np, pt);
	}
	jwt = fn_cache.jwt[ord_idx];
	e = fn_cache.e[ord_idx];

	sfn_t *u = get_fn(fu, ord_idx, ru, np, pt);
	sfn_t *v = get_fn(fv, ord_idx, rv, np, pt);

	user_data_t<scalar> ud;
	init_ext_fns(ud, mfv->ext, ord_idx, rv, np, pt);

	return mfv->fn(np, jwt, u_ext, u, v, &e, &ud);
}

scalar LinearProblem::eval_form(WeakForm::VectorFormVol *vfv, fn_t<scalar> *u_ext[], ShapeFunction *fv, RefMap *rv)
{
	_F_
	Element *elem = fv->get_active_element();

	// determine the integration order
	fn_t<ord_t> ov = init_fn(fv->get_fn_order());

	user_data_t<ord_t> fake_ud;
	init_ext_fns(fake_ud, vfv->ext);

	double fake_wt = 1.0;
	geom_t<ord_t> fake_e = init_geom(elem->marker);
	ord_t o = vfv->ord(1, &fake_wt, NULL, &ov, &fake_e, &fake_ud);
	order3_t order = rv->get_inv_ref_order();
	switch (order.type) {
		case MODE_TETRAHEDRON: order += order3_t(o.get_order()); break;
		case MODE_HEXAHEDRON: order += order3_t(o.get_order(), o.get_order(), o.get_order()); break;
	}
	order.limit();
	int ord_idx = order.get_idx();

	free_fn(&ov);

	// eval the form
	Quad3D *quad = get_quadrature(elem->get_mode());
	int np = quad->get_num_points(order);
	QuadPt3D *pt = quad->get_points(order);

	double *jwt = NULL;
	geom_t<double> e;
	if (!fn_cache.e.exists(ord_idx)) {
		fn_cache.jwt[ord_idx] = rv->get_jacobian(np, pt);
		fn_cache.e[ord_idx] = init_geom(elem->marker, rv, np, pt);
	}
	jwt = fn_cache.jwt[ord_idx];
	e = fn_cache.e[ord_idx];

	sfn_t *v = get_fn(fv, ord_idx, rv, np, pt);

	user_data_t<scalar> ud;
	init_ext_fns(ud, vfv->ext, ord_idx, rv, np, pt);

	return vfv->fn(np, jwt, u_ext, v, &e, &ud);
}

scalar LinearProblem::eval_form(WeakForm::MatrixFormSurf *mfs, fn_t<scalar> *u_ext[], ShapeFunction *fu, ShapeFunction *fv,
                             RefMap *ru, RefMap *rv, FacePos *fp)
{
	_F_

	// determine the integration order
	fn_t<ord_t> ou = init_fn(fu->get_fn_order());
	fn_t<ord_t> ov = init_fn(fv->get_fn_order());

	user_data_t<ord_t> fake_ud;
	init_ext_fns(fake_ud, mfs->ext);

	double fake_wt = 1.0;
	geom_t<ord_t> fake_e = init_geom(fp->marker);
	ord_t o = mfs->ord(1, &fake_wt, NULL, &ou, &ov, &fake_e, &fake_ud);
	order3_t order = ru->get_inv_ref_order();
	switch (order.type) {
		case MODE_TETRAHEDRON: order += order3_t(o.get_order()); break;
		case MODE_HEXAHEDRON: order += order3_t(o.get_order(), o.get_order(), o.get_order()); break;
	}
	order.limit();
	order2_t face_order = order.get_face_order(fp->face);
	int ord_idx = face_order.get_idx();

	free_fn(&ou);
	free_fn(&ov);

	// eval the form
	Quad3D *quad = get_quadrature(fu->get_active_element()->get_mode());
	int np = quad->get_face_num_points(fp->face, face_order);
	QuadPt3D *pt = quad->get_face_points(fp->face, face_order);

	double *jwt = NULL;
	geom_t<double> e;
	if (!fn_cache.e.exists(ord_idx)) {
		fn_cache.jwt[ord_idx] = ru->get_face_jacobian(fp->face, np, pt);
		fn_cache.e[ord_idx] = init_geom(fp->marker, ru, fp->face, np, pt);
	}
	jwt = fn_cache.jwt[ord_idx];
	e = fn_cache.e[ord_idx];

	sfn_t *u = get_fn(fu, ord_idx, ru, fp->face, np, pt);
	sfn_t *v = get_fn(fv, ord_idx, rv, fp->face, np, pt);

	user_data_t<scalar> ud;
	init_ext_fns(ud, mfs->ext, ord_idx, rv, np, pt);

	return mfs->fn(np, jwt, u_ext, u, v, &e, &ud);
}

scalar LinearProblem::eval_form(WeakForm::VectorFormSurf *vfs, fn_t<scalar> *u_ext[], ShapeFunction *fv, RefMap *rv, FacePos *fp)
{
	_F_

	// determine the integration order
	user_data_t<ord_t> fake_ud;
	init_ext_fns(fake_ud, vfs->ext);

	fn_t<ord_t> ov = init_fn(fv->get_fn_order());
	double fake_wt = 1.0;
	geom_t<ord_t> fake_e = init_geom(fp->marker);
	ord_t o = vfs->ord(1, &fake_wt, NULL, &ov, &fake_e, &fake_ud);
	order3_t order = rv->get_inv_ref_order();
	switch (order.type) {
		case MODE_TETRAHEDRON: order += order3_t(o.get_order()); break;
		case MODE_HEXAHEDRON: order += order3_t(o.get_order(), o.get_order(), o.get_order()); break;
	}
	order.limit();
	order2_t face_order = order.get_face_order(fp->face);
	int ord_idx = face_order.get_idx();

	free_fn(&ov);

	// eval the form
	Quad3D *quad = get_quadrature(fv->get_active_element()->get_mode());
	int np = quad->get_face_num_points(fp->face, face_order);
	QuadPt3D *pt = quad->get_face_points(fp->face, face_order);

	double *jwt = NULL;
	geom_t<double> e;
	if (!fn_cache.e.exists(ord_idx)) {
		fn_cache.jwt[ord_idx] = rv->get_face_jacobian(fp->face, np, pt);
		fn_cache.e[ord_idx] = init_geom(fp->marker, rv, fp->face, np, pt);
	}
	jwt = fn_cache.jwt[ord_idx];
	e = fn_cache.e[ord_idx];

	sfn_t *v = get_fn(fv, ord_idx, rv, fp->face, np, pt);

	user_data_t<scalar> ud;
	init_ext_fns(ud, vfs->ext, ord_idx, rv, np, pt);

	return vfs->fn(np, jwt, u_ext, v, &e, &ud);
}
