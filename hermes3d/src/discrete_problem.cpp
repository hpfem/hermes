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
#include "discrete_problem.h"
#include "matrix.h"
#include "traverse.h"
#include <common/error.h>
#include <common/callstack.h>


DiscreteProblem::FnCache::~FnCache()
{
	free();
}

void DiscreteProblem::FnCache::free()
{
	_F_
	for (Word_t i = jwt.first(); i != INVALID_IDX; i = jwt.next(i))
		delete [] jwt[i];
	jwt.remove_all();
	for (Word_t i = e.first(); i != INVALID_IDX; i = e.next(i))
		free_geom(&e[i]);
	e.remove_all();
	for (Word_t i = fn.first(); i != INVALID_IDX; i = fn.next(i))
		free_fn(fn[i]);
	fn.remove_all();
	for (Word_t i = ext.first(); i != INVALID_IDX; i = ext.next(i))
		delete ext[i];
	ext.remove_all();
	for (Word_t i = sln.first(); i != INVALID_IDX; i = sln.next(i))
		free_fn(sln[i]);
	sln.remove_all();
}

// DiscreteProblem ///////////////////////////////////////////////////////////////////////////////////////

DiscreteProblem::DiscreteProblem(WeakForm *wf)
{
	_F_
	this->wf = wf;

	this->ndof = -1;
	spaces = new Space *[wf->neq];
	have_spaces = false;
	sp_seq = new int[wf->neq];
	memset(sp_seq, -1, sizeof(int) * wf->neq);

	matrix_buffer = NULL;
	matrix_buffer_dim = 0;
}

DiscreteProblem::~DiscreteProblem()
{
	_F_
	delete [] spaces;
	delete [] sp_seq;
}

void DiscreteProblem::free()
{
	_F_
	memset(sp_seq, -1, sizeof(int) * wf->neq);
}

void DiscreteProblem::set_spaces(Tuple<Space *> sp)
{
	_F_
	  int n = sp.size();
	if (n <= 0 || n > wf->neq) error("Bad number of spaces.");
        if (n != this->wf->neq) 
          error("Number of spaces must match the number of equations in LinProblem::set_spaces()"); 
	for (int i = 0; i < wf->neq; i++) this->spaces[i] = sp[i];
	memset(sp_seq, -1, sizeof(int) * wf->neq);
	have_spaces = true;
}

int DiscreteProblem::get_num_dofs()
{
	_F_
	if (!is_up_to_date()) {
		this->ndof = 0;
		for (int i = 0; i < wf->neq; i++)
			this->ndof += spaces[i]->get_dof_count();
	}
	return this->ndof;
}

scalar **DiscreteProblem::get_matrix_buffer(int n)
{
	_F_
	if (n <= matrix_buffer_dim) return matrix_buffer;
	if (matrix_buffer != NULL) delete [] matrix_buffer;
	matrix_buffer_dim = n;
	return (matrix_buffer = new_matrix<scalar>(n, n));
}

void DiscreteProblem::assemble(Vector *rhs, Matrix *jac, Vector *x)
{
	_F_
        // Sanity checks.
	if (x->length() != this->ndof) error("Wrong vector length in assemble().");
	if (!have_spaces) error("You have to call set_spaces() before calling assemble().");
        for (int i=0; i<this->wf->neq; i++) {
          if (this->spaces[i] == NULL) error("A space is NULL in assemble().");
        }
 
        // Extract values from the vector 'x'.
	scalar *vv;
        if (x != NULL) {
          vv = new scalar[this->ndof]; MEM_CHECK(vv);
	  memset(vv, 0, this->ndof * sizeof(scalar));
	  x->extract(vv);
        }

        // Convert the coefficient vector 'x' into solutions Tuple 'u_ext'.
        Tuple<Solution*> u_ext;
	for (int i = 0; i < this->wf->neq; i++) {
	  if (x != NULL) {
            u_ext.push_back(new Solution(this->spaces[i]->get_mesh()));
	    u_ext[i]->set_fe_solution(this->spaces[i], vv);
          }
          else u_ext.push_back(NULL);
	}
	if (x != NULL) delete [] vv;

        // Perform assembling.
        this->assemble(rhs, jac, u_ext);

        // Delete temporary solutions.
	for (int i = 0; i < wf->neq; i++) {
	  if (u_ext[i] != NULL) {
            delete u_ext[i];
	    u_ext[i] = NULL;
          }
	}
}


void DiscreteProblem::assemble(Vector *rhs, Matrix *jac, Tuple<Solution*> u_ext)
{
	_F_
        // Sanity checks.
	if (!have_spaces) error("You have to call set_spaces() before calling assemble().");
        for (int i=0; i<this->wf->neq; i++) {
          if (this->spaces[i] == NULL) error("A space is NULL in assemble().");
        }

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
	wf->get_stages(spaces, stages, jac == NULL);

	// Loop through all assembling stages -- the purpose of this is increased performance
	// in multi-mesh calculations, where, e.g., only the right hand side uses two meshes.
	// In such a case, the matrix forms are assembled over one mesh, and only the rhs
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

				u_ext[j]->set_active_element(e[i]);
				u_ext[j]->force_transform(base_fn[j].get_transform(), base_fn[j].get_ctm());

				refmap[j].set_active_element(e[i]);
				refmap[j].force_transform(base_fn[j].get_transform(), base_fn[j].get_ctm());
			}
			int marker = e0->marker;

			fn_cache.free();
			if (jac != NULL) {
				// assemble volume matrix forms //////////////////////////////////////
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
								scalar bi = eval_form(mfv, u_ext, fu, fv, refmap + n, refmap + m)
									* an->coef[j] * am->coef[i];
								if (an->dof[j] != H3D_DIRICHLET_DOF) mat[i][j] = bi;
							}
						}
						else { // symmetric block
							for (int j = 0; j < an->cnt; j++) {
								if (j < i && an->dof[j] >= 0) continue;
								fu->set_active_shape(an->idx[j]);
								scalar bi = eval_form(mfv, u_ext, fu, fv, refmap + n, refmap + m)
									* an->coef[j] * am->coef[i];
								if (an->dof[j] != H3D_DIRICHLET_DOF) mat[i][j] = mat[j][i] = bi;
							}
						}
					}

					// insert the local stiffness matrix into the global one
					jac->add(am->cnt, an->cnt, mat, am->dof, an->dof);

					// insert also the off-diagonal (anti-)symmetric block, if required
					if (tra) {
						if (mfv->sym < 0) chsgn(mat, am->cnt, an->cnt);
						transpose(mat, am->cnt, an->cnt);
						jac->add(an->cnt, am->cnt, mat, an->dof, am->dof);
					}
				}
			}

			// assemble volume vector forms ////////////////////////////////////////
			if (rhs != NULL) {
				for (unsigned ww = 0; ww < s->vfvol.size(); ww++) {
					WeakForm::VectorFormVol *vfv = s->vfvol[ww];
					if (isempty[vfv->i]) continue;
					if (vfv->area != ANY && !wf->is_in_area(marker, vfv->area)) continue;
					int m = vfv->i;  fv = test_fn + m;  am = al + m;

					for (int i = 0; i < am->cnt; i++) {
						if (am->dof[i] == H3D_DIRICHLET_DOF) continue;
						fv->set_active_shape(am->idx[i]);
						rhs->add(am->dof[i], eval_form(vfv, u_ext, fv, refmap + m) * am->coef[i]);
					}
				}
			}

			/////////////////////////

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

				if (jac != NULL) {
					// assemble surface matrix forms ///////////////////////////////////
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
								scalar bi = eval_form(mfs, u_ext, fu, fv, refmap + n, refmap + m,
									fp + iface) * an->coef[j] * am->coef[i];
								if (an->dof[j] != H3D_DIRICHLET_DOF) mat[i][j] = bi;
							}
						}
						jac->add(am->cnt, an->cnt, mat, am->dof, an->dof);
					}
				}

				// assemble surface vector forms /////////////////////////////////////
				if (rhs != NULL) {
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
							rhs->add(am->dof[i],
								eval_form(vfs, u_ext, fv, refmap + m, fp + iface) * am->coef[i]);
						}
					}
				}
			}
		}
		trav.finish();
	}

	delete [] matrix_buffer;
	matrix_buffer = NULL;
	matrix_buffer_dim = 0;
}

bool DiscreteProblem::is_up_to_date()
{
	_F_
	// check if we can reuse the matrix structure
	bool up_to_date = true;
	for (int i = 0; i < wf->neq; i++)
		if (spaces[i]->get_seq() != sp_seq[i]) {
			up_to_date = false;
			break;
		}
	return up_to_date;
}

void DiscreteProblem::create(SparseMatrix *mat)
{
	_F_
	assert(mat != NULL);

	if (is_up_to_date()) {
		// reuse the matrix structure
		mat->zero();
		return;
	}
	mat->free();

	// int ndof = get_num_dofs();
	mat->prealloc(this->ndof);

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

void DiscreteProblem::init_ext_fns(user_data_t<scalar> &ud, std::vector<MeshFunction *> &ext, int order,
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

void DiscreteProblem::init_ext_fns(user_data_t<ord_t> &fake_ud, std::vector<MeshFunction *> &ext)
{
	_F_

	fake_ud.nf = ext.size();
	fn_t<ord_t> *fake_ext_fn = new fn_t<ord_t>[fake_ud.nf];
	for (int i = 0; i < fake_ud.nf; i++) {
		fake_ext_fn[i] = init_fn(ext[i]->get_fn_order());
	}
	fake_ud.ext = fake_ext_fn;
}

sfn_t *DiscreteProblem::get_fn(ShapeFunction *fu, int order, RefMap *rm, const int np, const QuadPt3D *pt)
{
	fn_key_t key(fu->get_active_shape(), order, fu->get_transform(), fu->get_shapeset()->id);
	sfn_t *u = NULL;
	if (!fn_cache.fn.lookup(key, u)) {
		u = init_fn(fu, rm, np, pt);
		fn_cache.fn.set(key, u);
	}
	return u;
}

sfn_t *DiscreteProblem::get_fn(ShapeFunction *fu, int order, RefMap *rm, int iface, const int np,
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

mfn_t *DiscreteProblem::get_fn(Solution *fu, int order, RefMap *rm, const int np, const QuadPt3D *pt)
{
	fn_key_t key(fu->seq, order, fu->get_transform());
	mfn_t *u = NULL;
	if (!fn_cache.sln.lookup(key, u)) {
		u = init_fn(fu, rm, np, pt);
		fn_cache.sln.set(key, u);
	}
	return u;
}

scalar DiscreteProblem::eval_form(WeakForm::MatrixFormVol *mfv, Tuple<Solution *> u_ext, ShapeFunction *fu,
                            ShapeFunction *fv, RefMap *ru, RefMap *rv)
{
	_F_
	Element *elem = fv->get_active_element();

	// determine the integration order
	fn_t<ord_t> *oi = new fn_t<ord_t>[wf->neq];
	for (int i = 0; i < wf->neq; i++) oi[i] = init_fn(u_ext[i]->get_fn_order());
	fn_t<ord_t> ou = init_fn(fu->get_fn_order());
	fn_t<ord_t> ov = init_fn(fv->get_fn_order());

	user_data_t<ord_t> fake_ud;
	init_ext_fns(fake_ud, mfv->ext);

	double fake_wt = 1.0;
	geom_t<ord_t> fake_e = init_geom(elem->marker);
	ord_t o = mfv->ord(1, &fake_wt, &oi, &ou, &ov, &fake_e, &fake_ud);
	order3_t order = ru->get_inv_ref_order();
	switch (order.type) {
		case MODE_TETRAHEDRON: order += order3_t(o.get_order()); break;
		case MODE_HEXAHEDRON: order += order3_t(o.get_order(), o.get_order(), o.get_order()); break;
	}
	order.limit();
	int ord_idx = order.get_idx();

	for (int i = 0; i < wf->neq; i++) free_fn(oi + i);
	delete [] oi;
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

	mfn_t *prev[wf->neq];
	for (int i = 0; i < wf->neq; i++) prev[i] = get_fn(u_ext[i], ord_idx, rv, np, pt);
	sfn_t *u = get_fn(fu, ord_idx, ru, np, pt);
        sfn_t *v = get_fn(fv, ord_idx, rv, np, pt);

	user_data_t<scalar> ud;
	init_ext_fns(ud, mfv->ext, ord_idx, rv, np, pt);

	return mfv->fn(np, jwt, prev, u, v, &e, &ud);
}

scalar DiscreteProblem::eval_form(WeakForm::VectorFormVol *vfv, Tuple<Solution *> u_ext, ShapeFunction *fv, RefMap *rv)
{
	_F_
	Element *elem = fv->get_active_element();

	// determine the integration order
	fn_t<ord_t> *oi = new fn_t<ord_t>[wf->neq];
	for (int i = 0; i < wf->neq; i++) oi[i] = init_fn(u_ext[i]->get_fn_order());
	fn_t<ord_t> ov = init_fn(fv->get_fn_order());

	user_data_t<ord_t> fake_ud;
	init_ext_fns(fake_ud, vfv->ext);

	double fake_wt = 1.0;
	geom_t<ord_t> fake_e = init_geom(elem->marker);
	ord_t o = vfv->ord(1, &fake_wt, &oi, &ov, &fake_e, &fake_ud);
	order3_t order = rv->get_inv_ref_order();
	switch (order.type) {
		case MODE_TETRAHEDRON: order += order3_t(o.get_order()); break;
		case MODE_HEXAHEDRON: order += order3_t(o.get_order(), o.get_order(), o.get_order()); break;
	}
	order.limit();
	int ord_idx = order.get_idx();

	for (int i = 0; i < wf->neq; i++) free_fn(oi + i);
	delete [] oi;
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

	mfn_t *prev[wf->neq];
	for (int i = 0; i < wf->neq; i++) prev[i] = get_fn(u_ext[i], ord_idx, rv, np, pt);
	sfn_t *v = get_fn(fv, ord_idx, rv, np, pt);

	user_data_t<scalar> ud;
	init_ext_fns(ud, vfv->ext, ord_idx, rv, np, pt);

	return vfv->fn(np, jwt, prev, v, &e, &ud);
}

scalar DiscreteProblem::eval_form(WeakForm::MatrixFormSurf *mfv, Tuple<Solution *> u_ext, ShapeFunction *fu,
                            ShapeFunction *fv, RefMap *ru, RefMap *rv, FacePos *fp)
{
	_F_

	// determine the integration order
	fn_t<ord_t> *oi = new fn_t<ord_t>[wf->neq];
	for (int i = 0; i < wf->neq; i++) oi[i] = init_fn(u_ext[i]->get_fn_order());
	fn_t<ord_t> ou = init_fn(fu->get_fn_order());
	fn_t<ord_t> ov = init_fn(fv->get_fn_order());

	user_data_t<ord_t> fake_ud;
	init_ext_fns(fake_ud, mfv->ext);

	double fake_wt = 1.0;
	geom_t<ord_t> fake_e = init_geom(fp->marker);
	ord_t o = mfv->ord(1, &fake_wt, &oi, &ou, &ov, &fake_e, &fake_ud);
	order3_t order = ru->get_inv_ref_order();
	switch (order.type) {
		case MODE_TETRAHEDRON: order += order3_t(o.get_order()); break;
		case MODE_HEXAHEDRON: order += order3_t(o.get_order(), o.get_order(), o.get_order()); break;
	}
	order.limit();
	order2_t face_order = order.get_face_order(fp->face);
	int ord_idx = face_order.get_idx();

	for (int i = 0; i < wf->neq; i++) free_fn(oi + i);
	delete [] oi;
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

	mfn_t *prev[wf->neq];
	for (int i = 0; i < wf->neq; i++) prev[i] = get_fn(u_ext[i], ord_idx, rv, np, pt);
	sfn_t *u = get_fn(fu, ord_idx, ru, fp->face, np, pt);
	sfn_t *v = get_fn(fv, ord_idx, rv, fp->face, np, pt);

	user_data_t<scalar> ud;
	init_ext_fns(ud, mfv->ext, ord_idx, rv, np, pt);

	return mfv->fn(np, jwt, prev, u, v, &e, &ud);
}

scalar DiscreteProblem::eval_form(WeakForm::VectorFormSurf *vfs, Tuple<Solution *> u_ext, 
                                  ShapeFunction *fv, RefMap *rv, FacePos *fp)
{
	_F_

	// determine the integration order
	user_data_t<ord_t> fake_ud;
	init_ext_fns(fake_ud, vfs->ext);

	fn_t<ord_t> *oi = new fn_t<ord_t>[wf->neq];
	for (int i = 0; i < wf->neq; i++) oi[i] = init_fn(u_ext[i]->get_fn_order());
	fn_t<ord_t> ov = init_fn(fv->get_fn_order());
	double fake_wt = 1.0;
	geom_t<ord_t> fake_e = init_geom(fp->marker);
	ord_t o = vfs->ord(1, &fake_wt, &oi, &ov, &fake_e, &fake_ud);
	order3_t order = rv->get_inv_ref_order();
	switch (order.type) {
		case MODE_TETRAHEDRON: order += order3_t(o.get_order()); break;
		case MODE_HEXAHEDRON: order += order3_t(o.get_order(), o.get_order(), o.get_order()); break;
	}
	order.limit();
	order2_t face_order = order.get_face_order(fp->face);
	int ord_idx = face_order.get_idx();

	for (int i = 0; i < wf->neq; i++) free_fn(oi + i);
	delete [] oi;
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

	mfn_t *prev[wf->neq];
	for (int i = 0; i < wf->neq; i++) prev[i] = get_fn(u_ext[i], ord_idx, rv, np, pt);
	sfn_t *v = get_fn(fv, ord_idx, rv, fp->face, np, pt);

	user_data_t<scalar> ud;
	init_ext_fns(ud, vfs->ext, ord_idx, rv, np, pt);

	return vfs->fn(np, jwt, prev, v, &e, &ud);
}
