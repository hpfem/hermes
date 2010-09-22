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

//  Solvers
#include "solver/amesos.h"
#include "solver/pardiso.h"
#include "solver/petsc.h"
#include "solver/mumps.h"
#include "solver/nox.h"
#include "solver/umfpack_solver.h"

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

DiscreteProblem::DiscreteProblem(WeakForm *wf, Tuple<Space *> spaces, bool is_linear)
{
  _F_
  // sanity checks
  int n = spaces.size();
  if (n != wf->neq) error("Bad number of spaces in DiscreteProblem::DiscreteProblem().");

  this->wf = wf;
  this->spaces = spaces;
  this->is_linear = is_linear;

  sp_seq = new int[wf->neq];
  memset(sp_seq, -1, sizeof(int) * wf->neq);
  wf_seq = -1;

  matrix_buffer = NULL;
  matrix_buffer_dim = 0;

  values_changed = true;
  struct_changed = true;

  have_matrix = false;

  this->spaces = Tuple<Space *>();
  for (int i = 0; i < wf->neq; i++) this->spaces.push_back(spaces[i]);
  have_spaces = true;

  // H2D is initializing precalc shapesets here.

  // Create global enumeration of dof and fill the ndof variable
  this->ndof = assign_dofs(this->spaces);
}

DiscreteProblem::~DiscreteProblem()
{
  _F_
  free();
  if (sp_seq != NULL) delete [] sp_seq;
  wf_seq = -1;
}

void DiscreteProblem::free()
{
  _F_
  struct_changed = values_changed = true;
  memset(sp_seq, -1, sizeof(int) * wf->neq);
  wf_seq = -1;
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

//// assembly //////////////////////////////////////////////////////////////////////////////////////

// Light version for linear problems.
void DiscreteProblem::assemble(SparseMatrix* mat, Vector* rhs, bool rhsonly) 
{
  _F_
  assemble(NULL, mat, rhs, rhsonly);
}

void DiscreteProblem::assemble(scalar* coeff_vec, SparseMatrix* mat, Vector* rhs, bool rhsonly)
{
  /* BEGIN IDENTICAL CODE WITH H2D */

  _F_
  // Sanity checks.
  if (coeff_vec == NULL && this->is_linear == false) error("coeff_vec is NULL in FeProblem::assemble().");
  if (rhs != NULL && rhs->length() != this->ndof) error("Wrong rhs_ext length in FeProblem::assemble().");
  if (!have_spaces) error("You have to call FeProblem::set_spaces() before calling assemble().");
  for (int i=0; i<this->wf->neq; i++) {
    if (this->spaces[i] == NULL) error("A space is NULL in assemble().");
  }
 
  this->create(mat, rhs, rhsonly);

  // Convert the coefficient vector 'coeff_vec' into solutions Tuple 'u_ext'.
  Tuple<Solution*> u_ext;
  for (int i = 0; i < this->wf->neq; i++) 
  {
    if (this->is_linear == false)
    {
      u_ext.push_back(new Solution(this->spaces[i]->get_mesh()));
      u_ext[i]->set_coeff_vector(this->spaces[i], coeff_vec);
    }
    else
      u_ext.push_back(NULL);
  }

  /* END IDENTICAL CODE WITH H3D */

  bool bnd[10];			    // FIXME: magic number - maximal possible number of element surfaces
  SurfPos surf_pos[10];
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

  // initialize matrix buffer
  matrix_buffer = NULL;
  matrix_buffer_dim = 0;
  get_matrix_buffer(10);

  // obtain a list of assembling stages
  std::vector<WeakForm::Stage> stages;
  wf->get_stages(spaces, this->is_linear ? NULL : u_ext, stages, rhsonly);

  // Loop through all assembling stages -- the purpose of this is increased performance
  // in multi-mesh calculations, where, e.g., only the right hand side uses two meshes.
  // In such a case, the matrix forms are assembled over one mesh, and only the rhs
  // traverses through the union mesh. On the other hand, if you don't use multi-mesh
  // at all, there will always be only one stage in which all forms are assembled as usual.
  Traverse trav;
  for (unsigned ss = 0; ss < stages.size(); ss++) {
    WeakForm::Stage *s = &stages[ss];
    for (unsigned i = 0; i < s->idx.size(); i++) s->fns[i] = &base_fn[s->idx[i]];
    trav.begin(s->meshes.size(), &(s->meshes.front()), &(s->fns.front()));

    // assemble one stage
    Element **e;
    while ((e = trav.get_next_state(bnd, surf_pos)) != NULL) 
    {
      // find a non-NULL e[i]
      Element *e0;
      for (unsigned int i = 0; i < s->idx.size(); i++)
	if ((e0 = e[i]) != NULL) break;
      if (e0 == NULL) continue;

      /* H2D has here:
      update_limit_table(e0->get_mode()); */

      // Obtain assembly lists for the element at all spaces of the stage, set appropriate mode for each pss.
      // NOTE: Active elements and transformations for external functions (including the solutions from previous
      // Newton's iteration) as well as basis functions (master PrecalcShapesets) have already been set in 
      // trav.get_next_state(...).
      memset(isempty, 0, sizeof(bool) * wf->neq);
      for (unsigned i = 0; i < s->idx.size(); i++) 
      {
        int j = s->idx[i];
        if (e[i] == NULL) { isempty[j] = true; continue; }

        // TODO: do not obtain again if the element was not changed.
        spaces[j]->get_element_assembly_list(e[i], al + j);

        // This is different in H2D (PrecalcShapeset is used).
        test_fn[j].set_active_element(e[i]);
        test_fn[j].set_transform(base_fn + j);

        // This is missing in H2D.
	u_ext[j]->set_active_element(e[i]);
	u_ext[j]->force_transform(base_fn[j].get_transform(), base_fn[j].get_ctm());

        // This is different in H2D (PrecalcShapeset is used).
	refmap[j].set_active_element(e[i]);
	refmap[j].force_transform(base_fn[j].get_transform(), base_fn[j].get_ctm());
      }
      int marker = e0->marker;

      fn_cache.free();  // This is different in H2D.

      if (mat != NULL) {
	// assemble volume matrix forms //////////////////////////////////////
	for (unsigned ww = 0; ww < s->mfvol.size(); ww++) {
	  WeakForm::MatrixFormVol *mfv = s->mfvol[ww];
	  if (isempty[mfv->i] || isempty[mfv->j]) continue;
	  if (mfv->area != HERMES_ANY && !wf->is_in_area(marker, mfv->area)) continue;
	  int m = mfv->i; fv = test_fn + m; am = al + m;
	  int n = mfv->j; fu = base_fn + n; an = al + n;
	  bool tra = (m != n) && (mfv->sym != UNSYM);
	  bool sym = (m == n) && (mfv->sym == SYM);

	  /* BEGIN IDENTICAL CODE WITH H2D */

	  // assemble the local stiffness matrix for the form mfv
          scalar **local_stiffness_matrix = get_matrix_buffer(std::max(am->cnt, an->cnt));
	  for (int i = 0; i < am->cnt; i++) 
          {
            if (!tra && am->dof[i] < 0) continue;
            fv->set_active_shape(am->idx[i]);

	    if (!sym) // unsymmetric block
            { 
	      for (int j = 0; j < an->cnt; j++) 
              {
	        fu->set_active_shape(an->idx[j]);
                if (an->dof[j] < 0) 
                {
                  // Linear problems only: Subtracting Dirichlet lift contribution from the RHS:
                  if (rhs != NULL && this->is_linear) 
                  {
                    scalar val = eval_form(mfv, u_ext, fu, fv, refmap + n, refmap + m) * an->coef[j] * am->coef[i];
                    rhs->add(am->dof[i], -val);
                  } 
                }
                else if (rhsonly == false) {
	          scalar val = eval_form(mfv, u_ext, fu, fv, refmap + n, refmap + m) * an->coef[j] * am->coef[i];
                  local_stiffness_matrix[i][j] = val;
                }
              }
	    }
	    else // symmetric block
            {
	      for (int j = 0; j < an->cnt; j++) 
              {
	        if (j < i && an->dof[j] >= 0) continue;
	        fu->set_active_shape(an->idx[j]);
                if (an->dof[j] < 0) 
                {
                  // Linear problems only: Subtracting Dirichlet lift contribution from the RHS:
                  if (rhs != NULL && this->is_linear) 
                  {
                    scalar val = eval_form(mfv, u_ext, fu, fv, refmap + n, refmap + m) * an->coef[j] * am->coef[i];
                    rhs->add(am->dof[i], -val);
                  }
                } 
                else if (rhsonly == false) 
                {
                  scalar val = eval_form(mfv, u_ext, fu, fv, refmap + n, refmap + m) * an->coef[j] * am->coef[i];
                  local_stiffness_matrix[i][j] = local_stiffness_matrix[j][i] = val;
                } 
	      }
	    }
	  }

	  // insert the local stiffness matrix into the global one
          if (rhsonly == false)
            mat->add(am->cnt, an->cnt, local_stiffness_matrix, am->dof, an->dof);

          // insert also the off-diagonal (anti-)symmetric block, if required
          if (tra) 
          {
            if (mfv->sym < 0) chsgn(local_stiffness_matrix, am->cnt, an->cnt);
            transpose(local_stiffness_matrix, am->cnt, an->cnt);

            if (rhsonly == false) 
              mat->add(an->cnt, am->cnt, local_stiffness_matrix, an->dof, am->dof);

            // Linear problems only: Subtracting Dirichlet lift contribution from the RHS:
            if (rhs != NULL && this->is_linear) 
            {
              for (int j = 0; j < am->cnt; j++) 
              {
                if (am->dof[j] < 0) 
                {
                  for (int i = 0; i < an->cnt; i++) 
                  {
                    if (an->dof[i] >= 0) 
                    {
                      rhs->add(an->dof[i], -local_stiffness_matrix[i][j]);
                    }
                  }
                }
              }
            }
          }
        }
      }

      /* END IDENTICAL CODE WITH H2D
         Assembling of volume vector forms below is almost identical, there
         is only one line of difference that is highlighted below */

      //// assemble volume vector forms ////////////////////////////////////////
      if (rhs != NULL) 
      {
        for (unsigned int ww = 0; ww < s->vfvol.size(); ww++) 
        {
          WeakForm::VectorFormVol* vfv = s->vfvol[ww];
          if (isempty[vfv->i]) continue;
          if (vfv->area != HERMES_ANY && !wf->is_in_area(marker, vfv->area)) continue;
          int m = vfv->i;  
          fv = test_fn + m;      // H2D uses fv = spss[m]
          am = al + m;

          for (int i = 0; i < am->cnt; i++) 
          {
            if (am->dof[i] < 0) continue;
            fv->set_active_shape(am->idx[i]);
            scalar val = eval_form(vfv, u_ext, fv, refmap + m) * am->coef[i];
            rhs->add(am->dof[i], val);
          }
        }
      }

      // assemble surface integrals now: loop through surfaces of the element
      for (unsigned int isurf = 0; isurf < e0->get_num_surf(); isurf++) 
      {
        fn_cache.free();  // This is not in H2D.

	if (!bnd[isurf]) continue;
        int marker = surf_pos[isurf].marker;

        // obtain the list of shape functions which are nonzero on this surface
        for (unsigned int i = 0; i < s->idx.size(); i++) 
        {
          if (e[i] == NULL) continue;
          int j = s->idx[i];
          if ((nat[j] = (spaces[j]->bc_type_callback(marker) == BC_NATURAL)))
             spaces[j]->get_boundary_assembly_list(e[i], isurf, al + j);
        }

        // assemble surface matrix forms ///////////////////////////////////
        if (mat != NULL) 
        {
          for (unsigned int ww = 0; ww < s->mfsurf.size(); ww++) 
          {
            WeakForm::MatrixFormSurf *mfs = s->mfsurf[ww];
            if (isempty[mfs->i] || isempty[mfs->j]) continue;
            if (mfs->area != HERMES_ANY && !wf->is_in_area(marker, mfs->area)) continue;
            int m = mfs->i; 
            int n = mfs->j; 
            fu = base_fn + n;    // This is different in H2D.
            fv = test_fn + m;    // This is different in H2D.
            am = al + m;
            an = al + n;

            if (!nat[m] || !nat[n]) continue;
            surf_pos[isurf].base = trav.get_base();
            surf_pos[isurf].space_v = spaces[m];
            surf_pos[isurf].space_u = spaces[n];

            scalar **local_stiffness_matrix = get_matrix_buffer(std::max(am->cnt, an->cnt));
            for (int i = 0; i < am->cnt; i++) 
            {
	      if (am->dof[i] < 0) continue;
	      fv->set_active_shape(am->idx[i]);
	      for (int j = 0; j < an->cnt; j++) 
                {
	        fu->set_active_shape(an->idx[j]);
                if (an->dof[j] < 0) {
                  // Linear problems only: Subtracting Dirichlet lift contribution from the RHS:
                  if (rhs != NULL && this->is_linear) 
                  {
                    scalar val = eval_form(mfs, u_ext, fu, fv, refmap + n, refmap + m, 
                                           surf_pos + isurf) * an->coef[j] * am->coef[i];
                    rhs->add(am->dof[i], -val);
                  }
                }
                else if (rhsonly == false) {
                  scalar val = eval_form(mfs, u_ext, fu, fv, refmap + n, refmap + m, 
                                         surf_pos + isurf) * an->coef[j] * am->coef[i];
                  local_stiffness_matrix[i][j] = val;
                } 
              }
            }
            if (rhsonly == false) {
              mat->add(am->cnt, an->cnt, local_stiffness_matrix, an->dof, am->dof);
            }
	  }
        }

        // assemble surface vector forms /////////////////////////////////////
        if (rhs != NULL) 
        {
          for (unsigned int ww = 0; ww < s->vfsurf.size(); ww++) 
          {
            WeakForm::VectorFormSurf* vfs = s->vfsurf[ww];
            if (isempty[vfs->i]) continue;
            if (vfs->area != HERMES_ANY && !wf->is_in_area(marker, vfs->area)) continue;
            int m = vfs->i; 
            fv = test_fn + m;      // This is different from H2D.  
            am = al + m;

            if (!nat[m]) continue;
            surf_pos[isurf].base = trav.get_base();
            surf_pos[isurf].space_v = spaces[m];

            for (int i = 0; i < am->cnt; i++) 
            {
      	      if (am->dof[i] < 0) continue;
	      fv->set_active_shape(am->idx[i]);
              scalar val = eval_form(vfs, u_ext, fv, refmap + m, surf_pos + isurf) * am->coef[i];
              rhs->add(am->dof[i], val);
	    }
	  }
	}
      }
  
      // H2D is deleting cache here.
    }
    trav.finish(); 
  }
 
  // Cleaning up.
  delete [] matrix_buffer;
  matrix_buffer = NULL;
  matrix_buffer_dim = 0;

  // Delete temporary solutions.
  for (int i = 0; i < wf->neq; i++) 
  {
    if (u_ext[i] != NULL) 
    {
      delete u_ext[i];
      u_ext[i] = NULL;
    }
  }
}

//// matrix structure precalculation ///////////////////////////////////////////////////////////////

// This functions is identical in H2D and H3D.
bool DiscreteProblem::is_up_to_date()
{
  _F_
  // check if we can reuse the matrix structure
  bool up_to_date = true;
  if (!have_matrix) up_to_date = false;
  for (int i = 0; i < wf->neq; i++)
    if (spaces[i]->get_seq() != sp_seq[i])
      { up_to_date = false; break; }
  if (wf->get_seq() != wf_seq)
    up_to_date = false;

  return up_to_date;
}

//// matrix creation ///////////////////////////////////////////////////////////////////////////////

// This function is identical to H2D.
void DiscreteProblem::create(SparseMatrix *mat, Vector* rhs, bool rhsonly)
{
  _F_
  assert(mat != NULL);

  if (is_up_to_date())
  {
    printf("Reusing matrix sparse structure.\n");
    if (!rhsonly)
      mat->zero();
    rhs->zero();
    return;
  }

  // spaces have changed: create the matrix from scratch
  mat->free();

  int ndof = get_num_dofs();
  mat->prealloc(this->ndof);

  // This is different in H2D.
  AsmList al[wf->neq];
  Mesh *meshes[wf->neq];
  bool **blocks = wf->get_blocks();

  // init multi-mesh traversal.
  for (int i = 0; i < wf->neq; i++)
    meshes[i] = spaces[i]->get_mesh();

  Traverse trav;
  trav.begin(wf->neq, meshes);

  // Loop through all elements.
  Element **e;
  while ((e = trav.get_next_state(NULL, NULL)) != NULL) 
  {
  // obtain assembly lists for the element at all spaces
  for (int i = 0; i < wf->neq; i++)
    // TODO: do not get the assembly list again if the element was not changed
    if (e[i] != NULL) spaces[i]->get_element_assembly_list(e[i], al + i);

    // go through all equation-blocks of the local stiffness matrix
    for (int m = 0; m < wf->neq; m++)
      for (int n = 0; n < wf->neq; n++)
	if (blocks[m][n] && e[m] != NULL && e[n] != NULL) 
        {
	  AsmList *am = al + m;
	  AsmList *an = al + n;

	  // pretend assembling of the element stiffness matrix
	  // register nonzero elements
	  for (int i = 0; i < am->cnt; i++)
	    if (am->dof[i] >= 0)
	      for (int j = 0; j < an->cnt; j++)
	        if (an->dof[j] >= 0)
	          mat->pre_add_ij(am->dof[i], an->dof[j]);
	}
  }

  trav.finish();
  delete [] blocks;

  mat->alloc();
  if (rhs != NULL) rhs->alloc(ndof);

  // save space seq numbers and weakform seq number, so we can detect their changes
  for (int i = 0; i < wf->neq; i++)
    sp_seq[i] = spaces[i]->get_seq();
  wf_seq = wf->get_seq();

  struct_changed = true;
  have_matrix = true;
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

sfn_t *DiscreteProblem::get_fn(ShapeFunction *fu, int order, RefMap *rm, int isurf, const int np,
                         const QuadPt3D *pt)
{
	fn_key_t key(fu->get_active_shape(), order, fu->get_transform(), fu->get_shapeset()->id);
	sfn_t *u = NULL;
	if (!fn_cache.fn.lookup(key, u)) {
		u = init_fn(fu, rm, isurf, np, pt);
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
                            ShapeFunction *fv, RefMap *ru, RefMap *rv, SurfPos *surf_pos)
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
	geom_t<ord_t> fake_e = init_geom(surf_pos->marker);
	ord_t o = mfv->ord(1, &fake_wt, &oi, &ou, &ov, &fake_e, &fake_ud);
	order3_t order = ru->get_inv_ref_order();
	switch (order.type) {
		case MODE_TETRAHEDRON: order += order3_t(o.get_order()); break;
		case MODE_HEXAHEDRON: order += order3_t(o.get_order(), o.get_order(), o.get_order()); break;
	}
	order.limit();
	order2_t face_order = order.get_face_order(surf_pos->surf_num);
	int ord_idx = face_order.get_idx();

	for (int i = 0; i < wf->neq; i++) free_fn(oi + i);
	delete [] oi;
	free_fn(&ou);
	free_fn(&ov);

	// eval the form
	Quad3D *quad = get_quadrature(fu->get_active_element()->get_mode());
	int np = quad->get_face_num_points(surf_pos->surf_num, face_order);
	QuadPt3D *pt = quad->get_face_points(surf_pos->surf_num, face_order);

	double *jwt = NULL;
	geom_t<double> e;
	if (!fn_cache.e.exists(ord_idx)) {
		fn_cache.jwt[ord_idx] = ru->get_face_jacobian(surf_pos->surf_num, np, pt);
		fn_cache.e[ord_idx] = init_geom(surf_pos->marker, ru, surf_pos->surf_num, np, pt);
	}
	jwt = fn_cache.jwt[ord_idx];
	e = fn_cache.e[ord_idx];

	mfn_t *prev[wf->neq];
	for (int i = 0; i < wf->neq; i++) prev[i] = get_fn(u_ext[i], ord_idx, rv, np, pt);
	sfn_t *u = get_fn(fu, ord_idx, ru, surf_pos->surf_num, np, pt);
	sfn_t *v = get_fn(fv, ord_idx, rv, surf_pos->surf_num, np, pt);

	user_data_t<scalar> ud;
	init_ext_fns(ud, mfv->ext, ord_idx, rv, np, pt);

	return mfv->fn(np, jwt, prev, u, v, &e, &ud);
}

scalar DiscreteProblem::eval_form(WeakForm::VectorFormSurf *vfs, Tuple<Solution *> u_ext, 
                                  ShapeFunction *fv, RefMap *rv, SurfPos *surf_pos)
{
	_F_

	// determine the integration order
	user_data_t<ord_t> fake_ud;
	init_ext_fns(fake_ud, vfs->ext);

	fn_t<ord_t> *oi = new fn_t<ord_t>[wf->neq];
	for (int i = 0; i < wf->neq; i++) oi[i] = init_fn(u_ext[i]->get_fn_order());
	fn_t<ord_t> ov = init_fn(fv->get_fn_order());
	double fake_wt = 1.0;
	geom_t<ord_t> fake_e = init_geom(surf_pos->marker);
	ord_t o = vfs->ord(1, &fake_wt, &oi, &ov, &fake_e, &fake_ud);
	order3_t order = rv->get_inv_ref_order();
	switch (order.type) {
		case MODE_TETRAHEDRON: order += order3_t(o.get_order()); break;
		case MODE_HEXAHEDRON: order += order3_t(o.get_order(), o.get_order(), o.get_order()); break;
	}
	order.limit();
	order2_t face_order = order.get_face_order(surf_pos->surf_num);
	int ord_idx = face_order.get_idx();

	for (int i = 0; i < wf->neq; i++) free_fn(oi + i);
	delete [] oi;
	free_fn(&ov);

	// eval the form
	Quad3D *quad = get_quadrature(fv->get_active_element()->get_mode());
	int np = quad->get_face_num_points(surf_pos->surf_num, face_order);
	QuadPt3D *pt = quad->get_face_points(surf_pos->surf_num, face_order);

	double *jwt = NULL;
	geom_t<double> e;
	if (!fn_cache.e.exists(ord_idx)) {
		fn_cache.jwt[ord_idx] = rv->get_face_jacobian(surf_pos->surf_num, np, pt);
		fn_cache.e[ord_idx] = init_geom(surf_pos->marker, rv, surf_pos->surf_num, np, pt);
	}
	jwt = fn_cache.jwt[ord_idx];
	e = fn_cache.e[ord_idx];

	mfn_t *prev[wf->neq];
	for (int i = 0; i < wf->neq; i++) prev[i] = get_fn(u_ext[i], ord_idx, rv, np, pt);
	sfn_t *v = get_fn(fv, ord_idx, rv, surf_pos->surf_num, np, pt);

	user_data_t<scalar> ud;
	init_ext_fns(ud, vfs->ext, ord_idx, rv, np, pt);

	return vfs->fn(np, jwt, prev, v, &e, &ud);
}

////////////////////////////////////////////////////////////////////////////////////////

// This function is identical in H2D and H3D.
Vector* create_vector(MatrixSolverType matrix_solver)
{
  switch (matrix_solver) 
  {
    case SOLVER_AMESOS:
      {
        return new EpetraVector;
        break;
      }
    case SOLVER_MUMPS: 
      {
        return new MumpsVector;
        break;
      }
    case SOLVER_PARDISO: 
      {
        return new PardisoVector;
        break;
      }
    case SOLVER_PETSC: 
      {
        return new PetscVector;
        break;
      }
    case SOLVER_UMFPACK: 
      {
        return new UMFPackVector;
        break;
      }
    default: 
      error("Unknown matrix solver requested.");
  }
}

// This function is identical in H2D and H3D.
SparseMatrix* create_matrix(MatrixSolverType matrix_solver)
{
  switch (matrix_solver) 
  {
    case SOLVER_AMESOS:
      {
        return new EpetraMatrix;
        break;
      }
    case SOLVER_MUMPS: 
      {
        return new MumpsMatrix;
        break;
      }
    case SOLVER_PARDISO: 
      {
        return new PardisoMatrix;
        break;
      }
    case SOLVER_PETSC: 
      {
        return new PetscMatrix;
        break;
      }
    case SOLVER_UMFPACK: 
      {
        return new UMFPackMatrix;
        break;
      }
    default: 
      error("Unknown matrix solver requested.");
  }
}

// This function is identical in H2D and H3D.
Solver* create_solver(MatrixSolverType matrix_solver, Matrix* matrix, Vector* rhs)
{
  switch (matrix_solver) 
  {
    case SOLVER_AMESOS:
      {
        return new AmesosSolver("Amesos_Klu", static_cast<EpetraMatrix*>(matrix), static_cast<EpetraVector*>(rhs));
        printf("Using Amesos.\n"); 
        break;
      }
    case SOLVER_MUMPS: 
      {
        return new MumpsSolver(static_cast<MumpsMatrix*>(matrix), static_cast<MumpsVector*>(rhs)); 
        printf("Using Mumps.\n"); 
        break;
      }
    case SOLVER_PARDISO: 
      {
        return new PardisoLinearSolver(static_cast<PardisoMatrix*>(matrix), static_cast<PardisoVector*>(rhs));
        printf("Using Pardiso.\n"); 
        break;
      }
    case SOLVER_PETSC: 
      {
        return new PetscLinearSolver(static_cast<PetscMatrix*>(matrix), static_cast<PetscVector*>(rhs)); 
        printf("Using PETSc.\n");
        break;
      }
    case SOLVER_UMFPACK: 
      {
        return new UMFPackLinearSolver(static_cast<UMFPackMatrix*>(matrix), static_cast<UMFPackVector*>(rhs)); 
        printf("Using UMFPack.\n"); 
        break;
      }
    default: 
      error("Unknown matrix solver requested.");
  }
}
////////////////////////////////////////////////////////////////////////////////////////

// The projection functionality below is identical in H2D and H3D.
template<typename Real, typename Scalar>
Scalar H1projection_biform(int n, double *wt, fn_t<Scalar> *u_ext[], fn_t<Real> *u, fn_t<Real> *v, geom_t<Real> *e, user_data_t<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->val[i] * v->val[i] + u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar H1projection_liform(int n, double *wt, fn_t<Scalar> *u_ext[], fn_t<Real> *v, geom_t<Real> *e, user_data_t<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (ext->fn[0]->val[i] * v->val[i] + ext->fn[0]->dx[i] * v->dx[i] + ext->fn[0]->dy[i] * v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar L2projection_biform(int n, double *wt, fn_t<Scalar> *u_ext[], fn_t<Real> *u, fn_t<Real> *v, geom_t<Real> *e, user_data_t<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->val[i] * v->val[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar L2projection_liform(int n, double *wt, fn_t<Scalar> *u_ext[], fn_t<Real> *v, geom_t<Real> *e, user_data_t<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (ext->fn[0]->val[i] * v->val[i]);
  return result;
}

// Hcurl projections
template<typename Real, typename Scalar>
Scalar Hcurlprojection_biform(int n, double *wt, fn_t<Scalar> *u_ext[], fn_t<Real> *u, 
                              fn_t<Real> *v, geom_t<Real> *e, user_data_t<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++) {
    result += wt[i] * (u->curl[i] * conj(v->curl[i]));
    result += wt[i] * (u->val0[i] * conj(v->val0[i]) + u->val1[i] * conj(v->val1[i]));
  }
  return result;
}

template<typename Real, typename Scalar>
Scalar Hcurlprojection_liform(int n, double *wt, fn_t<Scalar> *u_ext[], fn_t<Real> *v, 
                              geom_t<Real> *e, user_data_t<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++) {
    result += wt[i] * (ext->fn[0]->curl[i] * conj(v->curl[i]));
    result += wt[i] * (ext->fn[0]->val0[i] * conj(v->val0[i]) + ext->fn[0]->val1[i] * conj(v->val1[i]));
  }

  return result;
}

/*
double get_l2_norm(Vector* vec) 
{
  scalar val = 0;
  for (int i = 0; i < vec->length(); i++) {
    scalar inc = vec->get(i);
    val = val + inc*conj(inc);
  }
  return sqrt(std::abs(val));
}
*/
