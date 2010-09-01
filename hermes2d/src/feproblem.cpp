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

#include "common.h"
#include "limit_order.h"
#include "feproblem.h"
#include "traverse.h"
#include "space.h"
#include "precalc.h"
#include "refmap.h"
#include "solution.h"
#include "config.h"
#include "discrete_problem.h"

FeProblem::FeProblem(WeakForm* wf, Tuple<Space *> spaces)
{
  // sanity checks
  int n = spaces.size();
  this->wf = wf;
  if (n != wf->neq) error("Bad number of spaces in FeProblem.");

  this->wf = wf;
  this->spaces = spaces;

  sp_seq = new int[wf->neq];
  memset(sp_seq, -1, sizeof(int) * wf->neq);
  wf_seq = -1;
  pss = new PrecalcShapeset*[wf->neq];
  num_user_pss = 0;

  buffer = NULL;
  mat_size = 0;

  values_changed = true;
  struct_changed = true;
  have_matrix = false;

  memset(sp_seq, -1, sizeof(int) * wf->neq);
  have_spaces = true;

  // initialize precalc shapesets
  this->pss = new PrecalcShapeset*[this->wf->neq];
  for (int i=0; i < n; i++) this->pss[i] = NULL;
  this->num_user_pss = 0;
  for (int i = 0; i < n; i++){
    Shapeset *shapeset = spaces[i]->get_shapeset();
    if (shapeset == NULL) error("Internal in DiscreteProblem::init_spaces().");
    PrecalcShapeset *p = new PrecalcShapeset(shapeset);
    if (p == NULL) error("New PrecalcShapeset could not be allocated in DiscreteProblem::init_spaces().");
    this-> pss[i] = p;
    this->num_user_pss++;
  }  

  // Create global enumeration of dof.
  assign_dofs(this->spaces);
}

FeProblem::~FeProblem()
{
  free();
  if (sp_seq != NULL) delete [] sp_seq;
  if (pss != NULL) delete [] pss;
}

/*
void FeProblem::set_pss(Tuple<PrecalcShapeset*> pss)
{
  int n = pss.size();
  //printf("size = %d\n", n);
  if (n != this->wf->neq) error("Bad number of pss in FeProblem.");
  for (int i = 0; i < n; i++) this->pss[i] = pss[i];
}
*/

void FeProblem::free()
{
  struct_changed = values_changed = true;
  memset(sp_seq, -1, sizeof(int) * wf->neq);
  wf_seq = -1;
}

int FeProblem::get_num_dofs()
{
  ndof = 0;
  for (int i = 0; i < wf->neq; i++)
    ndof += spaces[i]->get_num_dofs();
  return ndof;
}

//// matrix structure precalculation ///////////////////////////////////////////////////////////////

bool FeProblem::is_up_to_date()
{
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

void FeProblem::create(SparseMatrix* mat)
{
  assert(mat != NULL);

  if (is_up_to_date())
  {
    verbose("Reusing matrix sparse structure.");
    mat->zero();
    return;
  }

  // spaces have changed: create the matrix from scratch
  mat->free();

  int ndof = get_num_dofs();
  mat->prealloc(ndof);

  AUTOLA_CL(AsmList, al, wf->neq);
  AUTOLA_OR(Mesh*, meshes, wf->neq);
  bool **blocks = wf->get_blocks();

  // init multi-mesh traversal
  for (int i = 0; i < wf->neq; i++)
    meshes[i] = spaces[i]->get_mesh();

  Traverse trav;
  trav.begin(wf->neq, meshes);

  // loop through all elements
  Element **e;
  while ((e = trav.get_next_state(NULL, NULL)) != NULL)
  {
    // obtain assembly lists for the element at all spaces
    for (int i = 0; i < wf->neq; i++)
      // TODO: do not get the assembly list again if the element was not changed
      if (e[i] != NULL)
        spaces[i]->get_element_assembly_list(e[i], al + i);

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

  // save space seq numbers and weakform seq number, so we can detect their changes
  for (int i = 0; i < wf->neq; i++)
    sp_seq[i] = spaces[i]->get_seq();
  wf_seq = wf->get_seq();

  struct_changed = true;
  have_matrix = true;
}

//// assembly //////////////////////////////////////////////////////////////////////////////////////

void FeProblem::assemble(_Vector *rhs, _Matrix *jac, _Vector *x)
{
        // Sanity checks.
	if (x->length() != this->ndof) error("Wrong vector length in assemble().");
	if (!have_spaces) error("You have to call set_spaces() before calling assemble().");
        for (int i=0; i<this->wf->neq; i++) {
          if (this->spaces[i] == NULL) error("A space is NULL in assemble().");
        }
 
        // Extract values from the vector 'x'.
	Vector* vv;
        if (x != NULL) {
          vv = new AVector(this->ndof);
	  memset(vv->get_c_array(), 0, this->ndof * sizeof(scalar));
	  for (int i=0; i<this->ndof; i++) vv->set(i, x->get(i));
        }

        // Convert the coefficient vector 'x' into solutions Tuple 'u_ext'.
        Tuple<Solution*> u_ext;
	for (int i = 0; i < this->wf->neq; i++) {
	  if (x != NULL) {
            u_ext.push_back(new Solution(this->spaces[i]->get_mesh()));
	    u_ext[i]->set_fe_solution(this->spaces[i], this->pss[i], vv);
          }
          else u_ext.push_back(NULL);
	}
	if (x != NULL) delete vv;

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

void FeProblem::assemble(_Vector *rhs, _Matrix *jac, Tuple<Solution*> u_ext)
{
  // Sanity checks.
  if (!have_spaces) error("You have to call set_spaces() before calling assemble().");
  for (int i=0; i<this->wf->neq; i++) {
    if (this->spaces[i] == NULL) error("A space is NULL in assemble().");
  }

  int k, m, n, marker;
  AUTOLA_CL(AsmList, al, wf->neq);
  AsmList *am, *an;
  bool bnd[4];
  AUTOLA_OR(bool, nat, wf->neq);
  AUTOLA_OR(bool, isempty, wf->neq);
  EdgePos ep[4];
  reset_warn_order();

  // create slave pss's for test functions, init quadrature points
  AUTOLA_OR(PrecalcShapeset*, spss, wf->neq);
  PrecalcShapeset *fu, *fv;
  AUTOLA_CL(RefMap, refmap, wf->neq);
  for (int i = 0; i < wf->neq; i++)
  {
    spss[i] = new PrecalcShapeset(pss[i]);
    pss [i]->set_quad_2d(&g_quad_2d_std);
    spss[i]->set_quad_2d(&g_quad_2d_std);
    refmap[i].set_quad_2d(&g_quad_2d_std);
  }

  // initialize buffer
  buffer = NULL;
  mat_size = 0;
  get_matrix_buffer(9);

  // obtain a list of assembling stages
  std::vector<WeakForm::Stage> stages;
  wf->get_stages(spaces, NULL, stages, jac == NULL);

  // Loop through all assembling stages -- the purpose of this is increased performance
  // in multi-mesh calculations, where, e.g., only the right hand side uses two meshes.
  // In such a case, the matrix forms are assembled over one mesh, and only the rhs
  // traverses through the union mesh. On the other hand, if you don't use multi-mesh
  // at all, there will always be only one stage in which all forms are assembled as usual.
  Traverse trav;
  for (unsigned ss = 0; ss < stages.size(); ss++)
  {
    WeakForm::Stage* s = &stages[ss];
    for (unsigned i = 0; i < s->idx.size(); i++)
      s->fns[i] = pss[s->idx[i]];
    for (unsigned i = 0; i < s->ext.size(); i++)
      s->ext[i]->set_quad_2d(&g_quad_2d_std);
    trav.begin(s->meshes.size(), &(s->meshes.front()), &(s->fns.front()));

    // assemble one stage
    Element** e;
    while ((e = trav.get_next_state(bnd, ep)) != NULL)
    {
      // find a non-NULL e[i]
      Element* e0;
      for (unsigned int i = 0; i < s->idx.size(); i++)
        if ((e0 = e[i]) != NULL) break;
      if (e0 == NULL) continue;

      // set maximum integration order for use in integrals, see limit_order()
      update_limit_table(e0->get_mode());

      // obtain assembly lists for the element at all spaces, set appropriate mode for each pss
      memset(isempty, 0, sizeof(bool) * wf->neq);
      for (unsigned int i = 0; i < s->idx.size(); i++)
      {
        int j = s->idx[i];
        if (e[i] == NULL) { isempty[j] = true; continue; }
        spaces[j]->get_element_assembly_list(e[i], al+j);

        spss[j]->set_active_element(e[i]);
        spss[j]->set_master_transform();
        refmap[j].set_active_element(e[i]);
        refmap[j].force_transform(pss[j]->get_transform(), pss[j]->get_ctm());

        u_ext[j]->set_active_element(e[i]);
        u_ext[j]->force_transform(pss[j]->get_transform(), pss[j]->get_ctm());
      }
      marker = e0->marker;

      init_cache();
      //// assemble volume matrix forms //////////////////////////////////////
      if (jac != NULL)
      {
        for (unsigned ww = 0; ww < s->mfvol.size(); ww++)
        {
          WeakForm::MatrixFormVol* mfv = s->mfvol[ww];
          if (isempty[mfv->i] || isempty[mfv->j]) continue;
          if (mfv->area != H2D_ANY && !wf->is_in_area(marker, mfv->area)) continue;
          m = mfv->i;  fv = spss[m];  am = &al[m];
          n = mfv->j;  fu = pss[n];   an = &al[n];
          bool tra = (m != n) && (mfv->sym != 0);
          bool sym = (m == n) && (mfv->sym == 1);

          // assemble the local stiffness matrix for the form mfv
          scalar bi, **mat = get_matrix_buffer(std::max(am->cnt, an->cnt));
          for (int i = 0; i < am->cnt; i++)
          {
            if (!tra && (k = am->dof[i]) < 0) continue;
            fv->set_active_shape(am->idx[i]);

            if (!sym) // unsymmetric block
            {
              for (int j = 0; j < an->cnt; j++) {
                fu->set_active_shape(an->idx[j]);
                bi = eval_form(mfv, u_ext, fu, fv, refmap+n, refmap+m) * an->coef[j] * am->coef[i];
                if (an->dof[j] >= 0) mat[i][j] = bi;
              }
            }
            else // symmetric block
            {
              for (int j = 0; j < an->cnt; j++) {
                if (j < i && an->dof[j] >= 0) continue;
                fu->set_active_shape(an->idx[j]);
                bi = eval_form(mfv, u_ext, fu, fv, refmap+n, refmap+m) * an->coef[j] * am->coef[i];
                if (an->dof[j] >= 0) mat[i][j] = mat[j][i] = bi;
              }
            }
          }
          // insert the local stiffness matrix into the global one
          jac->add(am->cnt, an->cnt, mat, am->dof, an->dof);

          // insert also the off-diagonal (anti-)symmetric block, if required
          if (tra)
          {
            if (mfv->sym < 0) chsgn(mat, am->cnt, an->cnt);
            transpose(mat, am->cnt, an->cnt);
            jac->add(am->cnt, an->cnt, mat, am->dof, an->dof);
          }
        }
      }

      //// assemble volume linear forms ////////////////////////////////////////
      if (rhs != NULL)
      {
        for (unsigned int ww = 0; ww < s->vfvol.size(); ww++)
        {
          WeakForm::VectorFormVol* vfv = s->vfvol[ww];
          if (isempty[vfv->i]) continue;
          if (vfv->area != H2D_ANY && !wf->is_in_area(marker, vfv->area)) continue;
          m = vfv->i;  fv = spss[m];  am = &al[m];

          for (int i = 0; i < am->cnt; i++)
          {
            if (am->dof[i] < 0) continue;
            fv->set_active_shape(am->idx[i]);
            rhs->add(am->dof[i], eval_form(vfv, u_ext, fv, refmap + m) * am->coef[i]);
          }
        }
      }

      // assemble surface integrals now: loop through boundary edges of the element
      for (unsigned int edge = 0; edge < e0->nvert; edge++)
      {
        if (!bnd[edge]) continue;
        marker = ep[edge].marker;

        // obtain the list of shape functions which are nonzero on this edge
        for (unsigned int i = 0; i < s->idx.size(); i++) {
          if (e[i] == NULL) continue;
          int j = s->idx[i];
          if ((nat[j] = (spaces[j]->bc_type_callback(marker) == BC_NATURAL)))
            spaces[j]->get_edge_assembly_list(e[i], edge, al + j);
        }

        // assemble surface matrix forms ///////////////////////////////////
        if (jac != NULL)
        {
          for (unsigned int ww = 0; ww < s->mfsurf.size(); ww++)
          {
            WeakForm::MatrixFormSurf* mfs = s->mfsurf[ww];
            if (isempty[mfs->i] || isempty[mfs->j]) continue;
            if (mfs->area != H2D_ANY && !wf->is_in_area(marker, mfs->area)) continue;
            m = mfs->i;  fv = spss[m];  am = &al[m];
            n = mfs->j;  fu = pss[n];   an = &al[n];

            if (!nat[m] || !nat[n]) continue;
            ep[edge].base = trav.get_base();
            ep[edge].space_v = spaces[m];
            ep[edge].space_u = spaces[n];

            scalar bi, **mat = get_matrix_buffer(std::max(am->cnt, an->cnt));
            for (int i = 0; i < am->cnt; i++)
            {
              if ((k = am->dof[i]) < 0) continue;
              fv->set_active_shape(am->idx[i]);
              for (int j = 0; j < an->cnt; j++)
              {
                fu->set_active_shape(an->idx[j]);
                bi = eval_form(mfs, u_ext, fu, fv, refmap+n, refmap+m, ep+edge) * an->coef[j] * am->coef[i];
                if (an->dof[j] >= 0) mat[i][j] = bi;
              }
            }
            jac->add(am->cnt, an->cnt, mat, am->dof, an->dof);
          }
        }
        // assemble surface linear forms /////////////////////////////////////
        if (rhs != NULL)
        {
          for (unsigned ww = 0; ww < s->vfsurf.size(); ww++)
          {
            WeakForm::VectorFormSurf* vfs = s->vfsurf[ww];
            if (isempty[vfs->i]) continue;
            if (vfs->area != H2D_ANY && !wf->is_in_area(marker, vfs->area)) continue;
            m = vfs->i;  fv = spss[m];  am = &al[m];

            if (!nat[m]) continue;
            ep[edge].base = trav.get_base();
            ep[edge].space_v = spaces[m];

            for (int i = 0; i < am->cnt; i++)
            {
              if (am->dof[i] < 0) continue;
              fv->set_active_shape(am->idx[i]);
              rhs->add(am->dof[i], eval_form(vfs, u_ext, fv, refmap+m, ep+edge) * am->coef[i]);
            }
          }
        }
      }
      delete_cache();
    }
    trav.finish();
  }

  for (int i = 0; i < wf->neq; i++) delete spss[i];
  delete [] buffer;
  buffer = NULL;
  mat_size = 0;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Initialize integration order for external functions
ExtData<Ord>* FeProblem::init_ext_fns_ord(std::vector<MeshFunction *> &ext)
{
  ExtData<Ord>* fake_ext = new ExtData<Ord>;
  fake_ext->nf = ext.size();
  Func<Ord>** fake_ext_fn = new Func<Ord>*[fake_ext->nf];
  for (int i = 0; i < fake_ext->nf; i++)
    fake_ext_fn[i] = init_fn_ord(ext[i]->get_fn_order());
  fake_ext->fn = fake_ext_fn;

  return fake_ext;
}

// Initialize external functions (obtain values, derivatives,...)
ExtData<scalar>* FeProblem::init_ext_fns(std::vector<MeshFunction *> &ext, RefMap *rm, const int order)
{
  ExtData<scalar>* ext_data = new ExtData<scalar>;
  Func<scalar>** ext_fn = new Func<scalar>*[ext.size()];
  for (unsigned i = 0; i < ext.size(); i++)
    ext_fn[i] = init_fn(ext[i], rm, order);
  ext_data->nf = ext.size();
  ext_data->fn = ext_fn;

  return ext_data;

}

// Initialize shape function values and derivatives (fill in the cache)
Func<double>* FeProblem::get_fn(PrecalcShapeset *fu, RefMap *rm, const int order)
{
  Key key(256 - fu->get_active_shape(), order, fu->get_transform(), fu->get_shapeset()->get_id());
  if (cache_fn[key] == NULL)
    cache_fn[key] = init_fn(fu, rm, order);

  return cache_fn[key];
}

// Caching transformed values
void FeProblem::init_cache()
{
  for (int i = 0; i < g_max_quad + 1 + 4; i++)
  {
    cache_e[i] = NULL;
    cache_jwt[i] = NULL;
  }
}

void FeProblem::delete_cache()
{
  for (int i = 0; i < g_max_quad + 1 + 4; i++)
  {
    if (cache_e[i] != NULL)
    {
      cache_e[i]->free(); delete cache_e[i];
      delete [] cache_jwt[i];
    }
  }
  for (std::map<Key, Func<double>*, Compare>::const_iterator it = cache_fn.begin(); it != cache_fn.end(); it++)
  {
    (it->second)->free_fn(); delete (it->second);
  }
  cache_fn.clear();
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Actual evaluation of volume matrix form (calculates integral)
scalar FeProblem::eval_form(WeakForm::MatrixFormVol *mfv, Tuple<Solution *> u_ext, PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv)
{
  // determine the integration order
  int inc = (fu->get_num_components() == 2) ? 1 : 0;
  AUTOLA_OR(Func<Ord>*, oi, wf->neq);
  for (int i = 0; i < wf->neq; i++) oi[i] = init_fn_ord(u_ext[i]->get_fn_order() + inc);
  Func<Ord>* ou = init_fn_ord(fu->get_fn_order() + inc);
  Func<Ord>* ov = init_fn_ord(fv->get_fn_order() + inc);
  ExtData<Ord>* fake_ext = init_ext_fns_ord(mfv->ext);

  double fake_wt = 1.0;
  Geom<Ord>* fake_e = init_geom_ord();
  Ord o = mfv->ord(1, &fake_wt, oi, ou, ov, fake_e, fake_ext);
  int order = ru->get_inv_ref_order();
  order += o.get_order();
  limit_order_nowarn(order);

  for (int i = 0; i < wf->neq; i++) {  oi[i]->free_ord(); delete oi[i]; }
  ou->free_ord(); delete ou;
  ov->free_ord(); delete ov;
  delete fake_e;
  fake_ext->free_ord(); delete fake_ext;

  // eval the form
  Quad2D* quad = fu->get_quad_2d();
  double3* pt = quad->get_points(order);
  int np = quad->get_num_points(order);

  // init geometry and jacobian*weights
  if (cache_e[order] == NULL)
  {
    cache_e[order] = init_geom_vol(ru, order);
    double* jac = ru->get_jacobian(order);
    cache_jwt[order] = new double[np];
    for(int i = 0; i < np; i++)
      cache_jwt[order][i] = pt[i][2] * jac[i];
  }
  Geom<double>* e = cache_e[order];
  double* jwt = cache_jwt[order];

  // function values and values of external functions
  AUTOLA_OR(Func<scalar>*, prev, wf->neq);
  for (int i = 0; i < wf->neq; i++) prev[i]  = init_fn(u_ext[i], rv, order);
  Func<double>* u = get_fn(fu, ru, order);
  Func<double>* v = get_fn(fv, rv, order);
  ExtData<scalar>* ext = init_ext_fns(mfv->ext, rv, order);

  scalar res = mfv->fn(np, jwt, prev, u, v, e, ext);

  for (int i = 0; i < wf->neq; i++) {  prev[i]->free_fn(); delete prev[i]; }
  ext->free(); delete ext;
  return res;
}


// Actual evaluation of volume linear form (calculates integral)
scalar FeProblem::eval_form(WeakForm::VectorFormVol *vfv, Tuple<Solution *> u_ext, PrecalcShapeset *fv, RefMap *rv)
{
  // determine the integration order
  int inc = (fv->get_num_components() == 2) ? 1 : 0;
  AUTOLA_OR(Func<Ord>*, oi, wf->neq);
  for (int i = 0; i < wf->neq; i++) oi[i] = init_fn_ord(u_ext[i]->get_fn_order() + inc);
  Func<Ord>* ov = init_fn_ord(fv->get_fn_order() + inc);
  ExtData<Ord>* fake_ext = init_ext_fns_ord(vfv->ext);

  double fake_wt = 1.0;
  Geom<Ord>* fake_e = init_geom_ord();
  Ord o = vfv->ord(1, &fake_wt, oi, ov, fake_e, fake_ext);
  int order = rv->get_inv_ref_order();
  order += o.get_order();
  limit_order_nowarn(order);

  for (int i = 0; i < wf->neq; i++) {  oi[i]->free_ord(); delete oi[i]; }
  ov->free_ord(); delete ov;
  delete fake_e;
  fake_ext->free_ord(); delete fake_ext;

  // eval the form
  Quad2D* quad = fv->get_quad_2d();
  double3* pt = quad->get_points(order);
  int np = quad->get_num_points(order);

  // init geometry and jacobian*weights
  if (cache_e[order] == NULL)
  {
    cache_e[order] = init_geom_vol(rv, order);
    double* jac = rv->get_jacobian(order);
    cache_jwt[order] = new double[np];
    for(int i = 0; i < np; i++)
      cache_jwt[order][i] = pt[i][2] * jac[i];
  }
  Geom<double>* e = cache_e[order];
  double* jwt = cache_jwt[order];

  // function values and values of external functions
  AUTOLA_OR(Func<scalar>*, prev, wf->neq);
  for (int i = 0; i < wf->neq; i++) prev[i]  = init_fn(u_ext[i], rv, order);
  Func<double>* v = get_fn(fv, rv, order);
  ExtData<scalar>* ext = init_ext_fns(vfv->ext, rv, order);

  scalar res = vfv->fn(np, jwt, prev, v, e, ext);

  for (int i = 0; i < wf->neq; i++) {  prev[i]->free_fn(); delete prev[i]; }
  ext->free(); delete ext;
  return res;

}


// Actual evaluation of surface matrix form (calculates integral)
scalar FeProblem::eval_form(WeakForm::MatrixFormSurf *mfs, Tuple<Solution *> u_ext, PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv, EdgePos* ep)
{
  // eval the form
  Quad2D* quad = fu->get_quad_2d();
  int eo = quad->get_edge_points(ep->edge);
  double3* pt = quad->get_points(eo);
  int np = quad->get_num_points(eo);

  // init geometry and jacobian*weights
  if (cache_e[eo] == NULL)
  {
    cache_e[eo] = init_geom_surf(ru, ep, eo);
    double3* tan = ru->get_tangent(ep->edge);
    cache_jwt[eo] = new double[np];
    for(int i = 0; i < np; i++)
      cache_jwt[eo][i] = pt[i][2] * tan[i][2];
  }
  Geom<double>* e = cache_e[eo];
  double* jwt = cache_jwt[eo];

  // function values and values of external functions
  AUTOLA_OR(Func<scalar>*, prev, wf->neq);
  for (int i = 0; i < wf->neq; i++) prev[i]  = init_fn(u_ext[i], rv, eo);
  Func<double>* u = get_fn(fu, ru, eo);
  Func<double>* v = get_fn(fv, rv, eo);
  ExtData<scalar>* ext = init_ext_fns(mfs->ext, rv, eo);

  scalar res = mfs->fn(np, jwt, prev, u, v, e, ext);

  for (int i = 0; i < wf->neq; i++) {  prev[i]->free_fn(); delete prev[i]; }
  ext->free(); delete ext;
  return 0.5 * res;
}


// Actual evaluation of surface linear form (calculates integral)
scalar FeProblem::eval_form(WeakForm::VectorFormSurf *vfs, Tuple<Solution *> u_ext, PrecalcShapeset *fv, RefMap *rv, EdgePos* ep)
{
  // eval the form
  Quad2D* quad = fv->get_quad_2d();
  int eo = quad->get_edge_points(ep->edge);
  double3* pt = quad->get_points(eo);
  int np = quad->get_num_points(eo);

  // init geometry and jacobian*weights
  if (cache_e[eo] == NULL)
  {
    cache_e[eo] = init_geom_surf(rv, ep, eo);
    double3* tan = rv->get_tangent(ep->edge);
    cache_jwt[eo] = new double[np];
    for(int i = 0; i < np; i++)
      cache_jwt[eo][i] = pt[i][2] * tan[i][2];
  }
  Geom<double>* e = cache_e[eo];
  double* jwt = cache_jwt[eo];

  // function values and values of external functions
  AUTOLA_OR(Func<scalar>*, prev, wf->neq);
  for (int i = 0; i < wf->neq; i++) prev[i]  = init_fn(u_ext[i], rv, eo);
  Func<double>* v = get_fn(fv, rv, eo);
  ExtData<scalar>* ext = init_ext_fns(vfs->ext, rv, eo);

  scalar res = vfs->fn(np, jwt, prev, v, e, ext);

  for (int i = 0; i < wf->neq; i++) {  prev[i]->free_fn(); delete prev[i]; }
  ext->free(); delete ext;
  return 0.5 * res;
}

////////////////////////////////////////////////////////////////////////////////////////

template<typename Real, typename Scalar>
Scalar H1projection_biform(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->val[i] * v->val[i] + u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar H1projection_liform(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (ext->fn[0]->val[i] * v->val[i] + ext->fn[0]->dx[i] * v->dx[i] + ext->fn[0]->dy[i] * v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar L2projection_biform(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->val[i] * v->val[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar L2projection_liform(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (ext->fn[0]->val[i] * v->val[i]);
  return result;
}


Projection::Projection(int n, ...)
{
  num = n;

  va_list ap;
  va_start(ap, n);
  for (int i = 0; i < num; i++)
    slns[i] = va_arg(ap, MeshFunction*);
  for (int i = 0; i < num; i++)
    spaces[i] = va_arg(ap, Space*);
  for (int i = 0; i < num; i++)
    pss[i] = va_arg(ap, PrecalcShapeset*);
  va_end(ap);
}

Projection::~Projection()
{
  delete [] vec;
}

void Projection::set_solver(Solver* solver)
{
  this->solver = solver;
}

scalar* Projection::project()
{
  error("project() in FeProblem does not work currently. Employ DiscreteProblem::project_global()");
  /*
  WeakForm wf(num);
  for (int i = 0; i < num; i++)
  {
    wf.add_matrix_form(i, i, callback(L2projection_biform));
    wf.add_vector_form(i, callback(L2projection_liform), H2D_ANY, 1, slns[i]);
  }

  DiscreteProblem ps(&wf, solver, // SPACES MISSING HERE ));
  ps.assemble();
  Solution temp;
  ps.solve(0);
  const scalar* sln_vec = ps.get_solution_vector();
  int ndof = ps.get_num_dofs();
  vec = new scalar[ndof];
  memcpy(vec, sln_vec, ndof * sizeof(scalar));
  return vec;
  */
// For Visual Studio compiler it is necessary to return a value.
#ifdef _MSC_VER
	return new scalar(0.0);
#endif
}


