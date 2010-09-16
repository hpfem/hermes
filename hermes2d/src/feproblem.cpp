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
#include "space/space.h"
#include "precalc.h"
#include "matrix.h"
#include "solver/solver.h"
#include "solver/umfpack_solver.h"
#include "refmap.h"
#include "solution.h"
#include "config.h"

//  Solvers
#include "solver/amesos.h"
#include "solver/pardiso.h"
#include "solver/petsc.h"
#include "solver/mumps.h"
#include "solver/nox.h"


FeProblem::FeProblem(WeakForm* wf, Tuple<Space *> spaces, bool is_linear)
{
  this->is_linear = is_linear;

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

  // Create global enumeration of dof and fill the ndof variable
  this->ndof = assign_dofs(this->spaces);
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

void FeProblem::create(SparseMatrix* mat, Vector* rhs)
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
  if (rhs != NULL) rhs->alloc(ndof);

  // save space seq numbers and weakform seq number, so we can detect their changes
  for (int i = 0; i < wf->neq; i++)
    sp_seq[i] = spaces[i]->get_seq();
  wf_seq = wf->get_seq();

  struct_changed = true;
  have_matrix = true;
}

//// assembly //////////////////////////////////////////////////////////////////////////////////////

void FeProblem::assemble(Vector* init_vec, Matrix* mat_ext, Vector* rhs_ext, bool rhsonly)
{
  // Sanity checks.
  if (init_vec == NULL && this->is_linear == false) error("init_vec is NULL in FeProblem::assemble().");
  if (!have_spaces) error("You have to call FeProblem::set_spaces() before calling assemble().");
  for (int i=0; i<this->wf->neq; i++) {
    if (this->spaces[i] == NULL) error("A space is NULL in assemble().");
  }
 
  // Convert the coefficient vector 'init_vec' into solutions Tuple 'u_ext'.
  Tuple<Solution*> u_ext;
  for (int i = 0; i < this->wf->neq; i++) 
  {
    if (this->is_linear == false)
    {
      u_ext.push_back(new Solution(this->spaces[i]->get_mesh()));
      u_ext[i]->set_coeff_vector(this->spaces[i], init_vec);
    }
    else
      u_ext.push_back(NULL);
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
  wf->get_stages(spaces, NULL, stages, mat_ext == NULL);

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

        if (u_ext[j] != NULL) {
          u_ext[j]->set_active_element(e[i]);
          u_ext[j]->force_transform(pss[j]->get_transform(), pss[j]->get_ctm());
        }
      }
      marker = e0->marker;

      init_cache();
      //// assemble volume matrix forms //////////////////////////////////////
      if (mat_ext != NULL)
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
          scalar **local_stiffness_matrix = get_matrix_buffer(std::max(am->cnt, an->cnt));
          for (int i = 0; i < am->cnt; i++)
          {
            if (!tra && am->dof[i] < 0) continue;
            fv->set_active_shape(am->idx[i]);

            if (!sym) // unsymmetric block
            {
              for (int j = 0; j < an->cnt; j++) {
                fu->set_active_shape(an->idx[j]);
                if (an->dof[j] < 0) {
                  // Linear problems only: Subtracting Dirichlet lift contribution from the RHS:
                  if (rhs_ext != NULL && this->is_linear) {
                    scalar val = eval_form(mfv, u_ext, fu, fv, &refmap[n], &refmap[m]) * an->coef[j] * am->coef[i];
                    //dir_ext[am->dof[i]] += val;
                    rhs_ext->add(am->dof[i], -val);
                  } 
                }
                else if (rhsonly == false) {
                  scalar val = eval_form(mfv, u_ext, fu, fv, &refmap[n], &refmap[m]) * an->coef[j] * am->coef[i];
                  local_stiffness_matrix[i][j] = val;
                }
              }
	    /* OLD CODE
              for (int j = 0; j < an->cnt; j++) {
                fu->set_active_shape(an->idx[j]);
                bi = eval_form(mfv, u_ext, fu, fv, refmap+n, refmap+m) * an->coef[j] * am->coef[i];
                if (an->dof[j] >= 0) local_stiffness_matrix[i][j] = bi;
              }
	    */
            }
            else {// symmetric block
              for (int j = 0; j < an->cnt; j++) {
                if (j < i && an->dof[j] >= 0) continue;
                fu->set_active_shape(an->idx[j]);
                if (an->dof[j] < 0) {
                  // Linear problems only: Subtracting Dirichlet lift contribution from the RHS:
                  if (rhs_ext != NULL && this->is_linear) {
                    scalar val = eval_form(mfv, u_ext, fu, fv, &refmap[n], &refmap[m]) * an->coef[j] * am->coef[i];
                    //dir_ext[am->dof[i]] += val;
                    rhs_ext->add(am->dof[i], -val);
                  }
                } 
                else if (rhsonly == false) {
                  scalar val = eval_form(mfv, u_ext, fu, fv, &refmap[n], &refmap[m]) * an->coef[j] * am->coef[i];
                  local_stiffness_matrix[i][j] = local_stiffness_matrix[j][i] = val;
                }
	      /* OLD CODE
                bi = eval_form(mfv, u_ext, fu, fv, refmap+n, refmap+m) * an->coef[j] * am->coef[i];
                if (an->dof[j] >= 0) local_stiffness_matrix[i][j] = local_stiffness_matrix[j][i] = bi;
	      */
              }
            }
          }
          // insert the local stiffness matrix into the global one
          mat_ext->add(am->cnt, an->cnt, local_stiffness_matrix, am->dof, an->dof);

          // insert also the off-diagonal (anti-)symmetric block, if required
          if (tra)
          {
            if (mfv->sym < 0) chsgn(local_stiffness_matrix, am->cnt, an->cnt);
            transpose(local_stiffness_matrix, am->cnt, an->cnt);
            if (rhsonly == false) {
              mat_ext->add(am->cnt, an->cnt, local_stiffness_matrix, am->dof, an->dof);
            }
  	    /* OLD CODE
            mat_ext->add(am->cnt, an->cnt, local_stiffness_matrix, am->dof, an->dof);
	    */

            // Linear problems only: Subtracting Dirichlet lift contribution from the RHS:
            if (rhs_ext != NULL && this->is_linear) {
              for (int j = 0; j < am->cnt; j++) {
                if (am->dof[j] < 0) {
                  for (int i = 0; i < an->cnt; i++) {
                    if (an->dof[i] >= 0) {
                      //dir_ext[an->dof[i]] += local_stiffness_matrix[i][j];
                      rhs_ext->add(an->dof[i], -local_stiffness_matrix[i][j]);
                    }
                  }
                }
              }
            }
          }
        }
      }

      //// assemble volume vector forms ////////////////////////////////////////
      if (rhs_ext != NULL)
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
            scalar val = eval_form(vfv, u_ext, fv, refmap + m) * am->coef[i];
            rhs_ext->add(am->dof[i], val);
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
        if (mat_ext != NULL)
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
                  if (rhs_ext != NULL && this->is_linear) {
                    scalar val = eval_form(mfs, u_ext, fu, fv, &refmap[n], &refmap[m], &(ep[edge])) 
                                 * an->coef[j] * am->coef[i];
                    //dir_ext[am->dof[i]] += val;
                    rhs_ext->add(am->dof[i], -val);
                  }
                }
                else if (rhsonly == false) {
                  scalar val = eval_form(mfs, u_ext, fu, fv, &refmap[n], &refmap[m], &(ep[edge])) 
                               * an->coef[j] * am->coef[i];
                  local_stiffness_matrix[i][j] = val;
                } 
              }
            }
            if (rhsonly == false) {
              mat_ext->add(am->cnt, an->cnt, local_stiffness_matrix, an->dof, am->dof);
            }
	    /* OLD CODE
                bi = eval_form(mfs, u_ext, fu, fv, refmap+n, refmap+m, ep+edge) * an->coef[j] * am->coef[i];
                if (an->dof[j] >= 0) local_stiffness_matrix[i][j] = bi;
              }
            }
            mat_ext->add(am->cnt, an->cnt, local_stiffness_matrix, am->dof, an->dof);
            */
          }
        }
        // assemble surface linear forms /////////////////////////////////////
        if (rhs_ext != NULL)
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
              scalar val = eval_form(vfs, u_ext, fv, refmap+m, ep+edge) * am->coef[i];
              rhs_ext->add(am->dof[i], val);
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

  // Delete temporary solutions.
  for (int i = 0; i < wf->neq; i++) {
    if (u_ext[i] != NULL) {
      delete u_ext[i];
      u_ext[i] = NULL;
    }
  }
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
  for (unsigned i = 0; i < ext.size(); i++) {
    if (ext[i] != NULL) ext_fn[i] = init_fn(ext[i], rm, order);
    else ext_fn[i] = NULL;
  }
  ext_data->nf = ext.size();
  ext_data->fn = ext_fn;

  return ext_data;

}

// Initialize integration order on a given edge for external functions
ExtData<Ord>* FeProblem::init_ext_fns_ord(std::vector<MeshFunction *> &ext, int edge)
{
  ExtData<Ord>* fake_ext = new ExtData<Ord>;
  fake_ext->nf = ext.size();
  Func<Ord>** fake_ext_fn = new Func<Ord>*[fake_ext->nf];
  for (int i = 0; i < fake_ext->nf; i++)
    fake_ext_fn[i] = init_fn_ord(ext[i]->get_edge_fn_order(edge));
  fake_ext->fn = fake_ext_fn;
  
  return fake_ext;
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
  for (int i = 0; i < g_max_quad + 1 + 4 * g_max_quad + 4; i++)
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

/* OLD CODE
// Actual evaluation of volume matrix form (calculates integral)
scalar FeProblem::eval_form(WeakForm::MatrixFormVol *mfv, Tuple<Solution *> u_ext, 
                  PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv)
{
  // Determine the integration order.
  int inc = (fu->get_num_components() == 2) ? 1 : 0;

  // Order of solutions from the previous Newton iteration.
  AUTOLA_OR(Func<Ord>*, oi, wf->neq);
  //for (int i = 0; i < wf->neq; i++) oi[i] = init_fn_ord(sln[i]->get_fn_order() + inc);
  if (u_ext != Tuple<Solution *>()) {
    for (int i = 0; i < wf->neq; i++) {
      if (u_ext[i] != NULL) oi[i] = init_fn_ord(u_ext[i]->get_fn_order() + inc);
      else oi[i] = init_fn_ord(0);
    }
  }
  else {
    for (int i = 0; i < wf->neq; i++) oi[i] = init_fn_ord(0);
  }
  
  // Order of shape functions.
  Func<Ord>* ou = init_fn_ord(fu->get_fn_order() + inc);
  Func<Ord>* ov = init_fn_ord(fv->get_fn_order() + inc);

  // Order of additional external functions.
  ExtData<Ord>* fake_ext = init_ext_fns_ord(mfv->ext);

  // Order of geometric attributes (eg. for multiplication of a solution with coordinates, normals, etc.).
  double fake_wt = 1.0;
  Geom<Ord>* fake_e = init_geom_ord();

  // Total order of the matrix form.
  Ord o = mfv->ord(1, &fake_wt, oi, ou, ov, fake_e, fake_ext);

  // Increase due to reference map.
  int order = ru->get_inv_ref_order();
  order += o.get_order();
  limit_order_nowarn(order);

  // Clean up.
  for (int i = 0; i < wf->neq; i++) {  
    if (oi[i] != NULL) { oi[i]->free_ord(); delete oi[i]; }
  }
  if (ou != NULL) {
    ou->free_ord(); delete ou;
  }
  if (ov != NULL) {
    ov->free_ord(); delete ov;
  }
  if (fake_e != NULL) delete fake_e;
  if (fake_ext != NULL) {fake_ext->free_ord(); delete fake_ext;}

  // Evaluate the form using the quadrature of the just calculated order.
  Quad2D* quad = fu->get_quad_2d();
  double3* pt = quad->get_points(order);
  int np = quad->get_num_points(order);

  // Init geometry and jacobian*weights.
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

  // Values of the previous Newton iteration, shape functions and external functions in quadrature points.
  AUTOLA_OR(Func<scalar>*, prev, wf->neq);
  //for (int i = 0; i < wf->neq; i++) prev[i]  = init_fn(sln[i], rv, order);
  if (u_ext != Tuple<Solution *>()) {
    for (int i = 0; i < wf->neq; i++) {
      if (u_ext[i] != NULL) prev[i] = init_fn(u_ext[i], rv, order);
      else prev[i] = NULL;
    }
  }
  else {
    for (int i = 0; i < wf->neq; i++) prev[i] = NULL;
  }

  Func<double>* u = get_fn(fu, ru, order);
  Func<double>* v = get_fn(fv, rv, order);
  ExtData<scalar>* ext = init_ext_fns(mfv->ext, rv, order);

  scalar res = mfv->fn(np, jwt, prev, u, v, e, ext);

  // Clean up.
  for (int i = 0; i < wf->neq; i++) {  
    if (prev[i] != NULL) prev[i]->free_fn(); delete prev[i]; 
  }
  if (ext != NULL) {ext->free(); delete ext;}

  return res;
}
*/

// Actual evaluation of volume matrix form (calculates integral)
scalar FeProblem::eval_form(WeakForm::MatrixFormVol *mfv, Tuple<Solution *> u_ext, 
                        PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv)
{
  // Determine the integration order.
  int inc = (fu->get_num_components() == 2) ? 1 : 0;
  
  // Order of solutions from the previous Newton iteration.
  AUTOLA_OR(Func<Ord>*, oi, wf->neq);
  //for (int i = 0; i < wf->neq; i++) oi[i] = init_fn_ord(sln[i]->get_fn_order() + inc);
  if (u_ext != Tuple<Solution *>()) {
    for (int i = 0; i < wf->neq; i++) {
      if (u_ext[i] != NULL) oi[i] = init_fn_ord(u_ext[i]->get_fn_order() + inc);
      else oi[i] = init_fn_ord(0);
    }
  }
  else {
    for (int i = 0; i < wf->neq; i++) oi[i] = init_fn_ord(0);
  }
  
  // Order of shape functions.
  Func<Ord>* ou = init_fn_ord(fu->get_fn_order() + inc);
  Func<Ord>* ov = init_fn_ord(fv->get_fn_order() + inc);
  
  // Order of additional external functions.
  ExtData<Ord>* fake_ext = init_ext_fns_ord(mfv->ext);
  
  // Order of geometric attributes (eg. for multiplication of a solution with coordinates, normals, etc.).
  double fake_wt = 1.0;
  Geom<Ord>* fake_e = init_geom_ord();
  
  // Total order of the matrix form.
  Ord o = mfv->ord(1, &fake_wt, oi, ou, ov, fake_e, fake_ext);
  
  // Increase due to reference map.
  int order = ru->get_inv_ref_order();
  order += o.get_order();
  limit_order_nowarn(order);
  
  // Clean up.
  for (int i = 0; i < wf->neq; i++) {  
    if (oi[i] != NULL) { oi[i]->free_ord(); delete oi[i]; }
  }
  if (ou != NULL) {
    ou->free_ord(); delete ou;
  }
  if (ov != NULL) {
    ov->free_ord(); delete ov;
  }
  if (fake_e != NULL) delete fake_e;
  if (fake_ext != NULL) {fake_ext->free_ord(); delete fake_ext;}
  
  // Evaluate the form using the quadrature of the just calculated order.
  Quad2D* quad = fu->get_quad_2d();
  double3* pt = quad->get_points(order);
  int np = quad->get_num_points(order);

  // Init geometry and jacobian*weights.
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

  // Values of the previous Newton iteration, shape functions and external functions in quadrature points.
  AUTOLA_OR(Func<scalar>*, prev, wf->neq);
  //for (int i = 0; i < wf->neq; i++) prev[i]  = init_fn(sln[i], rv, order);
  if (u_ext != Tuple<Solution *>()) {
    for (int i = 0; i < wf->neq; i++) {
      if (u_ext[i] != NULL) prev[i] = init_fn(u_ext[i], rv, order);
      else prev[i] = NULL;
    }
  }
  else {
    for (int i = 0; i < wf->neq; i++) prev[i] = NULL;
  }

  Func<double>* u = get_fn(fu, ru, order);
  Func<double>* v = get_fn(fv, rv, order);
  ExtData<scalar>* ext = init_ext_fns(mfv->ext, rv, order);

  scalar res = mfv->fn(np, jwt, prev, u, v, e, ext);

  // Clean up.
  for (int i = 0; i < wf->neq; i++) {  
    if (prev[i] != NULL) prev[i]->free_fn(); delete prev[i]; 
  }
  if (ext != NULL) {ext->free(); delete ext;}
  return res;
}

/* OLD CODE
// Actual evaluation of volume linear form (calculates integral)
scalar FeProblem::eval_form(WeakForm::VectorFormVol *vfv, Tuple<Solution *> u_ext, PrecalcShapeset *fv, RefMap *rv)
{
  // Determine the integration order.
  int inc = (fv->get_num_components() == 2) ? 1 : 0;

  // Order of solutions from the previous Newton iteration.
  AUTOLA_OR(Func<Ord>*, oi, wf->neq);
  //for (int i = 0; i < wf->neq; i++) oi[i] = init_fn_ord(sln[i]->get_fn_order() + inc);
  if (u_ext != Tuple<Solution *>()) {
    for (int i = 0; i < wf->neq; i++) {
      if (u_ext[i] != NULL) oi[i] = init_fn_ord(u_ext[i]->get_fn_order() + inc);
      else oi[i] = init_fn_ord(0);
    }
  }
  else {
    for (int i = 0; i < wf->neq; i++) oi[i] = init_fn_ord(0);
  }
 
  // Order of the shape function.
  Func<Ord>* ov = init_fn_ord(fv->get_fn_order() + inc);

  // Order of additional external functions.
  ExtData<Ord>* fake_ext = init_ext_fns_ord(vfv->ext);

  // Order of geometric attributes (eg. for multiplication of a solution with coordinates, normals, etc.).
  double fake_wt = 1.0;
  Geom<Ord>* fake_e = init_geom_ord();

  // Total order of the vector form.
  Ord o = vfv->ord(1, &fake_wt, oi, ov, fake_e, fake_ext);

  // Increase due to reference map.
  int order = rv->get_inv_ref_order();
  order += o.get_order();
  limit_order_nowarn(order);

  // Clean up.
  for (int i = 0; i < wf->neq; i++) { 
    if (oi[i] != NULL) {
      oi[i]->free_ord(); delete oi[i]; 
    }
  }
  if (ov != NULL) {ov->free_ord(); delete ov;}
  if (fake_e != NULL) delete fake_e;
  if (fake_ext != NULL) {fake_ext->free_ord(); delete fake_ext;}

  // Evaluate the form using the quadrature of the just calculated order.
  Quad2D* quad = fv->get_quad_2d();
  double3* pt = quad->get_points(order);
  int np = quad->get_num_points(order);

  // Init geometry and jacobian*weights.
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

  // Values of the previous Newton iteration, shape functions and external functions in quadrature points.
  AUTOLA_OR(Func<scalar>*, prev, wf->neq);
  //for (int i = 0; i < wf->neq; i++) prev[i]  = init_fn(sln[i], rv, order);
  if (u_ext != Tuple<Solution *>()) {
    for (int i = 0; i < wf->neq; i++) {
      if (u_ext[i] != NULL) prev[i]  = init_fn(u_ext[i], rv, order);
      else prev[i] = NULL;
    }
  }
  else {
    for (int i = 0; i < wf->neq; i++) prev[i] = NULL;
  }

  Func<double>* v = get_fn(fv, rv, order);
  ExtData<scalar>* ext = init_ext_fns(vfv->ext, rv, order);

  scalar res = vfv->fn(np, jwt, prev, v, e, ext);

  // Clean up.
  for (int i = 0; i < wf->neq; i++) { 
    if (prev[i] != NULL) {
      prev[i]->free_fn(); delete prev[i]; 
    }
  }
  if (ext != NULL) {ext->free(); delete ext;}

  return res;
}
*/

// Actual evaluation of volume vector form (calculates integral)
scalar FeProblem::eval_form(WeakForm::VectorFormVol *vfv, Tuple<Solution *> u_ext, PrecalcShapeset *fv, RefMap *rv)
{
  // Determine the integration order.
  int inc = (fv->get_num_components() == 2) ? 1 : 0;
  
  // Order of solutions from the previous Newton iteration.
  AUTOLA_OR(Func<Ord>*, oi, wf->neq);
  //for (int i = 0; i < wf->neq; i++) oi[i] = init_fn_ord(sln[i]->get_fn_order() + inc);
  if (u_ext != Tuple<Solution *>()) {
    for (int i = 0; i < wf->neq; i++) {
      if (u_ext[i] != NULL) oi[i] = init_fn_ord(u_ext[i]->get_fn_order() + inc);
      else oi[i] = init_fn_ord(0);
    }
  }
  else {
    for (int i = 0; i < wf->neq; i++) oi[i] = init_fn_ord(0);
  }
  
  // Order of the shape function.
  Func<Ord>* ov = init_fn_ord(fv->get_fn_order() + inc);
  
  // Order of additional external functions.
  ExtData<Ord>* fake_ext = init_ext_fns_ord(vfv->ext);
  
  // Order of geometric attributes (eg. for multiplication of a solution with coordinates, normals, etc.).
  double fake_wt = 1.0;
  Geom<Ord>* fake_e = init_geom_ord();
  
  // Total order of the vector form.
  Ord o = vfv->ord(1, &fake_wt, oi, ov, fake_e, fake_ext);
  
  // Increase due to reference map.
  int order = rv->get_inv_ref_order();
  order += o.get_order();
  limit_order_nowarn(order);

  // Clean up.
  for (int i = 0; i < wf->neq; i++) { 
    if (oi[i] != NULL) {
      oi[i]->free_ord(); delete oi[i]; 
    }
  }
  if (ov != NULL) {ov->free_ord(); delete ov;}
  if (fake_e != NULL) delete fake_e;
  if (fake_ext != NULL) {fake_ext->free_ord(); delete fake_ext;}

  // Evaluate the form using the quadrature of the just calculated order.
  Quad2D* quad = fv->get_quad_2d();
  double3* pt = quad->get_points(order);
  int np = quad->get_num_points(order);

  // Init geometry and jacobian*weights.
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

  // Values of the previous Newton iteration, shape functions and external functions in quadrature points.
  AUTOLA_OR(Func<scalar>*, prev, wf->neq);
  //for (int i = 0; i < wf->neq; i++) prev[i]  = init_fn(sln[i], rv, order);
  if (u_ext != Tuple<Solution *>()) {
    for (int i = 0; i < wf->neq; i++) {
      if (u_ext[i] != NULL) prev[i]  = init_fn(u_ext[i], rv, order);
      else prev[i] = NULL;
    }
  }
  else {
    for (int i = 0; i < wf->neq; i++) prev[i] = NULL;
  }

  Func<double>* v = get_fn(fv, rv, order);
  ExtData<scalar>* ext = init_ext_fns(vfv->ext, rv, order);

  scalar res = vfv->fn(np, jwt, prev, v, e, ext);

  // Clean up.
  for (int i = 0; i < wf->neq; i++) { 
    if (prev[i] != NULL) {
      prev[i]->free_fn(); delete prev[i]; 
    }
  }
  if (ext != NULL) {ext->free(); delete ext;}
  
  return res;
}

/* OLD CODE
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
  for (int i = 0; i < wf->neq; i++) {
    if (u_ext[i] != NULL) prev[i] = init_fn(u_ext[i], rv, eo);
    else prev[i] = NULL;
  }
  Func<double>* u = get_fn(fu, ru, eo);
  Func<double>* v = get_fn(fv, rv, eo);
  ExtData<scalar>* ext = init_ext_fns(mfs->ext, rv, eo);

  scalar res = mfs->fn(np, jwt, prev, u, v, e, ext);

  if (this->is_linear == false) {
    for (int i = 0; i < wf->neq; i++) {  
      prev[i]->free_fn(); 
      delete prev[i]; 
    }
  }
  ext->free(); delete ext;
  return 0.5 * res;
}
*/

// Actual evaluation of surface matrix forms (calculates integral)
scalar FeProblem::eval_form(WeakForm::MatrixFormSurf *mfs, Tuple<Solution *> u_ext, 
                        PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv, EdgePos* ep)
{
  // Determine the integration order.
  int inc = (fu->get_num_components() == 2) ? 1 : 0;
  
  // Order of solutions from the previous Newton iteration.
  AUTOLA_OR(Func<Ord>*, oi, wf->neq);
  //for (int i = 0; i < wf->neq; i++) oi[i] = init_fn_ord(sln[i]->get_fn_order() + inc);
  if (u_ext != Tuple<Solution *>()) {
    for (int i = 0; i < wf->neq; i++) {
      if (u_ext[i] != NULL) oi[i] = init_fn_ord(u_ext[i]->get_edge_fn_order(ep->edge) + inc);
      else oi[i] = init_fn_ord(0);
    }
  }
  else {
    for (int i = 0; i < wf->neq; i++) oi[i] = init_fn_ord(0);
  }
  
  // Order of shape functions.
  Func<Ord>* ou = init_fn_ord(fu->get_edge_fn_order(ep->edge) + inc);
  Func<Ord>* ov = init_fn_ord(fv->get_edge_fn_order(ep->edge) + inc);
  
  // Order of additional external functions.
  ExtData<Ord>* fake_ext = init_ext_fns_ord(mfs->ext, ep->edge);
  
  // Order of geometric attributes (eg. for multiplication of a solution with coordinates, normals, etc.).
  double fake_wt = 1.0;
  Geom<Ord>* fake_e = init_geom_ord();
  
  // Total order of the matrix form.
  Ord o = mfs->ord(1, &fake_wt, oi, ou, ov, fake_e, fake_ext);
  
  // Increase due to reference map.
  int order = ru->get_inv_ref_order();
  
  order += o.get_order();
  limit_order_nowarn(order);
  
  // Clean up.
  for (int i = 0; i < wf->neq; i++) {  
    if (oi[i] != NULL) { oi[i]->free_ord(); delete oi[i]; }
  }
  if (ou != NULL) {
    ou->free_ord(); delete ou;
  }
  if (ov != NULL) {
    ov->free_ord(); delete ov;
  }
  if (fake_e != NULL) delete fake_e;
  if (fake_ext != NULL) {fake_ext->free_ord(); delete fake_ext;}
  
  // Evaluate the form using the quadrature of the just calculated order.
  Quad2D* quad = fu->get_quad_2d();
  
  int eo = quad->get_edge_points(ep->edge, order);
  double3* pt = quad->get_points(eo);
  int np = quad->get_num_points(eo);

  // Init geometry and jacobian*weights.
  if (cache_e[eo] == NULL)
  {
    cache_e[eo] = init_geom_surf(ru, ep, eo);
    double3* tan = ru->get_tangent(ep->edge, eo);
    cache_jwt[eo] = new double[np];
    for(int i = 0; i < np; i++)
      cache_jwt[eo][i] = pt[i][2] * tan[i][2];
  }
  Geom<double>* e = cache_e[eo];
  double* jwt = cache_jwt[eo];

  // Values of the previous Newton iteration, shape functions and external functions in quadrature points.
  AUTOLA_OR(Func<scalar>*, prev, wf->neq);
  //for (int i = 0; i < wf->neq; i++) prev[i]  = init_fn(sln[i], rv, eo);
  if (u_ext != Tuple<Solution *>()) {
    for (int i = 0; i < wf->neq; i++) {
      if (u_ext[i] != NULL) prev[i]  = init_fn(u_ext[i], rv, eo);
      else prev[i] = NULL;
    }
  }
  else {
    for (int i = 0; i < wf->neq; i++) prev[i] = NULL;
  }

  Func<double>* u = get_fn(fu, ru, eo);
  Func<double>* v = get_fn(fv, rv, eo);
  ExtData<scalar>* ext = init_ext_fns(mfs->ext, rv, eo);

  scalar res = mfs->fn(np, jwt, prev, u, v, e, ext);

  // Clean up.
  for (int i = 0; i < wf->neq; i++) { 
    if (prev[i] != NULL) {
      prev[i]->free_fn(); delete prev[i]; 
    }
  }
  if (ext != NULL) {ext->free(); delete ext;}
  
  return 0.5 * res; // Edges are parameterized from 0 to 1 while integration weights
                    // are defined in (-1, 1). Thus multiplying with 0.5 to correct
                    // the weights.
}

/* OLD CODE
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
  for (int i = 0; i < wf->neq; i++) {
    if (u_ext[i] != NULL) prev[i] = init_fn(u_ext[i], rv, eo);
    else prev[i] = NULL;
  }
  Func<double>* v = get_fn(fv, rv, eo);
  ExtData<scalar>* ext = init_ext_fns(vfs->ext, rv, eo);

  scalar res = vfs->fn(np, jwt, prev, v, e, ext);

  if (this->is_linear == false) {
    for (int i = 0; i < wf->neq; i++) {  
      prev[i]->free_fn(); 
      delete prev[i]; 
    }
  }
  ext->free(); delete ext;
  return 0.5 * res;
}
*/

// Actual evaluation of surface vector form (calculates integral)
scalar FeProblem::eval_form(WeakForm::VectorFormSurf *vfs, Tuple<Solution *> u_ext, 
                        PrecalcShapeset *fv, RefMap *rv, EdgePos* ep)
{
  // Determine the integration order.
  int inc = (fv->get_num_components() == 2) ? 1 : 0;
  
  // Order of solutions from the previous Newton iteration.
  AUTOLA_OR(Func<Ord>*, oi, wf->neq);
  //for (int i = 0; i < wf->neq; i++) oi[i] = init_fn_ord(sln[i]->get_fn_order() + inc);
  if (u_ext != Tuple<Solution *>()) {
    for (int i = 0; i < wf->neq; i++) {
      if (u_ext[i] != NULL) oi[i] = init_fn_ord(u_ext[i]->get_edge_fn_order(ep->edge) + inc);
      else oi[i] = init_fn_ord(0);
    }
  }
  else {
    for (int i = 0; i < wf->neq; i++) oi[i] = init_fn_ord(0);
  }
  
  // Order of the shape function.
  Func<Ord>* ov = init_fn_ord(fv->get_edge_fn_order(ep->edge) + inc);
  
  // Order of additional external functions.
  ExtData<Ord>* fake_ext = init_ext_fns_ord(vfs->ext, ep->edge);
  
  // Order of geometric attributes (eg. for multiplication of a solution with coordinates, normals, etc.).
  double fake_wt = 1.0;
  Geom<Ord>* fake_e = init_geom_ord();
  
  // Total order of the vector form.
  Ord o = vfs->ord(1, &fake_wt, oi, ov, fake_e, fake_ext);
  
  // Increase due to reference map.
  int order = rv->get_inv_ref_order();
  
  order += o.get_order();
  limit_order_nowarn(order);
  
  // Clean up.
  for (int i = 0; i < wf->neq; i++) { 
    if (oi[i] != NULL) {
      oi[i]->free_ord(); delete oi[i]; 
    }
  }
  if (ov != NULL) {ov->free_ord(); delete ov;}
  if (fake_e != NULL) delete fake_e;
  if (fake_ext != NULL) {fake_ext->free_ord(); delete fake_ext;}
  
  // Evaluate the form using the quadrature of the just calculated order.
  Quad2D* quad = fv->get_quad_2d();
  
  int eo = quad->get_edge_points(ep->edge, order);
  double3* pt = quad->get_points(eo);
  int np = quad->get_num_points(eo);

  // Init geometry and jacobian*weights.
  if (cache_e[eo] == NULL)
  {
    cache_e[eo] = init_geom_surf(rv, ep, eo);
    double3* tan = rv->get_tangent(ep->edge, eo);
    cache_jwt[eo] = new double[np];
    for(int i = 0; i < np; i++)
      cache_jwt[eo][i] = pt[i][2] * tan[i][2];
  }
  Geom<double>* e = cache_e[eo];
  double* jwt = cache_jwt[eo];

  // Values of the previous Newton iteration, shape functions and external functions in quadrature points.
  AUTOLA_OR(Func<scalar>*, prev, wf->neq);
  //for (int i = 0; i < wf->neq; i++) prev[i]  = init_fn(sln[i], rv, eo);
  if (u_ext != Tuple<Solution *>()) {
    for (int i = 0; i < wf->neq; i++) {
      if (u_ext[i] != NULL) prev[i]  = init_fn(u_ext[i], rv, eo);
      else prev[i] = NULL;
    }
  }
  else {
    for (int i = 0; i < wf->neq; i++) prev[i] = NULL;
  }

  Func<double>* v = get_fn(fv, rv, eo);
  ExtData<scalar>* ext = init_ext_fns(vfs->ext, rv, eo);

  scalar res = vfs->fn(np, jwt, prev, v, e, ext);

  for (int i = 0; i < wf->neq; i++) {  
    if (prev[i] != NULL) {prev[i]->free_fn(); delete prev[i]; }
  }
  if (ext != NULL) {ext->free(); delete ext;}
  
  // Clean up.
  return 0.5 * res; // Edges are parameterized from 0 to 1 while integration weights
                    // are defined in (-1, 1). Thus multiplying with 0.5 to correct
                    // the weights.
}


////////////////////////////////////////////////////////////////////////////////////////

template<typename Real, typename Scalar>
Scalar H1projection_biform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->val[i] * v->val[i] + u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar H1projection_liform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (ext->fn[0]->val[i] * v->val[i] + ext->fn[0]->dx[i] * v->dx[i] + ext->fn[0]->dy[i] * v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar L2projection_biform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->val[i] * v->val[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar L2projection_liform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (ext->fn[0]->val[i] * v->val[i]);
  return result;
}

// Hcurl projections
template<typename Real, typename Scalar>
Scalar Hcurlprojection_biform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                              Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++) {
    result += wt[i] * (u->curl[i] * conj(v->curl[i]));
    result += wt[i] * (u->val0[i] * conj(v->val0[i]) + u->val1[i] * conj(v->val1[i]));
  }
  return result;
}

template<typename Real, typename Scalar>
Scalar Hcurlprojection_liform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                              Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++) {
    result += wt[i] * (ext->fn[0]->curl[i] * conj(v->curl[i]));
    result += wt[i] * (ext->fn[0]->val0[i] * conj(v->val0[i]) + ext->fn[0]->val1[i] * conj(v->val1[i]));
  }

  return result;
}

double get_l2_norm(Vector* vec) 
{
  scalar val = 0;
  for (int i = 0; i < vec->length(); i++) {
    scalar inc = vec->get(i);
    val = val + inc*conj(inc);
  }
  return sqrt(std::abs(val));
}

// Basic Newton's method, takes a coefficient vector and returns a coefficient vector. 
// Assumes that the matrix and vector weak forms are Jacobian and residual forms. 
bool solve_newton(Tuple<Space *> spaces, WeakForm* wf, Vector* init_vec, 
                  Matrix* mat, Solver* solver, Vector* rhs, double newton_tol, 
                  int newton_max_iter, bool verbose) 
{
  int ndof = get_num_dofs(spaces);
  
  // sanity checks
  if (init_vec == NULL) error("init_vec == NULL in solve_newton().");
  if (rhs == NULL) error("rhs == NULL in solve_newton().");
  int n = spaces.size();
  if (spaces.size() != wf->neq) 
    error("The number of spaces in newton_solve() must match the number of equation in the PDE system.");
  for (int i=0; i < n; i++) {
    if (spaces[i] == NULL) error("spaces[%d] is NULL in solve_newton().", i);
  }

  // Initialize the discrete problem.
  bool is_linear = false;
  FeProblem fep(wf, spaces, is_linear);
  //info("ndof = %d", ndof);


  int it = 1;
  while (1)
  {
    // Assemble the Jacobian matrix and residual vector.
    fep.assemble(init_vec, mat, rhs);

    // Multiply the residual vector with -1 since the matrix 
    // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
    for (int i = 0; i < ndof; i++) rhs->set(i, -rhs->get(i));
    
    // Calculate the l2-norm of residual vector.
    double res_l2_norm; 
    res_l2_norm = get_l2_norm(rhs);
    if (verbose) info("---- Newton iter %d, ndof %d, res. l2 norm %g", 
                        it, get_num_dofs(spaces), res_l2_norm);

    // If l2 norm of the residual vector is in tolerance, quit.
    if (res_l2_norm < newton_tol|| it > newton_max_iter) break;

    // Solve the matrix problem.
    if (!solver->solve()) error ("Matrix solver failed.\n");

    // Add \deltaY^{n+1} to Y^n.
    for (int i = 0; i < ndof; i++) init_vec->add(i, rhs->get(i));  
    
    it++;
  };

  return (it <= newton_max_iter);
}

int get_num_dofs(Tuple<Space *> spaces)
{
  int ndof = 0;
  for (int i=0; i<spaces.size(); i++) {
    ndof += spaces[i]->get_num_dofs();
  }
  return ndof;
}

// Underlying function for global orthogonal projection.
// Not intended for the user. NOTE: the weak form here must be 
// a special projection weak form, which is different from 
// the weak form of the PDE. If you supply a weak form of the 
// PDE, the PDE will just be solved. 
void project_internal(Tuple<Space *> spaces, WeakForm *wf, scalar* target_vec)
{
  int n = spaces.size();

  // sanity checks
  if (n <= 0 || n > 10) error("Wrong number of projected functions in project_internal().");
  for (int i = 0; i < n; i++) if(spaces[i] == NULL) error("this->spaces[%d] == NULL in project_internal().", i);
  if (spaces.size() != n) error("Number of spaces must matchnumber of projected functions in project_internal().");

  // this is needed since spaces may have their DOFs enumerated only locally.
  int ndof = assign_dofs(spaces);

  // Initialize FeProblem.
  bool is_linear = true;
  FeProblem* fep = new FeProblem(wf, spaces, is_linear);

  // Initialize the UMFpack solver.
  // FIXME: More solvers need to be enabled here.
  UMFPackLinearSolver* solver = new UMFPackLinearSolver(fep);

  // Calculate the coefficient vector.
  bool solved = solver->solve();
  scalar* coeffs;
  if (solved) coeffs = solver->get_solution();

  if (target_vec != NULL) 
    for (int i=0; i<ndof; i++) target_vec[i] = coeffs[i];
}

// global orthogonal projection
void project_global(Tuple<Space *> spaces, Tuple<int> proj_norms, Tuple<MeshFunction*> source_meshfns, 
                    scalar* target_vec)
{
  int n = spaces.size();  

  // define temporary projection weak form
  WeakForm* proj_wf = new WeakForm(n);
  int found[100];
  for (int i = 0; i < 100; i++) found[i] = 0;
  for (int i = 0; i < n; i++) {
    int norm;
    if (proj_norms == Tuple<int>()) norm = 1;
    else norm = proj_norms[i];
    if (norm == 0) {
      found[i] = 1;
      proj_wf->add_matrix_form(i, i, L2projection_biform<double, scalar>, L2projection_biform<Ord, Ord>);
      proj_wf->add_vector_form(i, L2projection_liform<double, scalar>, L2projection_liform<Ord, Ord>,
                     H2D_ANY, source_meshfns[i]);
    }
    if (norm == 1) {
      found[i] = 1;
      proj_wf->add_matrix_form(i, i, H1projection_biform<double, scalar>, H1projection_biform<Ord, Ord>);
      proj_wf->add_vector_form(i, H1projection_liform<double, scalar>, H1projection_liform<Ord, Ord>,
                     H2D_ANY, source_meshfns[i]);
    }
    if (norm == 2) {
      found[i] = 1;
      proj_wf->add_matrix_form(i, i, Hcurlprojection_biform<double, scalar>, Hcurlprojection_biform<Ord, Ord>);
      proj_wf->add_vector_form(i, Hcurlprojection_liform<double, scalar>, Hcurlprojection_liform<Ord, Ord>,
                     H2D_ANY, source_meshfns[i]);
    }
  }
  for (int i=0; i < n; i++) {
    if (found[i] == 0) {
      warn("index of component: %d\n", i);
      error("Wrong projection norm in project_global().");
    }
  }

  project_internal(spaces, proj_wf, target_vec);
}

void project_global(Tuple<Space *> spaces, matrix_forms_tuple_t proj_biforms, 
                    vector_forms_tuple_t proj_liforms, Tuple<MeshFunction*> source_meshfns, 
                    scalar* target_vec)
{
  int n = spaces.size();
  matrix_forms_tuple_t::size_type n_biforms = proj_biforms.size();
  if (n_biforms == 0)
    error("Please use the simpler version of project_global with the argument Tuple<int> proj_norms if you do not provide your own projection norm.");
  if (n_biforms != proj_liforms.size())
    error("Mismatched numbers of projection forms in project_global().");
  if (n != n_biforms)
    error("Mismatched numbers of projected functions and projection forms in project_global().");

  // This is needed since spaces may have their DOFs enumerated only locally
  // when they come here.
  int ndof = assign_dofs(spaces);

  // Define projection weak form.
  WeakForm proj_wf(n);
  for (int i = 0; i < n; i++) {
    proj_wf.add_matrix_form(i, i, proj_biforms[i].first, proj_biforms[i].second);
    proj_wf.add_vector_form(i, proj_liforms[i].first, proj_liforms[i].second,
                    H2D_ANY, source_meshfns[i]);
  }

  project_internal(spaces, &proj_wf, target_vec);
}

/// Global orthogonal projection of one vector-valued ExactFunction.
void project_global(Space *space, ExactFunction2 source_fn, scalar* target_vec)
{
  int proj_norm = 2; // Hcurl
  Mesh *mesh = space->get_mesh();
  if (mesh == NULL) error("Mesh is NULL in project_global().");
  Solution source_sln;
  source_sln.set_exact(mesh, source_fn);
  project_global(space, proj_norm, (MeshFunction*)&source_sln, target_vec);
};

/// Projection-based interpolation of an exact function. This is faster than the
/// global projection since no global matrix problem is solved.
void project_local(Space *space, int proj_norm, ExactFunction source_fn, Mesh* mesh,
                   scalar* target_vec)
{
  /// TODO
}

// Solve a typical linear problem (without automatic adaptivity).
// Feel free to adjust this function for more advanced applications.
bool solve_linear(Tuple<Space *> spaces, WeakForm* wf, MatrixSolverType matrix_solver, 
                  Tuple<Solution *> solutions, scalar* coeff_vec) 
{
  // Initialize the linear problem.
  bool is_linear = true;
  FeProblem fep(wf, spaces, is_linear);
  int ndof = get_num_dofs(spaces);

  // Initialize matrix solver.
  Solver* solver;
  switch (matrix_solver) {
    case SOLVER_AMESOS: solver = new AmesosSolver("Amesos_Klu", &fep); info("Using Amesos"); break;
    case SOLVER_MUMPS: solver = new MumpsSolver(&fep); info("Using Mumps"); break;
    case SOLVER_NOX: solver = new NoxSolver(&fep); info("Using Nox"); break;
    case SOLVER_PARDISO: solver = new PardisoLinearSolver(&fep); info("Using Pardiso"); break;
    case SOLVER_PETSC: solver = new PetscLinearSolver(&fep); info("Using PETSc"); break;
    case SOLVER_UMFPACK: solver = new UMFPackLinearSolver(&fep); info("Using UMFPack"); break;
    default: error("Unknown matrix solver requested.");
  }

  // Solve the matrix problem.
  if (!solver->solve()) error ("Matrix solver failed.\n");

  // Extract solution vector.
  scalar* coeffs = solver->get_solution();

  // Convert coefficient vector into a Solution.
  for (int i=0; i < solutions.size(); i++) {
    solutions[i]->set_coeff_vector(spaces[i], coeffs);
  }
	
  delete solver;

  return true;
}

// Solve a typical linear problem using automatic adaptivity.
// Feel free to adjust this function for more advanced applications.
bool solve_linear_adapt(Tuple<Space *> spaces, WeakForm* wf, scalar* coeff_vec, 
                        MatrixSolverType matrix_solver, Tuple<int> proj_norms, 
                        Tuple<Solution *> slns, Tuple<Solution *> ref_slns, 
                        Tuple<WinGeom *> sln_win_geom, Tuple<WinGeom *> mesh_win_geom, 
                        Tuple<RefinementSelectors::Selector *> selectors, AdaptivityParamType* apt,
                        bool verbose, Tuple<ExactSolution *> exact_slns) 
{
  // sanity checks
  if (spaces.size() != selectors.size()) 
    error("There must be a refinement selector for each solution component in solve_linear_adapt().");
  if (spaces.size() != proj_norms.size()) 
    error("There must be a projection norm for each solution component in solve_linear_adapt().");

  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Adaptivity parameters.
  double err_stop = apt->err_stop; 
  int ndof_stop = apt->ndof_stop;
  double threshold = apt->threshold;
  int strategy = apt->strategy; 
  int mesh_regularity = apt->mesh_regularity;
  double to_be_processed = apt->to_be_processed;
  int total_error_flag = apt->total_error_flag;
  int elem_error_flag = apt->elem_error_flag;

  // Number of physical fields in the problem.
  int num_comps = spaces.size();

  // Number of degreeso of freedom 
  int ndof = get_num_dofs(spaces);

  // Number of exact solutions given.
  if (exact_slns.size() != 0 && exact_slns.size() != num_comps) 
    error("Number of exact solutions does not match number of equations.");
  bool is_exact_solution;
  if (exact_slns.size() == num_comps) is_exact_solution = true;
  else is_exact_solution = false;

  // Initialize views.
  ScalarView* s_view[H2D_MAX_COMPONENTS];
  VectorView* v_view[H2D_MAX_COMPONENTS];
  OrderView*  o_view[H2D_MAX_COMPONENTS];
  for (int i = 0; i < num_comps; i++) {
    char* title = (char*)malloc(100*sizeof(char));
    if (sln_win_geom != Tuple<WinGeom *>() && sln_win_geom[i] != NULL) {
      if (num_comps == 1) sprintf(title, "Solution", i); 
      else sprintf(title, "Solution[%d]", i); 
      switch (proj_norms[i]) {
        case H2D_L2_NORM:    s_view[i] = new ScalarView(title, sln_win_geom[i]);
                             s_view[i]->show_mesh(false);
                             v_view[i] = NULL;
                             break;
        case H2D_H1_NORM:    s_view[i] = new ScalarView(title, sln_win_geom[i]);
                             s_view[i]->show_mesh(false);
                             v_view[i] = NULL;
                             break;
        case H2D_HCURL_NORM: s_view[i] = NULL;
                             v_view[i] = new VectorView(title, sln_win_geom[i]);
                             break;
        case H2D_HDIV_NORM:  s_view[i] = NULL;
		             v_view[i] = new VectorView(title, sln_win_geom[i]);
                             break;
      default: error("Unknown norm in solve_linear_adapt().");
      }
    }
    else {
      s_view[i] = NULL;
      v_view[i] = NULL;
    }
    if (mesh_win_geom != Tuple<WinGeom *>() && mesh_win_geom[i] != NULL) {
      if (num_comps == 1) sprintf(title, "Mesh", i); 
      else sprintf(title, "Mesh[%d]", i); 
      o_view[i] = new OrderView(title, mesh_win_geom[i]);
    }
    else o_view[i] = NULL;
  }

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est, graph_dof_exact, graph_cpu_exact;

  // FIXME: Conversion from Tuple<Solution *> to Tuple<MeshFunction *>
  // temporary so that project_global() below compiles. 
  Tuple<MeshFunction *> ref_slns_mf;
  for (int i = 0; i < num_comps; i++) {
    ref_slns_mf.push_back((MeshFunction*)ref_slns[i]);
  }

  int as = 1; bool done = false;
  do
  {
    if (verbose) {
      info("---- Adaptivity step %d:", as);
      info("Solving on reference mesh.");
    }

    // Construct globally refined reference mesh(es)
    // and setup reference space(s).
    Tuple<Space *> ref_spaces;
    Tuple<Mesh *> ref_meshes;
    for (int i = 0; i < num_comps; i++) {
      ref_meshes.push_back(new Mesh());
      Mesh *ref_mesh = ref_meshes.back();
      ref_mesh->copy(spaces[i]->get_mesh());
      ref_mesh->refine_all_elements();
      ref_spaces.push_back(spaces[i]->dup(ref_mesh));
      int order_increase = 1;
      ref_spaces[i]->copy_orders(spaces[i], order_increase);
    }

    // Solve the reference problem.
    // The NULL pointer means that we do not want the resulting coefficient vector.
    solve_linear(ref_spaces, wf, matrix_solver, ref_slns, NULL);

    // Project the reference solution on the coarse mesh.
    if (verbose) info("Projecting reference solution on coarse mesh.");
    scalar* coeff_vec = new scalar[ndof];
    project_global(spaces, proj_norms, ref_slns_mf, coeff_vec); 

    // Set projected Solutions on the coarse mesh.
    for (int i = 0; i < num_comps; i++)
      slns[i]->set_coeff_vector(spaces[i], coeff_vec);
    delete [] coeff_vec;

    // Time measurement.
    cpu_time.tick();

    // View the coarse mesh solution(s).
    for (int i = 0; i < num_comps; i++) {
      if (proj_norms[i] == H2D_H1_NORM || proj_norms[i] == H2D_L2_NORM) {
        if (s_view[i] != NULL) s_view[i]->show(slns[i]);
      }
      if (proj_norms[i] == H2D_HCURL_NORM || proj_norms[i] == H2D_HDIV_NORM) {
        if (v_view[i] != NULL) v_view[i]->show(slns[i]);
      }
      if (o_view[i] != NULL) o_view[i]->show(spaces[i]);
    }

    // Skip visualization time.
    cpu_time.tick(HERMES_SKIP);

    // Calculate element errors.
    if (verbose) info("Calculating error (est).");
    Adapt hp(spaces, proj_norms);
    // Pass special error forms if any.
    for (int k = 0; k < apt->error_form_val.size(); k++) {
      hp.set_error_form(apt->error_form_i[k], apt->error_form_j[k], 
                        apt->error_form_val[k], apt->error_form_ord[k]);
    }
    hp.set_solutions(slns, ref_slns);
    // Below, apt->total_error_flag = H2D_TOTAL_ERROR_REL or H2D_TOTAL_ERROR_ABS
    // and apt->elem_error_flag = H2D_ELEMENT_ERROR_REL or H2D_ELEMENT_ERROR_ABS
    hp.calc_elem_errors(total_error_flag | elem_error_flag);
 
    // Calculate error estimate for each solution component. Note, these can 
    // have different norms, such as L2, H1, Hcurl and Hdiv. 
    double err_est_abs[H2D_MAX_COMPONENTS];
    double norm_est[H2D_MAX_COMPONENTS];
    double err_est_abs_total = 0;
    double norm_est_total = 0;
    double err_est_rel_total;
    for (int i = 0; i < num_comps; i++) {
      err_est_abs[i] = calc_abs_error(slns[i], ref_slns[i], proj_norms[i]);
      norm_est[i] = calc_norm(ref_slns[i], proj_norms[i]);
      err_est_abs_total += err_est_abs[i] * err_est_abs[i];
      norm_est_total += norm_est[i] * norm_est[i];
    }
    err_est_abs_total = sqrt(err_est_abs_total);
    norm_est_total = sqrt(norm_est_total);
    err_est_rel_total = err_est_abs_total / norm_est_total * 100.;

    // Calculate exact error for each solution component.   
    double err_exact_abs[H2D_MAX_COMPONENTS];
    double norm_exact[H2D_MAX_COMPONENTS];
    double err_exact_abs_total = 0;
    double norm_exact_total = 0;
    double err_exact_rel_total;
    if (is_exact_solution == true) {
      for (int i = 0; i < num_comps; i++) {
        err_exact_abs[i] = calc_abs_error(slns[i], exact_slns[i], proj_norms[i]);
        norm_exact[i] = calc_norm(exact_slns[i], proj_norms[i]);
        err_exact_abs_total += err_exact_abs[i] * err_exact_abs[i];
        norm_exact_total += norm_exact[i] * norm_exact[i];
      }
      err_exact_abs_total = sqrt(err_exact_abs_total);
      norm_exact_total = sqrt(norm_exact_total);
      err_exact_rel_total = err_exact_abs_total / norm_exact_total * 100.;
    }

    // Report results.
    if (verbose) {
      if (num_comps == 1) {
        info("ndof: %d, ref_ndof: %d, err_est_rel_total: %g%%", 
             get_num_dofs(spaces), get_num_dofs(ref_spaces), err_est_rel_total);
        if (is_exact_solution == true) info("err_exact_rel_total: %g%%", err_exact_rel_total);
      }
      else {
        for (int i = 0; i < num_comps; i++) {
          info("ndof[%d]: %d, ref_ndof[%d]: %d, err_est_rel[%d]: %g%%", 
               i, spaces[i]->get_num_dofs(), i, ref_spaces[i]->get_num_dofs(),
               i, err_est_abs[i]/norm_est[i]*100);
          if (is_exact_solution == true) info("err_exact_rel[%d]: %g%%", 
                                              i, err_exact_abs[i]/norm_exact[i]*100);
        }
        info("ndof: %d, ref_ndof: %d, err_est_rel_total: %g%%", 
             get_num_dofs(spaces), get_num_dofs(ref_spaces), err_est_rel_total);
      }
    }

    // Add entry to DOF and CPU convergence graphs.
    graph_dof_est.add_values(get_num_dofs(spaces), err_est_rel_total);
    graph_dof_est.save("conv_dof_est.dat");
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est_rel_total);
    graph_cpu_est.save("conv_cpu_est.dat");
    if (is_exact_solution == true) {
      graph_dof_exact.add_values(get_num_dofs(spaces), err_exact_rel_total);
      graph_dof_exact.save("conv_dof_exact.dat");
      graph_cpu_exact.add_values(cpu_time.accumulated(), err_exact_rel_total);
      graph_cpu_exact.save("conv_cpu_exact.dat");
    }

    // If err_est too large, adapt the mesh.
    if (err_est_rel_total < err_stop) done = true;
    else {
      if (verbose) info("Adapting the coarse mesh.");
      done = hp.adapt(selectors, threshold, strategy, mesh_regularity, to_be_processed);

      if (get_num_dofs(spaces) >= ndof_stop) done = true;
    }

    // Free reference meshes and spaces.
    for (int i = 0; i < num_comps; i++) {
      ref_spaces[i]->free(); // This does not free the associated mesh, we must do it separately.
      ref_meshes[i]->free();
    }

    as++;
  }
  while (done == false);

  if (verbose) info("Total running time: %g s", cpu_time.accumulated());
	
  return true;
}






