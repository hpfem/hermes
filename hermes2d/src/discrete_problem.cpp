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

#include "h2d_common.h"
#include "limit_order.h"
#include "discrete_problem.h"
#include "traverse.h"
#include "space/space.h"
#include "precalc.h"
#include "../../hermes_common/matrix.h"
#include "refmap.h"
#include "solution.h"
#include "config.h"

DiscreteProblem::DiscreteProblem(WeakForm* wf, Tuple<Space *> spaces, bool is_linear)
{
  _F_
  // sanity checks
  int n = spaces.size();
  if (n != wf->neq) error("Bad number of spaces in DiscreteProblem.");

  this->wf = wf;
  this->spaces = spaces;
  this->is_linear = is_linear;

  sp_seq = new int[wf->neq];
  memset(sp_seq, -1, sizeof(int) * wf->neq);
  wf_seq = -1;

  // This is different from H3D.
  pss = new PrecalcShapeset*[wf->neq];
  num_user_pss = 0;

  matrix_buffer = NULL;
  matrix_buffer_dim = 0;

  values_changed = true;
  struct_changed = true;

  have_matrix = false;

  this->spaces = Tuple<Space *>();
  for (int i = 0; i < wf->neq; i++) this->spaces.push_back(spaces[i]);
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
    this->pss[i] = p;
    this->num_user_pss++;
  }  

  // Create global enumeration of dof and fill the ndof variable
  this->ndof = Space::assign_dofs(this->spaces);
}

DiscreteProblem::~DiscreteProblem()
{
  _F_
  free();
  if (sp_seq != NULL) delete [] sp_seq;
  if (pss != NULL) delete [] pss;
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
  ndof = 0;
  for (int i = 0; i < wf->neq; i++)
    ndof += Space::get_num_dofs(spaces[i]);
  return ndof;
}

scalar** DiscreteProblem::get_matrix_buffer(int n)
{
  _F_
  if (n <= matrix_buffer_dim) return matrix_buffer;
  if (matrix_buffer != NULL) delete [] matrix_buffer;
  matrix_buffer_dim = n;
  return (matrix_buffer = new_matrix<scalar>(n, n));
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
  {
    if (spaces[i]->get_seq() != sp_seq[i])
    { 
      up_to_date = false; 
      break; 
    }
  }
  
  if (wf->get_seq() != wf_seq)
    up_to_date = false;

  return up_to_date;
}

//// matrix creation ///////////////////////////////////////////////////////////////////////////////

// This functions is identical in H2D and H3D.
void DiscreteProblem::create(SparseMatrix* mat, Vector* rhs, bool rhsonly)
{
  _F_

  if (is_up_to_date())
  {
    if (!rhsonly && mat != NULL) 
    {
      verbose("Reusing matrix sparse structure.");
      mat->zero();
    }
    if (rhs != NULL) rhs->zero();
    return;
  }
  
  int ndof = get_num_dofs();
  
  if (mat != NULL)  // mat may be NULL when assembling the rhs for NOX
  {
    // spaces have changed: create the matrix from scratch
    mat->free();
    mat->prealloc(ndof);

    AUTOLA_CL(AsmList, al, wf->neq);
    AUTOLA_OR(Mesh*, meshes, wf->neq);
    bool **blocks = wf->get_blocks();

    // Init multi-mesh traversal.
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
      {
        // TODO: do not get the assembly list again if the element was not changed
        if (e[i] != NULL) spaces[i]->get_element_assembly_list(e[i], &(al[i]));
      }

      // go through all equation-blocks of the local stiffness matrix
      for (int m = 0; m < wf->neq; m++)
      {
        for (int n = 0; n < wf->neq; n++)
        {
          if (blocks[m][n] && e[m] != NULL && e[n] != NULL) 
          {
            AsmList *am = &(al[m]);
            AsmList *an = &(al[n]);

            // pretend assembling of the element stiffness matrix
            // register nonzero elements
            for (int i = 0; i < am->cnt; i++)
              if (am->dof[i] >= 0)
                for (int j = 0; j < an->cnt; j++)
                  if (an->dof[j] >= 0)
                    mat->pre_add_ij(am->dof[i], an->dof[j]);
          }
        }
      }
    }

    trav.finish();
    delete [] blocks;

    mat->alloc();
  }
  
  // WARNING: unlike Matrix::alloc(), Vector::alloc(ndof) frees the memory occupied 
  // by previous vector before allocating
  if (rhs != NULL) rhs->alloc(ndof);    

  // save space seq numbers and weakform seq number, so we can detect their changes
  for (int i = 0; i < wf->neq; i++)
    sp_seq[i] = spaces[i]->get_seq();
  
  wf_seq = wf->get_seq();

  struct_changed = true;
  have_matrix = true;
}

//// assembly //////////////////////////////////////////////////////////////////////////////////////

// Light version for linear problems.
void DiscreteProblem::assemble(SparseMatrix* mat, Vector* rhs, bool rhsonly) 
{
  _F_
  assemble(NULL, mat, rhs, rhsonly);
}

// General assembling function for nonlinear problem. For linear problems use the 
// light version above.
void DiscreteProblem::assemble(scalar* coeff_vec, SparseMatrix* mat, Vector* rhs, bool rhsonly)
{
  /* BEGIN IDENTICAL CODE WITH H3D */

	_F_
  // Sanity checks.
  if (coeff_vec == NULL && this->is_linear == false) error("coeff_vec is NULL in DiscreteProblem::assemble().");
  if (!have_spaces) error("You have to call DiscreteProblem::set_spaces() before calling assemble().");
  for (int i=0; i<this->wf->neq; i++)
  {
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
      Solution::vector_to_solution(coeff_vec, this->spaces[i], u_ext[i]);
    }
    else
      u_ext.push_back(NULL);
  }
 
  /* END IDENTICAL CODE WITH H3D */

  bool bnd[4];			    // FIXME: magic number - maximal possible number of element surfaces
  SurfPos surf_pos[4];
  AUTOLA_CL(AsmList, al, wf->neq);
  AUTOLA_OR(bool, nat, wf->neq);
  AUTOLA_OR(bool, isempty, wf->neq);
  AsmList *am, *an;
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

  // initialize matrix buffer
  matrix_buffer = NULL;
  matrix_buffer_dim = 0;
  if (mat != NULL) get_matrix_buffer(9);

  // obtain a list of assembling stages
  std::vector<WeakForm::Stage> stages;
  wf->get_stages(spaces, this->is_linear ? NULL : u_ext, stages, rhsonly);

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
    while ((e = trav.get_next_state(bnd, surf_pos)) != NULL)
    {
      // find a non-NULL e[i]
      Element* e0;
      for (unsigned int i = 0; i < s->idx.size(); i++)
        if ((e0 = e[i]) != NULL) break;
      if (e0 == NULL) continue;

      // set maximum integration order for use in integrals, see limit_order()
      update_limit_table(e0->get_mode());

      // Obtain assembly lists for the element at all spaces of the stage, set appropriate mode for each pss.
      // NOTE: Active elements and transformations for external functions (including the solutions from previous
      // Newton's iteration) as well as basis functions (master PrecalcShapesets) have already been set in 
      // trav.get_next_state(...).
      memset(isempty, 0, sizeof(bool) * wf->neq);
      for (unsigned int i = 0; i < s->idx.size(); i++)
      {
        int j = s->idx[i];
        if (e[i] == NULL) 
        { 
          isempty[j] = true; 
          continue; 
        }

        // TODO: do not obtain again if the element was not changed.
        spaces[j]->get_element_assembly_list(e[i], &(al[j]));

        // This is different in H3D (PrecalcShapeset is not used)
        spss[j]->set_active_element(e[i]);
        spss[j]->set_master_transform();

        // This is different in H2D (PrecalcShapeset is not used).
        refmap[j].set_active_element(e[i]);
        refmap[j].force_transform(pss[j]->get_transform(), pss[j]->get_ctm());
      }
      int marker = e0->marker;

      init_cache();     // This is different in H2D.

      //// assemble volume matrix forms //////////////////////////////////////
      if (mat != NULL)
      {
        for (unsigned ww = 0; ww < s->mfvol.size(); ww++)
        {
          WeakForm::MatrixFormVol* mfv = s->mfvol[ww];
          if (isempty[mfv->i] || isempty[mfv->j]) continue;
          if (mfv->area != HERMES_ANY && !wf->is_in_area(marker, mfv->area)) continue;
          int m = mfv->i;  
          int n = mfv->j;  
          fu = pss[n]; 
          fv = spss[m];  
          am = &al[m];  
          an = &al[n];
          bool tra = (m != n) && (mfv->sym != 0);
          bool sym = (m == n) && (mfv->sym == 1);

	  /* BEGIN IDENTICAL CODE WITH H3D */

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
                    scalar val = eval_form(mfv, u_ext, fu, fv, &(refmap[n]),
                            &(refmap[m])) * an->coef[j] * am->coef[i];
                    rhs->add(am->dof[i], -val);
                  } 
                }
                else if (rhsonly == false) 
                {
                  scalar val = eval_form(mfv, u_ext, fu, fv, &(refmap[n]),
                          &(refmap[m])) * an->coef[j] * am->coef[i];
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
                    scalar val = eval_form(mfv, u_ext, fu, fv, &(refmap[n]),
                            &(refmap[m])) * an->coef[j] * am->coef[i];
                    rhs->add(am->dof[i], -val);
                  }
                } 
                else if (rhsonly == false) 
                {
                  scalar val = eval_form(mfv, u_ext, fu, fv, &(refmap[n]),
                          &(refmap[m])) * an->coef[j] * am->coef[i];
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
            if (mfv->sym < 0) 
              chsgn(local_stiffness_matrix, am->cnt, an->cnt);
            
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

      /* END IDENTICAL CODE WITH H3D
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
          fv = spss[m];    // H3D uses fv = test_fn + m;
          am = &(al[m]);

          for (int i = 0; i < am->cnt; i++)
          {
            if (am->dof[i] < 0) continue;
            fv->set_active_shape(am->idx[i]);
            scalar val = eval_form(vfv, u_ext, fv, &(refmap[m])) * am->coef[i];
            rhs->add(am->dof[i], val);
          }
        }
      }

      // assemble surface integrals now: loop through surfaces of the element
      for (unsigned int isurf = 0; isurf < e0->get_num_surf(); isurf++)
      {
        // H3D is freeing a fn_cache at this point

        if (!bnd[isurf]) continue;
        
        int marker = surf_pos[isurf].marker;

        // obtain the list of shape functions which are nonzero on this surface
        for (unsigned int i = 0; i < s->idx.size(); i++) 
        {
          if (e[i] == NULL) continue;
          int j = s->idx[i];
          if ((nat[j] = (spaces[j]->bc_type_callback(marker) == BC_NATURAL)))
            spaces[j]->get_boundary_assembly_list(e[i], isurf, &(al[j]));
        }

        // assemble surface matrix forms ///////////////////////////////////
        if (mat != NULL)
        {
          for (unsigned int ww = 0; ww < s->mfsurf.size(); ww++)
          {
            WeakForm::MatrixFormSurf* mfs = s->mfsurf[ww];
            if (isempty[mfs->i] || isempty[mfs->j]) continue;
            if (mfs->area != HERMES_ANY && !wf->is_in_area(marker, mfs->area)) continue;
            int m = mfs->i;  
            int n = mfs->j;  
            fu = pss[n];      // This is different in H3D.
            fv = spss[m];     // This is different in H3D.
            am = &(al[m]);
            an = &(al[n]);

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
                if (an->dof[j] < 0) 
                {
                  // Linear problems only: Subtracting Dirichlet lift contribution from the RHS:
                  if (rhs != NULL && this->is_linear) 
                  {
                    scalar val = eval_form(mfs, u_ext, fu, fv, &(refmap[n]),
                            &(refmap[m]), surf_pos + isurf) * an->coef[j] * am->coef[i];
                    rhs->add(am->dof[i], -val);
                  }
                }
                else if (rhsonly == false) 
                {
                  scalar val = eval_form(mfs, u_ext, fu, fv, &(refmap[n]),
                          &(refmap[m]), surf_pos + isurf) * an->coef[j] * am->coef[i];
                  local_stiffness_matrix[i][j] = val;
                } 
              }
            }
            if (rhsonly == false) 
              mat->add(am->cnt, an->cnt, local_stiffness_matrix, am->dof, an->dof);
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
            fv = spss[m];        // This is different from H3D.  
            am = &(al[m]);

            if (!nat[m]) continue;
            surf_pos[isurf].base = trav.get_base();
            surf_pos[isurf].space_v = spaces[m];

            for (int i = 0; i < am->cnt; i++)
            {
              if (am->dof[i] < 0) continue;
              fv->set_active_shape(am->idx[i]);
              scalar val = eval_form(vfs, u_ext, fv, &(refmap[m]), surf_pos + isurf) * am->coef[i];
              rhs->add(am->dof[i], val);
            }
          }
        }
      }

      delete_cache();   // This is different in H3D.
    }

    if (mat != NULL) mat->finish();
    if (rhs != NULL) rhs->finish();
    trav.finish();
  }

  for (int i = 0; i < wf->neq; i++) delete spss[i];  // This is different from H3D.

  // Cleaning up.
  if (matrix_buffer != NULL) delete [] matrix_buffer;
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

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Initialize integration order for external functions
ExtData<Ord>* DiscreteProblem::init_ext_fns_ord(std::vector<MeshFunction *> &ext)
{
  _F_
  ExtData<Ord>* fake_ext = new ExtData<Ord>;
  fake_ext->nf = ext.size();
  Func<Ord>** fake_ext_fn = new Func<Ord>*[fake_ext->nf];
  for (int i = 0; i < fake_ext->nf; i++)
    fake_ext_fn[i] = init_fn_ord(ext[i]->get_fn_order());
  fake_ext->fn = fake_ext_fn;

  return fake_ext;
}

// Initialize external functions (obtain values, derivatives,...)
ExtData<scalar>* DiscreteProblem::init_ext_fns(std::vector<MeshFunction *> &ext, RefMap *rm, const int order)
{
  _F_
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
ExtData<Ord>* DiscreteProblem::init_ext_fns_ord(std::vector<MeshFunction *> &ext, int edge)
{
  _F_
  ExtData<Ord>* fake_ext = new ExtData<Ord>;
  fake_ext->nf = ext.size();
  Func<Ord>** fake_ext_fn = new Func<Ord>*[fake_ext->nf];
  for (int i = 0; i < fake_ext->nf; i++)
    fake_ext_fn[i] = init_fn_ord(ext[i]->get_edge_fn_order(edge));
  fake_ext->fn = fake_ext_fn;
  
  return fake_ext;
}

// Initialize shape function values and derivatives (fill in the cache)
Func<double>* DiscreteProblem::get_fn(PrecalcShapeset *fu, RefMap *rm, const int order)
{
  _F_
  PrecalcShapeset::Key key(256 - fu->get_active_shape(), order, fu->get_transform(), fu->get_shapeset()->get_id());
  if (cache_fn[key] == NULL)
    cache_fn[key] = init_fn(fu, rm, order);

  return cache_fn[key];
}

// Caching transformed values
void DiscreteProblem::init_cache()
{
  _F_
  for (int i = 0; i < g_max_quad + 1 + 4 * g_max_quad + 4; i++)
  {
    cache_e[i] = NULL;
    cache_jwt[i] = NULL;
  }
}

void DiscreteProblem::delete_cache()
{
  _F_
  for (int i = 0; i < g_max_quad + 1 + 4 * g_max_quad + 4; i++)
  {
    if (cache_e[i] != NULL)
    {
      cache_e[i]->free(); delete cache_e[i];
      delete [] cache_jwt[i];
    }
  }
  for (std::map<PrecalcShapeset::Key, Func<double>*, PrecalcShapeset::Compare>::const_iterator it = cache_fn.begin(); it != cache_fn.end(); it++)
  {
    (it->second)->free_fn(); delete (it->second);
  }
  cache_fn.clear();
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Actual evaluation of volume matrix form (calculates integral)
scalar DiscreteProblem::eval_form(WeakForm::MatrixFormVol *mfv, Tuple<Solution *> u_ext, 
                        PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv)
{
  _F_
  // Determine the integration order.
  int inc = (fu->get_num_components() == 2) ? 1 : 0;
  
  // Order of solutions from the previous Newton iteration.
  AUTOLA_OR(Func<Ord>*, oi, wf->neq);
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

// Actual evaluation of volume vector form (calculates integral)
scalar DiscreteProblem::eval_form(WeakForm::VectorFormVol *vfv, Tuple<Solution *> u_ext, PrecalcShapeset *fv, RefMap *rv)
{
  _F_
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

// Actual evaluation of surface matrix forms (calculates integral)
scalar DiscreteProblem::eval_form(WeakForm::MatrixFormSurf *mfs, Tuple<Solution *> u_ext, 
                        PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv, SurfPos* surf_pos)
{
  _F_
  // Determine the integration order.
  int inc = (fu->get_num_components() == 2) ? 1 : 0;
  
  // Order of solutions from the previous Newton iteration.
  AUTOLA_OR(Func<Ord>*, oi, wf->neq);
  //for (int i = 0; i < wf->neq; i++) oi[i] = init_fn_ord(sln[i]->get_fn_order() + inc);
  if (u_ext != Tuple<Solution *>()) {
    for (int i = 0; i < wf->neq; i++) {
      if (u_ext[i] != NULL) oi[i] = init_fn_ord(u_ext[i]->get_edge_fn_order(surf_pos->surf_num) + inc);
      else oi[i] = init_fn_ord(0);
    }
  }
  else {
    for (int i = 0; i < wf->neq; i++) oi[i] = init_fn_ord(0);
  }
  
  // Order of shape functions.
  Func<Ord>* ou = init_fn_ord(fu->get_edge_fn_order(surf_pos->surf_num) + inc);
  Func<Ord>* ov = init_fn_ord(fv->get_edge_fn_order(surf_pos->surf_num) + inc);
  
  // Order of additional external functions.
  ExtData<Ord>* fake_ext = init_ext_fns_ord(mfs->ext, surf_pos->surf_num);
  
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
  
  int eo = quad->get_edge_points(surf_pos->surf_num, order);
  double3* pt = quad->get_points(eo);
  int np = quad->get_num_points(eo);

  // Init geometry and jacobian*weights.
  if (cache_e[eo] == NULL)
  {
    cache_e[eo] = init_geom_surf(ru, surf_pos, eo);
    double3* tan = ru->get_tangent(surf_pos->surf_num, eo);
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

// Actual evaluation of surface vector form (calculates integral)
scalar DiscreteProblem::eval_form(WeakForm::VectorFormSurf *vfs, Tuple<Solution *> u_ext, 
                        PrecalcShapeset *fv, RefMap *rv, SurfPos* surf_pos)
{
  _F_
  // Determine the integration order.
  int inc = (fv->get_num_components() == 2) ? 1 : 0;
  
  // Order of solutions from the previous Newton iteration.
  AUTOLA_OR(Func<Ord>*, oi, wf->neq);
  //for (int i = 0; i < wf->neq; i++) oi[i] = init_fn_ord(sln[i]->get_fn_order() + inc);
  if (u_ext != Tuple<Solution *>()) {
    for (int i = 0; i < wf->neq; i++) {
      if (u_ext[i] != NULL) oi[i] = init_fn_ord(u_ext[i]->get_edge_fn_order(surf_pos->surf_num) + inc);
      else oi[i] = init_fn_ord(0);
    }
  }
  else {
    for (int i = 0; i < wf->neq; i++) oi[i] = init_fn_ord(0);
  }
  
  // Order of the shape function.
  Func<Ord>* ov = init_fn_ord(fv->get_edge_fn_order(surf_pos->surf_num) + inc);
  
  // Order of additional external functions.
  ExtData<Ord>* fake_ext = init_ext_fns_ord(vfs->ext, surf_pos->surf_num);
  
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
  
  int eo = quad->get_edge_points(surf_pos->surf_num, order);
  double3* pt = quad->get_points(eo);
  int np = quad->get_num_points(eo);

  // Init geometry and jacobian*weights.
  if (cache_e[eo] == NULL)
  {
    cache_e[eo] = init_geom_surf(rv, surf_pos, eo);
    double3* tan = rv->get_tangent(surf_pos->surf_num, eo);
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

double get_l2_norm(Vector* vec) 
{
  _F_
  scalar val = 0;
  for (int i = 0; i < vec->length(); i++) {
    scalar inc = vec->get(i);
    val = val + inc*conj(inc);
  }
  return sqrt(std::abs(val));
}

// Performs uniform global refinement of a FE space. 
Tuple<Space *> * construct_refined_spaces(Tuple<Space *> coarse, int order_increase)
{
  _F_
  Tuple<Space *> * ref_spaces = new Tuple<Space *>;
  for (int i = 0; i < coarse.size(); i++) 
  {
    Mesh* ref_mesh = new Mesh;
    ref_mesh->copy(coarse[i]->get_mesh());
    ref_mesh->refine_all_elements();
    ref_spaces->push_back(coarse[i]->dup(ref_mesh));
    (*ref_spaces)[i]->copy_orders(coarse[i], order_increase);
  }
  return ref_spaces;
}

// Light version for a single space.
Space* construct_refined_space(Space* coarse, int order_increase)
{
  _F_
  Mesh* ref_mesh = new Mesh;
  ref_mesh->copy(coarse->get_mesh());
  ref_mesh->refine_all_elements();
  Space* ref_space = coarse->dup(ref_mesh);
  ref_space->copy_orders(coarse, order_increase);
  return ref_space;
}
