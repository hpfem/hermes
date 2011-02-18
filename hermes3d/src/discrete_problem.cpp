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

#define HERMES_REPORT_INFO
#define HERMES_REPORT_WARN


#include "h3d_common.h"
#include "discrete_problem.h"
#include "traverse.h"
#include "../../hermes_common/error.h"
#include "../../hermes_common/callstack.h"



DiscreteProblem::FnCache::~FnCache()
{
  _F_
  free();
}

void DiscreteProblem::FnCache::free()
{
  _F_
  for(std::map<unsigned int, double *>::iterator it = jwt.begin(); it != jwt.end(); it++)
    delete [] it->second;
  jwt.clear();
  for(std::map<unsigned int, Geom<double> >::iterator it = e.begin(); it != e.end(); it++)
    free_geom(&it->second);
  e.clear();
  for (std::map<fn_key_t, sFunc*>::iterator it = fn.begin(); it != fn.end(); it++)
    free_fn(it->second);
  fn.clear();
  for (std::map<fn_key_t, mFunc*>::iterator it = ext.begin(); it != ext.end(); it++)
    free_fn(it->second);
  ext.clear();
  for (std::map<fn_key_t, mFunc*>::iterator it = sln.begin(); it != sln.end(); it++)
    free_fn(it->second);
  sln.clear();
}

// DiscreteProblem ///////////////////////////////////////////////////////////////////////////////////////

DiscreteProblem::DiscreteProblem(WeakForm *wf, Hermes::vector<Space *> spaces, bool is_linear)
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

  this->spaces = Hermes::vector<Space *>();
  for (int i = 0; i < wf->neq; i++) this->spaces.push_back(spaces[i]);
  have_spaces = true;

  // H2D is initializing precalc shapesets here.

  // Create global enumeration of dof and fill the ndof variable
  this->ndof = Space::assign_dofs(this->spaces);
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
  if (!is_up_to_date()) 
  {
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
void DiscreteProblem::create_sparse_structure(SparseMatrix *mat, Vector* rhs, bool rhsonly)
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

    AsmList *al = new AsmList[wf->neq];
    Mesh **meshes = new Mesh*[wf->neq];
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
        if (e[i] != NULL) spaces[i]->get_element_assembly_list(e[i], al + i);
      }

      // go through all equation-blocks of the local stiffness matrix
      for (int m = 0; m < wf->neq; m++)
      {
        for (int n = 0; n < wf->neq; n++)
        {
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
      }
    }

    trav.finish();
    delete [] al;
    delete [] meshes;
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

void DiscreteProblem::assemble(scalar* coeff_vec, SparseMatrix* mat, Vector* rhs, bool rhsonly)
{
  /* BEGIN IDENTICAL CODE WITH H2D */

  _F_
  // Sanity checks.
  if (coeff_vec == NULL && this->is_linear == false) error("coeff_vec is NULL in FeProblem::assemble().");
  if (!have_spaces) error("You have to call FeProblem::set_spaces() before calling assemble().");
  for (int i=0; i<this->wf->neq; i++)
  {
    if (this->spaces[i] == NULL) error("A space is NULL in assemble().");
  }
 
  this->create_sparse_structure(mat, rhs, rhsonly);

  // Convert the coefficient vector 'coeff_vec' into solutions Hermes::vector 'u_ext'.
  Hermes::vector<Solution*> u_ext;
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

  /* END IDENTICAL CODE WITH H2D */

  bool bnd[10];         // FIXME: magic number - maximal possible number of element surfaces
  SurfPos surf_pos[10];
  AsmList *al = new AsmList[wf->neq];
  bool *nat = new bool[wf->neq];
  bool *isempty = new bool[wf->neq];
  AsmList *am, *an;

  ShapeFunction *base_fn = new ShapeFunction[wf->neq];
  ShapeFunction *test_fn = new ShapeFunction[wf->neq];
  ShapeFunction *fu, *fv;
  RefMap * refmap = new RefMap[wf->neq];
  for (int i = 0; i < wf->neq; i++) 
  {
    base_fn[i].set_shapeset(spaces[i]->get_shapeset());
    test_fn[i].set_shapeset(spaces[i]->get_shapeset());
    refmap[i].set_mesh(spaces[i]->get_mesh());
  }

  // initialize matrix buffer
  matrix_buffer = NULL;
  matrix_buffer_dim = 0;
  if (mat != NULL) get_matrix_buffer(10);

  // obtain a list of assembling stages
  std::vector<WeakForm::Stage> stages;
  wf->get_stages(spaces, u_ext, stages, rhsonly);

  // Loop through all assembling stages -- the purpose of this is increased performance
  // in multi-mesh calculations, where, e.g., only the right hand side uses two meshes.
  // In such a case, the matrix forms are assembled over one mesh, and only the rhs
  // traverses through the union mesh. On the other hand, if you don't use multi-mesh
  // at all, there will always be only one stage in which all forms are assembled as usual.
  Traverse trav;
  for (unsigned ss = 0; ss < stages.size(); ss++)
  {
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

      // H2D has here:
      /* update_limit_table(e0->get_mode()); */

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
        spaces[j]->get_element_assembly_list(e[i], al + j);

        // This is different in H2D (PrecalcShapeset is used).
        test_fn[j].set_active_element(e[i]);
        test_fn[j].set_transform(base_fn + j);

        // This is different in H2D (PrecalcShapeset is used).
        refmap[j].set_active_element(e[i]);
        refmap[j].force_transform(base_fn[j].get_transform(), base_fn[j].get_ctm());
      }
      int marker = e0->marker;

      fn_cache.free();  // This is different in H2D.

      if (mat != NULL) 
      {
        // assemble volume matrix forms //////////////////////////////////////
        for (unsigned ww = 0; ww < s->mfvol.size(); ww++) 
        {
          WeakForm::MatrixFormVol *mfv = s->mfvol[ww];
          if (isempty[mfv->i] || isempty[mfv->j]) continue;
          if (mfv->area != HERMES_ANY && !wf->is_in_area(marker, mfv->area)) continue;
          int m = mfv->i; fv = test_fn + m; am = al + m;
          int n = mfv->j; fu = base_fn + n; an = al + n;
          bool tra = (m != n) && (mfv->sym != HERMES_NONSYM);
          bool sym = (m == n) && (mfv->sym == HERMES_SYM);

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
                else if (rhsonly == false) 
                {
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
      for (int isurf = 0; isurf < e0->get_num_surf(); isurf++)
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
            WeakForm::MatrixFormSurf* mfs = s->mfsurf[ww];
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
                if (an->dof[j] < 0) 
                {
                  // Linear problems only: Subtracting Dirichlet lift contribution from the RHS:
                  if (rhs != NULL && this->is_linear) 
                  {
                    scalar val = eval_form(mfs, u_ext, fu, fv, refmap + n, refmap + m, 
                                           surf_pos + isurf) * an->coef[j] * am->coef[i];
                    rhs->add(am->dof[i], -val);
                  }
                }
                else if (rhsonly == false) 
                {
                  scalar val = eval_form(mfs, u_ext, fu, fv, refmap + n, refmap + m, 
                                         surf_pos + isurf) * an->coef[j] * am->coef[i];
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

    if (mat != NULL) mat->finish();
    if (rhs != NULL) rhs->finish();
    trav.finish();
  }
 
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

  // Clean up.
  delete [] isempty;
  delete [] nat;
  delete [] al;
  delete [] base_fn;
  delete [] test_fn;
  delete [] refmap;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void DiscreteProblem::init_ext_fns(ExtData<scalar> &ext_data, std::vector<MeshFunction *> &ext, int order,
                             RefMap *rm, const int np, const QuadPt3D *pt)
{
  _F_

  ext_data.nf = ext.size();
  mFunc **ext_fn = new mFunc * [ext_data.nf];
  for (int i = 0; i < ext_data.nf; i++) 
  {
    fn_key_t key(ext[i]->seq, order, ext[i]->get_transform());
    if (fn_cache.ext.find(key) == fn_cache.ext.end()) 
    {
      fn_cache.ext[key] = init_fn(ext[i], rm, np, pt);
    }
    assert(fn_cache.ext[key] != NULL);
    ext_fn[i] = fn_cache.ext[key];
  }
  ext_data.fn = ext_fn;
}

void DiscreteProblem::init_ext_fns(ExtData<Ord> &fake_ext_data, std::vector<MeshFunction *> &ext)
{
  _F_

  fake_ext_data.nf = ext.size();
  Func<Ord> **fake_ext_fn = new Func<Ord> *[fake_ext_data.nf];
  
  for (int i = 0; i < fake_ext_data.nf; i++) 
    fake_ext_fn[i] = init_fn_ord(ext[i]->get_fn_order());
  
  fake_ext_data.fn = fake_ext_fn;
}

sFunc *DiscreteProblem::get_fn(ShapeFunction *fu, int order, RefMap *rm, const int np, const QuadPt3D *pt)
{
  fn_key_t key(fu->get_active_shape(), order, fu->get_transform(), fu->get_shapeset()->id);
  if (fn_cache.fn.find(key) == fn_cache.fn.end())
    fn_cache.fn[key] = init_fn(fu, rm, np, pt);
  return fn_cache.fn[key];
}

sFunc *DiscreteProblem::get_fn(ShapeFunction *fu, int order, RefMap *rm, int isurf, const int np,
                         const QuadPt3D *pt)
{
  fn_key_t key(fu->get_active_shape(), order, fu->get_transform(), fu->get_shapeset()->id);
  
  if (fn_cache.fn.find(key) == fn_cache.fn.end()) 
    fn_cache.fn[key] = init_fn(fu, rm, isurf, np, pt);
  return fn_cache.fn[key];
}

mFunc *DiscreteProblem::get_fn(Solution *fu, int order, RefMap *rm, const int np, const QuadPt3D *pt)
{
  fn_key_t key(fu->seq, order, fu->get_transform());

  if (fn_cache.sln.find(key) == fn_cache.sln.end()) 
    fn_cache.sln[key] = init_fn(fu, rm, np, pt);
  return fn_cache.sln[key];
}

scalar DiscreteProblem::eval_form(WeakForm::MatrixFormVol *mfv, Hermes::vector<Solution *> u_ext, ShapeFunction *fu,
                            ShapeFunction *fv, RefMap *ru, RefMap *rv)
{
  _F_
  // At this point H2D sets an increase of one if 
  // fu->get_num_components() == 2.

  // This is missing in H2D:
  Element *elem = fv->get_active_element();

  // Determine the integration order
  Func<Ord> **oi = new Func<Ord> *[wf->neq];

  // Order of solutions from the previous Newton iteration.
  if (u_ext != Hermes::vector<Solution *>()) 
  {
    for (int i = 0; i < wf->neq; i++) 
    {
      if (u_ext[i] != NULL) oi[i] = init_fn_ord(u_ext[i]->get_fn_order());
      else oi[i] = init_fn_ord(0);
    }
  } 
  else 
  {
    for (int i = 0; i < wf->neq; i++) oi[i] = init_fn_ord(0);
  }

  // Order of shape functions.
  Func<Ord> *ou = init_fn_ord(fu->get_fn_order());
  Func<Ord> *ov = init_fn_ord(fv->get_fn_order());

  // Order of additional external functions.
  ExtData<Ord> fake_ext;
  init_ext_fns(fake_ext, mfv->ext);

  // Order of geometric attributes (eg. for multiplication of a solution with coordinates, normals, etc.).
  double fake_wt = 1.0;
  Geom<Ord> fake_e = init_geom(elem->marker);

  // Total order of the matrix form.
  Ord o = mfv->ord(1, &fake_wt, oi, ou, ov, &fake_e, &fake_ext);

  // Increase due to reference map.
  Ord3 order = ru->get_inv_ref_order();
  switch (order.type) {
    case HERMES_MODE_TET: order += Ord3(o.get_order()); break;
    case HERMES_MODE_HEX: order += Ord3(o.get_order(), o.get_order(), o.get_order()); break;
  }
  order.limit();
  int ord_idx = order.get_idx();

  // Clean up.
  free_ext_fns_ord(&fake_ext);
  for (int i = 0; i < wf->neq; i++) {
    free_fn(oi[i]);
    delete oi[i];
  }
  delete [] oi;
  free_fn(ou);
  free_fn(ov);
  delete ou;
  delete ov;

  // Evaluate the form using the quadrature of the just calculated order.
  Quad3D *quad = get_quadrature(elem->get_mode());
  int np = quad->get_num_points(order);
  QuadPt3D *pt = quad->get_points(order);

  // Init geometry and jacobian*weights.
  double *jwt = NULL;
  Geom<double> e;
  if (fn_cache.e.find(ord_idx) == fn_cache.e.end()) 
  {
    fn_cache.jwt[ord_idx] = ru->get_jacobian(np, pt);
    fn_cache.e[ord_idx] = init_geom(elem->marker, ru, np, pt);
  }
  jwt = fn_cache.jwt[ord_idx];
  e = fn_cache.e[ord_idx];

  // Values of the previous Newton iteration, shape functions and external functions in quadrature points.
  mFunc **prev = new mFunc *[wf->neq];
  // OLD CODE: for (int i = 0; i < wf->neq; i++) prev[i] = get_fn(u_ext[i], ord_idx, rv, np, pt);
  if (u_ext != Hermes::vector<Solution *>()) 
  {
    for (int i = 0; i < wf->neq; i++) 
    {
      if (u_ext[i] != NULL) prev[i] = get_fn(u_ext[i], ord_idx, rv, np, pt);
      else prev[i] = NULL;
    }
  }
  else 
  {
    for (int i = 0; i < wf->neq; i++) prev[i] = NULL;
  }

  sFunc *u = get_fn(fu, ord_idx, ru, np, pt);
  sFunc *v = get_fn(fv, ord_idx, rv, np, pt);
  ExtData<scalar> ext;
  init_ext_fns(ext, mfv->ext, ord_idx, rv, np, pt);

  scalar res = mfv->fn(np, jwt, prev, u, v, &e, &ext);

  // Clean up.
  delete [] prev;
  //ext.free(); // FIXME: this needs to be unified with H2D.
  return res;
}

scalar DiscreteProblem::eval_form(WeakForm::VectorFormVol *vfv, Hermes::vector<Solution *> u_ext, ShapeFunction *fv, RefMap *rv)
{
  _F_
  // At this point H2D sets an increase of one if 
  // fu->get_num_components() == 2.

  // This is missing in H2D:
  Element *elem = fv->get_active_element();

  // Determine the integration order.
  Func<Ord> **oi = new Func<Ord> * [wf->neq];

  // Order of solutions from the previous Newton iteration.
  if (u_ext != Hermes::vector<Solution *>()) 
  {
    for (int i = 0; i < wf->neq; i++) 
    {
      if (u_ext[i] != NULL) oi[i] = init_fn_ord(u_ext[i]->get_fn_order());
      else oi[i] = init_fn_ord(0);
    }
  } 
  else 
  {
    for (int i = 0; i < wf->neq; i++) oi[i] = init_fn_ord(0);
  }

  // Order of the shape function.
  Func<Ord> *ov = init_fn_ord(fv->get_fn_order());

  // Order of additional external functions.
  ExtData<Ord> fake_ext;
  init_ext_fns(fake_ext, vfv->ext);

  // Order of geometric attributes (eg. for multiplication of a solution with coordinates, normals, etc.).
  double fake_wt = 1.0;
  Geom<Ord> fake_e = init_geom(elem->marker);

  // Total order of the vector form.
  Ord o = vfv->ord(1, &fake_wt, oi, ov, &fake_e, &fake_ext);

  // Increase due to reference map.
  Ord3 order = rv->get_inv_ref_order();
  switch (order.type) 
  {
    case HERMES_MODE_TET: order += Ord3(o.get_order()); break;
    case HERMES_MODE_HEX: order += Ord3(o.get_order(), o.get_order(), o.get_order()); break;
  }
  order.limit();
  int ord_idx = order.get_idx();

  // Clean up.
  free_ext_fns_ord(&fake_ext);
  for (int i = 0; i < wf->neq; i++) {
    free_fn(oi[i]);
    delete oi[i];
  }
  delete [] oi;
  free_fn(ov);
  delete ov;

  // Evaluate the form using the quadrature of the just calculated order.
  Quad3D *quad = get_quadrature(elem->get_mode());
  int np = quad->get_num_points(order);
  QuadPt3D *pt = quad->get_points(order);

        // Init geometry and jacobian*weights.
  double *jwt = NULL;
  Geom<double> e;
  if (fn_cache.e.find(ord_idx) == fn_cache.e.end()) 
  {
    fn_cache.jwt[ord_idx] = rv->get_jacobian(np, pt);
    fn_cache.e[ord_idx] = init_geom(elem->marker, rv, np, pt);
  }
  jwt = fn_cache.jwt[ord_idx];
  e = fn_cache.e[ord_idx];

  // Values of the previous Newton iteration, shape functions and external functions in quadrature points.
  mFunc ** prev = new mFunc *[wf->neq];
  // OLD CODE: for (int i = 0; i < wf->neq; i++) prev[i] = get_fn(u_ext[i], ord_idx, rv, np, pt);
  if (u_ext != Hermes::vector<Solution *>()) 
  {
    for (int i = 0; i < wf->neq; i++) 
    {
      if (u_ext[i] != NULL) prev[i] = get_fn(u_ext[i], ord_idx, rv, np, pt);
      else prev[i] = NULL;
    }
  }
  else 
  {
    for (int i = 0; i < wf->neq; i++) prev[i] = NULL;
  }
  sFunc *v = get_fn(fv, ord_idx, rv, np, pt);

  ExtData<scalar> ext;
  init_ext_fns(ext, vfv->ext, ord_idx, rv, np, pt);

  scalar res = vfv->fn(np, jwt, prev, v, &e, &ext);

  // Clean up.
  delete [] prev;
  //ext.free();// FIXME: this needs to be unified with H2D.
  
  return res;
}

scalar DiscreteProblem::eval_form(WeakForm::MatrixFormSurf *mfs, Hermes::vector<Solution *> u_ext, ShapeFunction *fu,
                                  ShapeFunction *fv, RefMap *ru, RefMap *rv, SurfPos *surf_pos)
{
  _F_
  // At this point H2D sets an increase of one if 
  // fu->get_num_components() == 2.

  // Determine the integration order.
  Func<Ord> **oi = new Func<Ord> *[wf->neq];
  
  // Order of solutions from the previous Newton iteration.
  if (u_ext != Hermes::vector<Solution *>()) 
  {
    for (int i = 0; i < wf->neq; i++) 
    {
      if (u_ext[i] != NULL) oi[i] = init_fn_ord(u_ext[i]->get_fn_order());
      else oi[i] = init_fn_ord(0);
    }
  } 
  else 
  {
    for (int i = 0; i < wf->neq; i++) oi[i] = init_fn_ord(0);
  }

  // Order of the shape functions.
  Func<Ord> *ou = init_fn_ord(fu->get_fn_order());
  Func<Ord> *ov = init_fn_ord(fv->get_fn_order());

  // Order of additional external functions.
  ExtData<Ord> fake_ext;
  init_ext_fns(fake_ext, mfs->ext);

  // Order of geometric attributes (eg. for multiplication of a solution with coordinates, normals, etc.).
  double fake_wt = 1.0;
  Geom<Ord> fake_e = init_geom(surf_pos->marker);

  // Total order of the surface matrix form.
  Ord o = mfs->ord(1, &fake_wt, oi, ou, ov, &fake_e, &fake_ext);

  // Increase due to reference map.
  Ord3 order = ru->get_inv_ref_order();
  switch (order.type) 
  {
    case HERMES_MODE_TET: order += Ord3(o.get_order()); break;
    case HERMES_MODE_HEX: order += Ord3(o.get_order(), o.get_order(), o.get_order()); break;
  }
  order.limit();
  Ord2 face_order = order.get_face_order(surf_pos->surf_num);
  int ord_idx = face_order.get_idx();

  // Clean up.
  free_ext_fns_ord(&fake_ext);
  for (int i = 0; i < wf->neq; i++) {
    free_fn(oi[i]);
    delete oi[i];
  }
  delete [] oi;
  free_fn(ou);
  free_fn(ov);
  delete ou;
  delete ov;

  // Evaluate the form using the quadrature of the just calculated order.
  Quad3D *quad = get_quadrature(fu->get_active_element()->get_mode());
  int np = quad->get_face_num_points(surf_pos->surf_num, face_order);
  QuadPt3D *pt = quad->get_face_points(surf_pos->surf_num, face_order);

        // Init geometry and jacobian*weights.
  double *jwt = NULL;
  Geom<double> e;
  if (fn_cache.e.find(ord_idx) == fn_cache.e.end()) 
  {
    fn_cache.jwt[ord_idx] = ru->get_face_jacobian(surf_pos->surf_num, np, pt);
    fn_cache.e[ord_idx] = init_geom(surf_pos->marker, ru, surf_pos->surf_num, np, pt);
  }
  jwt = fn_cache.jwt[ord_idx];
  e = fn_cache.e[ord_idx];

  // Values of the previous Newton iteration, shape functions and external functions in quadrature points.
  mFunc **prev = new mFunc *[wf->neq];
  // OLD CODE: for (int i = 0; i < wf->neq; i++) prev[i] = get_fn(u_ext[i], ord_idx, rv, np, pt);
  if (u_ext != Hermes::vector<Solution *>()) 
  {
    for (int i = 0; i < wf->neq; i++) 
    {
      if (u_ext[i] != NULL) prev[i] = get_fn(u_ext[i], ord_idx, rv, np, pt);
      else prev[i] = NULL;
    }
  }
  else 
  {
    for (int i = 0; i < wf->neq; i++) prev[i] = NULL;
  }

  sFunc *u = get_fn(fu, ord_idx, ru, surf_pos->surf_num, np, pt);
  sFunc *v = get_fn(fv, ord_idx, rv, surf_pos->surf_num, np, pt);

  ExtData<scalar> ext;
  init_ext_fns(ext, mfs->ext, ord_idx, rv, np, pt);

  scalar res = mfs->fn(np, jwt, prev, u, v, &e, &ext);

  // Clean up.
  delete [] prev;  
  //ext.free();// FIXME: this needs to be unified with H2D.
  
  return res;
}

scalar DiscreteProblem::eval_form(WeakForm::VectorFormSurf *vfs, Hermes::vector<Solution *> u_ext, 
                                  ShapeFunction *fv, RefMap *rv, SurfPos *surf_pos)
{
  _F_

  // Determine the integration order.
  Func<Ord> **oi = new Func<Ord> *[wf->neq];

  // Order of solutions from the previous Newton iteration.
  if (u_ext != Hermes::vector<Solution *>()) 
  {
    for (int i = 0; i < wf->neq; i++) 
    {
      if (u_ext[i] != NULL) oi[i] = init_fn_ord(u_ext[i]->get_fn_order());
      else oi[i] = init_fn_ord(0);
    }
  } 
  else 
  {
    for (int i = 0; i < wf->neq; i++) oi[i] = init_fn_ord(0);
  }

  // Order of the shape function.
  Func<Ord> *ov = init_fn_ord(fv->get_fn_order());

  // Order of additional external functions.
  ExtData<Ord> fake_ext;
  init_ext_fns(fake_ext, vfs->ext);

  // Order of geometric attributes (eg. for multiplication of a solution with coordinates, normals, etc.).
  double fake_wt = 1.0;
  Geom<Ord> fake_e = init_geom(surf_pos->marker);

  // Total order of the surface vector form.
  Ord o = vfs->ord(1, &fake_wt, oi, ov, &fake_e, &fake_ext);

  // Increase due to reference map.
  Ord3 order = rv->get_inv_ref_order();
  switch (order.type) 
  {
    case HERMES_MODE_TET: order += Ord3(o.get_order()); break;
    case HERMES_MODE_HEX: order += Ord3(o.get_order(), o.get_order(), o.get_order()); break;
  }
  order.limit();
  Ord2 face_order = order.get_face_order(surf_pos->surf_num);
  int ord_idx = face_order.get_idx();
 
  // Clean up.
  free_ext_fns_ord(&fake_ext);
  for (int i = 0; i < wf->neq; i++) {
    free_fn(oi[i]);
    delete oi[i];
  }
  delete [] oi;
  free_fn(ov);
  delete ov;

  // Evaluate the form using the quadrature of the just calculated order.
  Quad3D *quad = get_quadrature(fv->get_active_element()->get_mode());
  int np = quad->get_face_num_points(surf_pos->surf_num, face_order);
  QuadPt3D *pt = quad->get_face_points(surf_pos->surf_num, face_order);

  // Init geometry and jacobian*weights.
  double *jwt = NULL;
  Geom<double> e;
  if (fn_cache.e.find(ord_idx) == fn_cache.e.end()) {
    fn_cache.jwt[ord_idx] = rv->get_face_jacobian(surf_pos->surf_num, np, pt);
    fn_cache.e[ord_idx] = init_geom(surf_pos->marker, rv, surf_pos->surf_num, np, pt);
  }
  jwt = fn_cache.jwt[ord_idx];
  e = fn_cache.e[ord_idx];

  // Values of the previous Newton iteration, shape functions and external functions in quadrature points.
  mFunc **prev = new mFunc *[wf->neq];
  // OLD CODE: for (int i = 0; i < wf->neq; i++) prev[i] = get_fn(u_ext[i], ord_idx, rv, np, pt);
  if (u_ext != Hermes::vector<Solution *>()) 
  {
    for (int i = 0; i < wf->neq; i++) 
    {
      if (u_ext[i] != NULL) prev[i] = get_fn(u_ext[i], ord_idx, rv, np, pt);
      else prev[i] = NULL;
    }
  }
  else 
  {
    for (int i = 0; i < wf->neq; i++) prev[i] = NULL;
  }
  sFunc *v = get_fn(fv, ord_idx, rv, surf_pos->surf_num, np, pt);

  ExtData<scalar> ext;
  init_ext_fns(ext, vfs->ext, ord_idx, rv, np, pt);

  scalar res = vfs->fn(np, jwt, prev, v, &e, &ext);

  // Clean up.
  delete [] prev;
  //ext.free();// FIXME: this needs to be unified with H2D.
  
  return res;
}

////////////////////////////////////////////////////////////////////////////////////////

Hermes::vector<Space *> * construct_refined_spaces(Hermes::vector<Space *> coarse, int order_increase)
{
  _F_
  Hermes::vector<Space *> * ref_spaces = new Hermes::vector<Space *>;
  for (unsigned int i = 0; i < coarse.size(); i++) 
  {
    Mesh* ref_mesh = new Mesh;
    ref_mesh->copy(*coarse[i]->get_mesh());
    ref_mesh->refine_all_elements(H3D_H3D_H3D_REFT_HEX_XYZ);
    ref_spaces->push_back(coarse[i]->dup(ref_mesh));
    (*ref_spaces)[i]->copy_orders(*coarse[i], order_increase);
  }
  return ref_spaces;
}

// Light version for a single space.
Space* construct_refined_space(Space* coarse, int order_increase)
{
  _F_
  Mesh* ref_mesh = new Mesh;
  ref_mesh->copy(*coarse->get_mesh());
  ref_mesh->refine_all_elements(H3D_H3D_H3D_REFT_HEX_XYZ);
  Space* ref_space = coarse->dup(ref_mesh);
  ref_space->copy_orders(*coarse, order_increase);
  return ref_space;
}

// Perform Newton's iteration.
bool HERMES_RESIDUAL_AS_VECTOR = false;
bool solve_newton(scalar* coeff_vec, DiscreteProblem* dp, Solver* solver, SparseMatrix* matrix,
                  Vector* rhs, double newton_tol, int newton_max_iter, bool verbose,
                  double damping_coeff, double max_allowed_residual_norm)
{
  info("Solve newton");
  // Prepare solutions for measuring residual norm.
  int num_spaces = dp->get_num_spaces();

  Hermes::vector<Solution*> solutions;
  Hermes::vector<double> dir_lift_false;
  for (int i=0; i < num_spaces; i++) {
    solutions.push_back(new Solution(dp->get_space(i)->get_mesh()));
    dir_lift_false.push_back(0.0);      // No Dirichlet lifts will be considered.
  }

  // The Newton's loop.
  double residual_norm;
  int it = 1;
  while (1)
  {
    // Obtain the number of degrees of freedom.
    int ndof = dp->get_num_dofs();
    info("Newton iteration beginning, NDOFs %i", ndof);
    // Assemble the Jacobian matrix and residual vector.
    dp->assemble(coeff_vec, matrix, rhs, false);
    // Multiply the residual vector with -1 since the matrix
    // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
    for (int i = 0; i < ndof; i++) rhs->set(i, -rhs->get(i));

    // Measure the residual norm.


    scalar* rhs_coeffs = new scalar[rhs->length()];
    rhs->extract(rhs_coeffs);
    Solution::vector_to_solutions(rhs_coeffs, dp->get_spaces(), solutions, dir_lift_false);
    if (HERMES_RESIDUAL_AS_VECTOR) {
          // Calculate the l2-norm of residual vector.
          // residual_norm = get_l2_norm(rhs);
          // FIXME what's that?
    } else {

    info("vector_to_solutions completed");
      double norm_squares;
      for (unsigned int i = 0; i < solutions.size(); i++) {
        double norm = h1_norm(solutions[i]);
        norm_squares += norm * norm;
        residual_norm = sqrt(norm_squares);
      }
    }

    // Info for the user.
    if (verbose) {
      info("---- Newton iter %d, ndof %d, residual norm %g", it, ndof, residual_norm);
    }

    // If maximum allowed residual norm is exceeded, fail.
    if (residual_norm > max_allowed_residual_norm) {
      if (verbose) {
        info("Current residual norm: %g", residual_norm);
        info("Maximum allowed residual norm: %g", max_allowed_residual_norm);
        info("Newton solve not successful, returning false.");
      }
      return false;
    }

    // If residual norm is within tolerance, or the maximum number
    // of iteration has been reached, then quit.
    if ((residual_norm < newton_tol || it > newton_max_iter) && it > 1) break;

    // Solve the linear system.
    if(!solver->solve()) error ("Matrix solver failed.\n");

    // Add \deltaY^{n+1} to Y^n.
    for (int i = 0; i < ndof; i++) coeff_vec[i] += damping_coeff * solver->get_solution()[i];

    it++;
  }

  if (it >= newton_max_iter) {
    if (verbose) {
      info("Maximum allowed number of Newton iterations exceeded, returning false.");
    }
    return false;
  }

  return true;
}

