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
#include "discrete_problem.h"
#include "solver.h"
#include "traverse.h"
#include "space/space.h"
#include "precalc.h"
#include "shapeset/shapeset_h1_all.h"
#include "integrals_h1.h"
#include "refmap.h"
#include "solution.h"
#include "config.h"
#include "limit_order.h"
#include <algorithm>
#include "views/view.h"
#include "views/scalar_view.h"
#include "views/vector_view.h"
#include "tuple.h"
#include "norm.h"

int H2D_DEFAULT_PROJ_NORM = 1;

void qsort_int(int* pbase, size_t total_elems); // defined in qsort.cpp

//// interface /////////////////////////////////////////////////////////////////////////////////////

void DiscreteProblem::init(WeakForm* wf_)
{
  if (wf_ == NULL) error("DiscreteProblem: a weak form must be given.");
  this->wf = wf_;
  this->solver_default = NULL;
  this->solver = NULL;
  this->wf_seq = -1;

  this->mat_sym = false;

  this->spaces = NULL;
  this->pss = NULL;

  this->values_changed = true;
  this->struct_changed = true;
  this->have_spaces = false;
}

// this is needed because of a constructor in NonlinSystem
DiscreteProblem::DiscreteProblem() {}

DiscreteProblem::DiscreteProblem(WeakForm* wf_)
{
  this->init(wf_);
}

DiscreteProblem::DiscreteProblem(WeakForm* wf_, Hermes::Tuple<Space*> sp)
{
  this->init(wf_);
  this->init_spaces(sp);
}

DiscreteProblem::DiscreteProblem(WeakForm* wf_, Space* s_)
{
  this->init(wf_);
  this->init_space(s_);
}

DiscreteProblem::~DiscreteProblem()
{
  free();

  if (this->sp_seq != NULL) delete [] this->sp_seq;
  if (this->pss != NULL) {
    for (int i = 0; i < this->wf->neq; i++)
      if (this->pss[i] != NULL) delete this->pss[i];
    delete [] this->pss;
  }
  if (this->solver_default != NULL) delete this->solver_default;
}

// NOTE: This should not be called in the destructor to DiscreteProblem
// since Discrete Problems are used for intermediary operations such as 
// projections, and thustheir destruction should not destruct their spaces. 
void DiscreteProblem::free_spaces()
{
  // free spaces, making sure that duplicated ones do not get deleted twice
  if (this->spaces != Hermes::Tuple<Space *>())
  {
    for (int i = 0; i < this->wf->neq; i++) {
      // this loop skipped if there is only one space
      for (int j = i+1; j < this->wf->neq; j++) {
        if (this->spaces[j] == this->spaces[i]) this->spaces[j] = NULL;
      }
    }
    for (int i = 0; i < this->wf->neq; i++) {
      if (this->spaces[i] != NULL) {
        this->spaces[i]->free();
        //delete this->spaces[i];
        this->spaces[i] = NULL;
      }
    }
    this->spaces = Hermes::Tuple<Space *>();
  }
}

// Should not be called by the user.
void DiscreteProblem::init_spaces(Hermes::Tuple<Space*> spaces)
{
  int n = spaces.size();
  if (n != this->wf->neq)
    error("Number of spaces does not match number of equations in DiscreteProblem::init_spaces().");

  // initialize spaces
  this->spaces = spaces;
  this->sp_seq = new int[this->wf->neq];
  memset(sp_seq, -1, sizeof(int) * this->wf->neq);
  this->assign_dofs(); // Create global enumeration of DOF in all spaces in the system. NOTE: this
                       // overwrites possible existing local enumeration of DOF in the spaces.
  this->have_spaces = true;

  // initialize precalc shapesets
  this->pss = new PrecalcShapeset*[this->wf->neq];
  for (int i=0; i<this->wf->neq; i++) this->pss[i] = NULL;
  this->num_user_pss = 0;
  for (int i = 0; i < n; i++){
    Shapeset *shapeset = spaces[i]->get_shapeset();
    if (shapeset == NULL) error("Internal in DiscreteProblem::init_spaces().");
    PrecalcShapeset *p = new PrecalcShapeset(shapeset);
    if (p == NULL) error("New PrecalcShapeset could not be allocated in DiscreteProblem::init_spaces().");
    this-> pss[i] = p;
    this->num_user_pss++;
  }
}

// Should not be called by the user.
void DiscreteProblem::init_space(Space* s)
{
  if (this->wf->neq != 1)
    error("Do not call init_space() for PDE systems, call init_spaces() instead.");
  this->init_spaces(Hermes::Tuple<Space*>(s));
}

// Obsolete. Should be removed after FeProblem is removed.
void DiscreteProblem::set_spaces(Hermes::Tuple<Space*>spaces)
{
  this->init_spaces(spaces);
}

void DiscreteProblem::set_pss(Hermes::Tuple<PrecalcShapeset*> pss)
{
  warn("Call to deprecated function DiscreteProblem::set_pss().");
  int n = pss.size();
  if (n != this->wf->neq)
    error("The number of precalculated shapesets must match the number of equations.");

  for (int i = 0; i < n; i++) this->pss[i] = pss[i];
  num_user_pss = n;
}

void DiscreteProblem::set_pss(PrecalcShapeset* pss)
{
  this->set_pss(Hermes::Tuple<PrecalcShapeset*>(pss));
}

void DiscreteProblem::copy(DiscreteProblem* sys)
{
  error("Not implemented yet.");
}

void DiscreteProblem::free()
{
  this->struct_changed = this->values_changed = true;
  memset(this->sp_seq, -1, sizeof(int) * this->wf->neq);
  this->wf_seq = -1;
}

//// assembly //////////////////////////////////////////////////////////////////////////////////////

void DiscreteProblem::insert_block(Matrix *mat_ext, 
     scalar** mat, int* iidx, int* jidx, int ilen, int jlen)
{
  mat_ext->add(ilen, jlen, mat, iidx, jidx);
}

void DiscreteProblem::assemble(Vector* init_vec, Matrix* mat_ext, Vector* dir_ext, 
                               Vector* rhs_ext, bool rhsonly, bool is_complex)
{
  // sanity checks
  if (this->have_spaces == false)
    error("Before assemble(), you need to initialize spaces.");
  if (this->wf == NULL) 
		error("this->wf = NULL in DiscreteProblem::assemble().");
  if (this->spaces == Hermes::Tuple<Space*>()) 
		error("this->spaces is empty in DiscreteProblem::assemble().");
  int n = this->wf->neq;
  for (int i = 0; i < n; i++) 
	{
    if (this->spaces[i] == NULL) 
			error("this->spaces[%d] is NULL in DiscreteProblem::assemble().", i);
  }
  if (rhs_ext == NULL) 
		error("rhs_ext == NULL in DiscreteProblem::assemble().");
  if (rhsonly == false) 
	{
    if (mat_ext == NULL) 
			error("mat_ext == NULL in DiscreteProblem::assemble().");
    if (mat_ext->get_size() != rhs_ext->length()) 
		{
			printf("mat_ext matrix size = %d\n", mat_ext->get_size());
			printf("rhs_ext vector size = %d\n", rhs_ext->length());
      error("Mismatched mat_ext and rhs_ext vector sizes in DiscreteProblem::assemble().");
    }
  }

  // Assign dof in all spaces. 
  int ndof = this->assign_dofs();
  if (ndof == 0) error("ndof = 0 in DiscreteProblem::assemble().");
  //printf("ndof = %d\n", ndof);
  
  // Realloc mat_ext, dir_ext and rhs_ext if ndof changed, 
  // and clear dir_ext and rhs_ext. 
  // Do not touch the matrix if rhsonly == true. 
  if (rhsonly == false) {
    if (mat_ext->get_size() != ndof) mat_ext->init(is_complex);
  }
  if (dir_ext != NULL) {
    if (dir_ext->get_size() != ndof) {
      dir_ext->free_data();
      dir_ext->init(ndof, is_complex);
    }
    else dir_ext->set_zero();
  }
  if (rhs_ext->get_size() != ndof) {
    rhs_ext->free_data();
    rhs_ext->init(ndof, is_complex);
  }
  else rhs_ext->set_zero();

  // If init_vec != NULL, convert it to a Hermes::Tuple of solutions u_ext.
  Hermes::Tuple<Solution*> u_ext = Hermes::Tuple<Solution*>();
  for (int i = 0; i < wf->neq; i++) {
    if (init_vec != NULL) {
      u_ext.push_back(new Solution(spaces[i]->get_mesh()));
      if (!init_vec->is_complex())
        u_ext[i]->set_coeff_vector(spaces[i], this->pss[i], (scalar*)init_vec->get_c_array(), ndof);
      else
        u_ext[i]->set_coeff_vector(spaces[i], this->pss[i], (scalar*)init_vec->get_c_array_cplx(), ndof);
    }
    else u_ext.push_back(NULL);
  }

  int k, m, marker;
  std::vector<AsmList> al(wf->neq);
  AsmList* am, * an;
  bool bnd[4];
  std::vector<bool> nat(wf->neq), isempty(wf->neq);
  EdgePos ep[4];
  reset_warn_order();

  if (rhsonly == false) {
    trace("Creating matrix sparse structure...");
    mat_ext->free_data();
  }
  else trace("Reusing matrix sparse structure...");

  if (rhsonly == false) trace("Assembling stiffness matrix and rhs...");
  else trace("Assembling rhs only...");
  TimePeriod cpu_time;

  // create slave pss's for test functions, init quadrature points
  std::vector<PrecalcShapeset*> spss(wf->neq, static_cast<PrecalcShapeset*>(NULL));
  PrecalcShapeset *fu, *fv;
  std::vector<RefMap> refmap(wf->neq);
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
  // Returns assembling stages with correct meshes, ext_functions that are needed in a particular stage.
  wf->get_stages(spaces, u_ext, stages, rhsonly);

  // Loop through all assembling stages -- the purpose of this is increased performance
  // in multi-mesh calculations, where, e.g., only the right hand side uses two meshes.
  // In such a case, the matrix forms are assembled over one mesh, and only the rhs
  // traverses through the union mesh. On the other hand, if you don't use multi-mesh
  // at all, there will always be only one stage in which all forms are assembled as usual.
  Traverse trav;
  for (unsigned int ss = 0; ss < stages.size(); ss++)
  {
    WeakForm::Stage* s = &stages[ss];
    // Fills the solution and test functions into the stage s.
    for (unsigned int i = 0; i < s->idx.size(); i++)
			s->fns[i] = pss[s->idx[i]];
    for (unsigned int i = 0; i < s->ext.size(); i++)
			s->ext[i]->set_quad_2d(&g_quad_2d_std);
    // Tests whether the meshes in this stage are compatible and initializes the traverse process.
    trav.begin(s->meshes.size(), &(s->meshes.front()), &(s->fns.front()));

    // Assemble one stage.
    Element** e;
    // See Traverse::get_next_state for explanation.
    while ((e = trav.get_next_state(bnd, ep)) != NULL)
    {
      // Checking if at least on one mesh in this stage the element over which we are assembling is used.
      // If not, we continue with another state.
      Element* e0 = NULL;
      for (unsigned int i = 0; i < s->idx.size(); i++)
      if ((e0 = e[i]) != NULL) break;
      if (e0 == NULL) continue;

      // Set maximum integration order for use in integrals, see limit_order().
      update_limit_table(e0->get_mode());

      // Obtain assembly lists for the element at all spaces of the stage, set appropriate mode for each pss.
      // NOTE: Active elements and transformations for external functions (including the solutions from previous
      // Newton's iteration) as well as basis functions (master PrecalcShapesets) have already been set in 
      // trav.get_next_state(...).
      std::fill(isempty.begin(), isempty.end(), false);
      for (unsigned int i = 0; i < s->idx.size(); i++)
      {
        int j = s->idx[i];
        if (e[i] == NULL) { isempty[j] = true; continue; }
        // FIXME: Do not retrieve assembly list again if the element has not changed.
        spaces[j]->get_element_assembly_list(e[i], &al[j]);

	// Set active element to all test function PrecalcShapesets.
        spss[j]->set_active_element(e[i]);
	// Set the subelement transformation as it is set on all the appropriate solution function PrecalcShapeset.
        spss[j]->set_master_transform();
	// Set the active element to the reference mapping.
        refmap[j].set_active_element(e[i]);
	// Important : the reference mapping gets the same subelement transformation as the
	// appropriate PrecalcShapeset (~test function). This is used in eval_form functions.
        refmap[j].force_transform(pss[j]->get_transform(), pss[j]->get_ctm());
      }
      // Boundary marker.
      marker = e0->marker;

      init_cache();
      //// assemble volume matrix forms //////////////////////////////////////
      for (unsigned int ww = 0; ww < s->mfvol.size(); ww++)
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
                if (dir_ext != NULL) {
                  scalar val = eval_form(mfv, u_ext, fu, fv, &refmap[n], &refmap[m]) * an->coef[j] * am->coef[i];
                  dir_ext->add(am->dof[i], val);
                } 
              }
              else if (rhsonly == false) {
                scalar val = eval_form(mfv, u_ext, fu, fv, &refmap[n], &refmap[m]) * an->coef[j] * am->coef[i];
                local_stiffness_matrix[i][j] = val;
              }
            }
          }
          else // symmetric block
          {
            for (int j = 0; j < an->cnt; j++) {
              if (j < i && an->dof[j] >= 0) continue;
              fu->set_active_shape(an->idx[j]);
              if (an->dof[j] < 0) {
                if (dir_ext != NULL) {
                  scalar val = eval_form(mfv, u_ext, fu, fv, &refmap[n], &refmap[m]) * an->coef[j] * am->coef[i];
                  dir_ext->add(am->dof[i], val);
                }
              } 
              else if (rhsonly == false) {
                scalar val = eval_form(mfv, u_ext, fu, fv, &refmap[n], &refmap[m]) * an->coef[j] * am->coef[i];
                local_stiffness_matrix[i][j] = local_stiffness_matrix[j][i] = val;
              }
            }
          }
        }

        // insert the local stiffness matrix into the global one
        if (rhsonly == false) {
          insert_block(mat_ext, local_stiffness_matrix, am->dof, an->dof, am->cnt, an->cnt);
        }

        // insert also the off-diagonal (anti-)symmetric block, if required
        if (tra)
        {
          if (mfv->sym < 0) chsgn(local_stiffness_matrix, am->cnt, an->cnt);
          transpose(local_stiffness_matrix, am->cnt, an->cnt);
          if (rhsonly == false) {
            insert_block(mat_ext, local_stiffness_matrix, an->dof, am->dof, an->cnt, am->cnt);
          }

          // we also need to take care of the RHS...
          for (int j = 0; j < am->cnt; j++) {
            if (am->dof[j] < 0) {
              for (int i = 0; i < an->cnt; i++) {
                if (an->dof[i] >= 0) {
                  if (dir_ext != NULL) dir_ext->add(an->dof[i], local_stiffness_matrix[i][j]);
                }
              }
            }
          }
        }
      }

      //// assemble volume linear forms ////////////////////////////////////////
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
          scalar val = eval_form(vfv, u_ext, fv, &refmap[m]) * am->coef[i];
          rhs_ext->add(am->dof[i], val);
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
            spaces[j]->get_edge_assembly_list(e[i], edge, &al[j]);
        }

        // assemble surface matrix forms ///////////////////////////////////
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
                if (dir_ext != NULL) {
                  scalar val = eval_form(mfs, u_ext, fu, fv, &refmap[n], &refmap[m], &(ep[edge])) 
                               * an->coef[j] * am->coef[i];
                  dir_ext->add(am->dof[i], val);
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
            insert_block(mat_ext, local_stiffness_matrix, am->dof, an->dof, am->cnt, an->cnt);
          }
        }

        // assemble surface linear forms /////////////////////////////////////
        for (unsigned int ww = 0; ww < s->vfsurf.size(); ww++)
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
            scalar val = eval_form(vfs, u_ext, fv, &refmap[m], &(ep[edge])) * am->coef[i];
            rhs_ext->add(am->dof[i], val);
          }
        }
      }
      delete_cache();
    }
    trav.finish();
  }

  verbose("Stiffness matrix assembled (stages: %d)", stages.size());
  report_time("Stiffness matrix assembled in %g s", cpu_time.tick().last());
  for (int i = 0; i < wf->neq; i++) { 
    delete spss[i];
    if (u_ext[i] != NULL) delete u_ext[i];
  }
  delete [] buffer;

  if (rhsonly == false) values_changed = true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Initialize integration order for external functions
ExtData<Ord>* DiscreteProblem::init_ext_fns_ord(std::vector<MeshFunction *> &ext)
{
  ExtData<Ord>* fake_ext = new ExtData<Ord>;
  fake_ext->nf = ext.size();
  Func<Ord>** fake_ext_fn = new Func<Ord>*[fake_ext->nf];
  for (int i = 0; i < fake_ext->nf; i++)
    fake_ext_fn[i] = init_fn_ord(ext[i]->get_fn_order());
  fake_ext->fn = fake_ext_fn;

  return fake_ext;
}

// Initialize integration order on a given edge for external functions
ExtData<Ord>* DiscreteProblem::init_ext_fns_ord(std::vector<MeshFunction *> &ext, int edge)
{
  ExtData<Ord>* fake_ext = new ExtData<Ord>;
  fake_ext->nf = ext.size();
  Func<Ord>** fake_ext_fn = new Func<Ord>*[fake_ext->nf];
  for (int i = 0; i < fake_ext->nf; i++)
    fake_ext_fn[i] = init_fn_ord(ext[i]->get_edge_fn_order(edge));
  fake_ext->fn = fake_ext_fn;
  
  return fake_ext;
}

// Initialize external functions (obtain values, derivatives,...)
ExtData<scalar>* DiscreteProblem::init_ext_fns(std::vector<MeshFunction *> &ext, RefMap *rm, const int order)
{
  ExtData<scalar>* ext_data = new ExtData<scalar>;
  Func<scalar>** ext_fn = new Func<scalar>*[ext.size()];
  for (unsigned int i = 0; i < ext.size(); i++) {
    ext_fn[i] = init_fn(ext[i], rm, order);
  }
  ext_data->nf = ext.size();
  ext_data->fn = ext_fn;

  return ext_data;
}

// Initialize shape function values and derivatives (fill in the cache)
Func<double>* DiscreteProblem::get_fn(PrecalcShapeset *fu, RefMap *rm, const int order)
{
  Key key(256 - fu->get_active_shape(), order, fu->get_transform(), fu->get_shapeset()->get_id());
  if (cache_fn[key] == NULL)
    cache_fn[key] = init_fn(fu, rm, order);

  return cache_fn[key];
}

// Caching transformed values
void DiscreteProblem::init_cache()
{
  for (int i = 0; i < g_max_quad + 1 + 4 * g_max_quad + 4; i++)
  {
    cache_e[i] = NULL;
    cache_jwt[i] = NULL;
  }
}

void DiscreteProblem::delete_cache()
{
  for (int i = 0; i < g_max_quad + 1 + 4 * g_max_quad + 4; i++)
  {
    if (cache_e[i] != NULL)
    {
      cache_e[i]->free(); delete cache_e[i];
      delete [] cache_jwt[i];
    }
  }
  for (std::map<Key, Func<double>*, Compare>::iterator it = cache_fn.begin(); it != cache_fn.end(); it++)
  {
    (it->second)->free_fn(); delete (it->second);
  }
  cache_fn.clear();
}

//// evaluation of forms, general case ///////////////////////////////////////////////////////////

// Actual evaluation of volume Jacobian form (calculates integral)
scalar DiscreteProblem::eval_form(WeakForm::MatrixFormVol *mfv, Hermes::Tuple<Solution *> sln, 
                        PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv)
{
  // Determine the integration order.
  int inc = (fu->get_num_components() == 2) ? 1 : 0;
  
  // Order of solutions from the previous iteration level.
  AUTOLA_OR(Func<Ord>*, oi, wf->neq);
  //for (int i = 0; i < wf->neq; i++) oi[i] = init_fn_ord(sln[i]->get_fn_order() + inc);
  if (sln != Hermes::Tuple<Solution *>()) {
    for (int i = 0; i < wf->neq; i++) {
      if (sln[i] != NULL) oi[i] = init_fn_ord(sln[i]->get_fn_order() + inc);
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
  
  // Eval the form using the quadrature of the just calculated order.
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

  // Function values and values of external functions.
  AUTOLA_OR(Func<scalar>*, prev, wf->neq);
  //for (int i = 0; i < wf->neq; i++) prev[i]  = init_fn(sln[i], rv, order);
  if (sln != Hermes::Tuple<Solution *>()) {
    for (int i = 0; i < wf->neq; i++) {
      if (sln[i] != NULL) prev[i] = init_fn(sln[i], rv, order);
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
scalar DiscreteProblem::eval_form(WeakForm::VectorFormVol *vfv, Hermes::Tuple<Solution *> sln, PrecalcShapeset *fv, RefMap *rv)
{
  // Determine the integration order.
  int inc = (fv->get_num_components() == 2) ? 1 : 0;
  
  // Order of solutions from the previous iteration level.
  AUTOLA_OR(Func<Ord>*, oi, wf->neq);
  //for (int i = 0; i < wf->neq; i++) oi[i] = init_fn_ord(sln[i]->get_fn_order() + inc);
  if (sln != Hermes::Tuple<Solution *>()) {
    for (int i = 0; i < wf->neq; i++) {
      if (sln[i] != NULL) oi[i] = init_fn_ord(sln[i]->get_fn_order() + inc);
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

  // Eval the form using the quadrature of the just calculated order.
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

  // Function values and values of external functions.
  AUTOLA_OR(Func<scalar>*, prev, wf->neq);
  //for (int i = 0; i < wf->neq; i++) prev[i]  = init_fn(sln[i], rv, order);
  if (sln != Hermes::Tuple<Solution *>()) {
    for (int i = 0; i < wf->neq; i++) {
      if (sln[i] != NULL) prev[i]  = init_fn(sln[i], rv, order);
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

// Actual evaluation of surface Jacobian form (calculates integral)
scalar DiscreteProblem::eval_form(WeakForm::MatrixFormSurf *mfs, Hermes::Tuple<Solution *> sln, 
                        PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv, EdgePos* ep)
{
  // Determine the integration order.
  int inc = (fu->get_num_components() == 2) ? 1 : 0;
  
  // Order of solutions from the previous iteration level.
  AUTOLA_OR(Func<Ord>*, oi, wf->neq);
  //for (int i = 0; i < wf->neq; i++) oi[i] = init_fn_ord(sln[i]->get_fn_order() + inc);
  if (sln != Hermes::Tuple<Solution *>()) {
    for (int i = 0; i < wf->neq; i++) {
      if (sln[i] != NULL) oi[i] = init_fn_ord(sln[i]->get_edge_fn_order(ep->edge) + inc);
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
  
  // Eval the form using the quadrature of the just calculated order.
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

  // Function values and values of external functions.
  AUTOLA_OR(Func<scalar>*, prev, wf->neq);
  //for (int i = 0; i < wf->neq; i++) prev[i]  = init_fn(sln[i], rv, eo);
  if (sln != Hermes::Tuple<Solution *>()) {
    for (int i = 0; i < wf->neq; i++) {
      if (sln[i] != NULL) prev[i]  = init_fn(sln[i], rv, eo);
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
scalar DiscreteProblem::eval_form(WeakForm::VectorFormSurf *vfs, Hermes::Tuple<Solution *> sln, 
                        PrecalcShapeset *fv, RefMap *rv, EdgePos* ep)
{
  // Determine the integration order.
  int inc = (fv->get_num_components() == 2) ? 1 : 0;
  
  // Order of solutions from the previous iteration level.
  AUTOLA_OR(Func<Ord>*, oi, wf->neq);
  //for (int i = 0; i < wf->neq; i++) oi[i] = init_fn_ord(sln[i]->get_fn_order() + inc);
  if (sln != Hermes::Tuple<Solution *>()) {
    for (int i = 0; i < wf->neq; i++) {
      if (sln[i] != NULL) oi[i] = init_fn_ord(sln[i]->get_edge_fn_order(ep->edge) + inc);
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
  
  // Eval the form using the quadrature of the just calculated order.
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

  // Function values and values of external functions.
  AUTOLA_OR(Func<scalar>*, prev, wf->neq);
  //for (int i = 0; i < wf->neq; i++) prev[i]  = init_fn(sln[i], rv, eo);
  if (sln != Hermes::Tuple<Solution *>()) {
    for (int i = 0; i < wf->neq; i++) {
      if (sln[i] != NULL) prev[i]  = init_fn(sln[i], rv, eo);
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



//// solve /////////////////////////////////////////////////////////////////////////////////////////

bool DiscreteProblem::solve_matrix_problem(Matrix* mat, Vector* vec) 
{
  // check matrix size
  int ndof = this->get_num_dofs();
  if (ndof == 0) error("ndof = 0 in DiscreteProblem::solve().");
  if (ndof != mat->get_size())
    error("Matrix size does not match ndof in DiscreteProblem:solve().");
  if (ndof != vec->get_size())
    error("Vector size does not match ndof in DiscreteProblem:solve().");

  // FIXME: similar test should be done for the vector "vec" also, but we need
  // to access the information about its length.
  // ...

  // solve the matrix problem (and report time)
  TimePeriod cpu_time;
  bool flag = this->solver->solve(mat, vec);
  report_time("Matrix problem solved in %g s", cpu_time.tick().last());

  return flag;
}

bool DiscreteProblem::solve(Matrix* mat, Vector* rhs, Vector* vec)
{
  int ndof = this->get_num_dofs();

  // sanity checks
  if (mat == NULL) error("matrix is NULL in DiscreteProblem::solve().");
  if (rhs == NULL) error("rhs is NULL in DiscreteProblem::solve().");
  if (vec == NULL) error("vec is NULL in DiscreteProblem::solve().");
  if (ndof == 0) error("ndof = 0 in DiscreteProblem::solve().");
  if (ndof != mat->get_size())
    error("Matrix size does not match ndof in in DiscreteProblem:solve().");

  // copy "vec" into "delta" and solve the matrix problem with "mat", "delta"
  Vector* delta = new AVector(ndof);
  memcpy(delta->get_c_array(), rhs, sizeof(scalar) * ndof);
  bool flag = this->solve_matrix_problem(mat, delta);
  if (flag == false) return false;

  // add the result which is in "delta" to the previous 
  // solution vector which is in "vec"
  for (int i = 0; i < ndof; i++) vec->add(i, delta->get_c_array()[i]);
  delete delta;

  return true;
}

// L2 projections
template<typename Real, typename Scalar>
Scalar L2projection_biform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                           Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->val[i] * v->val[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar L2projection_liform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                           Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (ext->fn[0]->val[i] * v->val[i]);
  return result;
}

// H1 projections
template<typename Real, typename Scalar>
Scalar H1projection_biform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                           Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->val[i] * v->val[i] + u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar H1projection_liform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                           Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (ext->fn[0]->val[i] * v->val[i] + ext->fn[0]->dx[i] * 
                       v->dx[i] + ext->fn[0]->dy[i] * v->dy[i]);
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

int DiscreteProblem::assign_dofs()
{
  // sanity checks
  if (this->wf == NULL) error("this->wf = NULL in DiscreteProblem::assign_dofs().");

  // assigning dofs to each space
  if (this->spaces == Hermes::Tuple<Space *>()) error("this->spaces is empty in DiscreteProblem::assign_dofs().");
  int ndof = 0;
  for (int i = 0; i < this->wf->neq; i++) {
    if (this->spaces[i] == NULL) error("this->spaces[%d] is NULL in assign_dofs().", i);
    int inc = this->spaces[i]->assign_dofs(ndof);
    ndof += inc;
  }

  return ndof;
}

// Underlying function for global orthogonal projection.
// Not intended for the user. NOTE: the weak form here must be 
// a special projection weak form, which is different from 
// the weak form of the PDE. If you supply a weak form of the 
// PDE, the PDE will just be solved. 
void project_internal(Hermes::Tuple<Space *> spaces, WeakForm *wf, 
                      Hermes::Tuple<Solution*> target_slns, Vector* target_vec, bool is_complex)
{
  int n = spaces.size();

  // sanity checks
  if (n <= 0 || n > 10) error("Wrong number of projected functions in project_global().");
  for (int i = 0; i < n; i++) if(spaces[i] == NULL) error("this->spaces[%d] == NULL in project_global().", i);
  if (spaces.size() != n) error("Number of spaces must matchnumber of projected functions in project_global().");
  if (target_slns.size() != n && target_slns.size() != 0) 
    error("Mismatched numbers of projected functions and solutions in project_global().");

  // this is needed since spaces may have their DOFs enumerated only locally.
  int ndof = assign_dofs(spaces);

  // FIXME: enable other types of matrices and vectors.
  CooMatrix mat(ndof, is_complex);
  SolverSciPyUmfpack solver;
  Vector* dir = new AVector(ndof, is_complex);
  
  // Allocate the resulting coefficient vector.  
  Vector* rhs = new AVector(ndof, is_complex);

  //assembling the projection matrix, dir vector and rhs  
  DiscreteProblem dp(wf, spaces);
  bool rhsonly = false;
  // the NULL is there since no external solution vector is needed
  dp.assemble(NULL, &mat, dir, rhs, rhsonly, is_complex);
  // since this is a linear problem, subtract the dir vector from the right-hand side:
  if (is_complex)
      for (int i=0; i < ndof; i++) rhs->add(i, -dir->get_cplx(i));
  else
      for (int i=0; i < ndof; i++) rhs->add(i, -dir->get(i));

  // Calculate the Newton coefficient vector.
  solver.solve(&mat, rhs);

  // If the user wants the resulting Solutions.
  if (target_slns != Hermes::Tuple<Solution *>()) {
    for (int i=0; i < target_slns.size(); i++) {
      if (target_slns[i] != NULL) target_slns[i]->set_coeff_vector(spaces[i], rhs);
    }
  }

  // If the user wants the resulting coefficient vector
  // NOTE: this may change target_vector length.
  if (target_vec != NULL) {
    target_vec->free_data();
    if (target_vec->is_complex())
        target_vec->set_c_array_cplx(rhs->get_c_array_cplx(), rhs->get_size());
    else
        target_vec->set_c_array(rhs->get_c_array(), rhs->get_size());
    // NOTE: rhs must not be deleted.
  } else {
    delete rhs;
  }

  delete dir;
}

// global orthogonal projection
void project_global(Hermes::Tuple<Space *> spaces, Hermes::Tuple<int> proj_norms, Hermes::Tuple<MeshFunction*> source_meshfns, 
                    Hermes::Tuple<Solution*> target_slns, Vector* target_vec, bool is_complex)
{
  int n = spaces.size();  

  // define temporary projection weak form
  WeakForm* proj_wf = new WeakForm(n);
  int found[100];
  for (int i = 0; i < 100; i++) found[i] = 0;
  for (int i = 0; i < n; i++) {
    int norm;
    if (proj_norms == Hermes::Tuple<int>()) norm = 1;
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

  project_internal(spaces, proj_wf, target_slns, target_vec, is_complex);
}

void project_global(Hermes::Tuple<Space *> spaces, matrix_forms_tuple_t proj_biforms, 
                    vector_forms_tuple_t proj_liforms, Hermes::Tuple<MeshFunction*> source_meshfns, 
                    Hermes::Tuple<Solution*> target_slns, Vector* target_vec, bool is_complex)
{
  int n = spaces.size();
  matrix_forms_tuple_t::size_type n_biforms = proj_biforms.size();
  if (n_biforms == 0)
    error("Please use the simpler version of project_global with the argument Hermes::Tuple<int> proj_norms if you do not provide your own projection norm.");
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

  project_internal(spaces, &proj_wf, target_slns, target_vec, is_complex);
}

/// Global orthogonal projection of one vector-valued ExactFunction.
void project_global(Space *space, ExactFunction2 source_fn, Solution* target_sln, 
                    Vector* target_vec, bool is_complex)
{
  int proj_norm = 2; // Hcurl
  Mesh *mesh = space->get_mesh();
  if (mesh == NULL) error("Mesh is NULL in project_global().");
  Solution source_sln;
  source_sln.set_exact(mesh, source_fn);
  project_global(space, proj_norm, (MeshFunction*)&source_sln, target_sln, target_vec, is_complex);
};

/// Projection-based interpolation of an exact function. This is faster than the
/// global projection since no global matrix problem is solved.
void project_local(Space *space, int proj_norm, ExactFunction source_fn, Mesh* mesh,
                   Solution* target_sln, bool is_complex)
{
  /// TODO
}

int get_num_dofs(Hermes::Tuple<Space *> spaces)
{
  int ndof = 0;
  for (int i=0; i<spaces.size(); i++) {
    ndof += spaces[i]->get_num_dofs();
  }
  return ndof;
}

int DiscreteProblem::get_num_dofs()
{
  // sanity checks
  if (this->wf == NULL) error("this->wf is NULL in DiscreteProblem::get_num_dofs().");
  if (this->wf->neq == 0) error("this->wf->neq is 0 in DiscreteProblem::get_num_dofs().");
  if (this->spaces == Hermes::Tuple<Space *>()) error("this->spaces is empty in DiscreteProblem::get_num_dofs().");

  int ndof = 0;
  for (int i = 0; i < this->wf->neq; i++) {
    if (this->spaces[i] ==  NULL) error("this->spaces[%d] is NULL in DiscreteProblem::get_num_dofs().", i);
    ndof += this->get_num_dofs(i);
  }
  return ndof;
}

void DiscreteProblem::update_essential_bc_values()
{
  int n = this->wf->neq;
  for (int i=0; i<n; i++) this->spaces[i]->update_essential_bc_values();
}

void init_matrix_solver(MatrixSolverType matrix_solver, int ndof, 
                        Matrix* &mat, Vector* &rhs, 
                        Solver* &solver, bool is_complex) 
{
  // Initialize stiffness matrix, load vector, and matrix solver.
  // UMFpack.
  Matrix* mat_umfpack = new UMFPackMatrix();
  Vector* rhs_umfpack = new UMFPackVector();
  Solver* solver_umfpack = new UMFPackLinearSolver(mat_umfpack, rhs_umfpack);
  
  switch (matrix_solver) {
    case SOLVER_UMFPACK: 
      mat = mat_umfpack;
      rhs = rhs_umfpack;
      solver = solver_umfpack;
      break;
    case SOLVER_PETSC:  
      error("Petsc solver not implemented yet.");
      /*
      mat = &mat_petsc;
      rhs = &rhs_petsc;
      solver = &solver_petsc;
      */
      break;
    case SOLVER_MUMPS:  
      error("MUMPS solver not implemented yet.");
      /*
      mat = &mat_mumps;
      rhs = &rhs_mumps;
      solver = &solver_mumps;
      */
      break;
    default: error("Bad matrix solver in init_matrix_solver().");
  }
}

double get_l2_norm_real(Vector* vec) 
{
  double val = 0;
  for (int i = 0; i < vec->get_size(); i++) val += vec->get(i)*vec->get(i);
  val = sqrt(val);
  return val;
}

double get_l2_norm_cplx(Vector* vec) 
{
  double val = 0;
  for (int i = 0; i < vec->get_size(); i++) val += std::norm(vec->get_cplx(i));
  return sqrt(val);
}

// Basic Newton's method, takes a coefficient vector and returns a coefficient vector. 
bool solve_newton(Hermes::Tuple<Space *> spaces, WeakForm* wf, Vector* coeff_vec, 
                  MatrixSolverType matrix_solver, double newton_tol, 
                  int newton_max_iter, bool verbose, bool is_complex) 
{
  int ndof = get_num_dofs(spaces);
  
  // sanity checks
  if (coeff_vec == NULL) error("coeff_vec == NULL in solve_newton().");
  int n = spaces.size();
  if (spaces.size() != wf->neq) 
    error("The number of spaces in newton_solve() must match the number of equation in the PDE system.");
  for (int i=0; i < n; i++) {
    if (spaces[i] == NULL) error("spaces[%d] is NULL in solve_newton().", i);
  }
  if (coeff_vec->get_size() != ndof) error("Bad vector length in solve_newton().");

  // Initialize the discrete problem.
  DiscreteProblem dp(wf, spaces);
  //info("ndof = %d", dp.get_num_dofs());

  // Select matrix solver.
  Matrix* mat; Vector* rhs; Solver* solver;
  init_matrix_solver(matrix_solver, ndof, mat, rhs, solver, is_complex);

  int it = 1;
  while (1)
  {
    // Assemble the Jacobian matrix and residual vector.
    bool rhsonly = false;
    // the NULL stands for the dir vector which is not needed here
    dp.assemble(coeff_vec, mat, NULL, rhs, rhsonly, is_complex);

    // Multiply the residual vector with -1 since the matrix 
    // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
    if (is_complex)
      for (int i = 0; i < ndof; i++) rhs->set(i, -rhs->get_cplx(i));
    else
      for (int i = 0; i < ndof; i++) rhs->set(i, -rhs->get(i));

    // Calculate the l2-norm of residual vector.
    double res_l2_norm;
    if (!is_complex) res_l2_norm = get_l2_norm_real(rhs);
    else res_l2_norm = get_l2_norm_cplx(rhs);
    if (verbose) info("---- Newton iter %d, ndof %d, res. l2 norm %g", 
                        it, get_num_dofs(spaces), res_l2_norm);

    // If l2 norm of the residual vector is in tolerance, quit.
    if (res_l2_norm < newton_tol|| it > newton_max_iter) break;

    // Solve the matrix problem.
    if (!solver->solve(mat, rhs)) error ("Matrix solver failed.\n");

    // Add \deltaY^{n+1} to Y^n.
    if (is_complex)
      for (int i = 0; i < ndof; i++) coeff_vec->add(i, rhs->get_cplx(i));
    else
      for (int i = 0; i < ndof; i++) coeff_vec->add(i, rhs->get(i));
      
    
    it++;
  };

  delete rhs;
  delete mat;
  delete solver; // TODO: Create destructors for solvers.
  return (it <= newton_max_iter);
}

// Solves a typical nonlinear problem using the Newton's method and 
// automatic adaptivity. This function projects the init_meshfns
// to the coarse mesh(es), then solves the coarse mesh problem 
// using Newton, then projects the coarse mesh solution to the fine 
// mesh, and finally solves the fine mesh problem using Newton.
// So, this is not suitable for time-dependent problems.
// Feel free to adjust this function for more advanced applications.
HERMES_API bool solve_newton_adapt(Hermes::Tuple<Space *> spaces, WeakForm* wf, Vector *coeff_vec, 
                        MatrixSolverType matrix_solver, Hermes::Tuple<int>proj_norms, 
                        Hermes::Tuple<Solution *> target_slns, Hermes::Tuple<Solution *> target_ref_slns, 
                        Hermes::Tuple<WinGeom *> sln_win_geom, Hermes::Tuple<WinGeom *> mesh_win_geom, 
                        Hermes::Tuple<RefinementSelectors::Selector *> selectors, AdaptivityParamType* apt,
                        double newton_tol_coarse, double newton_tol_fine, int newton_max_iter, 
                        bool verbose, Hermes::Tuple<ExactSolution *> exact_slns, 
                        bool is_complex)
{
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

  // Number of degrees of freedom 
  int ndof = get_num_dofs(spaces);

  // Number of exact solutions given.
  if (exact_slns.size() != 0 && exact_slns.size() != num_comps) 
    error("Number of exact solutions does not match number of equations.");
  bool is_exact_solution;
  if (exact_slns.size() == num_comps) is_exact_solution = true;
  else is_exact_solution = false;

  // Initialize matrix solver.
  Matrix* mat; Vector *coeff_vec_dummy; Solver* solver;  
  init_matrix_solver(matrix_solver, ndof, mat, coeff_vec_dummy, solver, is_complex);
  if (coeff_vec->get_size() != coeff_vec_dummy->get_size()) 
    error("Mismatched matrix and vector size in solve_newton_adapt().");
  coeff_vec_dummy->free_data();

  // Initialize views.
  ScalarView* s_view[H2D_MAX_COMPONENTS];
  VectorView* v_view[H2D_MAX_COMPONENTS];
  OrderView*  o_view[H2D_MAX_COMPONENTS];
  for (int i = 0; i < num_comps; i++) {
    char* title = (char*)malloc(100*sizeof(char));
    if (sln_win_geom != Hermes::Tuple<WinGeom *>() && sln_win_geom[i] != NULL) {
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
    if (mesh_win_geom != Hermes::Tuple<WinGeom *>() && mesh_win_geom[i] != NULL) {
      if (num_comps == 1) sprintf(title, "Mesh", i); 
      else sprintf(title, "Mesh[%d]", i); 
      o_view[i] = new OrderView(title, mesh_win_geom[i]);
    }
    else o_view[i] = NULL;
  }

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est, graph_dof_exact, graph_cpu_exact;

  // Newton's loop on the coarse mesh.
  info("Solving on coarse mesh:");
  if (!solve_newton(spaces, wf, coeff_vec, matrix_solver, 
                    newton_tol_coarse, newton_max_iter, verbose, is_complex)) {
    error("Newton's method did not converge.");
  }

  // Temporary coarse and reference mesh solutions, and reference spaces.
  Hermes::Tuple<Solution *> slns = Hermes::Tuple<Solution *>();
  Hermes::Tuple<Solution *> ref_slns = Hermes::Tuple<Solution *>();
  for (int i = 0; i < num_comps; i++) {
    slns.push_back(new Solution);
    ref_slns.push_back(new Solution);
  }

  // Store the result in sln.
  for (int i = 0; i < num_comps; i++) slns[i]->set_coeff_vector(spaces[i], coeff_vec);
    
  // FIXME: this needs to be solved more elegantly.
  Hermes::Tuple<MeshFunction*> slns_mf = Hermes::Tuple<MeshFunction*>();
  Hermes::Tuple<MeshFunction*> ref_slns_mf = Hermes::Tuple<MeshFunction*>();
  for (int i = 0; i < num_comps; i++) {
    slns_mf.push_back((MeshFunction*)slns[i]);
    ref_slns_mf.push_back((MeshFunction*)ref_slns[i]);
  }

  // Adaptivity loop:
  bool done = false; int as = 1;
  double err_est;
  do {
    info("---- Adaptivity step %d:", as);

    // Time measurement..
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

    // Construct globally refined reference mesh(es)
    // and setup reference space(s).
    Hermes::Tuple<Space *> ref_spaces = Hermes::Tuple<Space *>();
    Hermes::Tuple<Mesh *> ref_meshes = Hermes::Tuple<Mesh *>();
    for (int i = 0; i < num_comps; i++) {
      ref_meshes.push_back(new Mesh());
      Mesh *ref_mesh = ref_meshes.back();
      ref_mesh->copy(spaces[i]->get_mesh());
      ref_mesh->refine_all_elements();
      ref_spaces.push_back(spaces[i]->dup(ref_mesh));
      int order_increase = 1;
      ref_spaces[i]->copy_orders(spaces[i], order_increase);
    }

    // Calculate initial coefficient vector for Newton on the fine mesh.
    if (as == 1) {
      info("Projecting coarse mesh solution to obtain initial vector on new fine mesh.");
      project_global(ref_spaces, proj_norms, slns_mf, Hermes::Tuple<Solution*>(), coeff_vec, is_complex);
    }
    else {
      info("Projecting previous fine mesh solution to obtain initial vector on new fine mesh.");
      project_global(ref_spaces, proj_norms, ref_slns_mf, Hermes::Tuple<Solution*>(), coeff_vec, is_complex);
    }

    // Newton's method on fine mesh
    info("Solving on fine mesh.");
    if (!solve_newton(ref_spaces, wf, coeff_vec, matrix_solver, newton_tol_fine, newton_max_iter, verbose, is_complex))
      error("Newton's method did not converge.");

    // Store the result in ref_sln.
    for (int i = 0; i < num_comps; i++) ref_slns[i]->set_coeff_vector(ref_spaces[i], coeff_vec);

    // Calculate element errors.
    if (verbose) info("Calculating error.");
    Adapt hp(spaces, proj_norms);
    // Pass special error forms if any.
    for (int k = 0; k < apt->error_form_val.size(); k++) {
      hp.set_error_form(apt->error_form_i[k], apt->error_form_j[k], 
                        apt->error_form_val[k], apt->error_form_ord[k]);
    }
    // Pass coarse mesh and reference solutions for error estimation.
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
          if (is_exact_solution == true) info("err_exact_rel[%d]: %g%%", i, err_exact_abs[i]/norm_exact[i]*100);
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

      if (get_num_dofs(spaces) >= ndof_stop) {
        done = true;
        break;
      }
      
      // Project last fine mesh solution on the new coarse mesh
      // to obtain new coarse mesh solution.
      if (verbose) info("Projecting reference solution on new coarse mesh.");
      // The NULL pointer means that we do not want the resulting coefficient vector.
      project_global(spaces, proj_norms, ref_slns_mf, slns, NULL, is_complex); 
    }

    for (int i = 0; i < num_comps; i++) {
      ref_meshes[i]->free();
      ref_spaces[i]->free();
    }

    as++;
  }
  while (done == false);

  // Obtain the coefficient vector on the coarse mesh.
  project_global(spaces, proj_norms, ref_slns_mf, NULL, coeff_vec, is_complex);

  // If the user wants, give him the coarse mesh solutions.
  if (target_slns != Hermes::Tuple<Solution *>()) {
    if (slns.size() != target_slns.size()) error("Mismatched solution Hermes::Tuple length in solve_newton_adapt().");
    for (int i = 0; i < num_comps; i++) {
      target_slns[i]->copy(slns[i]);
    }
  }

  // If the user wants, give him the reference mesh solutions.
  if (target_ref_slns != Hermes::Tuple<Solution *>()) {
    if (ref_slns.size() != target_ref_slns.size()) 
      error("Mismatched reference solution Hermes::Tuple length in solve_newton_adapt().");
    for (int i = 0; i < num_comps; i++) {
      target_ref_slns[i]->copy(ref_slns[i]);
    }
  }

  // Free temporary solutions.
  for (int i = 0; i < num_comps; i++) {
    slns[i]->free();
    ref_slns[i]->free();
  }

  if (verbose) info("Total running time: %g s", cpu_time.accumulated());
  return true;
}
