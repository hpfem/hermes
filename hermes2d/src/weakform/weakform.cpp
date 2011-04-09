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

#include "../h2d_common.h"
#include "weakform.h"
#include "../../../hermes_common/matrix.h"
#include "../function/forms.h"

//// interface /////////////////////////////////////////////////////////////////////////////////////

WeakForm::Form::Form(std::string area, Hermes::vector<MeshFunction *> ext, Hermes::vector<scalar> param, 
                     double scaling_factor, int u_ext_offset) :
  area(area), ext(ext), param(param), scaling_factor(scaling_factor), u_ext_offset(u_ext_offset)
{
  adapt_eval = false;
  stage_time = 0.0;
}

void WeakForm::Form::set_current_stage_time(double time)
{
  stage_time = time;
}

double WeakForm::Form::get_current_stage_time() const
{
  return stage_time;
}

scalar WeakForm::MatrixFormVol::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                                      Geom<double> *e, ExtData<scalar> *ext) const
{
  error("WeakForm::MatrixFormVol::value must be overridden.");
  return 0.0;
}

Ord WeakForm::MatrixFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                 Geom<Ord> *e, ExtData<Ord> *ext) const
{
  error("WeakForm::MatrixFormVol::ord must be overridden.");
  return Ord();
}

WeakForm::MatrixFormVol* WeakForm::MatrixFormVol::clone()
{
  error("WeakForm::MatrixFormVol::clone() must be overridden.");
  return NULL;
}

WeakForm::MatrixFormVol::MatrixFormVol(unsigned int i, unsigned int j, SymFlag sym,
                                       std::string area, Hermes::vector<MeshFunction *> ext, 
                                       Hermes::vector<scalar> param, double scaling_factor, int u_ext_offset) : 
  Form(area, ext, param, scaling_factor, u_ext_offset), i(i), j(j), sym(sym)
{
}

scalar WeakForm::MatrixFormSurf::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                                       Geom<double> *e, ExtData<scalar> *ext) const
{
  error("WeakForm::MatrixFormSurf::value must be overridden.");
  return 0.0;
}

Ord WeakForm::MatrixFormSurf::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                  Geom<Ord> *e, ExtData<Ord> *ext) const
{
  error("WeakForm::MatrixFormSurf::ord must be overridden.");
  return Ord();
}

WeakForm::MatrixFormSurf* WeakForm::MatrixFormSurf::clone()
{
  error("WeakForm::MatrixFormSurf::clone() must be overridden.");
  return NULL;
}

WeakForm::MatrixFormSurf::MatrixFormSurf(unsigned int i, unsigned int j, std::string area,
                                         Hermes::vector<MeshFunction *> ext, 
                                         Hermes::vector<scalar> param, double scaling_factor, int u_ext_offset) : 
  Form(area, ext, param, scaling_factor, u_ext_offset), i(i), j(j)
{
}

scalar WeakForm::VectorFormVol::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                                      Geom<double> *e, ExtData<scalar> *ext) const
{
  error("WeakForm::VectorFormVol::value must be overridden.");
  return 0.0;
}

Ord WeakForm::VectorFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                                 Geom<Ord> *e, ExtData<Ord> *ext) const
{
  error("WeakForm::VectorFormVol::ord must be overridden.");
  return Ord();
}

WeakForm::VectorFormVol* WeakForm::VectorFormVol::clone()
{
  error("WeakForm::VectorFormVol::clone() must be overridden.");
  return NULL;
}

WeakForm::VectorFormVol::VectorFormVol(unsigned int i, std::string area,
                                       Hermes::vector<MeshFunction *> ext, 
                                       Hermes::vector<scalar> param, double scaling_factor, int u_ext_offset) : 
  Form(area, ext, param, scaling_factor, u_ext_offset), i(i)
{
}

WeakForm::VectorFormSurf::VectorFormSurf(unsigned int i, std::string area,
                                         Hermes::vector<MeshFunction *> ext, 
                                         Hermes::vector<scalar> param, double scaling_factor, int u_ext_offset) : 
  Form(area, ext, param, scaling_factor, u_ext_offset), i(i)
{
}

scalar WeakForm::VectorFormSurf::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                                       Geom<double> *e, ExtData<scalar> *ext) const
{
  error("WeakForm::VectorFormSurf::value must be overridden.");
  return 0.0;
}

Ord WeakForm::VectorFormSurf::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                                  Geom<Ord> *e, ExtData<Ord> *ext) const
{
  error("WeakForm::VectorFormSurf::ord must be overridden.");
  return Ord();
}

WeakForm::VectorFormSurf* WeakForm::VectorFormSurf::clone()
{
  error("WeakForm::VectorFormSurf::clone() must be overridden.");
  return NULL;
}

// Multi component.
WeakForm::MultiComponentMatrixFormVol::MultiComponentMatrixFormVol(Hermes::vector<std::pair<unsigned int, unsigned int> >coordinates, SymFlag sym,
                                       std::string area, Hermes::vector<MeshFunction *> ext, 
                                       Hermes::vector<scalar> param, double scaling_factor, int u_ext_offset) : 
  Form(area, ext, param, scaling_factor, u_ext_offset), coordinates(coordinates), sym(sym)
{
}

WeakForm::MultiComponentMatrixFormVol* WeakForm::MultiComponentMatrixFormVol::clone()
{
  error("WeakForm::MultiComponentMatrixFormVol::clone() must be overridden.");
  return NULL;
}

WeakForm::MultiComponentMatrixFormSurf::MultiComponentMatrixFormSurf(Hermes::vector<std::pair<unsigned int, unsigned int> >coordinates, std::string area,
                                         Hermes::vector<MeshFunction *> ext, 
                                         Hermes::vector<scalar> param, double scaling_factor, int u_ext_offset) : 
  Form(area, ext, param, scaling_factor, u_ext_offset), coordinates(coordinates)
{
}

WeakForm::MultiComponentMatrixFormSurf* WeakForm::MultiComponentMatrixFormSurf::clone()
{
  error("WeakForm::MatrixFormSurf::clone() must be overridden.");
  return NULL;
}

WeakForm::MultiComponentVectorFormVol::MultiComponentVectorFormVol(Hermes::vector<unsigned int> coordinates, std::string area,
                                       Hermes::vector<MeshFunction *> ext, 
                                       Hermes::vector<scalar> param, double scaling_factor, int u_ext_offset) : 
  Form(area, ext, param, scaling_factor, u_ext_offset), coordinates(coordinates)
{
}
  
WeakForm::MultiComponentVectorFormVol* WeakForm::MultiComponentVectorFormVol::clone()
{
  error("WeakForm::VectorFormVol::clone() must be overridden.");
  return NULL;
}

WeakForm::MultiComponentVectorFormSurf::MultiComponentVectorFormSurf(Hermes::vector<unsigned int> coordinates, std::string area,
                                         Hermes::vector<MeshFunction *> ext, 
                                         Hermes::vector<scalar> param, double scaling_factor, int u_ext_offset) : 
  Form(area, ext, param, scaling_factor, u_ext_offset), coordinates(coordinates)
{
}

WeakForm::MultiComponentVectorFormSurf* WeakForm::MultiComponentVectorFormSurf::clone()
{
  error("WeakForm::VectorFormVol::clone() must be overridden.");
  return NULL;
}

WeakForm::WeakForm(unsigned int neq, bool mat_free)
{
  _F_

  this->neq = neq;
  this->seq = 0;
  this->is_matfree = mat_free;
}

WeakForm::~WeakForm()
{
}

void WeakForm::add_matrix_form(MatrixFormVol* form)
{
  _F_

  if (form->i >= neq || form->j >= neq)
    error("Invalid equation number.");
  if (form->sym < -1 || form->sym > 1)
    error("\"sym\" must be -1, 0 or 1.");
  if (form->sym < 0 && form->i == form->j)
    error("Only off-diagonal forms can be antisymmetric.");
  if (mfvol.size() > 100) {
    warn("Large number of forms (> 100). Is this the intent?");
  }

  form->set_weakform(this);
  mfvol.push_back(form);
  seq++;
}

void WeakForm::add_matrix_form_surf(MatrixFormSurf* form)
{
  _F_
  if (form->i >= neq || form->j >= neq)
    error("Invalid equation number.");

  form->set_weakform(this);
  mfsurf.push_back(form);
  seq++;
}

void WeakForm::add_vector_form(VectorFormVol* form)
{
  _F_
  if (form->i >= neq)
    error("Invalid equation number.");
  form->set_weakform(this);
  vfvol.push_back(form);
  seq++;
}

void WeakForm::add_vector_form_surf(VectorFormSurf* form)
{
  _F_
  if (form->i >= neq)
    error("Invalid equation number.");

  form->set_weakform(this);
  vfsurf.push_back(form);
  seq++;
}

void WeakForm::add_multicomponent_matrix_form(MultiComponentMatrixFormVol* form)
{
  _F_

  for(unsigned int form_i = 0; form_i < form->coordinates.size(); form_i++) {
    if(form->coordinates.at(form_i).first >= neq || form->coordinates.at(form_i).second >= neq)
      error("Invalid equation number.");
    if (form->sym < 0 && form->coordinates.at(form_i).first == form->coordinates.at(form_i).second)
      error("Only off-diagonal forms can be antisymmetric.");
  }
  if (form->sym < -1 || form->sym > 1)
    error("\"sym\" must be -1, 0 or 1.");

  if (mfvol_mc.size() > 100)
    warn("Large number of forms (> 100). Is this the intent?");

  form->set_weakform(this);
  
  mfvol_mc.push_back(form);
  seq++;
}

void WeakForm::add_multicomponent_matrix_form_surf(MultiComponentMatrixFormSurf* form)
{
  _F_
  for(unsigned int form_i = 0; form_i < form->coordinates.size(); form_i++)
    if(form->coordinates.at(form_i).first >= neq || form->coordinates.at(form_i).second >= neq)
      error("Invalid equation number.");

  form->set_weakform(this);
  mfsurf_mc.push_back(form);
  seq++;
}

void WeakForm::add_multicomponent_vector_form(MultiComponentVectorFormVol* form)
{
  _F_
  for(unsigned int form_i = 0; form_i < form->coordinates.size(); form_i++)
    if(form->coordinates.at(form_i) >= neq)
      error("Invalid equation number.");
  form->set_weakform(this);
  vfvol_mc.push_back(form);
  seq++;
}

void WeakForm::add_multicomponent_vector_form_surf(MultiComponentVectorFormSurf* form)
{
  _F_
  for(unsigned int form_i = 0; form_i < form->coordinates.size(); form_i++)
    if(form->coordinates.at(form_i) >= neq)
      error("Invalid equation number.");

  form->set_weakform(this);
  vfsurf_mc.push_back(form);
  seq++;
}

void WeakForm::set_ext_fns(void* fn, Hermes::vector<MeshFunction*> ext)
{
  _F_
  error("Not implemented yet.");
}


//// stages ////////////////////////////////////////////////////////////////////////////////////////

/// Constructs a list of assembling stages. Each stage contains a list of forms
/// that share the same meshes. Each stage is then assembled separately. This
/// improves the performance of multi-mesh assembling.
/// This function is identical in H2D and H3D.
///
void WeakForm::get_stages(Hermes::vector<Space *> spaces, Hermes::vector<Solution *>& u_ext,
                          std::vector<WeakForm::Stage>& stages, bool want_matrix, bool want_vector)
{
  _F_

  if (!want_matrix && !want_vector) return;

  unsigned int i;
  stages.clear();

  if (want_matrix || want_vector) {    // This is because of linear problems where 
                                       // matrix terms with the Dirichlet lift go to rhs.
    // Process volume matrix forms.
    for (i = 0; i < mfvol.size(); i++)
    {
      unsigned int ii = mfvol[i]->i, jj = mfvol[i]->j;
      Mesh* m1 = spaces[ii]->get_mesh();
      Mesh* m2 = spaces[jj]->get_mesh();
      Stage* s = find_stage(stages, ii, jj, m1, m2, mfvol[i]->ext, u_ext);
      s->mfvol.push_back(mfvol[i]);
    }

    // Process surface matrix forms.
    for (i = 0; i < mfsurf.size(); i++)
    {
      unsigned int ii = mfsurf[i]->i, jj = mfsurf[i]->j;
      Mesh* m1 = spaces[ii]->get_mesh();
      Mesh* m2 = spaces[jj]->get_mesh();
      Stage* s = find_stage(stages, ii, jj, m1, m2, mfsurf[i]->ext, u_ext);
      s->mfsurf.push_back(mfsurf[i]);
    }

    // Multi component forms.
    for (unsigned i = 0; i < mfvol_mc.size(); i++) {
      Mesh* the_one_mesh = spaces[mfvol_mc.at(i)->coordinates.at(0).first]->get_mesh();
      for(unsigned int form_i = 0; form_i < mfvol_mc.at(i)->coordinates.size(); form_i++) {
        if(spaces[mfvol_mc.at(i)->coordinates.at(form_i).first]->get_mesh()->get_seq() != the_one_mesh->get_seq())
          error("When using multi-component forms, the Meshes have to be identical.");
        if(spaces[mfvol_mc.at(i)->coordinates.at(form_i).second]->get_mesh()->get_seq() != the_one_mesh->get_seq())
          error("When using multi-component forms, the Meshes have to be identical.");
      }
    
      Stage* s = find_stage(stages, mfvol_mc.at(i)->coordinates, the_one_mesh, the_one_mesh, mfvol_mc[i]->ext, u_ext);
      s->mfvol_mc.push_back(mfvol_mc[i]);
    }
    for (unsigned i = 0; i < mfsurf_mc.size(); i++) {
      Mesh* the_one_mesh = spaces[mfsurf_mc.at(i)->coordinates.at(0).first]->get_mesh();
      for(unsigned int form_i = 0; form_i < mfsurf_mc.at(i)->coordinates.size(); form_i++) {
        if(spaces[mfsurf_mc.at(i)->coordinates.at(form_i).first]->get_mesh()->get_seq() != the_one_mesh->get_seq())
          error("When using multi-component forms, the Meshes have to be identical.");
        if(spaces[mfsurf_mc.at(i)->coordinates.at(form_i).second]->get_mesh()->get_seq() != the_one_mesh->get_seq())
          error("When using multi-component forms, the Meshes have to be identical.");
      }
    
      Stage* s = find_stage(stages, mfsurf_mc.at(i)->coordinates, the_one_mesh, the_one_mesh, mfsurf_mc[i]->ext, u_ext);
      s->mfsurf_mc.push_back(mfsurf_mc[i]);
    }
  }

  if (want_vector) {
    // Process volume vector forms.
    for (unsigned i = 0; i < vfvol.size(); i++) {
      unsigned int ii = vfvol[i]->i;
      Mesh *m = spaces[ii]->get_mesh();
      Stage *s = find_stage(stages, ii, ii, m, m, vfvol[i]->ext, u_ext);
      s->vfvol.push_back(vfvol[i]);
    }

    // Process surface vector forms.
    for (unsigned i = 0; i < vfsurf.size(); i++) {
      unsigned int ii = vfsurf[i]->i;
      Mesh *m = spaces[ii]->get_mesh();
      Stage *s = find_stage(stages, ii, ii, m, m, vfsurf[i]->ext, u_ext);
      s->vfsurf.push_back(vfsurf[i]);
    }

    // Multi component forms.
    for (unsigned i = 0; i < vfvol_mc.size(); i++) {
      Mesh* the_one_mesh = spaces[vfvol_mc.at(i)->coordinates.at(0)]->get_mesh();
      for(unsigned int form_i = 0; form_i < vfvol_mc.at(i)->coordinates.size(); form_i++)
        if(spaces[vfvol_mc.at(i)->coordinates.at(form_i)]->get_mesh()->get_seq() != the_one_mesh->get_seq())
          error("When using multi-component forms, the Meshes have to be identical.");
    
      Stage *s = find_stage(stages, vfvol_mc.at(i)->coordinates, the_one_mesh, the_one_mesh, vfvol_mc[i]->ext, u_ext);
      s->vfvol_mc.push_back(vfvol_mc[i]);
    }
    for (unsigned i = 0; i < vfsurf_mc.size(); i++) {
      Mesh* the_one_mesh = spaces[vfsurf_mc.at(i)->coordinates.at(0)]->get_mesh();
      for(unsigned int form_i = 0; form_i < vfsurf_mc.at(i)->coordinates.size(); form_i++)
        if(spaces[vfsurf_mc.at(i)->coordinates.at(form_i)]->get_mesh()->get_seq() != the_one_mesh->get_seq())
          error("When using multi-component forms, the Meshes have to be identical.");
    
      Stage *s = find_stage(stages, vfsurf_mc.at(i)->coordinates, the_one_mesh, the_one_mesh, vfsurf_mc[i]->ext, u_ext);
      s->vfsurf_mc.push_back(vfsurf_mc[i]);
    }
  }

  // Helper macro for iterating in a set,
#define set_for_each(myset, type) \
  for (std::set<type>::iterator it = (myset).begin(); it != (myset).end(); it++)

  // Initialize the arrays meshes and fns needed by Traverse for each stage.
  for (i = 0; i < stages.size(); i++)
  {
    Stage* s = &stages[i];

    // First, initialize arrays for the test functions. A pointer to the PrecalcShapeset
    // corresponding to each space will be assigned to s->fns later during assembling.
    set_for_each(s->idx_set, int)
    {
      s->idx.push_back(*it);
      s->meshes.push_back(spaces[*it]->get_mesh());
      s->fns.push_back(NULL);
    }

    // Next, append to the existing arrays the external functions (including the solutions
    // from previous Newton iteration) and their meshes. Also fill in a special array with
    // these external functions only.
    set_for_each(s->ext_set, MeshFunction*)
    {
      s->ext.push_back(*it);
      s->meshes.push_back((*it)->get_mesh());
      s->fns.push_back(*it);
    }

    s->idx_set.clear();
    s->seq_set.clear();
    s->ext_set.clear();
  }
}


/// Finds an assembling stage with the same set of meshes as [m1, m2, ext, u_ext]. If no such
/// stage can be found, a new one is created and returned.
/// This function is the same in H2D and H3D.
///
WeakForm::Stage* WeakForm::find_stage(std::vector<WeakForm::Stage>& stages, int ii, int jj,
                                      Mesh* m1, Mesh* m2,
                                      Hermes::vector<MeshFunction*>& ext, Hermes::vector<Solution*>& u_ext)
{
  _F_
  // first create a list of meshes the form uses
  std::set<unsigned> seq;
  seq.insert(m1->get_seq());
  seq.insert(m2->get_seq());
  Mesh *mmm;
  for (unsigned i = 0; i < ext.size(); i++) {
    mmm = ext[i]->get_mesh();
    if (mmm == NULL) error("NULL Mesh pointer detected in ExtData during assembling.\n  Have you initialized all external functions?");
    seq.insert(mmm->get_seq());
  }
  for (unsigned i = 0; i < u_ext.size(); i++) {
    if (u_ext[i] != NULL) {
      mmm = u_ext[i]->get_mesh();
      if (mmm == NULL) error("NULL Mesh pointer detected in u_ext during assembling.");
      seq.insert(mmm->get_seq());
    }
  }

  // find a suitable existing stage for the form
  Stage* s = NULL;
  for (unsigned i = 0; i < stages.size(); i++)
    if (seq.size() == stages[i].seq_set.size() &&
        equal(seq.begin(), seq.end(), stages[i].seq_set.begin())) {
      s = &stages[i]; break;
    }

  // create a new stage if not found
  if (s == NULL)
  {
    Stage newstage;
    stages.push_back(newstage);
    s = &stages.back();
    s->seq_set = seq;
  }

  // update and return the stage
  for (unsigned int i = 0; i < ext.size(); i++)
    s->ext_set.insert(ext[i]);
  for (unsigned int i = 0; i < u_ext.size(); i++)
    if (u_ext[i] != NULL)
      s->ext_set.insert(u_ext[i]);

  s->idx_set.insert(ii);
  s->idx_set.insert(jj);
  return s;
}

WeakForm::Stage* WeakForm::find_stage(std::vector<WeakForm::Stage>& stages, Hermes::vector<std::pair<unsigned int, unsigned int> > coordinates,
                                      Mesh* m1, Mesh* m2,
                                      Hermes::vector<MeshFunction*>& ext, Hermes::vector<Solution*>& u_ext)
{
  _F_
  // first create a list of meshes the form uses
  std::set<unsigned> seq;
  seq.insert(m1->get_seq());
  seq.insert(m2->get_seq());
  Mesh *mmm;
  for (unsigned i = 0; i < ext.size(); i++) {
    mmm = ext[i]->get_mesh();
    if (mmm == NULL) error("NULL Mesh pointer detected in ExtData during assembling.\n  Have you initialized all external functions?");
    seq.insert(mmm->get_seq());
  }
  for (unsigned i = 0; i < u_ext.size(); i++) {
    if (u_ext[i] != NULL) {
      mmm = u_ext[i]->get_mesh();
      if (mmm == NULL) error("NULL Mesh pointer detected in u_ext during assembling.");
      seq.insert(mmm->get_seq());
    }
  }

  // find a suitable existing stage for the form
  Stage* s = NULL;
  for (unsigned i = 0; i < stages.size(); i++)
    if (seq.size() == stages[i].seq_set.size() &&
        equal(seq.begin(), seq.end(), stages[i].seq_set.begin())) {
      s = &stages[i]; break;
    }

  // create a new stage if not found
  if (s == NULL)
  {
    Stage newstage;
    stages.push_back(newstage);
    s = &stages.back();
    s->seq_set = seq;
  }

  // update and return the stage
  for (unsigned int i = 0; i < ext.size(); i++)
    s->ext_set.insert(ext[i]);
  for (unsigned int i = 0; i < u_ext.size(); i++)
    if (u_ext[i] != NULL)
      s->ext_set.insert(u_ext[i]);

  for(unsigned int ii = 0; ii < coordinates.size(); ii++) {
    s->idx_set.insert(coordinates.at(ii).first);
    s->idx_set.insert(coordinates.at(ii).second);
  }
  return s;
}

WeakForm::Stage* WeakForm::find_stage(std::vector<WeakForm::Stage>& stages, Hermes::vector<unsigned int> coordinates,
                                      Mesh* m1, Mesh* m2,
                                      Hermes::vector<MeshFunction*>& ext, Hermes::vector<Solution*>& u_ext)
{
  _F_
  // first create a list of meshes the form uses
  std::set<unsigned> seq;
  seq.insert(m1->get_seq());
  seq.insert(m2->get_seq());
  Mesh *mmm;
  for (unsigned i = 0; i < ext.size(); i++) {
    mmm = ext[i]->get_mesh();
    if (mmm == NULL) error("NULL Mesh pointer detected in ExtData during assembling.\n  Have you initialized all external functions?");
    seq.insert(mmm->get_seq());
  }
  for (unsigned i = 0; i < u_ext.size(); i++) {
    if (u_ext[i] != NULL) {
      mmm = u_ext[i]->get_mesh();
      if (mmm == NULL) error("NULL Mesh pointer detected in u_ext during assembling.");
      seq.insert(mmm->get_seq());
    }
  }

  // find a suitable existing stage for the form
  Stage* s = NULL;
  for (unsigned i = 0; i < stages.size(); i++)
    if (seq.size() == stages[i].seq_set.size() &&
        equal(seq.begin(), seq.end(), stages[i].seq_set.begin())) {
      s = &stages[i]; break;
    }

  // create a new stage if not found
  if (s == NULL)
  {
    Stage newstage;
    stages.push_back(newstage);
    s = &stages.back();
    s->seq_set = seq;
  }

  // update and return the stage
  for (unsigned int i = 0; i < ext.size(); i++)
    s->ext_set.insert(ext[i]);
  for (unsigned int i = 0; i < u_ext.size(); i++)
    if (u_ext[i] != NULL)
      s->ext_set.insert(u_ext[i]);

  for(unsigned int ii = 0; ii < coordinates.size(); ii++)
    s->idx_set.insert(coordinates.at(ii));
  
  return s;
}

/// Returns a (neq x neq) array containing true in each element, if the corresponding
/// block of weak forms is used, and false otherwise.
/// This function is the same in H2D and H3D.
///
bool** WeakForm::get_blocks(bool force_diagonal_blocks)
{
  _F_
  bool** blocks = new_matrix<bool>(neq, neq);
  for (unsigned int i = 0; i < neq; i++) {
    for (unsigned int j = 0; j < neq; j++)
      blocks[i][j] = false;
    if (force_diagonal_blocks == true)
      blocks[i][i] = true;
  }
  for (unsigned i = 0; i < mfvol.size(); i++) {
    if (fabs(mfvol[i]->scaling_factor) > 1e-12)
      blocks[mfvol[i]->i][mfvol[i]->j] = true;
    if (mfvol[i]->sym)
      if (fabs(mfvol[i]->scaling_factor) > 1e-12)
        blocks[mfvol[i]->j][mfvol[i]->i] = true;
  }

  for (unsigned i = 0; i < mfvol_mc.size(); i++) {
    if (fabs(mfvol_mc[i]->scaling_factor) > 1e-12)
      for(unsigned int component_i = 0; component_i < mfvol_mc[i]->coordinates.size(); component_i++)
        blocks[mfvol_mc[i]->coordinates[component_i].first][mfvol_mc[i]->coordinates[component_i].second] = true;
    if (mfvol_mc[i]->sym)
      if (fabs(mfvol_mc[i]->scaling_factor) > 1e-12)
        for(unsigned int component_i = 0; component_i < mfvol_mc[i]->coordinates.size(); component_i++)
          blocks[mfvol_mc[i]->coordinates[component_i].second][mfvol_mc[i]->coordinates[component_i].first] = true;
  }

  for (unsigned i = 0; i < mfsurf.size(); i++)
    if (fabs(mfsurf[i]->scaling_factor) > 1e-12)
      blocks[mfsurf[i]->i][mfsurf[i]->j] = true;

  for (unsigned i = 0; i < mfsurf_mc.size(); i++)
    if (fabs(mfsurf_mc[i]->scaling_factor) > 1e-12)
      for(unsigned int component_i = 0; component_i < mfsurf_mc[i]->coordinates.size(); component_i++)
        blocks[mfsurf_mc[i]->coordinates[component_i].first][mfsurf_mc[i]->coordinates[component_i].second] = true;
  
  return blocks;
}

void WeakForm::set_current_time(double time)
{
  current_time = time;
}

double WeakForm::get_current_time() const
{
  return current_time;
}
