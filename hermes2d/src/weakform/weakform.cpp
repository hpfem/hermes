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
#include "forms.h"

//// interface /////////////////////////////////////////////////////////////////////////////////////

WeakForm::Form::Form(int area, Hermes::vector<MeshFunction *> ext, double scaling_factor, int u_ext_offset) :
  area(area), ext(ext), scaling_factor(scaling_factor), u_ext_offset(u_ext_offset)
{
}

scalar WeakForm::MatrixFormVol::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                                      Geom<double> *e, ExtData<scalar> *ext)
{
  error("WeakForm::MatrixFormVol::value must be overrided.");
  return 0.0;
}

Ord WeakForm::MatrixFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                 Geom<Ord> *e, ExtData<Ord> *ext)
{
  error("WeakForm::MatrixFormVol::ord must be overrided.");
  return Ord();
}

WeakForm::MatrixFormVol::MatrixFormVol(unsigned int i, unsigned int j, SymFlag sym,
                                       int area, Hermes::vector<MeshFunction *> ext, double scaling_factor, int u_ext_offset) : Form(area, ext, scaling_factor, u_ext_offset), i(i), j(j), sym(sym)
{
}

scalar WeakForm::MatrixFormSurf::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                                       Geom<double> *e, ExtData<scalar> *ext)
{
  error("WeakForm::MatrixFormSurf::value must be overrided.");
  return 0.0;
}

Ord WeakForm::MatrixFormSurf::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                  Geom<Ord> *e, ExtData<Ord> *ext)
{
  error("WeakForm::MatrixFormSurf::ord must be overrided.");
  return Ord();
}

WeakForm::MatrixFormSurf::MatrixFormSurf(unsigned int i, unsigned int j, int area,
                                         Hermes::vector<MeshFunction *> ext, double scaling_factor, int u_ext_offset) : Form(area, ext, scaling_factor, u_ext_offset), i(i), j(j)
{
}

scalar WeakForm::VectorFormVol::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext)
{
  error("WeakForm::VectorFormVol::value must be overrided.");
  return 0.0;
}

Ord WeakForm::VectorFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  error("WeakForm::VectorFormVol::ord must be overrided.");
  return Ord();
}

WeakForm::VectorFormVol::VectorFormVol(unsigned int i, int area,
                                       Hermes::vector<MeshFunction *> ext, double scaling_factor, int u_ext_offset) : Form(area, ext, scaling_factor, u_ext_offset), i(i)
{
}

WeakForm::VectorFormSurf::VectorFormSurf(unsigned int i, int area,
                                         Hermes::vector<MeshFunction *> ext, double scaling_factor, int u_ext_offset) : Form(area, ext, scaling_factor, u_ext_offset), i(i)
{
}

scalar WeakForm::VectorFormSurf::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext)
{
  error("WeakForm::VectorFormSurf::value must be overrided.");
  return 0.0;
}

Ord WeakForm::VectorFormSurf::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  error("WeakForm::VectorFormSurf::ord must be overrided.");
  return Ord();
}

<<<<<<< HEAD

=======
>>>>>>> Fixes & comments & old code deletion.
WeakForm::WeakForm(unsigned int neq, bool mat_free)
{
  _F_
  this->neq = neq;
  seq = 0;
  this->is_matfree = mat_free;
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
  if (form->area != HERMES_ANY && form->area < 0 && (unsigned) (-form->area) > areas.size())
    error("Invalid area number.");
  if (mfvol.size() > 100) {
    warn("Large number of forms (> 100). Is this the intent?");
  }

  form->set_weakform(this);
  mfvol.push_back(form);
  seq++;
}

void WeakForm::add_matrix_form_surf(MatrixFormSurf* form)
{
<<<<<<< HEAD
  _F_
      if (form->i >= neq || form->j >= neq)
      error("Invalid equation number.");
  if (form->area != HERMES_ANY && form->area != H2D_DG_BOUNDARY_EDGE && form->area !=
      H2D_DG_INNER_EDGE && form->area < 0 && (unsigned) (-form->area) > areas.size())
    error("Invalid area number.");

  form->set_weakform(this);
  mfsurf.push_back(form);
  seq++;
=======
    _F_
    if (form->i >= neq || form->j >= neq)
      error("Invalid equation number.");
    if (form->area != HERMES_ANY && form->area != H2D_DG_BOUNDARY_EDGE && form->area !=
    H2D_DG_INNER_EDGE && form->area < 0 && (unsigned) (-form->area) > areas.size())
      error("Invalid area number.");

    form->set_weakform(this);
    mfsurf.push_back(form);
    seq++;
>>>>>>> Fixes & comments & old code deletion.
}

void WeakForm::add_vector_form(VectorFormVol* form)
{
  _F_
  if (form->i >= neq)
    error("Invalid equation number.");
  if (form->area != HERMES_ANY && form->area < 0 && (unsigned) (-form->area) > areas.size())
    error("Invalid area number.");

  form->set_weakform(this);
  vfvol.push_back(form);
  seq++;
}

void WeakForm::add_vector_form_surf(VectorFormSurf* form)
{
  _F_
  if (form->i >= neq)
    error("Invalid equation number.");
  if (form->area != HERMES_ANY && form->area != H2D_DG_BOUNDARY_EDGE && form->area !=
  H2D_DG_INNER_EDGE && form->area < 0 && (unsigned) (-form->area) > areas.size())
    error("Invalid area number.");

  form->set_weakform(this);
  vfsurf.push_back(form);
  seq++;
}

void WeakForm::set_ext_fns(void* fn, Hermes::vector<MeshFunction*>ext)
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
          std::vector<WeakForm::Stage>& stages, bool rhsonly)
{
  _F_
  unsigned int i;
  stages.clear();

  // process volume matrix forms
  for (i = 0; i < mfvol.size(); i++)
  {
    unsigned int ii = mfvol[i]->i, jj = mfvol[i]->j;
    Mesh* m1 = spaces[ii]->get_mesh();
    Mesh* m2 = spaces[jj]->get_mesh();
    Stage* s = find_stage(stages, ii, jj, m1, m2, mfvol[i]->ext, u_ext);
    s->mfvol.push_back(mfvol[i]);
  }

  // process surface matrix forms
  for (i = 0; i < mfsurf.size(); i++)
  {
    unsigned int ii = mfsurf[i]->i, jj = mfsurf[i]->j;
    Mesh* m1 = spaces[ii]->get_mesh();
    Mesh* m2 = spaces[jj]->get_mesh();
    Stage* s = find_stage(stages, ii, jj, m1, m2, mfsurf[i]->ext, u_ext);
    s->mfsurf.push_back(mfsurf[i]);
  }

  // process volume vector forms
  for (unsigned i = 0; i < vfvol.size(); i++) {
    unsigned int ii = vfvol[i]->i;
    Mesh *m = spaces[ii]->get_mesh();
    Stage *s = find_stage(stages, ii, ii, m, m, vfvol[i]->ext, u_ext);
    s->vfvol.push_back(vfvol[i]);
  }

  // process surface vector forms
  for (unsigned i = 0; i < vfsurf.size(); i++) {
    unsigned int ii = vfsurf[i]->i;
    Mesh *m = spaces[ii]->get_mesh();
    Stage *s = find_stage(stages, ii, ii, m, m, vfsurf[i]->ext, u_ext);
    s->vfsurf.push_back(vfsurf[i]);
  }

  // helper macro for iterating in a set
#define set_for_each(myset, type) \
  for (std::set<type>::iterator it = (myset).begin(); it != (myset).end(); it++)

  // initialize the arrays meshes and fns needed by Traverse for each stage
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
<<<<<<< HEAD
        equal(seq.begin(), seq.end(), stages[i].seq_set.begin())) {
      s = &stages[i]; break;
=======
    equal(seq.begin(), seq.end(), stages[i].seq_set.begin())) { 
      s = &stages[i]; break; 
>>>>>>> Fixes & comments & old code deletion.
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


/// Returns a (neq x neq) array containing true in each element, if the corresponding
/// block of weak forms is used, and false otherwise.
/// This function is the same in H2D and H3D.
///
bool** WeakForm::get_blocks(bool force_diagonal_blocks)
{
  _F_
  bool** blocks = new_matrix<bool>(neq, neq);
  for (unsigned int i = 0; i < neq; i++) {
    for (unsigned int j = 0; j < neq; j++) {
      blocks[i][j] = false;
    }
    if (force_diagonal_blocks == true) blocks[i][i] = true;
  }
  for (unsigned i = 0; i < mfvol.size(); i++) {
    if (fabs(mfvol[i]->scaling_factor) > 1e-12) blocks[mfvol[i]->i][mfvol[i]->j] = true;
    if (mfvol[i]->sym) {
      if (fabs(mfvol[i]->scaling_factor) > 1e-12) blocks[mfvol[i]->j][mfvol[i]->i] = true;
    }
  }

  for (unsigned i = 0; i < mfsurf.size(); i++) {
    if (fabs(mfsurf[i]->scaling_factor) > 1e-12) blocks[mfsurf[i]->i][mfsurf[i]->j] = true;
  }
  return blocks;
}


//// areas /////////////////////////////////////////////////////////////////////////////////////////

bool WeakForm::is_in_area_2(int marker, int area) const
{
  _F_
  if (-area > (int)(areas.size())) error("Invalid area number.");
  const Area* a = &areas[-area-1];

  for (unsigned int i = 0; i < a->markers.size(); i++)
    if (a->markers[i] == marker)
  return true;

  return false;
}

// Function which according to the conversion table provided, updates the above members.
void WeakForm::update_markers_acc_to_conversion(Mesh::MarkersConversion* markers_conversion)
{
  Hermes::vector<MeshFunction*> vector_to_pass;
  // FIXME
  /*
    std::map<std::string, MatrixFormVol>::iterator it_mfv;
    for(it_mfv = mfvol_string_temp.begin(); it_mfv != mfvol_string_temp.end(); it_mfv++) {
    vector_to_pass = it_mfv->second.ext;
    add_matrix_form(it_mfv->second.i, it_mfv->second.j, it_mfv->second.fn,
        it_mfv->second.ord, (SymFlag)it_mfv->second.sym,
        markers_conversion->get_internal_boundary_marker(it_mfv->first),
        vector_to_pass);
    }
    std::map<std::string, MatrixFormSurf>::iterator it_mfs;
    for(it_mfs = mfsurf_string_temp.begin(); it_mfs != mfsurf_string_temp.end(); it_mfs++) {
    vector_to_pass = it_mfs->second.ext;
    add_matrix_form_surf(it_mfs->second.i, it_mfs->second.j,
             it_mfs->second.fn, it_mfs->second.ord,
             markers_conversion->get_internal_boundary_marker(it_mfs->first),
             vector_to_pass);
    }

    std::map<std::string, VectorFormVol>::iterator it_vfv;
    for(it_vfv = vfvol_string_temp.begin(); it_vfv != vfvol_string_temp.end(); it_vfv++) {
    vector_to_pass = it_vfv->second.ext;
    add_vector_form(it_vfv->second.i, it_vfv->second.fn, it_vfv->second.ord,
        markers_conversion->get_internal_boundary_marker(it_vfv->first),
        vector_to_pass);
    }

    std::map<std::string, VectorFormSurf>::iterator it_vfs;
    for(it_vfs = vfsurf_string_temp.begin(); it_vfs != vfsurf_string_temp.end(); it_vfs++) {
    vector_to_pass = it_vfs->second.ext;
    add_vector_form_surf(it_vfs->second.i, it_vfs->second.fn, it_vfs->second.ord,
             markers_conversion->get_internal_boundary_marker(it_vfs->first),
             vector_to_pass);
    }
    */
}


WeakForm::Form::Form(int area, Hermes::vector<MeshFunction *> ext, double scaling_factor, int u_ext_offset) : 
    area(area), ext(ext), scaling_factor(scaling_factor), u_ext_offset(u_ext_offset)
{
}

WeakForm::MatrixFormVol::MatrixFormVol(unsigned int i, unsigned int j, SymFlag sym, 
    int area, Hermes::vector<MeshFunction *> ext, double scaling_factor, int u_ext_offset) : Form(area, ext, scaling_factor, u_ext_offset), i(i), j(j), sym(sym)
{
}

WeakForm::MatrixFormSurf::MatrixFormSurf(unsigned int i, unsigned int j, int area, 
    Hermes::vector<MeshFunction *> ext, double scaling_factor, int u_ext_offset) : Form(area, ext, scaling_factor, u_ext_offset), i(i), j(j)
{
}

WeakForm::VectorFormVol::VectorFormVol(unsigned int i, int area, 
      Hermes::vector<MeshFunction *> ext, double scaling_factor, int u_ext_offset) : Form(area, ext, scaling_factor, u_ext_offset), i(i)
{
}

WeakForm::VectorFormSurf::VectorFormSurf(unsigned int i, int area, 
      Hermes::vector<MeshFunction *> ext, double scaling_factor, int u_ext_offset) : Form(area, ext, scaling_factor, u_ext_offset), i(i)
{
}
