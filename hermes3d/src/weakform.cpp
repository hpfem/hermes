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
//
// This file was taken from hermes2d and adjusted for hermes3d
//

#include "common.h"
#include "weakform.h"
#include "matrix.h"
#include <common/trace.h>
#include <common/error.h>
#include <common/callstack.h>

WeakForm::WeakForm(int neq, bool mat_free)
{
	_F_
	this->neq = neq;
	this->is_matfree = mat_free;
}

// single equation case
WeakForm::WeakForm(bool mat_free)
{
	_F_
	this->neq = 1;
	this->is_matfree = mat_free;
}

WeakForm::~WeakForm()
{
	_F_
}

void WeakForm::add_matrix_form(int i, int j, matrix_form_val_t fn, matrix_form_ord_t ord, SymFlag sym, int area,
                               Tuple<MeshFunction*> ext)
{
	_F_
	if (i < 0 || i >= neq || j < 0 || j >= neq) error("Invalid equation number.");
	if (sym != ANTISYM && sym != UNSYM && sym != SYM) error("\"sym\" must be ANTISYM, UNSYM or SYM.");
	if (sym < 0 && i == j) error("Only off-diagonal forms can be antisymmetric.");
	if (area != ANY && area < 0 && -area > (signed) areas.size()) error("Invalid area number.");
	if (mfvol.size() > 100) warning("Large number of forms (> 100). Is this the intent?");

	MatrixFormVol form = { i, j, sym, area, fn, ord };
	int nx = ext.size();
	for (int i = 0; i < nx; i++) form.ext.push_back(ext[i]);
	mfvol.push_back(form);
}

void WeakForm::add_matrix_form_surf(int i, int j, matrix_form_val_t fn, matrix_form_ord_t ord, int area, 
                                    Tuple<MeshFunction*> ext)
{
	_F_
	if (i < 0 || i >= neq || j < 0 || j >= neq) error("Invalid equation number.");
	if (area != ANY && area < 0 && -area > (signed) areas.size()) error("Invalid area number.");

	MatrixFormSurf form = { i, j, area, fn, ord };
	int nx = ext.size();
	for (int i = 0; i < nx; i++) form.ext.push_back(ext[i]);
	mfsurf.push_back(form);
}

void WeakForm::add_vector_form(int i, vector_form_val_t fn, vector_form_ord_t ord, int area, 
                               Tuple<MeshFunction*> ext)
{
	_F_
	if (i < 0 || i >= neq) error("Invalid equation number.");
	if (area != ANY && area < 0 && -area > (signed) areas.size()) error("Invalid area number.");

	VectorFormVol form = { i, area, fn, ord };
	int nx = ext.size();
	for (int i = 0; i < nx; i++) form.ext.push_back(ext[i]);
	vfvol.push_back(form);
}

void WeakForm::add_vector_form_surf(int i, vector_form_val_t fn, vector_form_ord_t ord, int area, 
                                    Tuple<MeshFunction*> ext)
{
	_F_
	if (i < 0 || i >= neq) error("Invalid equation number.");
	if (area != ANY && area < 0 && -area > (signed) areas.size()) error("Invalid area number.");

	VectorFormSurf form = { i, area, fn, ord };
	int nx = ext.size();
	for (int i = 0; i < nx; i++) form.ext.push_back(ext[i]);
	vfsurf.push_back(form);
}

void WeakForm::set_ext_fns(void *fn, Tuple<MeshFunction*> ext)
{
	EXIT(H3D_ERR_NOT_IMPLEMENTED);
}


//// stages ////////////////////////////////////////////////////////////////////////////////////////

/// Constructs a list of assembling stages. Each stage contains a list of forms
/// that share the same meshes. Each stage is then assembled separately. This
/// improves the performance of multi-mesh assembling.
///
void WeakForm::get_stages(Space **spaces, std::vector<WeakForm::Stage> &stages, bool rhsonly)
{
	_F_
	unsigned i;
	stages.clear();

	if (!rhsonly) {
		if (is_linear()) {
			// process volume biforms
			for (i = 0; i < mfvol.size(); i++) {
				int ii = mfvol[i].i, jj = mfvol[i].j;
				Mesh *m1 = spaces[ii]->get_mesh();
				Mesh *m2 = spaces[jj]->get_mesh();
				Stage *s = find_stage(stages, ii, jj, m1, m2, mfvol[i].ext);
				s->mfvol.push_back(&mfvol[i]);
			}

			// process surface biforms
			for (i = 0; i < mfsurf.size(); i++) {
				int ii = mfsurf[i].i, jj = mfsurf[i].j;
				Mesh *m1 = spaces[ii]->get_mesh();
				Mesh *m2 = spaces[jj]->get_mesh();
				Stage *s = find_stage(stages, ii, jj, m1, m2, mfsurf[i].ext);
				s->mfsurf.push_back(&mfsurf[i]);
			}
		}
		else {
			// process volume jac forms
			for (unsigned i = 0; i < mfvol.size(); i++) {
				int ii = mfvol[i].i, jj = mfvol[i].j;
				Mesh *m1 = spaces[ii]->get_mesh();
				Mesh *m2 = spaces[jj]->get_mesh();
				Stage *s = find_stage(stages, ii, jj, m1, m2, mfvol[i].ext);
				s->mfvol.push_back(&mfvol[i]);
			}

			// process surface jac forms
			for (unsigned i = 0; i < mfsurf.size(); i++) {
				int ii = mfsurf[i].i, jj = mfsurf[i].j;
				Mesh *m1 = spaces[ii]->get_mesh();
				Mesh *m2 = spaces[jj]->get_mesh();
				Stage *s = find_stage(stages, ii, jj, m1, m2, mfsurf[i].ext);
				s->mfsurf.push_back(&mfsurf[i]);
			}

		}
	}

	if (is_linear()) {
		// process volume liforms
		for (i = 0; i < vfvol.size(); i++) {
			int ii = vfvol[i].i;
			Mesh *m = spaces[ii]->get_mesh();
			Stage *s = find_stage(stages, ii, ii, m, m, vfvol[i].ext);
			s->vfvol.push_back(&vfvol[i]);
		}

		// process surface liforms
		for (i = 0; i < vfsurf.size(); i++) {
			int ii = vfsurf[i].i;
			Mesh *m = spaces[ii]->get_mesh();
			Stage *s = find_stage(stages, ii, ii, m, m, vfsurf[i].ext);
			s->vfsurf.push_back(&vfsurf[i]);
		}
	}
	else {
		// process volume res forms
		for (unsigned i = 0; i < vfvol.size(); i++) {
			int ii = vfvol[i].i;
			Mesh *m = spaces[ii]->get_mesh();
			Stage *s = find_stage(stages, ii, ii, m, m, vfvol[i].ext);
			s->vfvol.push_back(&vfvol[i]);
		}

		// process surface res forms
		for (unsigned i = 0; i < vfsurf.size(); i++) {
			int ii = vfsurf[i].i;
			Mesh *m = spaces[ii]->get_mesh();
			Stage *s = find_stage(stages, ii, ii, m, m, vfsurf[i].ext);
			s->vfsurf.push_back(&vfsurf[i]);
		}
	}

	// helper macro for iterating in a set
	#define SET_FOR_EACH(myset, type) \
		for (std::set<type>::iterator it = (myset).begin(); it != (myset).end(); it++)

	// initialize the arrays meshes and fns needed by Traverse for each stage
	for (i = 0; i < stages.size(); i++) {
		Stage *s = &stages[i];
		SET_FOR_EACH(s->idx_set, int) {
			s->idx.push_back(*it);
			s->meshes.push_back(spaces[*it]->get_mesh());
			s->fns.push_back(NULL);
		}
		SET_FOR_EACH(s->ext_set, MeshFunction *) {
			s->ext.push_back(*it);
			s->meshes.push_back((*it)->get_mesh());
			s->fns.push_back(*it);
		}
		s->idx_set.clear();
		s->seq_set.clear();
		s->ext_set.clear();
	}
}

/// Finds an assembling stage with the same set of meshes as [m1, m2, ext]. If no such
/// stage can be found, a new one is created and returned.
///
WeakForm::Stage *WeakForm::find_stage(std::vector<WeakForm::Stage> &stages, int ii, int jj,
                                      Mesh *m1, Mesh *m2, std::vector<MeshFunction *> &ext)
{
	_F_
	// first create a list of meshes the form uses
	std::set<unsigned> seq;
	seq.insert(m1->get_seq());
	seq.insert(m2->get_seq());
	for (unsigned i = 0; i < ext.size(); i++)
		seq.insert(ext[i]->get_mesh()->get_seq());

	// find a suitable existing stage for the form
	Stage *s = NULL;
	for (unsigned i = 0; i < stages.size(); i++)
		if (seq.size() == stages[i].seq_set.size() &&
			equal(seq.begin(), seq.end(), stages[i].seq_set.begin()))
		{
			s = &stages[i];
			break;
		}

	// create a new stage if not found
	if (s == NULL) {
		Stage newstage;
		stages.push_back(newstage);
		s = &stages.back();
		s->seq_set = seq;
	}

	// update and return the stage
	for (unsigned i = 0; i < ext.size(); i++)
		s->ext_set.insert(ext[i]);
	s->idx_set.insert(ii);
	s->idx_set.insert(jj);

	return s;
}


/// Returns a (neq x neq) array containing true in each element, if the corresponding
/// block of weak forms is used, and false otherwise.
///
bool **WeakForm::get_blocks()
{
	_F_
	bool **blocks = new_matrix<bool>(neq, neq);
	for (int i = 0; i < neq; i++)
		for (int j = 0; j < neq; j++)
			blocks[i][j] = false;

	if (is_linear()) {
		for (unsigned i = 0; i < mfvol.size(); i++) {
			blocks[mfvol[i].i][mfvol[i].j] = true;
			if (mfvol[i].sym)
				blocks[mfvol[i].j][mfvol[i].i] = true;
		}

		for (unsigned i = 0; i < mfsurf.size(); i++)
			blocks[mfsurf[i].i][mfsurf[i].j] = true;
	}
	else {
		for (unsigned i = 0; i < mfvol.size(); i++) {
			blocks[mfvol[i].i][mfvol[i].j] = true;
			if (mfvol[i].sym)
				blocks[mfvol[i].j][mfvol[i].i] = true;
		}

		for (unsigned i = 0; i < mfsurf.size(); i++)
			blocks[mfsurf[i].i][mfsurf[i].j] = true;
	}

	return blocks;
}

//// areas /////////////////////////////////////////////////////////////////////////////////////////

int WeakForm::def_area(Tuple<int> area_markers)
{
	_F_
	Area newarea;
        int n = area_markers.size();
	for (int i = 0; i < n; i++) newarea.markers.push_back(area_markers[i]);

	areas.push_back(newarea);
	return -areas.size();
}


bool WeakForm::is_in_area_2(int marker, int area) const
{
	_F_
	if (-area > (signed) areas.size()) error("Invalid area number.");
	const Area *a = &areas[-area - 1];

	for (unsigned i = 0; i < a->markers.size(); i++)
		if (a->markers[i] == marker)
			return true;

	return false;
}
