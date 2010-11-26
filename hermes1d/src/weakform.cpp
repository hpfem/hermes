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

#include "hermes1d.h"
#include "weakform.h"
#include "../../hermes_common/matrix.h"


WeakForm::WeakForm(int neq, bool mat_free)
{
  _F_
  this->neq = neq;
  seq = 0;
  this->is_matfree = mat_free;
}

//// interface /////////////////////////////////////////////////////////////////////////////////////

void WeakForm::add_matrix_form(int i, int j, matrix_form fn, Space* space, int marker)
{
    if (marker != ANY && marker < 0) error("Invalid element marker.");
    MatrixFormVol form = {i, j, fn, marker, space};
    this->matrix_forms_vol.push_back(form);
}

void WeakForm::add_vector_form(int i, vector_form fn, Space* space, int marker)
{
    if (marker != ANY && marker < 0) error("Invalid element marker.");
	  VectorFormVol form = {i, fn, marker, space};
    this->vector_forms_vol.push_back(form);
}

void WeakForm::add_matrix_form_surf(int i, int j, matrix_form_surf fn, int bdy_index)
{
    MatrixFormSurf form = {i, j, bdy_index, fn};
    this->matrix_forms_surf.push_back(form);
}

void WeakForm::add_vector_form_surf(int i, vector_form_surf fn, int bdy_index)
{
    VectorFormSurf form = {i, bdy_index, fn};
    this->vector_forms_surf.push_back(form);
}

void WeakForm::add_matrix_form(matrix_form fn, Space* space, int marker)
{
    if (marker != ANY && marker < 0) error("Invalid element marker.");
    MatrixFormVol form = {0, 0, fn, marker, space};
    this->matrix_forms_vol.push_back(form);
}

void WeakForm::add_vector_form(vector_form fn, Space* space, int marker)
{
    if (marker != ANY && marker < 0) error("Invalid element marker.");
	  VectorFormVol form = {0, fn, marker, space};
    this->vector_forms_vol.push_back(form);
}

void WeakForm::add_matrix_form_surf(matrix_form_surf fn, int bdy_index)
{
    MatrixFormSurf form = {0, 0, bdy_index, fn};
    this->matrix_forms_surf.push_back(form);
}

void WeakForm::add_vector_form_surf(vector_form_surf fn, int bdy_index)
{
    VectorFormSurf form = {0, bdy_index, fn};
    this->vector_forms_surf.push_back(form);
}