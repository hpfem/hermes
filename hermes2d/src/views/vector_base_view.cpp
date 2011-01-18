// This file is part of Hermes2D.
//
// Copyright 2005-2008 Jakub Cerveny <jakub.cerveny@gmail.com>
// Copyright 2005-2008 Lenka Dubcova <dubcova@gmail.com>
// Copyright 2005-2008 Pavel Solin <solin@unr.edu>
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

// $Id: view4.cpp 1086 2008-10-21 09:05:44Z jakub $

#ifndef NOGLUT

#include <GL/freeglut.h>
#include "../h2d_common.h"
#include "vector_base_view.h"

void VectorBaseView::show(Space* space)
{
  free();
  pss = new PrecalcShapeset(space->get_shapeset());
  sln = new Solution();
  this->space = space;
  ndof = Space::get_num_dofs(space);
  base_index = 0;
  update_solution();
}


void VectorBaseView::free()
{
  if (pss != NULL) { delete pss; pss = NULL; }
  if (sln != NULL) { delete sln; sln = NULL; }
}


void VectorBaseView::update_solution()
{
  scalar* coeffs = new scalar[ndof + 1];
  memset(coeffs, 0, sizeof(scalar) * (ndof + 1));
  if (base_index >= -1 && base_index < ndof)
    coeffs[base_index + 1] = 1.0;
  Solution::vector_to_solution(coeffs, space, sln, pss);

  VectorView::show(sln,  sln, 0.001, H2D_FN_VAL_0, H2D_FN_VAL_1);
  update_title();

  delete [] coeffs;
}


void VectorBaseView::update_title()
{
  std::stringstream str;
  str << basic_title << " - dof = " << base_index;
  if (base_index < 0)
    str << " (Dirichlet lift)";
  View::set_title(str.str().c_str());
}


void VectorBaseView::on_special_key(int key, int x, int y)
{
  switch (key)
  {
    case GLUT_KEY_LEFT:
      if (base_index > -1) base_index--;
      update_solution();
      break;

    case GLUT_KEY_RIGHT:
      if (base_index < ndof-1) base_index++;
      update_solution();
      break;

    default:
      VectorView::on_special_key(key, x, y);
  }
}


const char* VectorBaseView::get_help_text() const
{
  return
  "VectorBaseView\n\n"
  "Controls:\n"
  "  Left mouse - pan\n"
  "  Right mouse - zoom\n"
  "  Left arrow - previous basis function\n"
  "  Right arrow - next basis function\n"
  "  C - center image\n"
  "  F - toggle smooth palette\n"
  "  X - toggle hexagonal grid\n"
  "  H - render high-quality frame\n"
  "  M - toggle mesh\n"
  "  P - cycle palettes\n"
  "  S - save screenshot\n"
  "  F1 - this help\n"
  "  Esc, Q - quit";
}

#endif // NOGLUT
