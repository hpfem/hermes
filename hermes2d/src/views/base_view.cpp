// This file is part of Hermes2D.
//
// Copyright 2005-2008 Jakub Cerveny <jakub.cerveny@gmail.com>
// Copyright 2005-2008 Lenka Dubcova <dubcova@gmail.com>
// Copyright 2005-2008 Pavel Solin <solin@unr.edu>
// Copyright 2009-2010 Ivo Hanak <hanak@byte.cz>
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
#include "../common.h"
#include "base_view.h"

BaseView::BaseView(const char* title, int x, int y, int width, int height)
        : ScalarView(title, x, y, width, height)
{
  pss = NULL;
  sln = NULL;
  show_edges = true;
}

#ifndef _MSC_VER
BaseView::BaseView(const char* title, WinGeom* wg)
        : ScalarView(title, wg)
{
  pss = NULL;
  sln = NULL;
  show_edges = true;
}
#endif

BaseView::BaseView(char* title, WinGeom* wg)
        : ScalarView(title, wg)
{
  pss = NULL;
  sln = NULL;
  show_edges = true;
}

void BaseView::show(Space* space, double eps, int item)
{
  free();
  // todo: copy space
  pss = new PrecalcShapeset(space->get_shapeset());
  sln = new Solution();
  ndof = space->get_num_dofs();
  base_index = 0;
  this->space = space;
  this->eps = eps;
  this->item = item;
  update_solution();
}


void BaseView::free()
{
  if (pss != NULL) { delete pss; pss = NULL; }
  if (sln != NULL) { delete sln; sln = NULL; }
}


void BaseView::update_solution()
{
  Vector* vec = new AVector(ndof);
  memset(vec->get_c_array(), 0, sizeof(scalar) * ndof);
  if (base_index >= 0)
  {
    if (base_index < ndof) vec->set(base_index, 1.0);
    sln->set_fe_solution(space, pss, vec, 0.0);
  }
  else
  {
    sln->set_fe_solution(space, pss, vec, 1.0);
  }

  ScalarView::show(sln, eps, item);
  update_title();
  delete vec;
}


void BaseView::update_title()
{
  std::stringstream str;
  str << title << " - dof = " << base_index;
  if (base_index < 0)
    str << " (Dirichlet lift)";
  set_title(str.str().c_str());
}


void BaseView::on_special_key(int key, int x, int y)
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
      ScalarView::on_special_key(key, x, y);
  }
}


const char* BaseView::get_help_text() const
{
  return
  "BaseView\n\n"
  "Controls:\n"
  "  Left mouse - pan\n"
  "  Right mouse - zoom\n"
  "  Left arrow - previous basis function\n"
  "  Right arrow - next basis function\n"
  "  3 - toggle 3D mode\n"
  "  C - center image\n"
  "  F - toggle smooth palette\n"
  "  H - render high-quality frame\n"
  "  M - toggle mesh\n"
  "  P - cycle palettes\n"
  "  S - save screenshot\n"
  "  F1 - this help\n"
  "  Esc, Q - quit\n\n"
  "3D mode:\n"
  "  Left mouse - rotate\n"
  "  Right mouse - zoom\n"
  "  Middle mouse - pan\n"
  "  * - increase Z scale\n"
  "  / - decrease Z scale";
}

#endif // NOGLUT
