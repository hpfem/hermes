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
// along with Hermes2D. If not, see <http://www.gnu.org/licenses/>.

#ifndef __H2D_DEFINITIONS_H
#define __H2D_DEFINITIONS_H

// Bilinear form symmetry flag, see WeakForm::add_matrix_form
enum SymFlag
{
  HERMES_ANTISYM = -1,
  HERMES_NONSYM = 0,
  HERMES_SYM = 1
};

// Geometrical type of weak forms.
enum GeomType
{
  HERMES_PLANAR = 0,         // Planar problem.
  HERMES_AXISYM_X = 1,       // Axisymmetric problem where x-axis is the axis of symmetry.
  HERMES_AXISYM_Y = 2        // Axisymmetric problem where y-axis is the axis of symmetry.
};

// Basic type of boundary conditions.
enum BCType {
  HERMES_NATURAL = 0,        // Natural boundary condition (does not eliminate DOF).
  HERMES_ESSENTIAL = 1       // Essential boundary condition (does eliminate DOF).
};

// Boundary conditions for H1 problems.
enum BCTypeH1 {
  HERMES_DIRICHLET = 0,      // Dirichlet boundary condition.
  HERMES_NEUMANN = 1,        // Neumann boundary condition.
  HERMES_NEWTON = 2          // Newton boundary condition.
};

// Analysis type.
enum AnalysisType {
  HERMES_STEADY_STATE = 0,
  HERMES_TRANSIENT = 1,
  HERMES_HARMONIC = 2
};

#endif
