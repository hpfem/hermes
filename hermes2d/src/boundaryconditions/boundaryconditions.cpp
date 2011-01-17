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

#include "../hermes2d.h"

// General BC.
BoundaryCondition::BoundaryCondition() {
  current_time = 0.0;
};

BoundaryCondition::~BoundaryCondition() {};

void BoundaryCondition::set_current_time(double time) {
  current_time = time;
};

// Dirichlet BC.
DirichletBoundaryCondition::DirichletBoundaryCondition() {
  value = 0.0;
};

DirichletBoundaryCondition::~DirichletBoundaryCondition() {};

scalar DirichletBoundaryCondition::function(double x, double y) const {
  error("DirichletBoundaryCondition::Function used either for a constant condition, or not redefined for nonconstant condition.");
  return 0.0;
};

BoundaryCondition::BoundaryConditionType DirichletBoundaryCondition::get_type() {
  return BoundaryCondition::BC_DIRICHLET;
};

// Neumann BC.
NeumannBoundaryCondition::NeumannBoundaryCondition() {
  value = 0.0;
};

NeumannBoundaryCondition::~NeumannBoundaryCondition() {};

scalar NeumannBoundaryCondition::function(double x, double y) const {
  error("NeumannBoundaryCondition::Function used either for a constant condition, or not redefined for nonconstant condition.");
  return 0.0;
};

BoundaryCondition::BoundaryConditionType NeumannBoundaryCondition::get_type() {
  return BoundaryCondition::BC_NEUMANN;
};

// Newton BC.
NewtonBoundaryCondition::NewtonBoundaryCondition() {
  value_u = 0.0;
  value_g = 0.0;
};

NewtonBoundaryCondition::~NewtonBoundaryCondition() {};

scalar NewtonBoundaryCondition::function_u(double x, double y) const {
  error("NewtonBoundaryCondition::Function used either for a constant condition, or not redefined for nonconstant condition.");
  return 0.0;
};

scalar NewtonBoundaryCondition::function_g(double x, double y) const {
  error("NeumannBoundaryCondition::Function used either for a constant condition, or not redefined for nonconstant condition.");
  return 0.0;
};

BoundaryCondition::BoundaryConditionType NewtonBoundaryCondition::get_type() {
  return BoundaryCondition::BC_NEWTON;
};

// BoundaryConditions.

BoundaryConditions::BoundaryConditions() {};

BoundaryConditions::BoundaryConditions(Hermes::vector<BoundaryCondition *> boundary_conditions) {
  all = boundary_conditions;
  for(all_iterator = all_begin(); all_iterator != all_end(); all_iterator++)
    switch((*all_iterator)->get_type()) {
      case BoundaryCondition::BC_DIRICHLET:
        dirichlet.push_back(static_cast<DirichletBoundaryCondition *>(*all_iterator));
        break;
      case BoundaryCondition::BC_NEUMANN:
        neumann.push_back(static_cast<NeumannBoundaryCondition *>(*all_iterator));
        break;
      case BoundaryCondition::BC_NEWTON:
        newton.push_back(static_cast<NewtonBoundaryCondition *>(*all_iterator));
        break;
    }
};

Hermes::vector<BoundaryCondition *>::const_iterator BoundaryConditions::all_begin() const {
  return all.begin();
}

Hermes::vector<BoundaryCondition *>::const_iterator BoundaryConditions::all_end() const {
  return all.end();
}

Hermes::vector<DirichletBoundaryCondition *>::const_iterator BoundaryConditions::dirichlet_begin() const {
  return dirichlet.begin();
}
Hermes::vector<DirichletBoundaryCondition *>::const_iterator BoundaryConditions::dirichlet_end() const {
  return dirichlet.end();
}

Hermes::vector<NeumannBoundaryCondition *>::const_iterator BoundaryConditions::neumann_begin() const {
  return neumann.begin();
}
Hermes::vector<NeumannBoundaryCondition *>::const_iterator BoundaryConditions::neumann_end() const {
  return neumann.end();
}

Hermes::vector<NewtonBoundaryCondition *>::const_iterator BoundaryConditions::newton_begin() const {
  return newton.begin();
}
Hermes::vector<NewtonBoundaryCondition *>::const_iterator BoundaryConditions::newton_end() const {
  return newton.end();
}
BoundaryConditions::~BoundaryConditions() {};

