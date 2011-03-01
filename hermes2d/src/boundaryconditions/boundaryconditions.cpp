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
BoundaryCondition::BoundaryCondition(Hermes::vector<std::string> markers) : markers(markers) {
  current_time = 0.0;
};

BoundaryCondition::~BoundaryCondition() {};

void BoundaryCondition::set_current_time(double time) {
  current_time = time;
};

// Dirichlet BC.
DirichletBoundaryCondition::DirichletBoundaryCondition(Hermes::vector<std::string> markers) : BoundaryCondition(markers) {
  value = 0.0;
};

DirichletBoundaryCondition::~DirichletBoundaryCondition() {};

scalar DirichletBoundaryCondition::function(double x, double y) const {
  error("DirichletBoundaryCondition::Function used either for a constant condition, or not redefined for nonconstant condition.");
  return 0.0;
};

// Dirichlet BC Value
DirichletValueBoundaryCondition::DirichletValueBoundaryCondition(Hermes::vector<std::string> markers, scalar value) : DirichletBoundaryCondition(markers) {
  this->value = value;
}

DirichletValueBoundaryCondition::DirichletValueBoundaryCondition(std::string marker, scalar value) : DirichletBoundaryCondition(Hermes::vector<std::string>()) {
  this->value = value;
  markers.push_back(marker);
}

// Neumann BC.
NeumannBoundaryCondition::NeumannBoundaryCondition(Hermes::vector<std::string> markers) : BoundaryCondition(markers) {
  value = 0.0;
};

NeumannBoundaryCondition::~NeumannBoundaryCondition() {};

scalar NeumannBoundaryCondition::function(double x, double y) const {
  error("NeumannBoundaryCondition::Function used either for a constant condition, or not redefined for nonconstant condition.");
  return 0.0;
};

// Neumann BC Value
NeumannValueBoundaryCondition::NeumannValueBoundaryCondition(Hermes::vector<std::string> markers, scalar value) : NeumannBoundaryCondition(markers) {
  this->value = value;
}

NeumannValueBoundaryCondition::NeumannValueBoundaryCondition(std::string marker, scalar value) : NeumannBoundaryCondition(Hermes::vector<std::string>()) {
  this->value = value;
  markers.push_back(marker);
}

// Newton BC.
NewtonBoundaryCondition::NewtonBoundaryCondition(Hermes::vector<std::string> markers) : BoundaryCondition(markers) {
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

// BoundaryConditions.

BoundaryConditions::BoundaryConditions() {
  std::ostringstream oss;
  oss << "Everywhere where no other is.";
  Hermes::vector<std::string> markers_to_pass;
  markers_to_pass.push_back(oss.str());
  this->empty_condition = new EmptyBoundaryCondition(markers_to_pass);
};

BoundaryConditions::BoundaryConditions(Hermes::vector<BoundaryCondition *> boundary_conditions) : all(boundary_conditions) {
  add_boundary_conditions(boundary_conditions);
};

BoundaryConditions::BoundaryConditions(BoundaryCondition * boundary_condition) {
  Hermes::vector<BoundaryCondition *> boundary_conditions;
  boundary_conditions.push_back(boundary_condition);
  add_boundary_conditions(boundary_conditions);
};

void BoundaryConditions::add_boundary_conditions(Hermes::vector<BoundaryCondition *> boundary_conditions) {
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

  create_marker_cache();

  std::ostringstream oss;
  oss << "Everywhere where no other is.";
  Hermes::vector<std::string> markers_to_pass;
  markers_to_pass.push_back(oss.str());
  this->empty_condition = new EmptyBoundaryCondition(markers_to_pass);
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
std::map<std::string, BoundaryCondition *>::const_iterator BoundaryConditions::markers_begin() const {
  return markers.begin();
}
std::map<std::string, BoundaryCondition *>::const_iterator BoundaryConditions::markers_end() const {
  return markers.end();
}
BoundaryConditions::~BoundaryConditions() {
  delete this->empty_condition;
};

void BoundaryConditions::create_marker_cache() {
  for(all_iterator = all_begin(); all_iterator != all_end(); all_iterator++)
    for(Hermes::vector<std::string>::const_iterator it = (*all_iterator)->markers.begin(); it != (*all_iterator)->markers.end(); it++)
    {
      // FIXME - check consistency
      if (markers[*it])
        error("Attempt to define more than one description of the BC on the same part of the boundary with marker '%d'.", *it);

      std::cout << "marker = " << *it
                << ", type = " << (*all_iterator)->get_type()
                << std::endl;

      markers[*it] = *all_iterator;
    }
}


BoundaryCondition* BoundaryConditions::get_boundary_condition(std::string marker) {
  if(markers.find(marker) == markers.end())
    return this->empty_condition;
  else
    return markers[marker];
}

void BoundaryConditions::set_current_time(double time) {
  for(all_iterator = all_begin(); all_iterator != all_end(); all_iterator++)
    (*all_iterator)->set_current_time(time);
};
