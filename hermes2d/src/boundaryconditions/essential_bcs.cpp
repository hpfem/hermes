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

#define HERMES_REPORT_WARN

#include "../hermes2d.h"
// Essential boundary condition.
template<typename Scalar>
EssentialBoundaryCondition<Scalar>::EssentialBoundaryCondition(Hermes::vector<std::string> markers) : markers(markers) {
  current_time = 0.0;
  value_const = 0.0;
};

template<typename Scalar>
EssentialBoundaryCondition<Scalar>::EssentialBoundaryCondition(std::string marker) {
  markers.push_back(marker);
  current_time = 0.0;
  value_const = 0.0;
};

template<typename Scalar>
EssentialBoundaryCondition<Scalar>::~EssentialBoundaryCondition() {};

template<typename Scalar>
Scalar EssentialBoundaryCondition<Scalar>::value(double x, double y) const {
  warn("EssentialBoundaryCondition<Scalar>::Function used either for a constant condition, or not redefined for nonconstant condition.");
  return 0.0;
};

template<typename Scalar>
void EssentialBoundaryCondition<Scalar>::set_current_time(double time) {
  this->current_time = time;
}

template<typename Scalar>
double EssentialBoundaryCondition<Scalar>::get_current_time() const {
  return current_time;
}

// Essential BoundaryCondition Constant.
template<typename Scalar>
DefaultEssentialBCConst<Scalar>::DefaultEssentialBCConst(Hermes::vector<std::string> markers, Scalar value_const) : EssentialBoundaryCondition(markers) {
  this->value_const = value_const;
}

template<typename Scalar>
DefaultEssentialBCConst<Scalar>::DefaultEssentialBCConst(std::string marker, Scalar value_const) : EssentialBoundaryCondition(Hermes::vector<std::string>()) {
  this->value_const = value_const;
  markers.push_back(marker);
}

// Essential BoundaryCondition NonConstant.
template<typename Scalar>
DefaultEssentialBCNonConst<Scalar>::DefaultEssentialBCNonConst(Hermes::vector<std::string> markers_, 
                                                       ExactSolutionScalar<Scalar>* exact_solution)  
       : EssentialBoundaryCondition(Hermes::vector<std::string>()), exact_solution(exact_solution) 
{
  for (unsigned int i=0; i < markers.size(); i++) markers.push_back(markers_[i]);
};

template<typename Scalar>
DefaultEssentialBCNonConst<Scalar>::DefaultEssentialBCNonConst(std::string marker, ExactSolutionScalar<Scalar>* exact_solution) 
       : EssentialBoundaryCondition(Hermes::vector<std::string>()), exact_solution(exact_solution) 
{
  markers.push_back(marker);
};

template<typename Scalar>
Scalar DefaultEssentialBCNonConst<Scalar>::value(double x, double y) const 
{
  return exact_solution->value(x, y);
};

// EssentialBCs.
template<typename Scalar>
EssentialBCs<Scalar>::EssentialBCs() {
};

template<typename Scalar>
EssentialBCs<Scalar>::EssentialBCs(Hermes::vector<EssentialBoundaryCondition<Scalar> *> essential_bcs) {
  add_boundary_conditions(essential_bcs);
};

template<typename Scalar>
EssentialBCs<Scalar>::EssentialBCs(EssentialBoundaryCondition<Scalar> * boundary_condition) {
  Hermes::vector<EssentialBoundaryCondition<Scalar> *> boundary_conditions;
  boundary_conditions.push_back(boundary_condition);
  add_boundary_conditions(boundary_conditions);
};

template<typename Scalar>
void EssentialBCs<Scalar>::add_boundary_conditions(Hermes::vector<EssentialBoundaryCondition<Scalar> *> boundary_conditions) {
  for(Hermes::vector<EssentialBoundaryCondition<Scalar> *>::iterator it = boundary_conditions.begin(); it != boundary_conditions.end(); it++)
        all.push_back(*it);

  markers.clear();
  create_marker_cache();
};

template<typename Scalar>
void EssentialBCs<Scalar>::add_boundary_condition(EssentialBoundaryCondition<Scalar> * boundary_condition) {
  Hermes::vector<EssentialBoundaryCondition<Scalar> *> boundary_conditions;
  boundary_conditions.push_back(boundary_condition);
  add_boundary_conditions(boundary_conditions);
};

template<typename Scalar>
typename Hermes::vector<EssentialBoundaryCondition<Scalar> *>::const_iterator EssentialBCs<Scalar>::begin() const {
  return all.begin();
}

template<typename Scalar>
typename Hermes::vector<EssentialBoundaryCondition<Scalar> *>::const_iterator EssentialBCs<Scalar>::end() const {
  return all.end();
}

template<typename Scalar>
EssentialBCs<Scalar>::~EssentialBCs() {
};

template<typename Scalar>
void EssentialBCs<Scalar>::create_marker_cache() {
  for(this->iterator = begin(); iterator != end(); iterator++)
    for(Hermes::vector<std::string>::const_iterator it = (*iterator)->markers.begin(); it != (*iterator)->markers.end(); it++) {
      if (markers[*it] != NULL)
        error("Attempt to define more than one description of the BC on the same part of the boundary with marker '%s'.", it->c_str());
      markers[*it] = *iterator;
    }
}


template<typename Scalar>
EssentialBoundaryCondition<Scalar>* EssentialBCs<Scalar>::get_boundary_condition(std::string marker) {
  if(markers.find(marker) == markers.end())
    return NULL;
  else
    return markers[marker];
}

template<typename Scalar>
void EssentialBCs<Scalar>::set_current_time(double time) {
  for(iterator = begin(); iterator != end(); iterator++)
    (*iterator)->set_current_time(time);
};

template HERMES_API class EssentialBoundaryCondition<double>;
template HERMES_API class EssentialBoundaryCondition<std::complex<double>>;
template HERMES_API class DefaultEssentialBCConst<double>;
template HERMES_API class DefaultEssentialBCConst<std::complex<double>>;
template HERMES_API class DefaultEssentialBCNonConst<double>;
template HERMES_API class DefaultEssentialBCNonConst<std::complex<double>>;
template HERMES_API class EssentialBCs<double>;
template HERMES_API class EssentialBCs<std::complex<double>>;