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

#ifndef __HERMES_COMMON_BOUNDARYCONDITIONS_H
#define __HERMES_COMMON_BOUNDARYCONDITIONS_H

#include "vector.h"
#include "map"

//#include "../function/solution.h"
class ExactSolutionScalar;
class ExactSolutionVector;

/// Abstract class representing Essential boundary condition of the form u|_{\Gamma_Essential} = u_Essential.
class HERMES_API EssentialBoundaryCondition
{
public:
  /// Default constructor.
  EssentialBoundaryCondition(Hermes::vector<std::string> markers);
  EssentialBoundaryCondition(std::string marker);

  /// Virtual destructor.
  virtual ~EssentialBoundaryCondition();

  /// Types of description of boundary values, either a function (callback), or a constant.
  enum EssentialBCValueType {
    BC_FUNCTION,
    BC_CONST
  };

  /// Pure virtual function reporting the type of the essential boundary condition.
  virtual EssentialBCValueType get_value_type() const = 0;

  /// Represents a function prescribed on the boundary. Gets the boundary point coordinate as well as the 
  /// normal and tangential vectors.
  virtual scalar value(double x, double y, double n_x, double n_y, double t_x, double t_y) const = 0;

  /// Special case of a constant function.
  scalar value_const;

  /// Sets the current time for time-dependent boundary conditions.
  void set_current_time(double time);
  double get_current_time() const;

protected:
  /// Current time.
  double current_time;
  
  // Markers.
  Hermes::vector<std::string> markers;

  // Friend class.
  friend class EssentialBCs;
  friend class Space;
};

/// Class representing constant essential boundary condition.
class HERMES_API DefaultEssentialBCConst : public EssentialBoundaryCondition {
public:
  /// Constructors.
  DefaultEssentialBCConst(Hermes::vector<std::string> markers, scalar value_const);
  DefaultEssentialBCConst(std::string marker, scalar value_const);

  /// Function reporting the type of the essential boundary condition.
  inline EssentialBCValueType get_value_type() const { return EssentialBoundaryCondition::BC_CONST; }
  virtual scalar value(double x, double y, double n_x, double n_y, double t_x, double t_y) const;
};

/// Class representing non-constant essential boundary condition for scalar approximation.
class HERMES_API DefaultEssentialBCNonConst : public EssentialBoundaryCondition
{
public:
  // Function values given by a scalar exact solution.
  DefaultEssentialBCNonConst(Hermes::vector<std::string> markers_, 
                             ExactSolutionScalar* exact_solution); 
  DefaultEssentialBCNonConst(std::string marker, ExactSolutionScalar* exact_solution); 
 
  ~DefaultEssentialBCNonConst() {};

  virtual scalar value(double x, double y, double n_x, double n_y, double t_x, double t_y) const;

  /// Function reporting the type of the essential boundary condition.
  inline EssentialBCValueType get_value_type() const { return EssentialBoundaryCondition::BC_FUNCTION; }

  ExactSolutionScalar* exact_solution;
};

/// Class representing non-constant essential boundary condition 
/// (tangential component for Hcurl approximations).
class HERMES_API DefaultEssentialBCNonConstHcurl : public EssentialBoundaryCondition
{
public:
  // Tangential values given by a vector-valued solution.
  DefaultEssentialBCNonConstHcurl(Hermes::vector<std::string> markers_, 
                                  ExactSolutionVector* exact_solution2); 
  DefaultEssentialBCNonConstHcurl(std::string marker, ExactSolutionVector* exact_solution2); 
 
  ~DefaultEssentialBCNonConstHcurl() {};

  virtual scalar value(double x, double y, double n_x, double n_y, double t_x, double t_y) const;

  /// Function giving info that u_Essential is a non-constant function.
  inline EssentialBCValueType get_value_type() const { return EssentialBoundaryCondition::BC_FUNCTION; }

  ExactSolutionVector* exact_solution2;
};

/// Class encapsulating all boundary conditions of one problem.
/// Using the class EssentialBCs and its descendants.
class HERMES_API EssentialBCs {
public:
  /// Default constructor.
  EssentialBCs();

  /// Constructor with all boundary conditions of a problem.
  EssentialBCs(Hermes::vector<EssentialBoundaryCondition *> essential_bcs);
  EssentialBCs(EssentialBoundaryCondition* boundary_condition);

  /// Default destructor.
  ~EssentialBCs();

  /// Initializes the class, fills the structures.
  void add_boundary_conditions(Hermes::vector<EssentialBoundaryCondition *> essential_bcs);
  void add_boundary_condition(EssentialBoundaryCondition* essential_bc);

  /// Public iterators for the private data structures.
  Hermes::vector<EssentialBoundaryCondition *>::const_iterator iterator;
  Hermes::vector<EssentialBoundaryCondition *>::const_iterator begin() const;
  Hermes::vector<EssentialBoundaryCondition *>::const_iterator end() const;
  
  EssentialBoundaryCondition* get_boundary_condition(std::string marker);

  /// Sets the current time for time-dependent boundary conditions.
  void set_current_time(double time);

private:
  /// All boundary conditions together.
  Hermes::vector<EssentialBoundaryCondition *> all;

  /// Boundary markers cache
  std::map<std::string, EssentialBoundaryCondition *> markers;

  /// Create boundary markers cache for assembling
  void create_marker_cache();
};

#endif

