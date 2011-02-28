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

/// Abstract class representing a general boundary condition.
class HERMES_API BoundaryCondition {
public: 
  /// Types of boundary conditions.
  /// There is no need to a special marker BC_NONE, as the default condition is zero Dirichlet.
  enum BoundaryConditionType {
    BC_DIRICHLET, ///< Dirichlet BC.
    BC_NEUMANN,   ///< Neumann BC.
    BC_NEWTON,    ///< Newton BC.
    BC_NONE       ///< Empty (none) BC.
  };

  /// Types of description of boundary values, either a function (callback), or a constant.
  enum BoundaryConditionValueType {
    BC_FUNCTION,
    BC_VALUE
  };

  /// Default constructor.
  BoundaryCondition(Hermes::vector<std::string> markers);

  /// Default destructor
  ~BoundaryCondition();

  /// Sets the current time for time-dependent boundary conditions.
  void set_current_time(double time);

  /// Gets the type of this boundary condition. Pure virtual function making this class abstract.
  virtual BoundaryConditionType get_type() const = 0;

protected:
  /// Current time.
  double current_time;

  /// Markers where this boundary condition is present.
  /// This facilitates the creation, that one condition can be imposed on multiple parts of the boundary.
  Hermes::vector<std::string> markers;

  /// Friend class.
  friend class BoundaryConditions;
};



/// Abstract class representing Dirichlet boundary condition of the form u|_{\Gamma_Dirichlet} = u_Dirichlet.
class HERMES_API DirichletBoundaryCondition : public BoundaryCondition {
public:
  /// Default constructor.
  DirichletBoundaryCondition(Hermes::vector<std::string> markers);

  /// Virtual destructor.
  virtual ~DirichletBoundaryCondition();

  /// Pure virtual function giving info whether u_Dirichlet is a constant or a function.
  virtual BoundaryConditionValueType get_value_type() const = 0;

  /// Gets the type of this boundary condition.
  inline BoundaryConditionType get_type() const { return BoundaryCondition::BC_DIRICHLET; }

  /// Represents a function prescribed on the boundary.
  virtual scalar function(double x, double y) const;

  /// Special case of a constant function.
  scalar value;
};

/// Class representing Dirichlet boundary condition of the form u|_{\Gamma_Dirichlet} = u_Dirichlet given by value.
class HERMES_API DirichletValueBoundaryCondition : public DirichletBoundaryCondition {
public:
  /// Constructor.
  DirichletValueBoundaryCondition(Hermes::vector<std::string> markers, scalar value);

  /// Function giving info that u_Dirichlet is a constant.
  inline BoundaryConditionValueType get_value_type() const { return BoundaryCondition::BC_VALUE; }
};



/// Abstract class representing Neumann boundary condition of the form \nabla u \cdot n |_{\Gamma_Neumann} = u_Neumann.
class HERMES_API NeumannBoundaryCondition : public BoundaryCondition {
public:
  /// Default constructor.
  NeumannBoundaryCondition(Hermes::vector<std::string> markers);

  /// Virtual destructor.
  virtual ~NeumannBoundaryCondition();

  /// Pure virtual function giving info whether u_Neumann is a constant or a function.
  virtual BoundaryConditionValueType get_value_type() const = 0;

  /// Gets the type of this boundary condition.
  inline BoundaryConditionType get_type() const { return BoundaryCondition::BC_NEUMANN; }

  /// Represents a function prescribed on the boundary.
  virtual scalar function(double x, double y) const;

  /// Special case of a constant function.
  scalar value;
};

/// Class representing Neumann boundary condition of the form \nabla u \cdot n |_{\Gamma_Neumann} = u_Neumann given by value.
class HERMES_API NeumannValueBoundaryCondition : public NeumannBoundaryCondition {
public:
  /// Constructor.
  NeumannValueBoundaryCondition(Hermes::vector<std::string> markers, scalar value);

  /// Function giving info that u_Neumann is a constant.
  inline BoundaryConditionValueType get_value_type() const { return BoundaryCondition::BC_VALUE; }
};


/// Abstract class representing Newton (mixed) boundary condition of the form \nabla u \cdot n |_{\Gamma_Neumann} + g * u = u_Newton.
class HERMES_API NewtonBoundaryCondition : public BoundaryCondition {
public:
  /// Default constructor.
  NewtonBoundaryCondition(Hermes::vector<std::string> markers);

  /// Virtual destructor.
  virtual ~NewtonBoundaryCondition();

  /// Pure virtual function giving info whether u_Newton is a constant or a function.
  virtual BoundaryConditionValueType get_value_type_u() const = 0;

  /// Pure virtual function giving info whether g is a constant or a function.
  virtual BoundaryConditionValueType get_value_type_g() const = 0;

  /// Gets the type of this boundary condition.
  inline BoundaryConditionType get_type() const { return BoundaryCondition::BC_NEWTON; }

protected:
  /// Represents a function prescribed on the boundary.
  virtual scalar function_u(double x, double y) const;

  /// Special case of a constant function.
  scalar value_u;

  /// Represents a function prescribed on the boundary.
  virtual scalar function_g(double x, double y) const;

  /// Special case of a constant function.
  scalar value_g;
};

/// Class representing empty boundary condition.
class HERMES_API EmptyBoundaryCondition : public BoundaryCondition {
public:
  /// Constructor.
  EmptyBoundaryCondition(Hermes::vector<std::string> markers) : BoundaryCondition(markers) {};

  /// Function giving info that u_Neumann is a constant.
  inline BoundaryConditionType get_type() const { return BoundaryCondition::BC_NONE; };
};


/// Class encapsulating all boundary conditions of one problem.
/// Using the class BoundaryCondition and its descendants.
class HERMES_API BoundaryConditions {
public:
  /// Default constructor.
  BoundaryConditions();

  /// Constructor with all boundary conditions of a problem.
  BoundaryConditions(Hermes::vector<BoundaryCondition *> boundary_conditions);

  /// Default destructor.
  ~BoundaryConditions();

  /// Initializes the class, fills the structures.
  void add_boundary_conditions(Hermes::vector<BoundaryCondition *> boundary_conditions);

  /// Public iterators for the private data structures.
  Hermes::vector<BoundaryCondition *>::const_iterator all_iterator;
  Hermes::vector<BoundaryCondition *>::const_iterator all_begin() const;
  Hermes::vector<BoundaryCondition *>::const_iterator all_end() const;

  Hermes::vector<DirichletBoundaryCondition *>::const_iterator dirichlet_iterator;
  Hermes::vector<DirichletBoundaryCondition *>::const_iterator dirichlet_begin() const;
  Hermes::vector<DirichletBoundaryCondition *>::const_iterator dirichlet_end() const;

  Hermes::vector<NeumannBoundaryCondition *>::const_iterator neumann_iterator;
  Hermes::vector<NeumannBoundaryCondition *>::const_iterator neumann_begin() const;
  Hermes::vector<NeumannBoundaryCondition *>::const_iterator neumann_end() const;

  Hermes::vector<NewtonBoundaryCondition *>::const_iterator newton_iterator;
  Hermes::vector<NewtonBoundaryCondition *>::const_iterator newton_begin() const;
  Hermes::vector<NewtonBoundaryCondition *>::const_iterator newton_end() const;

  std::map<std::string, BoundaryCondition *>::const_iterator markers_iterator;
  std::map<std::string, BoundaryCondition *>::const_iterator markers_begin() const;
  std::map<std::string, BoundaryCondition *>::const_iterator markers_end() const;

  BoundaryCondition* get_boundary_condition(std::string marker);

  /// Sets the current time for time-dependent boundary conditions.
  void set_current_time(double time);

private:
  /// All boundary conditions together.
  Hermes::vector<BoundaryCondition *> all;

  /// Dirichlet boundary conditions.
  Hermes::vector<DirichletBoundaryCondition *> dirichlet;

  /// Neumann boundary conditions.
  Hermes::vector<NeumannBoundaryCondition *> neumann;

  /// Newton boundary conditions.
  Hermes::vector<NewtonBoundaryCondition *> newton;

  /// Boundary markers cache
  std::map<std::string, BoundaryCondition *> markers;

  /// Create boundary markers cache for assembling
  void create_marker_cache();

  /// For purpose of filling non-user-provided BCs.
  EmptyBoundaryCondition* empty_condition;
};

#endif

