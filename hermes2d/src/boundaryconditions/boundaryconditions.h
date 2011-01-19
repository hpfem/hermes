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

# include "vector.h"

/// Abstract class representing a general boundary condition.
class HERMES_API BoundaryCondition {
public:
  /// Default constructor.
  BoundaryCondition();

  /// Default destructor
  ~BoundaryCondition();

  /// Sets the current time for time-dependent boundary conditions.
  void set_current_time(double time);

  /// Gets the type of this boundary condition. Pure virtual function making this class abstract.
  virtual BCType get_type() const = 0;

protected:
  /// Current time.
  double current_time;

  /// Types of boundary conditions.
  /// There is no need to a special marker BC_NONE, as the default condition is zero Dirichlet.
  enum BoundaryConditionType {
    BC_DIRICHLET, ///< Dirichlet BC.
    BC_NEUMANN,   ///< Neumann BC.
    BC_NEWTON     ///< Newton BC.
  };

  /// Types of description of boundary values, either a function (callback), or a constant.
  enum BoundaryConditionValueType {
    BC_FUNCTION,
    BC_VALUE
  };

  /// Markers where this boundary condition is present.
  /// This facilitates the creation, that one condition can be imposed on multiple parts of the boundary.
  Hermes::vector<int> markers;

  /// Friend class.
  friend class BoundaryConditions;
};



/// Abstract class representing Dirichlet boundary condition of the form u|_{\Gamma_Dirichlet} = u_Dirichlet.
class HERMES_API DirichletBoundaryCondition : public BoundaryCondition {
public:
  /// Default constructor.
  DirichletBoundaryCondition();

  /// Virtual destructor.
  virtual ~DirichletBoundaryCondition();

  /// Pure virtual function giving info whether u_Dirichlet is a constant or a function.
  virtual BoundaryConditionValueType get_value_type() const = 0;

  /// Gets the type of this boundary condition.
  BoundaryConditionType get_type();

protected:
  /// Represents a function prescribed on the boundary.
  virtual scalar function(double x, double y) const;

  /// Special case of a constant function.
  scalar value;
};



/// Abstract class representing Neumann boundary condition of the form \nabla u \cdot n |_{\Gamma_Neumann} = u_Neumann.
class HERMES_API NeumannBoundaryCondition : public BoundaryCondition {
public:
  /// Default constructor.
  NeumannBoundaryCondition();

  /// Virtual destructor.
  virtual ~NeumannBoundaryCondition();

  /// Pure virtual function giving info whether u_Neumann is a constant or a function.
  virtual BoundaryConditionValueType get_value_type() const = 0;

  /// Gets the type of this boundary condition.
  BoundaryConditionType get_type();

protected:
  /// Represents a function prescribed on the boundary.
  virtual scalar function(double x, double y) const;

  /// Special case of a constant function.
  scalar value;
};



/// Abstract class representing Newton (mixed) boundary condition of the form \nabla u \cdot n |_{\Gamma_Neumann} + g * u = u_Newton.
class HERMES_API NewtonBoundaryCondition : public BoundaryCondition {
public:
  /// Default constructor.
  NewtonBoundaryCondition();

  /// Virtual destructor.
  virtual ~NewtonBoundaryCondition();

  /// Pure virtual function giving info whether u_Newton is a constant or a function.
  virtual BoundaryConditionValueType get_value_type_u() const = 0;

  /// Pure virtual function giving info whether g is a constant or a function.
  virtual BoundaryConditionValueType get_value_type_g() const = 0;

  /// Gets the type of this boundary condition.
  BoundaryConditionType get_type();

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

private:
  /// All boundary conditions together.
  Hermes::vector<BoundaryCondition *> all;

  /// Dirichlet boundary conditions.
  Hermes::vector<DirichletBoundaryCondition *> dirichlet;

  /// Neumann boundary conditions.
  Hermes::vector<NeumannBoundaryCondition *> neumann;

  /// Newton boundary conditions.
  Hermes::vector<NewtonBoundaryCondition *> newton;
};

#endif

