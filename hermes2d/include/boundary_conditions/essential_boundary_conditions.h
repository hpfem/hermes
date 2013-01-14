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

#include "../global.h"
namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar> class ExactSolutionScalar;
    template<typename Scalar> class ExactSolutionVector;
    template<typename Scalar> class EssentialBCs;

    /// Abstract class representing Essential boundary condition of the form u|_{\Gamma_Essential} = u_Essential.
    /// @ingroup inner
    template<typename Scalar>
    class HERMES_API EssentialBoundaryCondition : public Hermes::Mixins::Loggable
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
      virtual Scalar value(double x, double y, double n_x, double n_y, double t_x, double t_y) const = 0;

      /// Sets the current time for time-dependent boundary conditions.
      void set_current_time(double time);
      double get_current_time() const;

    protected:
      /// Special case of a constant function.
      Scalar value_const;

      /// Current time.
      double current_time;

      /// Markers.
      Hermes::vector<std::string> markers;

      template<typename T> friend class EssentialBCs;
      template<typename T> friend class Space;
      template<typename T> friend class H1Space;
      template<typename T> friend class L2Space;
      template<typename T> friend class HcurlSpace;
      template<typename T> friend class HdivSpace;
    };

    /// Class representing constant essential boundary condition.
    template<typename Scalar>
    class HERMES_API DefaultEssentialBCConst : public EssentialBoundaryCondition<Scalar> {
    public:
      /// Constructors.
      DefaultEssentialBCConst(Hermes::vector<std::string> markers, Scalar value_const);
      DefaultEssentialBCConst(std::string marker, Scalar value_const);

      virtual Scalar value(double x, double y, double n_x, double n_y, double t_x, double t_y) const;

      /// Function giving info that u_Essential is a constant.
      inline typename EssentialBoundaryCondition<Scalar>::EssentialBCValueType get_value_type() const { return EssentialBoundaryCondition<Scalar>::BC_CONST; }
    };

    /// Class representing non-constant essential boundary condition for Scalar approximation.
    template<typename Scalar>
    class HERMES_API DefaultEssentialBCNonConst : public EssentialBoundaryCondition<Scalar>
    {
    public:
      DefaultEssentialBCNonConst(Hermes::vector<std::string> markers_,
        ExactSolutionScalar<Scalar>* exact_solution);

      DefaultEssentialBCNonConst(std::string marker, ExactSolutionScalar<Scalar>* exact_solution);

      ~DefaultEssentialBCNonConst() {};

      virtual Scalar value(double x, double y, double n_x, double n_y, double t_x, double t_y) const;

      /// Function giving info that u_Essential is a non-constant function.
      inline typename EssentialBoundaryCondition<Scalar>::EssentialBCValueType get_value_type() const { return EssentialBoundaryCondition<Scalar>::BC_FUNCTION; }

      ExactSolutionScalar<Scalar>* exact_solution;
    };

    /// Class representing non-constant essential boundary condition
    /// (tangential component for Hcurl approximations).
    template<typename Scalar>
    class HERMES_API DefaultEssentialBCNonConstHcurl : public EssentialBoundaryCondition<Scalar>
    {
    public:
      // Tangential values given by a vector-valued solution.
      DefaultEssentialBCNonConstHcurl(Hermes::vector<std::string> markers_,
        ExactSolutionVector<Scalar>* exact_solution2);
      DefaultEssentialBCNonConstHcurl(std::string marker, ExactSolutionVector<Scalar>* exact_solution2);

      ~DefaultEssentialBCNonConstHcurl() {};

      virtual Scalar value(double x, double y, double n_x, double n_y, double t_x, double t_y) const;

      /// Function giving info that u_Essential is a non-constant function.
      inline typename EssentialBoundaryCondition<Scalar>::EssentialBCValueType get_value_type() const { return EssentialBoundaryCondition<Scalar>::BC_FUNCTION; }

      ExactSolutionVector<Scalar>* exact_solution2;
    };

    /// Class encapsulating all boundary conditions of one problem.
    /// Using the class EssentialBCs and its descendants.
    template<typename Scalar>
    class HERMES_API EssentialBCs {
    public:
      /// Default constructor.
      EssentialBCs();

      /// Constructor with all boundary conditions of a problem.
      EssentialBCs(Hermes::vector<EssentialBoundaryCondition<Scalar> *> essential_bcs);
      EssentialBCs(EssentialBoundaryCondition<Scalar>* boundary_condition);

      /// Default destructor.
      ~EssentialBCs();

      /// Initializes the class, fills the structures.
      void add_boundary_conditions(Hermes::vector<EssentialBoundaryCondition<Scalar> *> essential_bcs);
      void add_boundary_condition(EssentialBoundaryCondition<Scalar>* essential_bc);

      /// Public iterators for the private data structures.
      typename Hermes::vector<EssentialBoundaryCondition<Scalar> *>::const_iterator iterator;
      typename Hermes::vector<EssentialBoundaryCondition<Scalar> *>::const_iterator begin() const;
      typename Hermes::vector<EssentialBoundaryCondition<Scalar> *>::const_iterator end() const;

      EssentialBoundaryCondition<Scalar>* get_boundary_condition(std::string marker);

      /// Sets the current time for time-dependent boundary conditions.
      void set_current_time(double time);

    private:
      /// All boundary conditions together.
      Hermes::vector<EssentialBoundaryCondition<Scalar> *> all;

      /// Boundary markers.
      Hermes::vector<std::string> markers;
      /// Boundary conditions with the same order.
      Hermes::vector<EssentialBoundaryCondition<Scalar> *> BCs;

      /// Special boundary condition when it is defined on all boundary markers.
      EssentialBoundaryCondition<Scalar> * HermesAnyBC;

      /// Create boundary markers cache for assembling
      void create_marker_cache();

      template<typename T> friend class EssentialBCs;
      template<typename T> friend class Space;
      template<typename T> friend class H1Space;
      template<typename T> friend class L2Space;
      template<typename T> friend class HcurlSpace;
      template<typename T> friend class HdivSpace;
    };
  }
}
#endif
