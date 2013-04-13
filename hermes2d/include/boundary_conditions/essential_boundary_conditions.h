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

    /// Base, abstract class representing Essential boundary condition of the form u|_{\Gamma_Essential} = u_Essential.
    /// Internal.
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
      /// \param[in] x x-coordinate of the point where the value is evaluated.
      /// \param[in] y y-coordinate of the point where the value is evaluated.
      /// \param[in] n_x the x-component of the unit outer normal.
      /// \param[in] n_y the y-component of the unit outer normal.
      /// \param[in] t_x the x-component of the tangent(perpendicular to normal).
      /// \param[in] t_y the y-component of the tangent(perpendicular to normal).
      virtual Scalar value(double x, double y, double n_x, double n_y, double t_x, double t_y) const = 0;

      /// Set the current time for time-dependent boundary conditions.
      void set_current_time(double time);

      /// Get the current time for time-dependent boundary conditions.
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

      Scalar value(double x, double y, double n_x, double n_y, double t_x, double t_y) const;

      /// Function giving info that u_Essential is a constant.
      inline typename EssentialBoundaryCondition<Scalar>::EssentialBCValueType get_value_type() const { return EssentialBoundaryCondition<Scalar>::BC_CONST; }
    };

    /// Class representing non-constant essential boundary condition for Scalar approximation.
    /// Typical usage - this example is a non-const Dirichlet boundary condition for the incompressible Navier-Stokes equations.
    /// class MyEssentialBCNonConst : public DefaultEssentialBCNonConst<double>
    /// {
    /// public:
    ///&nbsp;// Constructor with multiple markers, setting velocity at inlet(vel_inlet), domain height(H) and startup time(startup_time) - all parameters are of this example - DERIVED - class.
    ///&nbsp;MyEssentialBCNonConst(Hermes::vector<std::string> markers, double vel_inlet, double H, double startup_time) : 
    ///&nbsp;  EssentialBoundaryCondition<double>(markers), vel_inlet(vel_inlet), H(H), startup_time(startup_time) {};
    ///
    ///  // VERY IMPORTANT - overriding the method of the base class (DefaultEssentialBCNonConst::value) with a custom implementation.
    ///  // NOTE - one can use the top-level base class (EssentialBoundaryCondition)'s methods for handling the time variable for time-dependent problems: get_current_time().
    ///  // NOTE - the 'virtual' keyword is not here anymore - because we will not need to further derive from this class and override this method.
    ///  double value(double x, double y, double n_x, double n_y, double t_x, double t_y) const {
    ///&nbsp; double val_y = vel_inlet * y*(H-y) / (H/2.)/(H/2.);
    ///&nbsp; if (get_current_time() <= startup_time) 
    ///&nbsp;   return val_y * get_current_time()/startup_time;
    ///&nbsp; else 
    ///&nbsp;   return val_y;
    ///  };
    ///
    /// protected:
    ///&nbsp;// Members of MyEssentialBCNonConst.
    ///&nbsp;double vel_inlet;
    ///&nbsp;double H;
    ///&nbsp;double startup_time;
    ///}; 

    template<typename Scalar>
    class HERMES_API DefaultEssentialBCNonConst : public EssentialBoundaryCondition<Scalar>
    {
    public:
      DefaultEssentialBCNonConst(Hermes::vector<std::string> markers_,
        MeshFunctionSharedPtr<Scalar> exact_solution);

      DefaultEssentialBCNonConst(std::string marker, MeshFunctionSharedPtr<Scalar> exact_solution);

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
        MeshFunctionSharedPtr<Scalar> exact_solution2);
      DefaultEssentialBCNonConstHcurl(std::string marker, MeshFunctionSharedPtr<Scalar> exact_solution2);

      ~DefaultEssentialBCNonConstHcurl() {};

      virtual Scalar value(double x, double y, double n_x, double n_y, double t_x, double t_y) const;

      /// Function giving info that u_Essential is a non-constant function.
      inline typename EssentialBoundaryCondition<Scalar>::EssentialBCValueType get_value_type() const { return EssentialBoundaryCondition<Scalar>::BC_FUNCTION; }

      ExactSolutionVector<Scalar>* exact_solution2;
    };

    /// Class encapsulating all boundary conditions of one problem.
    /// Using the class EssentialBCs and its descendants.
    /// Usage: for passing to Hermes2D::Space in the constructor or set_essential_bcs() method.
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
