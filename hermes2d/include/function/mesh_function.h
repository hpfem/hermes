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
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __H2D_MESH_FUNCTION_H
#define __H2D_MESH_FUNCTION_H

#include "function.h"
#include "../mesh/refmap.h"
#include "../space/space.h"
#include "exceptions.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar> class MeshFunction;
  }
}

template<typename Scalar>
#ifdef _WINDOWS
class HERMES_API MeshFunctionSharedPtr : public std::shared_ptr<Hermes::Hermes2D::MeshFunction<Scalar> >
#else
class HERMES_API MeshFunctionSharedPtr : public std::tr1::shared_ptr<Hermes::Hermes2D::MeshFunction<Scalar> >
#endif
{
public:
  MeshFunctionSharedPtr(Hermes::Hermes2D::MeshFunction<Scalar>* ptr = NULL);

  MeshFunctionSharedPtr(const MeshFunctionSharedPtr<Scalar>& other);

  void operator=(const MeshFunctionSharedPtr<Scalar>& other);

  Hermes::Hermes2D::Solution<Scalar>* get_solution();
      
  ~MeshFunctionSharedPtr();
};


namespace Hermes
{
  namespace Hermes2D
  {
    /** \defgroup meshFunctions Mesh functions
    * \brief Collection of classes that represent various functions of the mesh coordinates, i.e. defined on the Mesh.
    * These comprise solutions, exact &amp; initial solutions, filters (functions of the solutions) etc.
    */

    /// @ingroup meshFunctions
    /// \brief Represents a function defined on a mesh.
    ///
    /// MeshFunction is a base class for all classes representing an arbitrary function
    /// superimposed on a mesh (ie., domain). These include the Solution, ExactSolution
    /// and Filter classes, which define the concrete behavior and the way the function
    /// is (pre)calculated. Any such function can later be visualized.
    ///
    /// (This is an abstract class and cannot be instantiated.)
    ///
    template<typename Scalar>
    class HERMES_API MeshFunction : public Function<Scalar>, public Hermes::Hermes2D::Mixins::StateQueryable
    {
    public:
      /// Empty constructor.
      MeshFunction();

      /// Constructor.
      MeshFunction(MeshSharedPtr mesh);

      /// Destructor.
      virtual ~MeshFunction();

      /// Return the mesh.
      MeshSharedPtr get_mesh() const;

      /// Copy from sln to this instance.
      virtual void copy(const MeshFunction<Scalar>* sln);
      virtual void copy(MeshFunctionSharedPtr<Scalar> sln);

      /// Return the reference mapping.
      RefMap* get_refmap(bool update = true);

      /// Return the value at the coordinates x,y.
      virtual Func<Scalar>* get_pt_value(double x, double y, bool use_MeshHashGrid = false, Element* e = NULL) = 0;

      /// Cloning function - for parallel OpenMP blocks.
      /// Designed to return an identical clone of this instance.
      virtual MeshFunction<Scalar>* clone() const
      {
        throw Hermes::Exceptions::Exception("You need to implement MeshFunctionSharedPtr::clone() to be able to use paralellization");
        return NULL;
      }
      
      /// Multiplies the function represented by this class by the given coefficient.
      virtual void multiply(Scalar coef);

      /// Adds another mesh function on the given space.
      /// ! Resulting mesh function is a solution.
      virtual void add(MeshFunctionSharedPtr<Scalar> other_mesh_function, SpaceSharedPtr<Scalar> target_space);

      /// Return the approximate maximum value of this instance.
      virtual Scalar get_approx_max_value(int item = H2D_FN_VAL_0);

      /// Return the approximate minimum value of this instance.
      virtual Scalar get_approx_min_value(int item = H2D_FN_VAL_0);

      /// State querying helpers.
      virtual bool isOkay() const;

      /// Internal.
      inline std::string getClassName() const { return "MeshFunction"; }

      /// Internal.
      virtual void init();

      virtual void free();

      /// Internal.
      virtual void reinit();

      /// Set the quadrature rule.
      /// Internal.
      virtual void set_quad_2d(Quad2D* quad_2d);

      /// Set the active element.
      /// Internal.
      virtual void set_active_element(Element* e);

      /// See Transformable::push_transform.
      /// Internal.
      virtual void push_transform(int son);

      /// See Transformable::pop_transform.
      /// Internal.
      virtual void pop_transform();

      /// Set the reference mapping.
      /// Internal.
      void set_refmap(RefMap* refmap_to_set);

      /// Returns the order of the edge number edge of the current active element.
      virtual int get_edge_fn_order(int edge);
    protected:
      ElementMode2D mode;
      MeshSharedPtr mesh;
      RefMap* refmap;

      void force_transform(MeshFunctionSharedPtr<Scalar> mf);

      void update_refmap();

      void force_transform(uint64_t sub_idx, Trf* ctm);

      friend class RefMap;
      template<typename T> friend class KellyTypeAdapt;
      template<typename T> friend class Adapt;

      template<typename T> friend class Func;
      template<typename T> friend class Geom;

      template<typename T> friend HERMES_API Func<T>* init_fn(MeshFunction<T>*fu, const int order);

      template<typename T> friend class DiscontinuousFunc;
      template<typename T> friend class DiscreteProblem;
      template<typename T> friend class NeighborSearch;
    };
  }
}

#endif
