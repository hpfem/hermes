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
#include "../mesh/mesh.h"
#include "exceptions.h"

namespace Hermes
{
  namespace Hermes2D
  {
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

      MeshFunction();
      MeshFunction(const Mesh *mesh);
      virtual ~MeshFunction();

      /// State querying helpers.
      virtual bool isOkay() const;
      inline std::string getClassName() const { return "MeshFunction"; }

      virtual void init();
      virtual void reinit();

      virtual void set_quad_2d(Quad2D* quad_2d);

      virtual void set_active_element(Element* e);

      const Mesh* get_mesh() const;
      RefMap* get_refmap(bool update = true);
      void set_refmap(RefMap* refmap_to_set);

      virtual int get_edge_fn_order(int edge);

      /// Returns the value at the coordinates x,y.
      virtual Func<Scalar>* get_pt_value(double x, double y) = 0;

      /// Virtual function handling overflows. Has to be virtual, because
      /// the necessary iterators in the templated class do not work with GCC.
      virtual void handle_overflow_idx();

      /// See Transformable::push_transform.
      virtual void push_transform(int son);

      virtual void pop_transform();

      virtual MeshFunction<Scalar>* clone() const
      {
        throw Hermes::Exceptions::Exception("You need to implement MeshFunction::clone() to be able to use paralellization");
        return NULL;
      }

    protected:
      ElementMode2D mode;
      const Mesh* mesh;
      RefMap* refmap;
    
      void force_transform(MeshFunction<Scalar>* mf);

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
      template<typename T> friend class DiscreteProblemLinear;
      template<typename T> friend class NeighborSearch;
    };
  }
}

#endif
