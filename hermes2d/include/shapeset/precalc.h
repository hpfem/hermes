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

#ifndef __H2D_PRECALC_H
#define __H2D_PRECALC_H

#include "../function/function.h"
#include "../shapeset/shapeset.h"

namespace Hermes
{
  namespace Hermes2D
  {
    enum SpaceType;
    /// @ingroup meshFunctions
    /// \brief Caches precalculated shape function values.
    ///
    /// PrecalcShapeset is a cache of precalculated shape function values.
    ///
    ///
    class HERMES_API PrecalcShapeset : public Function<double>
    {
    public:
      /// Returns type of space
      SpaceType get_space_type() const;

      /// \brief Constructs a standard (master) precalculated shapeset class.
      /// \param shapeset[in] Pointer to the shapeset to be precalculated.
      PrecalcShapeset(Shapeset* shapeset);

      /// \brief Constructs a slave precalculated shapeset class.
      /// \details The slave instance does not hold any precalculated tables.
      /// Instead, it refers to those contained in the master instance. However,
      /// the slave can have different shape function active, different transform
      /// selected, etc. Slave pss's are used for test functions when calling
      /// bilinear forms, inside Solution so as not to disrupt user's pss, etc.
      /// \param master_pss[in] Master precalculated shapeset pointer.
      PrecalcShapeset(PrecalcShapeset* master_pss);

      /// Destructor.
      virtual ~PrecalcShapeset();

      /// Ensures subsequent calls to get_active_element() will be returning 'e'.
      /// Switches the class to the appropriate mode (triangle, quad).
      virtual void set_active_element(Element* e);

      /// Activates a shape function given by its index. The values of the shape function
      /// can then be obtained by setting the required integration rule order by calling
      /// set_quad_order() and after that calling get_values(), get_dx_values(), etc.
      /// \param index[in] Shape index.
      void set_active_shape(int index);

    private:
      virtual void set_quad_2d(Quad2D* quad_2d);

      /// \brief Frees all precalculated tables.
      virtual void free();

      /// Virtual function handling overflows. Has to be virtual, because
      /// the necessary iterators in the templated class do not work with GCC.
      virtual void handle_overflow_idx();


      /// Returns the index of the active shape (can be negative if the shape is constrained).
      int get_active_shape() const;

      /// Returns a pointer to the shapeset which is being precalculated.
      Shapeset* get_shapeset() const;

      /// For internal use only.
      void set_master_transform();

      /// Returns the polynomial order of the active shape function on given edge.
      virtual int get_edge_fn_order(int edge);

      /// See Transformable::push_transform.
      virtual void push_transform(int son);

      virtual void pop_transform();

      Shapeset* shapeset;

      /// Main structure.
      /// There is a 3-layer structure of the precalculated tables.
      /// The first (the lowest) one is the layer where mapping of integral orders to
      /// Function::Node takes place. See function.h for details.
      /// The second one is the layer with mapping of sub-element transformation to
      /// a table from the lowest layer.
      /// The highest and most complicated one maps a key formed by
      /// quadrature table selector (0-7), mode of the shape function (triangle/quad),
      /// and shape function index to a table from the middle layer.
      LightArray<std::map<uint64_t, LightArray<Node*>*>*> tables;

      int index;

      int max_index[2];

      PrecalcShapeset* master_pss;

      /// Returns true iff this is a precalculated shapeset for test functions.
      bool is_slave() const;

      virtual void precalculate(int order, int mask);

      void update_max_index();

      /// Forces a transform without using push_transform() etc.
      /// Used by the Solution class.
      void force_transform(uint64_t sub_idx, Trf* ctm);

      friend class RefMap;
      template<typename T> friend class KellyTypeAdapt;
      template<typename T> friend class Adapt;
      template<typename T> friend class Func;
      template<typename T> friend class Solution;
      template<typename T> friend class DiscontinuousFunc;
      template<typename T> friend class DiscreteProblem;
      template<typename T> friend class NeighborSearch;
      friend class CurvMap;
    };
  }
}
#endif
