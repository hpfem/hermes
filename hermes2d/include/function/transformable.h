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

#ifndef __H2D_TRANSFORM_H
#define __H2D_TRANSFORM_H

#include "../global.h"
#include "exceptions.h"
namespace Hermes
{
  namespace Hermes2D
  {
    class Element;
    template<typename Scalar> class Func;
    namespace Views
    {
      class Vectorizer;
    }

    /// @defgroup meshFunctions Mesh functions

    /// 2D transformation.
    struct Trf
    {
      double2 m; /// The 2x2 diagonal transformation matrix.
      double2 t; /// Translation vector.
    };

    /// A table of triangle sub-subdomain transforms. Only first ::H2D_TRF_TRI_NUM transformations are valid, the rest are identity transformation.
    extern HERMES_API Trf tri_trf[H2D_TRF_NUM];
    /// A table of quad sub-subdomain transforms. Only first ::H2D_TRF_QUAD_NUM transformations are valid, the rest are identity transformation.
    extern HERMES_API Trf quad_trf[H2D_TRF_NUM];

    /// @ingroup meshFunctions
    /// Transformable is a base class for all classes that perform some kind of precalculation of
    /// function values on elements. These classes (PrecalcShapeset, Solution, RefMap) inherit
    /// from Transformable the ability to transform integration points to the sub-elements
    /// of an element.
    class HERMES_API Transformable : public Hermes::Mixins::Loggable
    {
    public:
      /// \return The element associated with the function being represented by the class.
      Element* get_active_element() const;

      /// Sets the current transform at once as if it was created by multiple calls to push_transform().
      /// \param idx[in] The number of the sub-element, as returned by get_transform().
      void set_transform(uint64_t idx);

      /// \return The current transform index.
      uint64_t get_transform() const;

      virtual ~Transformable();

      /// Multiplies the current transformation matrix on the right by a transformation to the
      /// specified son element and pushes it on top of the matrix stack. All integration
      /// points will then be transformed to this sub-element. This process can be repeated.
      /// \param son[in] Son element number in the range[0-3] for triangles and[0-7] for quads.
      virtual void push_transform(int son);

    protected:

      Transformable();

      /// Called by the assembling procedure and by other functions. In PrecalcShapeset it
      /// sets an internal variable that can be later retrieved by get_active_element().
      /// In Solution it selects the element to retrieve solution values for, etc.
      /// \param e[in] Element associated with the function being represented by the class.
      virtual void set_active_element(Element* e);

      /// Removes the current transformation matrix from the top of the stack. The new top becomes
      /// the current transformation matrix. This returns the transform to the state before the
      /// last push_transform() was performed.
      virtual void pop_transform();

      /// Empties the stack, loads identity transform.
      void reset_transform();

      /// \return The jacobian of the current transformation matrix.
      inline double get_transform_jacobian() const { return ctm->m[0] * ctm->m[1]; }

      /// \return The current transformation matrix.
      inline Trf* get_ctm() const { return ctm; }

      /// \return The depth of the current transformation.
      inline unsigned int get_depth() const { return top; }

      static void push_transforms(std::set<Transformable *>& transformables, int son);
      static void pop_transforms(std::set<Transformable *>& transformables);
      static const unsigned int H2D_MAX_TRN_LEVEL = 15;

      /// The active element.
      Element* element;

      /// Current sub-element transformation matrix.
      Trf* ctm;
      /// Sub-element transformation index.
      uint64_t sub_idx;

      /// The largest sub_idx for top <= 10.
      /// FIXME: Why it is only 0x4000?
      static const uint64_t H2D_MAX_IDX = (1ULL << 3 * H2D_MAX_TRN_LEVEL) - 1;

      /// Transformation matrix stack.
      Trf stack[21];
      /// Stack top.
      unsigned int top;

      template<typename T> friend class KellyTypeAdapt;
      template<typename T> friend class Adapt;
      template<typename T> friend class Func;
      template<typename T> friend class Function;
      template<typename T> friend class DiscontinuousFunc;
      template<typename T> friend class DiscreteProblem;
      template<typename T> friend class DiscreteProblemLinear;
      template<typename T> friend class NeighborSearch;
      friend class CurvMap;
      friend class Traverse;
      friend class Views::Vectorizer;
    };
  }
}
#endif