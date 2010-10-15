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

#include "mesh.h"

/// 2D transformatin.
struct Trf
{
  double2 m; ///< 2x2 diagonal transformation matrix
  double2 t; ///< translation vector
};

#define H2D_TRF_TRI_NUM 4 ///< A total number of valid transformation of a triangle to a sub-domain.
#define H2D_TRF_QUAD_NUM 8 ///< A total number of valid transformation of a quad to a sub-domain.
#define H2D_TRF_NUM (H2D_TRF_QUAD_NUM + 1) ///< A total number of transformations.
#define H2D_TRF_IDENTITY H2D_TRF_QUAD_NUM ///< An index of identity transformation.

extern HERMES_API Trf tri_trf[H2D_TRF_NUM];  ///< A table of triangle sub-subdomain transforms. Only first ::H2D_TRF_TRI_NUM transformations are valid, the rest are identity transformation.
extern HERMES_API Trf quad_trf[H2D_TRF_NUM]; ///< A table of quad sub-subdomain transforms. Only first ::H2D_TRF_QUAD_NUM transformations are valid, the rest are identity transformation.


/// Transformable is a base class for all classes that perform some kind of precalculation of
/// function values on elements. These classes (PrecalcShapeset, Solution, RefMap) inherit
/// from Transformable the ability to transform integration points to the sub-elements
/// of an element.
///
class HERMES_API Transformable
{
public:

  Transformable()
  {
    #ifndef NDEBUG
    memset(stack, 0, sizeof(stack));
    #endif
    reset_transform();
    element = NULL;
  }

  virtual ~Transformable() {}

  /// Called by the assembling procedure and by other functions. In PrecalcShapeset it
  /// sets an internal variable that can be later retrieved by get_active_element().
  /// In Solution it selects the element to retrieve solution values for, etc.
  /// \param e [in] Element associated with the function being represented by the class.
  virtual void set_active_element(Element* e) { element = e; }

  /// \return The element associated with the function being represented by the class.
  Element* get_active_element() const { return element; }

  /// Multiplies the current transformation matrix on the right by a transformation to the
  /// specified son element and pushes it on top of the matrix stack. All integration
  /// points will then be transformed to this sub-element. This process can be repeated.
  /// \param son [in] Son element number in the range [0-3] for triangles and [0-7] for quads.
  virtual void push_transform(int son)
  {
    assert(element != NULL);
    if (top >= 20) error("Too deep transform.");

    Trf* mat = stack + (++top);
    Trf* tr = (element->is_triangle() ? tri_trf + son : quad_trf + son);

    mat->m[0] = ctm->m[0] * tr->m[0];
    mat->m[1] = ctm->m[1] * tr->m[1];
    mat->t[0] = ctm->m[0] * tr->t[0] + ctm->t[0];
    mat->t[1] = ctm->m[1] * tr->t[1] + ctm->t[1];

    ctm = mat;
    sub_idx = (sub_idx << 3) + son + 1; // see traverse.cpp if this changes
  }

  /// Removes the current transformation matrix from the top of the stack. The new top becomes
  /// the current transformation matrix. This returns the transform to the state before the
  /// last push_transform() was performed.
  virtual void pop_transform()
  {
    assert(top > 0);
    ctm = stack + (--top);
    sub_idx = (sub_idx - 1) >> 3;
  }

  /// Sets the current transform at once as if it was created by multiple calls to push_transform().
  /// \param idx [in] The number of the sub-element, as returned by get_transform().
  void set_transform(uint64_t idx);

  /// \return The current transform index.
  uint64_t get_transform() const { return sub_idx; }

  /// Empties the stack, loads identity transform.
  void reset_transform()
  {
    stack[0].m[0] = stack[0].m[1] = 1.0;
    stack[0].t[0] = stack[0].t[1] = 0.0;
    ctm = stack;
    sub_idx = top = 0;
  }

  /// \return The jacobian of the current transformation matrix.
  double get_transform_jacobian() const { return ctm->m[0] * ctm->m[1]; }

  /// \return The current transformation matrix.
  Trf* get_ctm() const { return ctm; }

  /// \return The depth of the current transformation.
  int get_depth() const { return top; }


protected:

  Element* element; ///< the active element

  Trf* ctm;  ///< current sub-element transformation matrix
  uint64_t sub_idx; ///< sub-element transformation index. Data type is equal to the type used by Judy.
  //static const unsigned H2D_MAX_IDX = 0x49249248; ///< largest sub_idx for top <= 10
  static const uint64_t H2D_MAX_IDX = 0x4000; ///< largest sub_idx for top <= 10. Data type is equal to the type used by Judy.

  Trf stack[21]; ///< transformation matrix stack
  int top;       ///< stack top

};


#endif
