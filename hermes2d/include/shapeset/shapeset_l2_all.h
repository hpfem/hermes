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

#ifndef __H2D_SHAPESET_L2_ALL
#define __H2D_SHAPESET_L2_ALL

// This file is a common header for all L2 shapesets.

#include "shapeset.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// L2 shapeset - products of legendre polynomials
    /// @ingroup spaces
    class HERMES_API L2ShapesetLegendre : public Shapeset
    {
    public:
      L2ShapesetLegendre();
      virtual Shapeset* clone() { return new L2ShapesetLegendre(*this); };
      virtual SpaceType get_space_type() const { return HERMES_L2_SPACE; }
      virtual int get_max_index(ElementMode2D mode);
      virtual int get_id() const { return 30; }
      
      static const int max_index[H2D_NUM_MODES];
    };

    /// L2 Taylor shapeset - Taylor basis functions as proposed by Kuzmin, Luo
    /// @ingroup spaces
    class HERMES_API L2ShapesetTaylor : public Shapeset
    {
    public:
      L2ShapesetTaylor(bool contains_means = true);
      virtual Shapeset* clone() { return new L2ShapesetTaylor(*this); };
      virtual SpaceType get_space_type() const { return HERMES_L2_SPACE; }
      virtual int get_max_index(ElementMode2D mode);
      virtual int get_id() const { return 31; }

      /// Returns a complete set of indices of bubble functions for an element of the given order.
      /// Reimplemented because this shapeset uses linear (not bi-linear), quadratic (not bi-quadratic) etc. polynomials.
      int* get_bubble_indices(int order, ElementMode2D mode) const;

      /// Returns the number of bubble functions for an element of the given order.
      /// Reimplemented because this shapeset uses linear (not bi-linear), quadratic (not bi-quadratic) etc. polynomials.
      virtual int get_num_bubbles(int order, ElementMode2D mode) const;
      
      static const int max_index[H2D_NUM_MODES];
    };

    /// This is the default shapeset typedef
    typedef L2ShapesetLegendre L2Shapeset;
  }
}
#endif