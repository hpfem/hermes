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
#include "global.h"
#include "order_permutator.h"
namespace Hermes
{
  namespace Hermes2D
  {
    namespace RefinementSelectors
    {
      OrderPermutator::OrderPermutator(int start_quad_order, int end_quad_order, bool iso_p, int* tgt_quad_order)
        : start_order_h(H2D_GET_H_ORDER(start_quad_order)), start_order_v(H2D_GET_V_ORDER(start_quad_order))
        , end_order_h(H2D_GET_H_ORDER(end_quad_order)), end_order_v(H2D_GET_V_ORDER(end_quad_order))
        , iso_p(iso_p), tgt_quad_order(tgt_quad_order)
      {
        reset();
      }

      bool OrderPermutator::next()
      {
        if(iso_p)
        {
          if(order_h >= end_order_h || order_v >= end_order_v)
            return false;

          order_h++;
          order_v++;
        }
        else
        {
          if(order_h >= end_order_h && order_v >= end_order_v)
            return false;

          order_h++;
          if(order_h > end_order_h)
          {
            order_h = start_order_h;
            order_v++;
          }
        }

        if(tgt_quad_order != NULL)
          *tgt_quad_order = H2D_MAKE_QUAD_ORDER(order_h, order_v);
        return true;
      }

      void OrderPermutator::reset()
      {
        order_h = start_order_h;
        order_v = start_order_v;
        if(tgt_quad_order != NULL)
          *tgt_quad_order = H2D_MAKE_QUAD_ORDER(order_h, order_v);
      }

      int OrderPermutator::get_order_h() const
      {
        return order_h;
      }

      int OrderPermutator::get_order_v() const
      {
        return order_v;
      }

      int OrderPermutator::get_quad_order() const
      {
        return H2D_MAKE_QUAD_ORDER(order_h, order_v);
      }

      int OrderPermutator::get_start_quad_order() const
      {
        return H2D_MAKE_QUAD_ORDER(start_order_h, start_order_v);
      }

      int OrderPermutator::get_end_quad_order() const
      {
        return H2D_MAKE_QUAD_ORDER(end_order_h, end_order_v);
      }
    }
  }
}