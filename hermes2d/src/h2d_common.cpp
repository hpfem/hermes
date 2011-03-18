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

#include "h2d_common.h"

const std::string Hermes2D::get_quad_order_str(const int quad_order) {
  std::stringstream str;
  str << "(H:" << H2D_GET_H_ORDER(quad_order) << ";V:" << H2D_GET_V_ORDER(quad_order) << ")";
  return str.str();
}

int Hermes2D::make_edge_order(int mode, int edge, int encoded_order)
{
  assert(edge < 4);

  if (mode == HERMES_MODE_TRIANGLE || edge == 0 || edge == 2)
    return H2D_GET_H_ORDER(encoded_order);
  else
    return H2D_GET_V_ORDER(encoded_order);
}

//// python support //////////////////////////////////////////////////////////////////////////////////
HERMES_API void throw_exception(char *text)
{
  throw std::runtime_error(text);
}
