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

#include "function.h"
namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    Function<Scalar>::Function()
      : Transformable()
    {
      order = 0;
      max_mem = total_mem = 0;
      cur_node = NULL;
      sub_tables = NULL;
      nodes = NULL;
      overflow_nodes = NULL;
      memset(quads, 0, sizeof(quads));
    }

    template<typename Scalar>
    void Function<Scalar>::set_quad_2d(Quad2D* quad_2d)
    {
      int i;

      // check to see if we already have the quadrature
      for (i = 0; i < 8; i++)
        if (quads[i] == quad_2d) {
          cur_quad = i;
          return;
        }

        // if not, add the quadrature to a free slot
        for (i = 0; i < 8; i++)
          if (quads[i] == NULL) {
            quads[i] = quad_2d;
            cur_quad = i;
            return;
          }

          error("too many quadratures.");
    }

    template<typename Scalar>
    int Function<Scalar>::idx2mask[6][2] =
    {
      { H2D_FN_VAL_0, H2D_FN_VAL_1 }, { H2D_FN_DX_0,  H2D_FN_DX_1  }, { H2D_FN_DY_0,  H2D_FN_DY_1  },
      { H2D_FN_DXX_0, H2D_FN_DXX_1 }, { H2D_FN_DYY_0, H2D_FN_DYY_1 }, { H2D_FN_DXY_0, H2D_FN_DXY_1 }
    };


    template<typename Scalar>
    typename Function<Scalar>::Node* Function<Scalar>::new_node(int mask, int num_points)
    {
      // get the number of tables
      int nt = 0, m = mask;
      if (num_components < 2) m &= H2D_FN_VAL_0 | H2D_FN_DX_0 | H2D_FN_DY_0 | H2D_FN_DXX_0 | H2D_FN_DYY_0 | H2D_FN_DXY_0;
      while (m) { nt += m & 1; m >>= 1; }

      // allocate a node including its data part, init table pointers
      int size = H2D_Node_HDR_SIZE + sizeof(Scalar) * num_points * nt; //Due to impl. reasons, the structure Node has non-zero length of data even though they can be zero.
      Node* node = (Node*) malloc(size);
      node->mask = mask;
      node->size = size;
      memset(node->values, 0, sizeof(node->values));
      Scalar* data = node->data;
      for (int j = 0; j < num_components; j++) {
        for (int i = 0; i < 6; i++)
          if (mask & idx2mask[i][j]) {
            node->values[j][i] = data;
            data += num_points;
          }
      }

      total_mem += size;
      if (max_mem < total_mem) max_mem = total_mem;
      return node;
    }

    template class HERMES_API Function<double>;
    template class HERMES_API Function<std::complex<double> >;
  }
}