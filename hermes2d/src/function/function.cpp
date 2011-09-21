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
    int Function<Scalar>::get_fn_order() const
    { 
      return order;
    }

    template<typename Scalar>
    int Function<Scalar>::get_edge_fn_order(int edge) const
    {
      return order;
    }

    template<typename Scalar>
    int Function<Scalar>::get_num_components() const
    {
      return num_components;
    }

    template<typename Scalar>
    void Function<Scalar>::set_quad_order(unsigned int order, int mask)
    {
      if(nodes->present(order)) {
        cur_node = nodes->get(order);
        // If the mask has changed.
        if((cur_node->mask & mask) != mask) {
          precalculate(order, mask);
          nodes->add(cur_node, order);
        }
      }
      else {
        // The value had not existed.
        cur_node = NULL;
        precalculate(order, mask);
        nodes->add(cur_node, order);
      }
    }

    template<typename Scalar>
    Scalar* Function<Scalar>::get_fn_values(int component)
    {
      check_params(component, cur_node, num_components); check_table(component, cur_node, 0, "Function values");
      return cur_node->values[component][0];
    }

    template<typename Scalar>
    Scalar* Function<Scalar>::get_dx_values(int component)
    {
      check_params(component, cur_node, num_components); check_table(component, cur_node, 1, "DX values");
      return cur_node->values[component][1];
    }

    template<typename Scalar>
    Scalar* Function<Scalar>::get_dy_values(int component)
    {
      check_params(component, cur_node, num_components); check_table(component, cur_node, 2, "DY values");
      return cur_node->values[component][2];
    }

    template<typename Scalar>
    void Function<Scalar>::get_dx_dy_values(Scalar*& dx, Scalar*& dy, int component)
    {
      check_params(component, cur_node, num_components); check_table(component, cur_node, 1, "DX values"); check_table(component, cur_node, 2, "DY values");
      dx = cur_node->values[component][1];
      dy = cur_node->values[component][2];
    }

    template<typename Scalar>
    Scalar* Function<Scalar>::get_dxx_values(int component)
    {
      check_params(component, cur_node, num_components); check_table(component, cur_node, 3, "DXX values");
      return cur_node->values[component][3];
    }

    template<typename Scalar>
    Scalar* Function<Scalar>::get_dyy_values(int component)
    {
      check_params(component, cur_node, num_components); check_table(component, cur_node, 4, "DYY values");
      return cur_node->values[component][4];
    }

    template<typename Scalar>
    Scalar* Function<Scalar>::get_dxy_values(int component)
    {
      check_params(component, cur_node, num_components); check_table(component, cur_node, 5, "DXY values");
      return cur_node->values[component][5];
    }

    template<typename Scalar>
    Scalar* Function<Scalar>::get_values(int a, int b)
    {
      return cur_node->values[a][b];
    }

    template<typename Scalar>
    void Function<Scalar>::set_quad_2d(Quad2D* quad_2d)
    {
      int i;

      // check to see if we already have the quadrature
      for (i = 0; i < 4; i++)
        if (quads[i] == quad_2d)
        {
          cur_quad = i;
          return;
        }

      // if not, add the quadrature to a free slot
      for (i = 0; i < 4; i++)
        if (quads[i] == NULL)
        {
          quads[i] = quad_2d;
          cur_quad = i;
          return;
        }

      error("too many quadratures.");
    }

    template<typename Scalar>
    Quad2D* Function<Scalar>::get_quad_2d() const 
    {
      return quads[cur_quad]; 
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
      int size = (sizeof(Node) - sizeof(Scalar)) + sizeof(Scalar) * num_points * nt; //Due to impl. reasons, the structure Node has non-zero length of data even though they can be zero.
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

    template<typename Scalar>
    void Function<Scalar>::update_nodes_ptr()
    {
      if (sub_idx > H2D_MAX_IDX)
        handle_overflow_idx();
      else {
        if(sub_tables->find(sub_idx) == sub_tables->end())
          sub_tables->insert(std::pair<uint64_t, LightArray<Node*>*>(sub_idx, new LightArray<Node*>));
        nodes = sub_tables->find(sub_idx)->second;
      }
    }

    template<typename Scalar>
    void Function<Scalar>::force_transform(uint64_t sub_idx, Trf* ctm)
    {
      this->sub_idx = sub_idx;
      this->ctm = ctm;
      update_nodes_ptr();
    }

    template<typename Scalar>
    void Function<Scalar>::replace_cur_node(Node* node)
    {
      if (node == NULL) throw Exceptions::NullException(1);
      if (cur_node != NULL) {
        total_mem -= cur_node->size;
        ::free(cur_node);
      }
      cur_node = node;
    }

    template<typename Scalar>
    void Function<Scalar>::check_order(Quad2D* quad, int order)
    {
      if (order < 0 || order >= quad->get_num_tables())
        error("Hermes::Order out of range (%d, %d).", order, quad->get_num_tables());
    }

    template class HERMES_API Function<double>;
    template class HERMES_API Function<std::complex<double> >;
  }
}
