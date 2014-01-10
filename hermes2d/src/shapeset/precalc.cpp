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
#include "quad_all.h"
#include "precalc.h"
#include "mesh.h"
namespace Hermes
{
  namespace Hermes2D
  {
    PrecalcShapeset::PrecalcShapeset(Shapeset* shapeset) : Function<double>()
    {
      if (shapeset == nullptr)
        throw Exceptions::NullException(0);
      this->shapeset = shapeset;
      num_components = shapeset->get_num_components();
      assert(num_components == 1 || num_components == 2);
      update_max_index();
      set_quad_2d(&g_quad_2d_std);
    }

    void PrecalcShapeset::update_max_index()
    {
      max_index[0] = shapeset->get_max_index(HERMES_MODE_TRIANGLE);
      max_index[1] = shapeset->get_max_index(HERMES_MODE_QUAD);
    }

    void PrecalcShapeset::set_quad_2d(Quad2D* quad_2d)
    {
      Function<double>::set_quad_2d(quad_2d);
    }

    void PrecalcShapeset::set_active_shape(int index)
    {
      // Update the Node table.
      update_nodes_ptr();

      this->index = index;

      this->order = std::max(H2D_GET_H_ORDER(shapeset->get_order(index, element->get_mode())), H2D_GET_V_ORDER(shapeset->get_order(index, element->get_mode())));
    }

    void PrecalcShapeset::set_active_element(Element* e)
    {
      Transformable::set_active_element(e);
    }

    void PrecalcShapeset::precalculate(int order, int mask)
    {
      int i, j, k;

      int np = this->quads[cur_quad]->get_num_points(order, this->element->get_mode());
      double3* pt = this->quads[cur_quad]->get_points(order, this->element->get_mode());

      // precalculate all required tables
      for (j = 0; j < num_components; j++)
      {
        for (k = 0; k < 6; k++)
        {
          if (mask & idx2mask[k][j])
          {
            for (i = 0; i < np; i++)
              cur_node.values[j][k][i] = shapeset->get_value(k, index, ctm->m[0] * pt[i][0] + ctm->t[0],
              ctm->m[1] * pt[i][1] + ctm->t[1], j, element->get_mode());
          }
        }
      }
    }

    void PrecalcShapeset::free()
    {
    }

    extern PrecalcShapeset ref_map_pss;

    PrecalcShapeset::~PrecalcShapeset()
    {
      free();
    }

    int PrecalcShapeset::get_active_shape() const
    {
      return index;
    };

    Shapeset* PrecalcShapeset::get_shapeset() const
    {
      return shapeset;
    }

    SpaceType PrecalcShapeset::get_space_type() const
    {
      return shapeset->get_space_type();
    }

    int PrecalcShapeset::get_edge_fn_order(int edge)
    {
      return H2D_MAKE_EDGE_ORDER(element->get_mode(), edge, shapeset->get_order(index, element->get_mode()));
    }

    void PrecalcShapeset::force_transform(uint64_t sub_idx, Trf* ctm)
    {
      this->sub_idx = sub_idx;
      this->ctm = ctm;
    }
  }
}