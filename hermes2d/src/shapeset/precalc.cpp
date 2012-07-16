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
      if(shapeset == NULL)
        throw Exceptions::NullException(0);
      this->shapeset = shapeset;
      master_pss = NULL;
      num_components = shapeset->get_num_components();
      assert(num_components == 1 || num_components == 2);
      update_max_index();
      set_quad_2d(&g_quad_2d_std);
    }

    PrecalcShapeset::PrecalcShapeset(PrecalcShapeset* pss) : Function<double>()
    {
      while (pss->is_slave())
        pss = pss->master_pss;
      master_pss = pss;
      shapeset = pss->shapeset;
      num_components = pss->num_components;
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

    void PrecalcShapeset::handle_overflow_idx()
    {
      if(overflow_nodes != NULL)
      {
        for(unsigned int i = 0; i < overflow_nodes->get_size(); i++)
          if(overflow_nodes->present(i))
            ::free(overflow_nodes->get(i));
        delete overflow_nodes;
      }
      nodes = new LightArray<Node *>;
      overflow_nodes = nodes;
    }

    void PrecalcShapeset::set_active_shape(int index)
    {
      // Key creation.
      unsigned key = cur_quad | (element->get_mode() << 3) | ((unsigned) (max_index[element->get_mode()] - index) << 4);

      if(master_pss == NULL)
      {
        if(!tables.present(key))
          tables.add(new std::map<uint64_t, LightArray<Node*>*>, key);
        sub_tables = tables.get(key);
      }
      else
      {
        if(!master_pss->tables.present(key))
          master_pss->tables.add(new std::map<uint64_t, LightArray<Node*>*>, key);
        sub_tables = master_pss->tables.get(key);
      }

      // Update the Node table.
      update_nodes_ptr();

      this->index = index;
      order = std::max(H2D_GET_H_ORDER(shapeset->get_order(index, element->get_mode())), H2D_GET_V_ORDER(shapeset->get_order(index, element->get_mode())));
    }

    void PrecalcShapeset::set_active_element(Element* e)
    {
      Transformable::set_active_element(e);
    }

    void PrecalcShapeset::precalculate(int order, int mask)
    {
      int i, j, k;

      // initialization
      Quad2D* quad = get_quad_2d();
      int np = quad->get_num_points(order, this->element->get_mode());
      double3* pt = quad->get_points(order, this->element->get_mode());

      int oldmask = (cur_node != NULL) ? cur_node->mask : 0;
      int newmask = mask | oldmask;
      Node* node = new_node(newmask, np);

      // precalculate all required tables
      for (j = 0; j < num_components; j++)
      {
        for (k = 0; k < 6; k++)
        {
          if(newmask & idx2mask[k][j])
          {
            if(oldmask & idx2mask[k][j])
              memcpy(node->values[j][k], cur_node->values[j][k], np * sizeof(double));
            else
              for (i = 0; i < np; i++)
                node->values[j][k][i] = shapeset->get_value(k, index, ctm->m[0] * pt[i][0] + ctm->t[0],
                ctm->m[1] * pt[i][1] + ctm->t[1], j, element->get_mode());
          }
        }
      }
      if(nodes->present(order))
      {
        assert(nodes->get(order) == cur_node);
        ::free(nodes->get(order));
      }
      nodes->add(node, order);
      cur_node = node;
    }

    void PrecalcShapeset::free()
    {
      if(master_pss != NULL) return;

      for(unsigned int i = 0; i < tables.get_size(); i++)
        if(tables.present(i))
        {
          for(std::map<uint64_t, LightArray<Node*>*>::iterator it = tables.get(i)->begin(); it != tables.get(i)->end(); it++)
          {
            for(unsigned int k = 0; k < it->second->get_size(); k++)
              if(it->second->present(k))
                ::free(it->second->get(k));
            delete it->second;
          }
          delete tables.get(i);
        }

        if(overflow_nodes != NULL)
        {
          for(unsigned int i = 0; i < overflow_nodes->get_size(); i++)
            if(overflow_nodes->present(i))
              ::free(overflow_nodes->get(i));
          delete overflow_nodes;
        }
    }

    extern PrecalcShapeset ref_map_pss;

    PrecalcShapeset::~PrecalcShapeset()
    {
      free();
    }

    void PrecalcShapeset::push_transform(int son)
    {
      Transformable::push_transform(son);
      if(sub_tables != NULL)
        update_nodes_ptr();
    }

    void PrecalcShapeset::pop_transform()
    {
      Transformable::pop_transform();
      if(sub_tables != NULL)
        update_nodes_ptr();
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

    void PrecalcShapeset::set_master_transform()
    {
      assert(master_pss != NULL);
      sub_idx = master_pss->sub_idx;
      top = master_pss->top;
      stack[top] = *(master_pss->ctm);
      ctm = stack + top;
    }

    int PrecalcShapeset::get_edge_fn_order(int edge)
    {
      return H2D_MAKE_EDGE_ORDER(element->get_mode(), edge, shapeset->get_order(index, element->get_mode()));
    }

    bool PrecalcShapeset::is_slave() const
    {
      return master_pss != NULL;
    }

    void PrecalcShapeset::force_transform(uint64_t sub_idx, Trf* ctm)
    {
      this->sub_idx = sub_idx;
      this->ctm = ctm;
    }
  }
}