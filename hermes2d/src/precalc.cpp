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
#include "quad.h"
#include "precalc.h"



PrecalcShapeset::PrecalcShapeset(Shapeset* shapeset)
               : RealFunction()
{
  assert_msg(shapeset != NULL, "Shapeset cannot be NULL.");
  this->shapeset = shapeset;
  master_pss = NULL;
  num_components = shapeset->get_num_components();
  assert(num_components == 1 || num_components == 2);
  update_max_index();
  set_quad_2d(&g_quad_2d_std);
}


PrecalcShapeset::PrecalcShapeset(PrecalcShapeset* pss)
               : RealFunction()
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
  shapeset->set_mode(HERMES_MODE_TRIANGLE);
  max_index[0] = shapeset->get_max_index();
  shapeset->set_mode(HERMES_MODE_QUAD);
  max_index[1] = shapeset->get_max_index();
}


void PrecalcShapeset::set_quad_2d(Quad2D* quad_2d)
{
  RealFunction::set_quad_2d(quad_2d);
}

void PrecalcShapeset::handle_overflow_idx()
{
  if(overflow_nodes != NULL) {
    for(std::map<unsigned int, Node*>::iterator it = overflow_nodes->begin(); it != overflow_nodes->end(); it++)
      ::free(it->second);
    delete overflow_nodes;
  }
  nodes = new std::map<unsigned int, Node*>;
  overflow_nodes = nodes;
}

void PrecalcShapeset::set_active_shape(int index)
{
  // Key creation.
  unsigned key = cur_quad | (mode << 3) | ((unsigned) (max_index[mode] - index) << 4);
  
  // Blank value.
  std::map<uint64_t, std::map<unsigned int, Node*>*>* updated_nodes = new std::map<uint64_t, std::map<unsigned int, Node*>*>;

  // Decision based on this PrecalcShapeset being a slave or not.
  std::map<unsigned int, std::map<uint64_t, std::map<unsigned int, Node*>*>*>* tab = (master_pss == NULL) ? &tables : &(master_pss->tables);

  if(tab->insert(make_pair(key, updated_nodes)).second == false)
    // If the value had existed.
    delete updated_nodes;
  
  // Get the proper sub-element tables.
  sub_tables = (*tab)[key];
  
  // Update the Node table.
  update_nodes_ptr();

  this->index = index;
  order = std::max(H2D_GET_H_ORDER(shapeset->get_order(index)), H2D_GET_V_ORDER(shapeset->get_order(index)));
}


void PrecalcShapeset::set_active_element(Element* e)
{
  mode = e->get_mode();
  shapeset->set_mode(mode);
  get_quad_2d()->set_mode(mode);
  element = e;
}


void PrecalcShapeset::set_mode(int mode)  // used in curved.cpp
{
  this->mode = mode;
  shapeset->set_mode(mode);
  get_quad_2d()->set_mode(mode);
  element = NULL;
}


void PrecalcShapeset::precalculate(int order, int mask)
{
  int i, j, k;

  // initialization
  Quad2D* quad = get_quad_2d();
  quad->set_mode(mode);
  H2D_CHECK_ORDER(quad, order);
  int np = quad->get_num_points(order);
  double3* pt = quad->get_points(order);

  int oldmask = (cur_node != NULL) ? cur_node->mask : 0;
  int newmask = mask | oldmask;
  Node* node = new_node(newmask, np);

  // precalculate all required tables
  for (j = 0; j < num_components; j++)
  {
    for (k = 0; k < 6; k++)
    {
      if (newmask & idx2mask[k][j]) {
        if (oldmask & idx2mask[k][j])
          memcpy(node->values[j][k], cur_node->values[j][k], np * sizeof(double));
        else
          for (i = 0; i < np; i++)
            node->values[j][k][i] = shapeset->get_value(k, index, ctm->m[0] * pt[i][0] + ctm->t[0],
                                                                  ctm->m[1] * pt[i][1] + ctm->t[1], j);
      }
    }
  }
  if((*nodes)[order] != NULL) {
    assert((*nodes)[order] == cur_node);
    ::free((*nodes)[order]);
  }
  (*nodes)[order] = node;
  cur_node = node;
}


void PrecalcShapeset::free()
{
  if (master_pss != NULL) return;
  std::map<unsigned int, std::map<uint64_t, std::map<unsigned int, Node*>*>*>::iterator it;
  for(it = tables.begin(); it != tables.end(); it++) {  
    std::map<uint64_t, std::map<unsigned int, Node*>*>::iterator it_inner;
    for (it_inner = it->second->begin(); it_inner != it->second->end(); it_inner++) {
      std::map<unsigned int, Node*>::iterator it_inner_inner;
      for (it_inner_inner = it_inner->second->begin(); it_inner_inner != it_inner->second->end(); it_inner_inner++)
        ::free(it_inner_inner->second);
      it_inner->second->clear();
    }
    it->second->clear();
  }
  if(overflow_nodes != NULL) {
    for(std::map<unsigned int, Node*>::iterator it = overflow_nodes->begin(); it != overflow_nodes->end(); it++)
      ::free(it->second);
    delete overflow_nodes;
  }
}

extern PrecalcShapeset ref_map_pss;

PrecalcShapeset::~PrecalcShapeset()
{
  free();
  {
    verbose("~PrecalcShapeset(): peak size of precalculated tables: %d B (%0.1lf MB)%s", max_mem,
            (double) max_mem / (1024 * 1024), (this == &ref_map_pss) ? " (refmap)" : "");
  }
}