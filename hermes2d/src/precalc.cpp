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
  tables = NULL;
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
  tables = NULL;
  update_max_index();
  set_quad_2d(&g_quad_2d_std);
}


void PrecalcShapeset::update_max_index()
{
  shapeset->set_mode(H2D_MODE_TRIANGLE);
  max_index[0] = shapeset->get_max_index();
  shapeset->set_mode(H2D_MODE_QUAD);
  max_index[1] = shapeset->get_max_index();
}


void PrecalcShapeset::set_quad_2d(Quad2D* quad_2d)
{
  //((master_pss != NULL) ? master_pss : this)->set_quad_2d(quad_2d);
  RealFunction::set_quad_2d(quad_2d);

  mode = 0;
  set_active_shape(0); // FIXME: this depends a lot on the order of calling
                       // set_quad_2d, set_active_element, set_active_shape.
                       // This is because segfaults in push_transform called right
                       // after set_quad_2d and set_active_element, without
                       // set_active_shape.
}


void PrecalcShapeset::set_active_shape(int index)
{
  // Each precalculated table is accessed and uniquely identified by the
  // following seven items:
  //
  //   - cur_quad:  quadrature table selector (0-7)
  //   - mode:      mode of the shape function (triangle/quad)
  //   - index:     shape function index
  //   - sub_idx:   the index of the sub-element
  //   - order:     integration rule order
  //   - component: shape function component (0-1)
  //   - val/d/dd:  values, dx, dy, ddx, ddy (0-4)
  //
  // The table database is implemented as a three-way chained Judy array.
  // The key to the first Judy array ('tables') is formed by cur_quad,
  // mode and index. This gives a pointer to the second Judy array, which
  // is indexed solely by sub_idx. The last Judy array is the node table,
  // understood by the base class and indexed by order. The component and
  // val/d/dd indices are used directly in the Node structure.

  unsigned key = cur_quad | (mode << 3) | ((unsigned) (max_index[mode] - index) << 4);
  void** tab = (master_pss == NULL) ? &tables : &(master_pss->tables);
  sub_tables = (void**) JudyLIns(tab, key, NULL);
  update_nodes_ptr();

  this->index = index;
  order = shapeset->get_order(index);
  order = std::max(H2D_GET_H_ORDER(order), H2D_GET_V_ORDER(order));
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

  // remove the old node and attach the new one to the Judy array
  replace_cur_node(node);
}


void PrecalcShapeset::free()
{
  if (master_pss != NULL) return;

  // iterate through the primary Judy array
  unsigned long key = 0;
  void** sub = (void**) JudyLFirst(tables, &key, NULL);
  while (sub != NULL)
  {
    // free the secondary and tertiary Judy arrays
    free_sub_tables(sub);
    JudyLDel(&tables, key, NULL);
    sub = JudyLNext(tables, &key, NULL);
  }
}


void PrecalcShapeset::dump_info(int quad, const char* filename)
{
  FILE* f = fopen(filename, "w");
  if (f == NULL) error("Could not open %s for writing.", filename);

  unsigned long key = 0, n1 = 0, m1 = 0, n2 = 0, n3 = 0, size = 0;
  void** sub = (void**) JudyLFirst(tables, &key, NULL);
  while (sub != NULL)
  {
    if ((key & 7) == (unsigned) quad)
    {
      fprintf(f, "PRIMARY TABLE, mode=%ld, index=%ld\n", (key >> 3) & 1, max_index[mode] - (key >> 4));
      unsigned long idx = 0;
      void** nodes = (void**) JudyLFirst(*sub, &idx, NULL);
      while (nodes != NULL)
      {
        fprintf(f, "   SUB TABLE, sub_idx=%ld\n      NODES: ", idx); n2++;
        unsigned long order = 0;
        void** pp = (void**) JudyLFirst(*nodes, &order, NULL);
        while (pp != NULL)
        {
          fprintf(f, "%ld ", order); n3++;
          size += ((Node*) *pp)->size;
          pp = JudyLNext(*nodes, &order, NULL);
        }
        fprintf(f, "\n");
        nodes = JudyLNext(*sub, &idx, NULL);
      }
      fprintf(f, "\n\n"); n1++;
    }
    sub = JudyLNext(tables, &key, NULL); m1++;
  }

  fprintf(f, "Number of primary tables: %ld (%ld for all quadratures)\n"
             "Avg. size of sub table:   %g\n"
             "Avg. number of nodes:     %g\n"
             "Total number of nodes:    %ld\n"
             "Total size of all nodes:  %ld bytes\n",
              n1, m1, (double) n2 / n1, (double) n3 / n2, n3, size);
  fclose(f);
}


extern PrecalcShapeset ref_map_pss; // see below

PrecalcShapeset::~PrecalcShapeset()
{
  free();
  JudyLFreeArray(&tables, NULL);

  /*if (master_pss == NULL)
  {
    verbose("~PrecalcShapeset(): peak size of precalculated tables: %d B (%0.1lf MB)%s", max_mem,
            (double) max_mem / (1024 * 1024), (this == &ref_map_pss) ? " (refmap)" : "");
  }*/
}
