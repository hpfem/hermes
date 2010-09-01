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

#include "common.h"
#include "space_l2.h"
#include "matrix_old.h"
#include "quad_all.h"
#include "shapeset_l2_all.h"

L2Space::L2Space(Mesh* mesh, int p_init, Shapeset* shapeset)
  : Space(mesh, shapeset, NULL, NULL, p_init)
{
  if (shapeset == NULL) this->shapeset = new L2Shapeset;
  ldata = NULL;
  lsize = 0;

  // set uniform poly order in elements
  if (p_init < 0) error("P_INIT must be >= 0 in an Hcurl space.");
  else this->set_uniform_order_internal(p_init);

  // enumerate basis functions
  this->assign_dofs();
}


L2Space::~L2Space()
{
  ::free(ldata);
}


Space* L2Space::dup(Mesh* mesh) const
{
  L2Space* space = new L2Space(mesh, 0, shapeset);
  space->copy_callbacks(this);
  return space;
}


//// dof assignment ////////////////////////////////////////////////////////////////////////////////

void L2Space::resize_tables()
{
  if (lsize < mesh->get_max_element_id())
  {
    if (!lsize) lsize = 1000;
    while (lsize < mesh->get_max_element_id()) lsize = lsize * 3 / 2;
    ldata = (L2Data*) realloc(ldata, sizeof(L2Data) * lsize);
  }
  Space::resize_tables();
}


void L2Space::assign_bubble_dofs()
{
  Element* e;
  for_all_active_elements(e, mesh)
  {
    shapeset->set_mode(e->get_mode());
    ElementData* ed = &edata[e->id];
    ed->bdof = next_dof;
    ed->n = shapeset->get_num_bubbles(ed->order); //FIXME: this function might return invalid value because retrieved bubble functions for non-uniform orders might be invalid for the given order.
    next_dof += ed->n * stride;
  }
}


//// assembly lists ////////////////////////////////////////////////////////////////////////////////

void L2Space::get_element_assembly_list(Element* e, AsmList* al)
{
  int i;

  // some checks
  if (e->id >= esize || edata[e->id].order < 0)
    error("Uninitialized element order (id = #%d).", e->id);
  if (!is_up_to_date())
    error("The space is out of date. You need to update it with assign_dofs()"
          " any time the mesh changes.");

  // add vertex, edge and bubble functions to the assembly list
  al->clear();
  shapeset->set_mode(e->get_mode());
  /*
  for (i = 0; i < e->nvert; i++)
    get_vertex_assembly_list(e, i, al);
  for (i = 0; i < e->nvert; i++)
    get_edge_assembly_list_internal(e, i, al);
  */
  get_bubble_assembly_list(e, al);
}

void L2Space::get_bubble_assembly_list(Element* e, AsmList* al)
{
  ElementData* ed = &edata[e->id];
  if (!ed->n) return;

  int* indices = shapeset->get_bubble_indices(ed->order);
  for (int i = 0, dof = ed->bdof; i < ed->n; i++, dof += stride) {
    //printf("triplet: %d, %d, %f\n", *indices, dof, 1.0);
    al->add_triplet(*indices++, dof, 1.0);
  }
}


void L2Space::get_edge_assembly_list_internal(Element* e, int ie, AsmList* al)
{
    this->get_bubble_assembly_list(e, al);
}

scalar* L2Space::get_bc_projection(EdgePos* ep, int order)
{
  assert(order >= 1);
  scalar* proj = new scalar[order + 1];

  // obtain linear part of the projection
  ep->t = ep->lo;
  proj[0] = bc_value_callback_by_edge(ep);
  ep->t = ep->hi;
  proj[1] = bc_value_callback_by_edge(ep);

  if (order-- > 1)
  {
    Quad1DStd quad1d;
    scalar* rhs = proj + 2;
    int mo = quad1d.get_max_order();
    double2* pt = quad1d.get_points(mo);

    // get boundary values at integration points, construct rhs
    for (int i = 0; i < order; i++)
    {
      rhs[i] = 0.0;
      int ii = shapeset->get_edge_index(0, 0, i+2);
      for (int j = 0; j < quad1d.get_num_points(mo); j++)
      {
        double t = (pt[j][0] + 1) * 0.5, s = 1.0 - t;
        scalar l = proj[0] * s + proj[1] * t;
        ep->t = ep->lo * s + ep->hi * t;
        rhs[i] += pt[j][1] * shapeset->get_fn_value(ii, pt[j][0], -1.0, 0)
                           * (bc_value_callback_by_edge(ep) - l);
      }
    }

    // solve the system using a precalculated Cholesky decomposed projection matrix
    cholsl(proj_mat, order, chol_p, rhs, rhs);
  }

  return proj;
}
