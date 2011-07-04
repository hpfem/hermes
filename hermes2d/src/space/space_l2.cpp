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

#include "hermes2d_common_defs.h"
#include "space_l2.h"
#include "matrix.h"
#include "quad_all.h"
#include "shapeset/shapeset_l2_all.h"
#include "boundary_conditions/essential_boundary_conditions.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    void L2Space<Scalar>::init(Shapeset* shapeset, Ord2 p_init)
    {
      if (shapeset == NULL)
      {
        this->shapeset = new L2Shapeset;
        this->own_shapeset = true;
      }
      ldata = NULL;
      lsize = 0;

      // set uniform poly order in elements
      if (p_init.order_h < 0 || p_init.order_v < 0) error("P_INIT must be >= 0 in an L2 space.");
      else this->set_uniform_order_internal(p_init, HERMES_ANY_INT);

      // enumerate basis functions
      this->assign_dofs();
    }

    template<typename Scalar>
    L2Space<Scalar>::L2Space(Mesh* mesh, int p_init, Shapeset* shapeset)
      : Space<Scalar>(mesh, shapeset, NULL, Ord2(p_init, p_init))
    {
      _F_;
      init(shapeset, Ord2(p_init, p_init));
    }

    template<typename Scalar>
    L2Space<Scalar>::~L2Space()
    {
      ::free(ldata);
      if (this->own_shapeset)
        delete this->shapeset;
    }

    template<typename Scalar>
    Space<Scalar>* L2Space<Scalar>::dup(Mesh* mesh, int order_increase) const
    {
      L2Space<Scalar>* space = new L2Space(mesh, 0, this->shapeset);
      space->copy_orders(this, order_increase);
      return space;
    }

    template<typename Scalar>
    void L2Space<Scalar>::set_shapeset(Shapeset *shapeset)
    {
      if(shapeset->get_id() < 40 && shapeset->get_id() > 29)
      {
        this->shapeset = shapeset;
        this->own_shapeset = false;
      }
      else
        error("Wrong shapeset type in L2Space<Scalar>::set_shapeset()");
    }

    template<typename Scalar>
    void L2Space<Scalar>::resize_tables()
    {
      if (lsize < this->mesh->get_max_element_id())
      {
        if (!lsize) lsize = 1000;
        while (lsize < this->mesh->get_max_element_id()) lsize = lsize * 3 / 2;
        ldata = (L2Data*) realloc(ldata, sizeof(L2Data) * lsize);
      }
      Space<Scalar>::resize_tables();
    }

    template<typename Scalar>
    void L2Space<Scalar>::assign_bubble_dofs()
    {
      Element* e;
      for_all_active_elements(e, this->mesh)
      {
        this->shapeset->set_mode(e->get_mode());
        typename Space<Scalar>::ElementData* ed = &this->edata[e->id];
        ed->bdof = this->next_dof;
        ed->n = this->shapeset->get_num_bubbles(ed->order); //FIXME: this function might return invalid value because retrieved bubble functions for non-uniform orders might be invalid for the given order.
        this->next_dof += ed->n * this->stride;
      }
    }

    template<typename Scalar>
    void L2Space<Scalar>::get_vertex_assembly_list(Element* e, int iv, AsmList<Scalar>* al)
    {}

    template<typename Scalar>
    void L2Space<Scalar>::get_element_assembly_list(Element* e, AsmList<Scalar>* al)
    {
      // some checks
      if (e->id >= this->esize || this->edata[e->id].order < 0)
        error("Uninitialized element order (id = #%d).", e->id);
      if (!this->is_up_to_date())
        error("The space is out of date. You need to update it with assign_dofs()"
        " any time the mesh changes.");

      // add bubble functions to the assembly list
      this->shapeset->set_mode(e->get_mode());
      get_bubble_assembly_list(e, al);
    }

    template<typename Scalar>
    void L2Space<Scalar>::get_bubble_assembly_list(Element* e, AsmList<Scalar>* al)
    {
      typename Space<Scalar>::ElementData* ed = &this->edata[e->id];
      if (!ed->n) return;

      int* indices = this->shapeset->get_bubble_indices(ed->order);
      for (int i = 0, dof = ed->bdof; i < ed->n; i++, dof += this->stride)
      {
        //printf("triplet: %d, %d, %f\n", *indices, dof, 1.0);
        al->add_triplet(*indices++, dof, 1.0);
      }
    }

    template<typename Scalar>
    void L2Space<Scalar>::get_boundary_assembly_list_internal(Element* e, int surf_num, AsmList<Scalar>* al)
    {
      this->get_bubble_assembly_list(e, al);
    }

    template<typename Scalar>
    Scalar* L2Space<Scalar>::get_bc_projection(SurfPos* surf_pos, int order)
    {
      error("Method get_bc_projection() called from an L2Space.");
      return NULL;
    }

    template HERMES_API class L2Space<double>;
    template HERMES_API class L2Space<std::complex<double> >;
  }
}