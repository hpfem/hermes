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

#include "space_h1.h"
#include "space_l2.h"
#include "shapeset_hc_all.h"
#include "shapeset_hd_all.h"
#include "shapeset_h1_all.h"
#include "shapeset_l2_all.h"
#include "space_hcurl.h"
#include "space_hdiv.h"
#include "space_h2d_xml.h"
#include "api2d.h"
#include <iostream>

namespace Hermes
{
  namespace Hermes2D
  {
    unsigned g_space_seq = 0;

		template<>
		void Space<double>::init()
		{
			this->default_tri_order = -1;
			this->default_quad_order = -1;
			this->ndata = NULL;
			this->edata = NULL;
			this->nsize = esize = 0;
			this->ndata_allocated = 0;
			this->mesh_seq = -1;
			this->seq = g_space_seq;
			this->was_assigned = -1;
			this->ndof = 0;
      this->proj_mat = NULL;
      this->chol_p = NULL;
      this->vertex_functions_count = this->edge_functions_count = this->bubble_functions_count = 0;

			if(essential_bcs != NULL)
				for(Hermes::vector<EssentialBoundaryCondition<double>*>::const_iterator it = essential_bcs->begin(); it != essential_bcs->end(); it++)
					for(unsigned int i = 0; i < (*it)->markers.size(); i++)
					{
						if(mesh->boundary_markers_conversion.conversion_table_inverse.find((*it)->markers.at(i)) == mesh->boundary_markers_conversion.conversion_table_inverse.end() && (*it)->markers.at(i) != HERMES_ANY)
							throw Hermes::Exceptions::Exception("A boundary condition defined on a non-existent marker %s.", (*it)->markers.at(i).c_str());
					}

					own_shapeset = (shapeset == NULL);
		}

		template<>
		void Space<std::complex<double> >::init()
		{
			this->default_tri_order = -1;
			this->default_quad_order = -1;
			this->ndata = NULL;
			this->edata = NULL;
			this->nsize = esize = 0;
			this->ndata_allocated = 0;
			this->mesh_seq = -1;
			this->seq = g_space_seq;
			this->was_assigned = -1;
			this->ndof = 0;
      this->proj_mat = NULL;
      this->chol_p = NULL;
      this->vertex_functions_count = this->edge_functions_count = this->bubble_functions_count = 0;

			if(essential_bcs != NULL)
				for(Hermes::vector<EssentialBoundaryCondition<std::complex<double> >*>::const_iterator it = essential_bcs->begin(); it != essential_bcs->end(); it++)
					for(unsigned int i = 0; i < (*it)->markers.size(); i++)
					{
						if(mesh->boundary_markers_conversion.conversion_table_inverse.find((*it)->markers.at(i)) == mesh->boundary_markers_conversion.conversion_table_inverse.end() && (*it)->markers.at(i) != HERMES_ANY)
							throw Hermes::Exceptions::Exception("A boundary condition defined on a non-existent marker %s.", (*it)->markers.at(i).c_str());
					}

					own_shapeset = (shapeset == NULL);
		}

		template<>
		void Space<double>::free()
		{
			free_bc_data();
			if(nsize) { ::free(ndata); nsize = 0; ndata = NULL; }
			if(esize) { ::free(edata); edata = 0; edata = NULL; }
			this->seq = -1;
		}

		template<>
		void Space<std::complex<double> >::free()
		{
			free_bc_data();
			if(nsize) { ::free(ndata); nsize = 0; ndata = NULL; }
			if(esize) { ::free(edata); edata = 0; edata = NULL; }
			this->seq = -1;
		}

    template<>
    Space<double>::Space() : shapeset(NULL), essential_bcs(NULL), mesh(NULL)
    {
      this->init();
    }

    template<>
    Space<std::complex<double> >::Space() : shapeset(NULL), essential_bcs(NULL), mesh(NULL)
    {
      this->init();
    }

     template<>
    Space<double>::Space(const Mesh* mesh, Shapeset* shapeset, EssentialBCs<double>* essential_bcs)
      : shapeset(shapeset), essential_bcs(essential_bcs), mesh(mesh)
    {
      if(mesh == NULL)
        throw Hermes::Exceptions::NullException(0);
      this->init();
    }

    template<>
    Space<std::complex<double> >::Space(const Mesh* mesh, Shapeset* shapeset, EssentialBCs<std::complex<double> >* essential_bcs)
      : shapeset(shapeset), essential_bcs(essential_bcs), mesh(mesh)
    {
      if(mesh == NULL)
        throw Hermes::Exceptions::NullException(0);
      this->init();
    }
    
    template<typename Scalar>
    bool Space<Scalar>::isOkay() const
    {
      if(ndata == NULL || edata == NULL || !nsize || !esize)
        return false;
      if(seq < 0)
        return false;
      if(this->mesh == NULL)
        return false;	

      this->mesh->check();

      if(edata == NULL)
      {
          throw Hermes::Exceptions::Exception("NULL edata detected in Space<Scalar>::get_element_order().");
          return false;
      }

      if(!this->is_up_to_date())
      {
        throw Hermes::Exceptions::Exception("Space is not up to date.");
        return false;
      }
      else
        return true;
    }

    template<>
    Space<double>::~Space()
    {
      free();

      if(this->proj_mat != NULL)
        delete [] this->proj_mat;
      if(this->chol_p != NULL)
        delete [] this->chol_p;

    }

    template<>
    Space<std::complex<double> >::~Space()
    {
      free();

      if(this->proj_mat != NULL)
        delete [] this->proj_mat;
      if(this->chol_p != NULL)
        delete [] this->chol_p;

    }

		template<typename Scalar>
    Node* Space<Scalar>::get_mid_edge_vertex_node(Element* e, int i, int j)
    {
      if(e->is_triangle())
        return e->sons[3]->vn[e->prev_vert(i)];
      else if(e->sons[2] == NULL)
        return i == 1 ? e->sons[0]->vn[2] : i == 3 ? e->sons[0]->vn[3] : NULL;
      else if(e->sons[0] == NULL)
        return i == 0 ? e->sons[2]->vn[1] : i == 2 ? e->sons[2]->vn[2] : NULL;
      else return e->sons[i]->vn[j];
    }

    template<typename Scalar>
    void Space<Scalar>::resize_tables()
    {
      if((nsize < mesh->get_max_node_id()) || (ndata == NULL))
      {
        //HACK: definition of allocated size and the result number of elements
        nsize = mesh->get_max_node_id();
        if((nsize > ndata_allocated) || (ndata == NULL))
        {
          int prev_allocated = ndata_allocated;
          if(ndata_allocated == 0)
            ndata_allocated = 1024;
          while (ndata_allocated < nsize)
            ndata_allocated = ndata_allocated * 3 / 2;
          ndata = (NodeData*)realloc(ndata, ndata_allocated * sizeof(NodeData));
          for(int i = prev_allocated; i < ndata_allocated; i++)
            ndata[i].edge_bc_proj = NULL;
        }
      }

      if((esize < mesh->get_max_element_id()) || (edata == NULL))
      {
        int oldsize = esize;
        if(!esize)
          esize = 1024;
        while (esize < mesh->get_max_element_id())
          esize = esize * 3 / 2;
        edata = (ElementData*) realloc(edata, sizeof(ElementData) * esize);
        for (int i = oldsize; i < esize; i++)
          edata[i].order = -1;
        for (int i = oldsize; i < esize; i++)
          edata[i].changed_in_last_adaptation = true;
      }
    }

    template<typename Scalar>
    void Space<Scalar>::copy(const Space<Scalar>* space, Mesh* new_mesh)
    {
      this->free();
      
      this->vertex_functions_count = this->edge_functions_count = this->bubble_functions_count = 0;

      this->essential_bcs = space->essential_bcs;
      this->shapeset = space->shapeset->clone();

      new_mesh->copy(space->get_mesh());
      this->mesh = new_mesh;

      this->resize_tables();

      Element* e;
      for_all_active_elements(e, this->mesh)
      {
        this->set_element_order_internal(e->id, space->get_element_order(e->id));
      }

      this->seq = g_space_seq++;
      
      for_all_active_elements(e, this->mesh)
      {
        if(space->edata[e->id].changed_in_last_adaptation)
          this->edata[e->id].changed_in_last_adaptation = true;
        else
          this->edata[e->id].changed_in_last_adaptation = false;
      }
    }

    template<typename Scalar>
    Shapeset* Space<Scalar>::get_shapeset() const
    {
      return this->shapeset;
    }

    template<typename Scalar>
    int Space<Scalar>::get_num_dofs() const
    {
      check();
      return ndof;
    }

    template<typename Scalar>
    int Space<Scalar>::get_max_dof() const
    {
      check();
      return next_dof - stride;
    }

    template<typename Scalar>
    Mesh* Space<Scalar>::get_mesh() const
    {
      return const_cast<Mesh*>(mesh);
    }

    template<typename Scalar>
    bool Space<Scalar>::is_up_to_date() const
    {
      return was_assigned == this->seq && mesh_seq == (int) mesh->get_seq();
    }

    template<typename Scalar>
    EssentialBCs<Scalar>* Space<Scalar>::get_essential_bcs() const
    {
      check();
      return essential_bcs;
    }

    template<typename Scalar>
    void Space<Scalar>::set_element_order(int id, int order)
    {
      set_element_order_internal(id, order);
    }

    template<typename Scalar>
    void Space<Scalar>::set_element_order_internal(int id, int order)
    {
      if(id < 0 || id >= mesh->get_max_element_id())
        throw Hermes::Exceptions::Exception("Space<Scalar>::set_element_order_internal: Invalid element id.");

      resize_tables();

      if(mesh->get_element(id)->is_quad() && get_type() != HERMES_L2_SPACE && H2D_GET_V_ORDER(order) == 0)
        order = H2D_MAKE_QUAD_ORDER(order, order);

      edata[id].order = order;
      seq = g_space_seq++;
    }

    template<typename Scalar>
    void Space<Scalar>::update_essential_bc_values(Hermes::vector<Space<Scalar>*> spaces, double time)
    {
      int n = spaces.size();
      for (int i = 0; i < n; i++)
      {
        if(spaces[i]->get_essential_bcs() != NULL)
          spaces[i]->get_essential_bcs()->set_current_time(time);
        spaces[i]->update_essential_bc_values();
      }
    }

    template<typename Scalar>
    void Space<Scalar>::update_essential_bc_values(Space<Scalar>*s, double time)
    {
      s->get_essential_bcs()->set_current_time(time);
      s->update_essential_bc_values();
    }

    template<typename Scalar>
    int Space<Scalar>::get_num_dofs(Hermes::vector<const Space<Scalar>*> spaces)
    {
      int ndof = 0;
      for (unsigned int i = 0; i<spaces.size(); i++)
        ndof += spaces[i]->get_num_dofs();
      return ndof;
    }

    template<typename Scalar>
    int Space<Scalar>::get_num_dofs(Hermes::vector<Space<Scalar>*> spaces)
    {
      int ndof = 0;
      for (unsigned int i = 0; i<spaces.size(); i++)
        ndof += spaces[i]->get_num_dofs();
      return ndof;
    }

    template<typename Scalar>
    int Space<Scalar>::get_num_dofs(const Space<Scalar>* space)
    {
      return space->get_num_dofs();
    }

    template<typename Scalar>
    int Space<Scalar>::get_num_dofs(Space<Scalar>* space)
    {
      return space->get_num_dofs();
    }

    template<typename Scalar>
    int Space<Scalar>::assign_dofs(Hermes::vector<Space<Scalar>*> spaces)
    {
      int n = spaces.size();

      int ndof = 0;
      for (int i = 0; i < n; i++) {
        ndof += spaces[i]->assign_dofs(ndof);
      }

      return ndof;
    }

    template<typename Scalar>
    int Space<Scalar>::get_element_order(int id) const
    {
      if(id >= esize)
      {
        this->warn("Element index %d in Space<Scalar>::get_element_order() while maximum is %d.", id, esize);
        throw Hermes::Exceptions::Exception("Wrong element index in Space<Scalar>::get_element_order().");
      }

      return edata[id].order;
    }

    template<typename Scalar>
    void Space<Scalar>::set_uniform_order(int order, std::string marker)
    {
      if(marker == HERMES_ANY)
        set_uniform_order_internal(order, -1234);
      else
        set_uniform_order_internal(order, mesh->element_markers_conversion.get_internal_marker(marker).marker);
    }

    template<typename Scalar>
    void Space<Scalar>::set_uniform_order_internal(int order, int marker)
    {
      resize_tables();
      int quad_order = H2D_MAKE_QUAD_ORDER(order, order);

      Element* e;
      for_all_active_elements(e, mesh)
      {
        if(marker == HERMES_ANY_INT || e->marker == marker)
        {
          ElementData* ed = &edata[e->id];
          if(e->is_triangle())
            ed->order = order;
          else
            ed->order = quad_order;
        }
      }
      seq = g_space_seq++;
    }

    template<typename Scalar>
    void Space<Scalar>::set_element_orders(int* elem_orders_)
    {
      resize_tables();

      Element* e;
      int counter = 0;
      for_all_elements(e, mesh)
      {
        assert(elem_orders_[counter] >= 0 && elem_orders_[counter] <= shapeset->get_max_order());
        ElementData* ed = &edata[e->id];
        if(e->is_triangle())
          ed->order = elem_orders_[counter];
        else
          ed->order = H2D_MAKE_QUAD_ORDER(elem_orders_[counter], elem_orders_[counter]);
        counter++;
      }
    }

    template<typename Scalar>
    void Space<Scalar>::adjust_element_order(int order_change, int min_order)
    {
      Element* e;
      for_all_active_elements(e, this->get_mesh())
      {
        if(e->is_triangle())
          set_element_order_internal(e->id, std::max<int>(min_order, get_element_order(e->id) + order_change));
        else
        {
          int h_order, v_order;

          if(H2D_GET_H_ORDER(get_element_order(e->id)) + order_change < min_order)
            h_order = min_order;
          else
            h_order = H2D_GET_H_ORDER(get_element_order(e->id)) + order_change;

          if(H2D_GET_V_ORDER(get_element_order(e->id)) + order_change < min_order)
            v_order = min_order;
          else
            v_order = H2D_GET_V_ORDER(get_element_order(e->id)) + order_change;

          set_element_order_internal(e->id, H2D_MAKE_QUAD_ORDER(h_order, v_order));
        }
      }
    }

    template<typename Scalar>
    void Space<Scalar>::adjust_element_order(int horizontal_order_change, int vertical_order_change, unsigned int horizontal_min_order, unsigned int vertical_min_order)
    {
      Element* e;
      for_all_active_elements(e, this->get_mesh())
      {
        if(e->is_triangle())
        {
          this->warn("Using quad version of Space<Scalar>::adjust_element_order(), only horizontal orders will be used.");
          set_element_order_internal(e->id, std::max<int>(horizontal_min_order, get_element_order(e->id) + horizontal_order_change));
        }
        else
          if(get_element_order(e->id) == -1)
            set_element_order_internal(e->id, H2D_MAKE_QUAD_ORDER(horizontal_min_order, vertical_min_order));
          else
            set_element_order_internal(e->id, H2D_MAKE_QUAD_ORDER(std::max<int>(H2D_GET_H_ORDER(get_element_order(e->id)) + horizontal_order_change, horizontal_min_order), std::max<int>(H2D_GET_V_ORDER(get_element_order(e->id)) + vertical_order_change, vertical_min_order)));
      }
    }

    template<typename Scalar>
    void Space<Scalar>::update_element_orders_after_refinement()
    {
      Element* e;
      for_all_active_elements(e, this->mesh)
      {
        if(this->get_element_order(e->id) < 0)
          this->set_element_order_internal(e->id, this->get_element_order(e->parent->id));
      }
    }

    template<typename Scalar>
    void Space<Scalar>::unrefine_all_mesh_elements(bool keep_initial_refinements)
    {
      check();
      // find inactive elements with active sons
      Hermes::vector<int> list;
      Element* e;
      for_all_inactive_elements(e, this->mesh)
      {
        bool found = true;
        for (unsigned int i = 0; i < 4; i++)
          if(e->sons[i] != NULL &&
            (!e->sons[i]->active || (keep_initial_refinements && e->sons[i]->id < this->mesh->ninitial))
            )
          { found = false; break; }

          if(found) list.push_back(e->id);
      }

      // unrefine the found elements
      for (unsigned int i = 0; i < list.size(); i++)
      {
        unsigned int order = 0, h_order = 0, v_order = 0;
        unsigned int num_sons = 0;
        if(this->mesh->get_element_fast(list[i])->bsplit())
        {
          num_sons = 4;
          for (int sons_i = 0; sons_i < 4; sons_i++)
          {
            if(this->mesh->get_element_fast(list[i])->sons[sons_i]->active)
            {
              if(this->mesh->get_element_fast(list[i])->sons[sons_i]->is_triangle())
                order += this->get_element_order(this->mesh->get_element_fast(list[i])->sons[sons_i]->id);
              else
              {
                h_order += H2D_GET_H_ORDER(this->get_element_order(this->mesh->get_element_fast(list[i])->sons[sons_i]->id));
                v_order += H2D_GET_V_ORDER(this->get_element_order(this->mesh->get_element_fast(list[i])->sons[sons_i]->id));
              }
            }
          }
        }
        else
        {
          if(this->mesh->get_element_fast(list[i])->hsplit())
          {
            num_sons = 2;
            if(this->mesh->get_element_fast(list[i])->sons[0]->active)
            {
              if(this->mesh->get_element_fast(list[i])->sons[0]->is_triangle())
                order += this->get_element_order(this->mesh->get_element_fast(list[i])->sons[0]->id);
              else
              {
                h_order += H2D_GET_H_ORDER(this->get_element_order(this->mesh->get_element_fast(list[i])->sons[0]->id));
                v_order += H2D_GET_V_ORDER(this->get_element_order(this->mesh->get_element_fast(list[i])->sons[0]->id));
              }
            }
            if(this->mesh->get_element_fast(list[i])->sons[1]->active)
            {
              if(this->mesh->get_element_fast(list[i])->sons[1]->is_triangle())
                order += this->get_element_order(this->mesh->get_element_fast(list[i])->sons[1]->id);
              else
              {
                h_order += H2D_GET_H_ORDER(this->get_element_order(this->mesh->get_element_fast(list[i])->sons[1]->id));
                v_order += H2D_GET_V_ORDER(this->get_element_order(this->mesh->get_element_fast(list[i])->sons[1]->id));
              }
            }
          }
          else
          {
            num_sons = 2;
            if(this->mesh->get_element_fast(list[i])->sons[2]->active)
            {
              if(this->mesh->get_element_fast(list[i])->sons[2]->is_triangle())
                order += this->get_element_order(this->mesh->get_element_fast(list[i])->sons[2]->id);
              else
              {
                h_order += H2D_GET_H_ORDER(this->get_element_order(this->mesh->get_element_fast(list[i])->sons[2]->id));
                v_order += H2D_GET_V_ORDER(this->get_element_order(this->mesh->get_element_fast(list[i])->sons[2]->id));
              }
            }
            if(this->mesh->get_element_fast(list[i])->sons[3]->active)
            {
              if(this->mesh->get_element_fast(list[i])->sons[3]->is_triangle())
                order += this->get_element_order(this->mesh->get_element_fast(list[i])->sons[3]->id);
              else
              {
                h_order += H2D_GET_H_ORDER(this->get_element_order(this->mesh->get_element_fast(list[i])->sons[3]->id));
                v_order += H2D_GET_V_ORDER(this->get_element_order(this->mesh->get_element_fast(list[i])->sons[3]->id));
              }
            }
          }
        }
        order = (unsigned int)(order / num_sons);
        h_order = (unsigned int)(h_order / num_sons);
        v_order = (unsigned int)(v_order / num_sons);

        if(this->mesh->get_element_fast(list[i])->is_triangle())
          edata[list[i]].order = order;
        else
          edata[list[i]].order = H2D_MAKE_QUAD_ORDER(h_order, v_order);
        const_cast<Mesh*>(this->mesh)->unrefine_element_id(list[i]);
      }

      // Recalculate all integrals, do not use previous adaptivity step.
      for_all_active_elements(e, this->mesh)
        this->edata[e->id].changed_in_last_adaptation = true;

    }

    template<typename Scalar>
    Space<Scalar>::ReferenceSpaceCreator::ReferenceSpaceCreator(const Space<Scalar>* coarse_space, const Mesh* ref_mesh, unsigned int order_increase) : coarse_space(coarse_space), ref_mesh(ref_mesh), order_increase(order_increase)
    {
    }

    template<typename Scalar>
    void Space<Scalar>::ReferenceSpaceCreator::handle_orders(Space<Scalar>* ref_space)
    {
      Element* e;
      for_all_active_elements(e, coarse_space->get_mesh())
      {
        // This is the element id on the COARSE mesh. One can use it in the logic of increasing polynomial order (selectively).
        int coarse_element_id = e->id;

        // Get the current order (may be triangular or quadrilateral - in the latter case, it is both horizontal and vertical order encoded in one number).
        int current_order = coarse_space->get_element_order(coarse_element_id);
        // Sanity check.
        if(current_order < 0)
          throw Hermes::Exceptions::Exception("Source space has an uninitialized order (element id = %d)", coarse_element_id);

        // New order calculation.
        int new_order;
        if(e->is_triangle())
        {
          // The triangular order is just the current order.
          int current_order_triangular = current_order;
          /// This is the default setup, we ALWAYS increase by attribute of this functor class: 'order_increase'.
          /// This is OVERRIDABLE. Plus when overriding, one does not have to care about min / max possible values (due to shapeset, Space requirements).
          /// These are guarded internally.
          new_order = current_order_triangular + this->order_increase;
        }
        else
        {
          // We have to get the proper parts of the encoded order.
          // - horizontal.
          int current_order_quadrilateral_horizontal = H2D_GET_H_ORDER(current_order);
          // - vertical.
          int current_order_quadrilateral_vertical = H2D_GET_V_ORDER(current_order);

          // And now we have to create the new encoded order.
          int new_order_quadrilateral_horizontal = current_order_quadrilateral_horizontal + this->order_increase;
          int new_order_quadrilateral_vertical = current_order_quadrilateral_vertical + this->order_increase;
          new_order = H2D_MAKE_QUAD_ORDER(new_order_quadrilateral_horizontal, new_order_quadrilateral_vertical);
        }

        // And now call this method that does the magic for us (sets the new_order to .
        ref_space->update_orders_recurrent(ref_space->mesh->get_element(coarse_element_id), new_order);
      }
    }
    
    template<typename Scalar>
    Space<Scalar>* Space<Scalar>::ReferenceSpaceCreator::create_ref_space(bool assign_dofs)
    {
      /// Initialization.
      Space<Scalar>* ref_space = NULL;
      if(dynamic_cast<const L2Space<Scalar>*>(this->coarse_space) != NULL)
        ref_space = this->init_construction_l2();
      if(dynamic_cast<const H1Space<Scalar>*>(this->coarse_space) != NULL)
        ref_space = this->init_construction_h1();
      if(dynamic_cast<const HcurlSpace<Scalar>*>(this->coarse_space) != NULL)
        ref_space = this->init_construction_hcurl();
      if(dynamic_cast<const HdivSpace<Scalar>*>(this->coarse_space) != NULL)
        ref_space = this->init_construction_hdiv();

      if(ref_space == NULL)
        throw Exceptions::Exception("Something went wrong in ReferenceSpaceCreator::create_ref_space().");

      /// Call to the OVERRIDABLE handling method.
      this->handle_orders(ref_space);

      /// Finish - MUST BE CALLED BEFORE RETURN.
      this->finish_construction(ref_space);

      // Assign dofs?
      if(assign_dofs)
        ref_space->assign_dofs();

      /// Return.
      return ref_space;
    }

    template<typename Scalar>
    L2Space<Scalar>* Space<Scalar>::ReferenceSpaceCreator::init_construction_l2()
    {
      L2Space<Scalar>* ref_space;
      if(this->coarse_space->own_shapeset)
        ref_space = new L2Space<Scalar>(this->ref_mesh, 0);
      else
        ref_space = new L2Space<Scalar>(this->ref_mesh, 0, this->coarse_space->get_shapeset());
      return ref_space;
    }

    template<typename Scalar>
    H1Space<Scalar>* Space<Scalar>::ReferenceSpaceCreator::init_construction_h1()
    {
      H1Space<Scalar>* ref_space;
      if(this->coarse_space->own_shapeset)
        ref_space = new H1Space<Scalar>(this->ref_mesh, this->coarse_space->get_essential_bcs(), 1);
      else
        ref_space = new H1Space<Scalar>(this->ref_mesh, this->coarse_space->get_essential_bcs(), 1, this->coarse_space->get_shapeset());
      return ref_space;
    }

    template<typename Scalar>
    HcurlSpace<Scalar>* Space<Scalar>::ReferenceSpaceCreator::init_construction_hcurl()
    {
      HcurlSpace<Scalar>* ref_space;
      if(this->coarse_space->own_shapeset)
        ref_space = new HcurlSpace<Scalar>(this->ref_mesh, this->coarse_space->essential_bcs, 1);
      else
        ref_space = new HcurlSpace<Scalar>(this->ref_mesh, this->coarse_space->essential_bcs, 1, this->coarse_space->shapeset);
      return ref_space;
    }

    template<typename Scalar>
    HdivSpace<Scalar>* Space<Scalar>::ReferenceSpaceCreator::init_construction_hdiv()
    {
      HdivSpace<Scalar>* ref_space;
      if(this->coarse_space->own_shapeset)
        ref_space = new HdivSpace<Scalar>(this->ref_mesh, this->coarse_space->essential_bcs, 1);
      else
        ref_space = new HdivSpace<Scalar>(this->ref_mesh, this->coarse_space->essential_bcs, 1, this->coarse_space->shapeset);
      return ref_space;
    }

    template<typename Scalar>
    void Space<Scalar>::ReferenceSpaceCreator::finish_construction(Space<Scalar>* ref_space)
    {
      ref_space->seq = g_space_seq++;

      Element *e;
      for_all_active_elements(e, coarse_space->get_mesh())
      {
        if(this->coarse_space->edata[e->id].changed_in_last_adaptation)
        {
          if(ref_space->mesh->get_element(e->id)->active)
            ref_space->edata[e->id].changed_in_last_adaptation = true;
          else
            for(unsigned int i = 0; i < 4; i++)
              if(ref_space->mesh->get_element(e->id)->sons[i] != NULL)
                if(ref_space->mesh->get_element(e->id)->sons[i]->active)
                  ref_space->edata[ref_space->mesh->get_element(e->id)->sons[i]->id].changed_in_last_adaptation = true;
        }
      }
    }

    template<typename Scalar>
    void Space<Scalar>::update_orders_recurrent(Element* e, int order)
    {
      // Adjust wrt. max and min possible orders.
      int mo = shapeset->get_max_order();
      int lower_limit = (get_type() == HERMES_L2_SPACE || get_type() == HERMES_HCURL_SPACE) ? 0 : 1; // L2 and Hcurl may use zero orders.
      int ho = std::max(lower_limit, std::min(H2D_GET_H_ORDER(order), mo));
      int vo = std::max(lower_limit, std::min(H2D_GET_V_ORDER(order), mo));
      order = e->is_triangle() ? ho : H2D_MAKE_QUAD_ORDER(ho, vo);

      if(e->active)
        edata[e->id].order = order;
      else
        for (int i = 0; i < 4; i++)
          if(e->sons[i] != NULL)
            update_orders_recurrent(e->sons[i], order);
    }

    template<typename Scalar>
    int Space<Scalar>::get_edge_order(Element* e, int edge) const
    {
      Node* en = e->en[edge];
      if(en->id >= nsize || edge >= (int)e->get_nvert()) return 0;

      if(ndata[en->id].n == -1)
        return get_edge_order_internal(ndata[en->id].base); // constrained node
      else
        return get_edge_order_internal(en);
    }

    template<typename Scalar>
    int Space<Scalar>::get_edge_order_internal(Node* en) const
    {
      assert(en->type == HERMES_TYPE_EDGE);
      Element** e = en->elem;
      int o1 = 1000, o2 = 1000;
      assert(e[0] != NULL || e[1] != NULL);

      if(e[0] != NULL)
      {
        if(e[0]->is_triangle() || en == e[0]->en[0] || en == e[0]->en[2])
          o1 = H2D_GET_H_ORDER(edata[e[0]->id].order);
        else
          o1 = H2D_GET_V_ORDER(edata[e[0]->id].order);
      }

      if(e[1] != NULL)
      {
        if(e[1]->is_triangle() || en == e[1]->en[0] || en == e[1]->en[2])
          o2 = H2D_GET_H_ORDER(edata[e[1]->id].order);
        else
          o2 = H2D_GET_V_ORDER(edata[e[1]->id].order);
      }

      if(o1 == 0) return o2 == 1000 ? 0 : o2;
      if(o2 == 0) return o1 == 1000 ? 0 : o1;
      return std::min(o1, o2);
    }

    template<typename Scalar>
    void Space<Scalar>::set_mesh(Mesh* mesh)
    {
      if(this->mesh == mesh) 
        return;
      free();
      this->mesh = mesh;
      this->mesh_seq = mesh->get_seq();
      seq = g_space_seq++;
    }

    template<typename Scalar>
    void Space<Scalar>::set_mesh_seq(int seq)
    {
      this->mesh_seq = seq;
      const_cast<Mesh*>(this->mesh)->set_seq(seq);
    }

    template<typename Scalar>
    void Space<Scalar>::update_constraints()
    {
    }

    template<typename Scalar>
    void Space<Scalar>::post_assign()
    {
    }

    template<typename Scalar>
    int Space<Scalar>::get_seq() const
    {
      return seq;
    }

    template<typename Scalar>
    void Space<Scalar>::distribute_orders(Mesh* mesh, int* parents)
    {
      int num = mesh->get_max_element_id();
      int* orders = new int[num + 1];
      Element* e;
      for_all_active_elements(e, mesh)
      {
        int p = get_element_order(parents[e->id]);
        if(e->is_triangle() && (H2D_GET_V_ORDER(p) != 0))
          p = std::max(H2D_GET_H_ORDER(p), H2D_GET_V_ORDER(p));
        orders[e->id] = p;
      }
      for_all_active_elements(e, mesh)
        set_element_order_internal(e->id, orders[e->id]);
      delete [] orders;
    }

    template<typename Scalar>
    int Space<Scalar>::assign_dofs(int first_dof, int stride)
    {
      if(ndata == NULL || edata == NULL || !nsize || !esize)
        return false;
      if(seq < 0)
        return false;
      if(this->mesh == NULL)
        return false;

      this->mesh->check();

      if(edata == NULL)
      {
          throw Hermes::Exceptions::Exception("NULL edata detected in Space<Scalar>::get_element_order().");
          return false;
      }

      if(first_dof < 0)
        throw Hermes::Exceptions::ValueException("first_dof", first_dof, 0);
      if(stride < 1)
        throw Hermes::Exceptions::ValueException("stride", stride, 1);

      resize_tables();

      Element* e;

      //check validity of orders
      for_all_active_elements(e, mesh)
      {
        if(e->id >= esize || edata[e->id].order < 0)
        {
          printf("e->id = %d\n", e->id);
          printf("esize = %d\n", esize);
          printf("edata[%d].order = %d\n", e->id, edata[e->id].order);
          throw
            Hermes::Exceptions::Exception("Uninitialized element order in Space::assign_dofs().");
        }
        this->edata[e->id].changed_in_last_adaptation = true;
      }

      this->first_dof = next_dof = first_dof;
      this->stride = stride;

      reset_dof_assignment();
      assign_vertex_dofs();
      assign_edge_dofs();
      assign_bubble_dofs();

      free_bc_data();
      update_essential_bc_values();
      update_constraints();
      post_assign();

      mesh_seq = mesh->get_seq();
      was_assigned = this->seq;
      this->ndof = (next_dof - first_dof) / stride;

      return this->ndof;
      check();
    }

    template<typename Scalar>
    void Space<Scalar>::reset_dof_assignment()
    {
      // First assume that all vertex nodes are part of a natural BC. the member NodeData::n
      // is misused for this purpose, since it stores nothing at this point. Also assume
      // that all DOFs are unassigned.
      int i, j;
      for (i = 0; i < mesh->get_max_node_id(); i++)
      {
        ndata[i].n = 1; // Natural boundary condition. The point is that it is not (0 == Dirichlet).
        ndata[i].dof = H2D_UNASSIGNED_DOF;
      }

      // next go through all boundary edge nodes constituting an essential BC and mark their
      // neighboring vertex nodes also as essential
      Element* e;
      for_all_active_elements(e, mesh)
      {
        for (unsigned int i = 0; i < e->get_nvert(); i++)
        {
          if(e->en[i]->bnd)
            if(essential_bcs != NULL)
              if(essential_bcs->get_boundary_condition(mesh->boundary_markers_conversion.get_user_marker(e->en[i]->marker).marker) != NULL)
              {
                j = e->next_vert(i);
                ndata[e->vn[i]->id].n = 0;
                ndata[e->vn[j]->id].n = 0;
              }
        }
      }
    }

    template<typename Scalar>
    int Space<Scalar>::get_vertex_functions_count()
    {
      return this->vertex_functions_count;
    }

    template<typename Scalar>
    int Space<Scalar>::get_edge_functions_count()
    {
      return this->edge_functions_count;
    }

    template<typename Scalar>
    int Space<Scalar>::get_bubble_functions_count()
    {
      return this->bubble_functions_count;
    }

    template<typename Scalar>
    void Space<Scalar>::get_element_assembly_list(Element* e, AsmList<Scalar>* al, unsigned int first_dof) const
    {
      this->check();
      // some checks
      if(e->id >= esize || edata[e->id].order < 0)
        throw Hermes::Exceptions::Exception("Uninitialized element order in get_element_assembly_list(id = #%d).", e->id);
      if(!is_up_to_date())
        throw Hermes::Exceptions::Exception("The space in get_element_assembly_list() is out of date. You need to update it with assign_dofs()"
        " any time the mesh changes.");

      // add vertex, edge and bubble functions to the assembly list
      al->cnt = 0;
      for (unsigned int i = 0; i < e->get_nvert(); i++)
        get_vertex_assembly_list(e, i, al);
      for (unsigned int i = 0; i < e->get_nvert(); i++)
        get_boundary_assembly_list_internal(e, i, al);
      get_bubble_assembly_list(e, al);
      for(unsigned int i = 0; i < al->cnt; i++)
        if(al->dof[i] >= 0)
          al->dof[i] += first_dof;
    }

    template<typename Scalar>
    void Space<Scalar>::get_boundary_assembly_list(Element* e, int surf_num, AsmList<Scalar>* al, unsigned int first_dof) const
    {
      this->check();
      al->cnt = 0;
      get_vertex_assembly_list(e, surf_num, al);
      get_vertex_assembly_list(e, e->next_vert(surf_num), al);
      get_boundary_assembly_list_internal(e, surf_num, al);
      for(unsigned int i = 0; i < al->cnt; i++)
        if(al->dof[i] >= 0)
          al->dof[i] += first_dof;
    }

    template<typename Scalar>
    void Space<Scalar>::get_bubble_assembly_list(Element* e, AsmList<Scalar>* al) const
    {
      this->check();
      ElementData* ed = &edata[e->id];

      if(!ed->n) return;

      int* indices = shapeset->get_bubble_indices(ed->order, e->get_mode());
      for (int i = 0, dof = ed->bdof; i < ed->n; i++, dof += stride, indices++)
        al->add_triplet(*indices, dof, 1.0);
    }

    template<typename Scalar>
    void Space<Scalar>::set_essential_bcs(EssentialBCs<Scalar>* essential_bcs)
    {
      this->essential_bcs = essential_bcs;
    }

    template<typename Scalar>
    void Space<Scalar>::precalculate_projection_matrix(int nv, double**& mat, double*& p)
    {
      int n = shapeset->get_max_order() + 1 - nv;
      mat = new_matrix<double>(n, n);
      int component = (get_type() == HERMES_HDIV_SPACE) ? 1 : 0;

      Quad1DStd quad1d;
      for (int i = 0; i < n; i++)
      {
        for (int j = i; j < n; j++)
        {
          int o = i + j + 4;
          double2* pt = quad1d.get_points(o);
          int ii = shapeset->get_edge_index(0, 0, i + nv, HERMES_MODE_QUAD);
          int ij = shapeset->get_edge_index(0, 0, j + nv, HERMES_MODE_QUAD);
          double val = 0.0;
          for (int k = 0; k < quad1d.get_num_points(o); k++)
          {
            val += pt[k][1] * shapeset->get_fn_value(ii, pt[k][0], -1.0, component, HERMES_MODE_QUAD)
              * shapeset->get_fn_value(ij, pt[k][0], -1.0, component, HERMES_MODE_QUAD);
          }
          mat[i][j] = val;
        }
      }

      p = new double[n];
      choldc(mat, n, p);
    }

    template<typename Scalar>
    void Space<Scalar>::update_edge_bc(Element* e, SurfPos* surf_pos)
    {
      if(e->active)
      {
        Node* en = e->en[surf_pos->surf_num];
        NodeData* nd = &ndata[en->id];
        nd->edge_bc_proj = NULL;

        if(nd->dof != H2D_UNASSIGNED_DOF && en->bnd)
          if(essential_bcs != NULL)
          {
            EssentialBoundaryCondition<Scalar> *bc = this->essential_bcs->get_boundary_condition(this->mesh->boundary_markers_conversion.get_user_marker(en->marker).marker);
            if(bc != NULL)
            {
              int order = get_edge_order_internal(en);
              surf_pos->marker = en->marker;
              nd->edge_bc_proj = get_bc_projection(surf_pos, order, bc);
              bc_data.push_back(nd->edge_bc_proj);

              int i = surf_pos->surf_num, j = e->next_vert(i);
              ndata[e->vn[i]->id].vertex_bc_coef = nd->edge_bc_proj + 0;
              ndata[e->vn[j]->id].vertex_bc_coef = nd->edge_bc_proj + 1;
            }
          }
      }
      else
      {
        int son1, son2;
        if(mesh->get_edge_sons(e, surf_pos->surf_num, son1, son2) == 2)
        {
          double mid = (surf_pos->lo + surf_pos->hi) * 0.5, tmp = surf_pos->hi;
          surf_pos->hi = mid;
          update_edge_bc(e->sons[son1], surf_pos);
          surf_pos->lo = mid; surf_pos->hi = tmp;
          update_edge_bc(e->sons[son2], surf_pos);
        }
        else
          update_edge_bc(e->sons[son1], surf_pos);
      }
    }

    template<typename Scalar>
    void Space<Scalar>::update_essential_bc_values()
    {
      Element* e;
      for_all_base_elements(e, mesh)
      {
        for (unsigned int i = 0; i < e->get_nvert(); i++)
        {
          int j = e->next_vert(i);
          if(e->vn[i]->bnd && e->vn[j]->bnd)
          {
            SurfPos surf_pos = {0, i, e, e->vn[i]->id, e->vn[j]->id, 0.0, 0.0, 1.0};
            update_edge_bc(e, &surf_pos);
          }
        }
      }
    }

    template<typename Scalar>
    void Space<Scalar>::free_bc_data()
    {
      for (unsigned int i = 0; i < bc_data.size(); i++)
        delete (bc_data[i]);
      bc_data.clear();
    }

    template<typename Scalar>
    bool Space<Scalar>::save(const char *filename) const
    {
      this->check();
      XMLSpace::space xmlspace;

      switch(this->get_type())
      {
        case HERMES_H1_SPACE:
            xmlspace.spaceType().set("h1");
            break;
        case HERMES_HCURL_SPACE:
            xmlspace.spaceType().set("hcurl");
            break;
        case HERMES_HDIV_SPACE:
            xmlspace.spaceType().set("hdiv");
            break;
        case HERMES_L2_SPACE:
            xmlspace.spaceType().set("l2");
            break;
        default:
            return false;
      }

      // Utility pointer.
      Element *e;
      for_all_elements(e, this->get_mesh())
        xmlspace.element_data().push_back(XMLSpace::space::element_data_type(e->id, this->edata[e->id].order, this->edata[e->id].bdof, this->edata[e->id].n, this->edata[e->id].changed_in_last_adaptation));

      std::string space_schema_location(Hermes2DApi.get_text_param_value(xmlSchemasDirPath));
      space_schema_location.append("/space_h2d_xml.xsd");
      ::xml_schema::namespace_info namespace_info_space("XMLSpace", space_schema_location);

      ::xml_schema::namespace_infomap namespace_info_map;
      namespace_info_map.insert(std::pair<std::basic_string<char>, xml_schema::namespace_info>("space", namespace_info_space));

      std::ofstream out(filename);

      ::xml_schema::flags parsing_flags = ::xml_schema::flags::dont_pretty_print;

      XMLSpace::space_(out, xmlspace, namespace_info_map, "UTF-8", parsing_flags);
      out.close();

      return true;
    }

    template<typename Scalar>
    Space<Scalar>* Space<Scalar>::load(const char *filename, Mesh* mesh, bool validate, EssentialBCs<Scalar>* essential_bcs, Shapeset* shapeset)
    {
      try
      {
				Space<Scalar>* space;

        ::xml_schema::flags parsing_flags = 0;

        if(!validate)
          parsing_flags = xml_schema::flags::dont_validate;

        std::auto_ptr<XMLSpace::space> parsed_xml_space (XMLSpace::space_(filename, parsing_flags));

				if(!strcmp(parsed_xml_space->spaceType().get().c_str(),"h1"))
				{
					space = new H1Space<Scalar>();
					space->mesh = mesh;

					if(shapeset == NULL)
					{
						space->shapeset = new H1Shapeset;
						space->own_shapeset = true;
					}
					else
					{
						if(shapeset->get_space_type() != HERMES_H1_SPACE)
							throw Hermes::Exceptions::SpaceLoadFailureException("Wrong shapeset / Wrong spaceType in the Solution XML file %s in Space::load.", filename);
						else
							space->shapeset = shapeset;
					}

					space->precalculate_projection_matrix(2, space->proj_mat, space->chol_p);
				}
				else if (!strcmp(parsed_xml_space->spaceType().get().c_str(),"hcurl"))
				{
					space = new HcurlSpace<Scalar>();
					space->mesh = mesh;

					if(shapeset == NULL)
					{
						space->shapeset = new HcurlShapeset;
						space->own_shapeset = true;
					}
				else
					{
						if(shapeset->get_num_components() < 2)
							throw Hermes::Exceptions::Exception("HcurlSpace requires a vector shapeset in Space::load.");
						if(shapeset->get_space_type() != HERMES_HCURL_SPACE)
							throw Hermes::Exceptions::SpaceLoadFailureException("Wrong shapeset / Wrong spaceType in the Solution XML file %s in Space::load.", filename);
						else
							space->shapeset = shapeset;
					}

					space->precalculate_projection_matrix(0, space->proj_mat, space->chol_p);
				}
				else if(!!strcmp(parsed_xml_space->spaceType().get().c_str(),"hdiv"))
				{
					space = new HdivSpace<Scalar>();
					space->mesh = mesh;

					if(shapeset == NULL)
					{
						space->shapeset = new HdivShapeset;
						space->own_shapeset = true;
					}
					else
					{
						if(shapeset->get_num_components() < 2)
							throw Hermes::Exceptions::Exception("HdivSpace requires a vector shapeset in Space::load.");
						if(shapeset->get_space_type() != HERMES_HDIV_SPACE)
							throw Hermes::Exceptions::SpaceLoadFailureException("Wrong shapeset / Wrong spaceType in the Solution XML file %s in Space::load.", filename);
						else
							space->shapeset = shapeset;
					}

					space->precalculate_projection_matrix(0, space->proj_mat, space->chol_p);
				}
				else if(strcmp(parsed_xml_space->spaceType().get().c_str(),"l2"))
				{
					space = new L2Space<Scalar>();
					space->mesh = mesh;

					if(shapeset == NULL)
					{
						space->shapeset = new L2Shapeset;
						space->own_shapeset = true;
					}
					{
						if(shapeset->get_space_type() != HERMES_L2_SPACE)
							throw Hermes::Exceptions::SpaceLoadFailureException("Wrong shapeset / Wrong spaceType in the Solution XML file %s in Space::load.", filename);
						else
							space->shapeset = shapeset;
					}

					static_cast<L2Space<Scalar>*>(space)->ldata = NULL;
					static_cast<L2Space<Scalar>*>(space)->lsize = 0;
				}
				else
				{
					throw Exceptions::SpaceLoadFailureException("Wrong spaceType in the Solution XML file %s in Space::load.", filename);
					return NULL;
				}

				space->essential_bcs = essential_bcs;
				space->mesh_seq = space->mesh->get_seq();

				// L2 space does not have any (strong) essential BCs.
				if(essential_bcs != NULL && parsed_xml_space->spaceType().get().c_str() != "l2")
					for(typename Hermes::vector<EssentialBoundaryCondition<Scalar>*>::const_iterator it = essential_bcs->begin(); it != essential_bcs->end(); it++)
						for(unsigned int i = 0; i < (*it)->markers.size(); i++)
							if(space->get_mesh()->boundary_markers_conversion.conversion_table_inverse.find((*it)->markers.at(i)) == space->get_mesh()->boundary_markers_conversion.conversion_table_inverse.end())
								throw Hermes::Exceptions::Exception("A boundary condition defined on a non-existent marker.");

				space->resize_tables();

        // Element data //
        unsigned int elem_data_count = parsed_xml_space->element_data().size();
        for (unsigned int elem_data_i = 0; elem_data_i < elem_data_count; elem_data_i++)
        {
          space->edata[parsed_xml_space->element_data().at(elem_data_i).e_id()].order = parsed_xml_space->element_data().at(elem_data_i).ord();
          space->edata[parsed_xml_space->element_data().at(elem_data_i).e_id()].bdof = parsed_xml_space->element_data().at(elem_data_i).bd();
          space->edata[parsed_xml_space->element_data().at(elem_data_i).e_id()].n = parsed_xml_space->element_data().at(elem_data_i).n();
          space->edata[parsed_xml_space->element_data().at(elem_data_i).e_id()].changed_in_last_adaptation = parsed_xml_space->element_data().at(elem_data_i).chgd();
        }

        space->seq = g_space_seq++;

        space->assign_dofs();

        return space;
      }
      catch (const xml_schema::exception& e)
      {
        throw Hermes::Exceptions::SpaceLoadFailureException(e.what());
      }
    }

    template class HERMES_API Space<double>;
    template class HERMES_API Space<std::complex<double> >;
  }
}
