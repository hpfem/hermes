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

#include "space.h"
namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    Space<Scalar>::Space(Mesh* mesh, Shapeset* shapeset, EssentialBCs<Scalar>* essential_bcs, Ord2 p_init)
      : shapeset(shapeset), essential_bcs(essential_bcs), mesh(mesh)
    {
      _F_
        if (mesh == NULL) error("Space must be initialized with an existing mesh.");
      this->default_tri_order = -1;
      this->default_quad_order = -1;
      this->ndata = NULL;
      this->edata = NULL;
      this->nsize = esize = 0;
      this->ndata_allocated = 0;
      this->mesh_seq = -1;
      this->seq = 0;
      this->was_assigned = false;
      this->ndof = 0;

      if(essential_bcs != NULL)
        for(typename Hermes::vector<EssentialBoundaryCondition<Scalar>*>::const_iterator it = essential_bcs->begin(); it != essential_bcs->end(); it++)
          for(unsigned int i = 0; i < (*it)->markers.size(); i++)
            if(mesh->get_boundary_markers_conversion().conversion_table_inverse->find((*it)->markers.at(i)) == mesh->get_boundary_markers_conversion().conversion_table_inverse->end())
              error("A boundary condition defined on a non-existent marker.");

      own_shapeset = (shapeset == NULL);
    }

    template<typename Scalar>
    Space<Scalar>::~Space()
    {
      _F_
        free();
    }

    template<typename Scalar>
    void Space<Scalar>::free()
    {
      _F_
        free_extra_data();
      if (nsize) { ::free(ndata); ndata=NULL; }
      if (esize) { ::free(edata); edata=NULL; }
    }

    template<typename Scalar>
    void Space<Scalar>::resize_tables()
    {
      _F_
        if ((nsize < mesh->get_max_node_id()) || (ndata == NULL))
        {
          //HACK: definition of allocated size and the result number of elements
          nsize = mesh->get_max_node_id();
          if ((nsize > ndata_allocated) || (ndata == NULL))
          {
            int prev_allocated = ndata_allocated;
            if (ndata_allocated == 0)
              ndata_allocated = 1024;
            while (ndata_allocated < nsize)
              ndata_allocated = ndata_allocated * 3 / 2;
            ndata = (NodeData*)realloc(ndata, ndata_allocated * sizeof(NodeData));
            for(int i = prev_allocated; i < ndata_allocated; i++)
              ndata[i].edge_bc_proj = NULL;
          }
        }

        if ((esize < mesh->get_max_element_id()) || (edata == NULL))
        {
          int oldsize = esize;
          if (!esize) esize = 1024;
          while (esize < mesh->get_max_element_id()) esize = esize * 3 / 2;
          edata = (ElementData*) realloc(edata, sizeof(ElementData) * esize);
          for (int i = oldsize; i < esize; i++)
            edata[i].order = -1;
        }
    }

    template<typename Scalar>
    void Space<Scalar>::H2D_CHECK_ORDER(int order)
    {
      _F_
        if (H2D_GET_H_ORDER(order) < 0 || H2D_GET_V_ORDER(order) < 0)
          error("Hermes::Order cannot be negative.");
      if (H2D_GET_H_ORDER(order) > 10 || H2D_GET_V_ORDER(order) > 10)
        error("Hermes::Order = %d, maximum is 10.", order);
    }

    template<typename Scalar>
    void Space<Scalar>::set_element_order(int id, int order)
    {
      _F_
        set_element_order_internal(id, order);

      // since space changed, enumerate basis functions
      this->assign_dofs();
    }

    template<typename Scalar>
    void Space<Scalar>::set_element_order_internal(int id, int order)
    {
      _F_
        //NOTE: We need to take into account that L2 and Hcurl may use zero orders. The latter has its own version of this method, however.
        assert_msg(mesh->get_element(id)->is_triangle() || get_type() == HERMES_L2_SPACE || H2D_GET_V_ORDER(order) != 0, "Element #%d is quad but given vertical order is zero", id);
      assert_msg(mesh->get_element(id)->is_quad() || H2D_GET_V_ORDER(order) == 0, "Element #%d is triangle but vertical is not zero", id);
      if (id < 0 || id >= mesh->get_max_element_id())
        error("Invalid element id.");
      H2D_CHECK_ORDER(order);

      resize_tables();
      if (mesh->get_element(id)->is_quad() && get_type() != HERMES_L2_SPACE && H2D_GET_V_ORDER(order) == 0)
        order = H2D_MAKE_QUAD_ORDER(order, order);
      edata[id].order = order;
      seq++;
    }


    template<typename Scalar>
    int Space<Scalar>::get_element_order(int id) const
    {
      _F_
        // sanity checks (for internal purposes)
        if (this->mesh == NULL) error("NULL Mesh pointer detected in Space<Scalar>::get_element_order().");
      if(edata == NULL) error("NULL edata detected in Space<Scalar>::get_element_order().");
      if (id >= esize) 
      {
        warn("Element index %d in Space<Scalar>::get_element_order() while maximum is %d.", id, esize);
        error("Wring element index in Space<Scalar>::get_element_order().");
      }
      return edata[id].order;
    }


    template<typename Scalar>
    void Space<Scalar>::set_uniform_order(int order, std::string marker)
    {
      _F_
        if(marker == HERMES_ANY)
          set_uniform_order_internal(Ord2(order,order), -1234);
        else
          set_uniform_order_internal(Ord2(order,order), mesh->element_markers_conversion.get_internal_marker(marker));

      // since space changed, enumerate basis functions
      this->assign_dofs();
    }

    template<typename Scalar>
    void Space<Scalar>::set_uniform_order_internal(Ord2 order, int marker)
    {
      _F_
        resize_tables();
      if (order.order_h < 0 || order.order_v < 0)
        error("Hermes::Order cannot be negative.");
      if (order.order_h > 10 || order.order_v > 10)
        error("Hermes::Order = %d x %d, maximum is 10.", order.order_h, order.order_v);
      int quad_order = H2D_MAKE_QUAD_ORDER(order.order_h, order.order_v);

      Element* e;
      for_all_active_elements(e, mesh)
      {
        if (marker == HERMES_ANY_INT || e->marker == marker)
        {
          ElementData* ed = &edata[e->id];
          if (e->is_triangle())
            if(order.order_h != order.order_v)
              error("Hermes::Orders do not match and triangles are present in the mesh.");
            else
              ed->order = order.order_h;
          else
            ed->order = quad_order;
        }
      }
      seq++;
    }

    template<typename Scalar>
    void Space<Scalar>::set_element_orders(int* elem_orders_)
    {
      _F_
        resize_tables();

      Element* e;
      int counter = 0;
      for_all_elements(e, mesh)
      {
        H2D_CHECK_ORDER(elem_orders_[counter]);
        ElementData* ed = &edata[e->id];
        if (e->is_triangle())
          ed->order = elem_orders_[counter];
        else
          ed->order = H2D_MAKE_QUAD_ORDER(elem_orders_[counter], elem_orders_[counter]);
        counter++;
      }
    }

    template<typename Scalar>
    void Space<Scalar>::set_default_order(int tri_order, int quad_order)
    {
      _F_
        if (quad_order == -1) quad_order = H2D_MAKE_QUAD_ORDER(tri_order, tri_order);
      default_tri_order = tri_order;
      default_quad_order = quad_order;
    }

    template<typename Scalar>
    void Space<Scalar>::adjust_element_order(int order_change, int min_order)
    {
      _F_
        Element* e;
      for_all_active_elements(e, this->get_mesh()) 
      {
        if(e->is_triangle())
          set_element_order_internal(e->id, std::max<int>(min_order, get_element_order(e->id) + order_change));
        else
        {
          int h_order, v_order;
          // check that we are not imposing smaller than minimal orders.
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
      assign_dofs();
    }

    template<typename Scalar>
    void Space<Scalar>::adjust_element_order(int horizontal_order_change, int vertical_order_change, unsigned int horizontal_min_order, unsigned int vertical_min_order)
    {
      _F_
        Element* e;
      for_all_active_elements(e, this->get_mesh()) 
      {
        if(e->is_triangle()) 
        {
          warn("Using quad version of Space<Scalar>::adjust_element_order(), only horizontal orders will be used.");
          set_element_order_internal(e->id, std::max<int>(horizontal_min_order, get_element_order(e->id) + horizontal_order_change));
        }
        else
          set_element_order_internal(e->id, std::max<int>
          (H2D_MAKE_QUAD_ORDER(horizontal_min_order, vertical_min_order), 
          H2D_MAKE_QUAD_ORDER(H2D_GET_H_ORDER(get_element_order(e->id)) + horizontal_order_change, H2D_GET_V_ORDER(get_element_order(e->id)) + vertical_order_change)));
      }
      assign_dofs();
    }

    template<typename Scalar>
    void Space<Scalar>::unrefine_all_mesh_elements(bool keep_initial_refinements)
    {
      // find inactive elements with active sons
      Hermes::vector<int> list;
      Element* e;
      for_all_inactive_elements(e, this->mesh)
      {
        bool found = true;
        for (unsigned int i = 0; i < 4; i++)
          if (e->sons[i] != NULL && 
            (!e->sons[i]->active || (keep_initial_refinements && e->sons[i]->id < this->mesh->ninitial))  
            )
          { found = false; break; }

          if (found) list.push_back(e->id);
      }

      // unrefine the found elements
      for (unsigned int i = 0; i < list.size(); i++) 
      {
        unsigned int order = 0, h_order = 0, v_order = 0;
        unsigned int num_sons = 0;
        if (this->mesh->get_element_fast(list[i])->bsplit()) 
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
          if (this->mesh->get_element_fast(list[i])->hsplit()) 
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
        this->mesh->unrefine_element_id(list[i]);
      }

      this->assign_dofs();
    }

    template<typename Scalar>
    void Space<Scalar>::copy_orders_recurrent(Element* e, int order)
    {
      _F_
        if (e->active)
          edata[e->id].order = order;
        else
          for (int i = 0; i < 4; i++)
            if (e->sons[i] != NULL)
              copy_orders_recurrent(e->sons[i], order);
    }


    template<typename Scalar>
    void Space<Scalar>::copy_orders(const Space<Scalar>* space, int inc)
    {
      _F_
        Element* e;
      resize_tables();
      for_all_active_elements(e, space->get_mesh())
      {
        int o = space->get_element_order(e->id);
        if (o < 0) error("Source space has an uninitialized order (element id = %d)", e->id);

        int mo = shapeset->get_max_order();
        int lower_limit = (get_type() == HERMES_L2_SPACE || get_type() == HERMES_HCURL_SPACE) ? 0 : 1; // L2 and Hcurl may use zero orders.
        int ho = std::max(lower_limit, std::min(H2D_GET_H_ORDER(o) + inc, mo));
        int vo = std::max(lower_limit, std::min(H2D_GET_V_ORDER(o) + inc, mo));
        o = e->is_triangle() ? ho : H2D_MAKE_QUAD_ORDER(ho, vo);

        H2D_CHECK_ORDER(o);
        copy_orders_recurrent(mesh->get_element/*sic!*/(e->id), o);
      }
      seq++;

      // since space changed, enumerate basis functions
      this->assign_dofs();
    }


    template<typename Scalar>
    int Space<Scalar>::get_edge_order(Element* e, int edge)
    {
      _F_
        Node* en = e->en[edge];
      if (en->id >= nsize || edge >= (int)e->nvert) return 0;

      if (ndata[en->id].n == -1)
        return get_edge_order_internal(ndata[en->id].base); // constrained node
      else
        return get_edge_order_internal(en);
    }


    template<typename Scalar>
    int Space<Scalar>::get_edge_order_internal(Node* en)
    {
      _F_
        assert(en->type == HERMES_TYPE_EDGE);
      Element** e = en->elem;
      int o1 = 1000, o2 = 1000;
      assert(e[0] != NULL || e[1] != NULL);

      if (e[0] != NULL)
      {
        if (e[0]->is_triangle() || en == e[0]->en[0] || en == e[0]->en[2])
          o1 = H2D_GET_H_ORDER(edata[e[0]->id].order);
        else
          o1 = H2D_GET_V_ORDER(edata[e[0]->id].order);
      }

      if (e[1] != NULL)
      {
        if (e[1]->is_triangle() || en == e[1]->en[0] || en == e[1]->en[2])
          o2 = H2D_GET_H_ORDER(edata[e[1]->id].order);
        else
          o2 = H2D_GET_V_ORDER(edata[e[1]->id].order);
      }

      if (o1 == 0) return o2 == 1000 ? 0 : o2;
      if (o2 == 0) return o1 == 1000 ? 0 : o1;
      return std::min(o1, o2);
    }


    template<typename Scalar>
    void Space<Scalar>::set_mesh(Mesh* mesh)
    {
      _F_
        if (this->mesh == mesh) return;
      free();
      this->mesh = mesh;
      seq++;

      // since space changed, enumerate basis functions
      this->assign_dofs();
    }


    template<typename Scalar>
    void Space<Scalar>::propagate_zero_orders(Element* e)
    {
      _F_
        warn_if(get_element_order(e->id) != 0, "zeroing order of an element ID:%d, original order (H:%d; V:%d)", e->id, H2D_GET_H_ORDER(get_element_order(e->id)), H2D_GET_V_ORDER(get_element_order(e->id)));
      set_element_order_internal(e->id, 0);
      if (!e->active)
        for (int i = 0; i < 4; i++)
          if (e->sons[i] != NULL)
            propagate_zero_orders(e->sons[i]);
    }


    template<typename Scalar>
    void Space<Scalar>::distribute_orders(Mesh* mesh, int* parents)
    {
      _F_
        int num = mesh->get_max_element_id();
      int* orders = new int[num+1];
      Element* e;
      for_all_active_elements(e, mesh)
      {
        int p = get_element_order(parents[e->id]);
        if (e->is_triangle() && (H2D_GET_V_ORDER(p) != 0))
          p = std::max(H2D_GET_H_ORDER(p), H2D_GET_V_ORDER(p));
        orders[e->id] = p;
      }
      for_all_active_elements(e, mesh)
        set_element_order_internal(e->id, orders[e->id]);
      delete [] orders;
    }


    //// dof assignment ////////////////////////////////////////////////////////////////////////////////

    template<typename Scalar>
    int Space<Scalar>::assign_dofs(int first_dof, int stride)
    {
      _F_
        if (first_dof < 0) error("Invalid first_dof.");
      if (stride < 1)    error("Invalid stride.");

      resize_tables();

      Element* e;
      /** \todo Find out whether the following code this is crucial.
      *  If uncommented, this enforces 0 order for all sons if the base element has 0 order.
      *  In this case, an element with 0 order means an element which is left out from solution. */
      //for_all_base_elements(e, mesh)
      //  if (get_element_order(e->id) == 0)
      //    propagate_zero_orders(e);

      //check validity of orders
      for_all_active_elements(e, mesh) 
      {
        if (e->id >= esize || edata[e->id].order < 0) 
        {
          printf("e->id = %d\n", e->id);
          printf("esize = %d\n", esize);
          printf("edata[%d].order = %d\n", e->id, edata[e->id].order);
          error("Uninitialized element order.");
        }
      }

      this->first_dof = next_dof = first_dof;
      this->stride = stride;

      reset_dof_assignment();
      assign_vertex_dofs();
      assign_edge_dofs();
      assign_bubble_dofs();

      free_extra_data();
      update_essential_bc_values();
      update_constraints();
      post_assign();

      mesh_seq = mesh->get_seq();
      was_assigned = true;
      this->ndof = (next_dof - first_dof) / stride;

      return this->ndof;
    }

    template<typename Scalar>
    void Space<Scalar>::reset_dof_assignment()
    {
      _F_
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
        for (unsigned int i = 0; i < e->nvert; i++)
        {
          if (e->en[i]->bnd)
            if(essential_bcs != NULL)
              if(essential_bcs->get_boundary_condition(mesh->boundary_markers_conversion.get_user_marker(e->en[i]->marker)) != NULL) 
              {
                j = e->next_vert(i);
                ndata[e->vn[i]->id].n = 0;
                ndata[e->vn[j]->id].n = 0;
              }
        }
      }
    }

    template<typename Scalar>
    void Space<Scalar>::get_element_assembly_list(Element* e, AsmList<Scalar>* al)
    {
      _F_
        // some checks
        if (e->id >= esize || edata[e->id].order < 0)
          error("Uninitialized element order (id = #%d).", e->id);
      if (!is_up_to_date())
        error("The space is out of date. You need to update it with assign_dofs()"
        " any time the mesh changes.");

      // add vertex, edge and bubble functions to the assembly list
      al->clear();
      shapeset->set_mode(e->get_mode());
      for (unsigned int i = 0; i < e->nvert; i++)
        get_vertex_assembly_list(e, i, al);
      for (unsigned int i = 0; i < e->nvert; i++)
        get_boundary_assembly_list_internal(e, i, al);
      get_bubble_assembly_list(e, al);
    }


    template<typename Scalar>
    void Space<Scalar>::get_boundary_assembly_list(Element* e, int surf_num, AsmList<Scalar>* al)
    {
      _F_
        al->clear();
      shapeset->set_mode(e->get_mode());
      get_vertex_assembly_list(e, surf_num, al);
      get_vertex_assembly_list(e, e->next_vert(surf_num), al);
      get_boundary_assembly_list_internal(e, surf_num, al);
    }


    template<typename Scalar>
    void Space<Scalar>::get_bubble_assembly_list(Element* e, AsmList<Scalar>* al)
    {
      _F_
        ElementData* ed = &edata[e->id];

      if (!ed->n) return;

      int* indices = shapeset->get_bubble_indices(ed->order);
      for (int i = 0, dof = ed->bdof; i < ed->n; i++, dof += stride, indices++)
        al->add_triplet(*indices, dof, 1.0);
    }

    //// BC stuff /////////////////////////////////////////////////////////////////////////////////////
    template<typename Scalar>
    void Space<Scalar>::set_essential_bcs(EssentialBCs<Scalar>* essential_bcs)
    {
      _F_
        this->essential_bcs = essential_bcs;

      // since space changed, enumerate basis functions
      this->assign_dofs();
    }

    template<typename Scalar>
    void Space<Scalar>::precalculate_projection_matrix(int nv, double**& mat, double*& p)
    {
      _F_
        int n = shapeset->get_max_order() + 1 - nv;
      mat = new_matrix<double>(n, n);
      int component = (get_type() == HERMES_HDIV_SPACE) ? 1 : 0;

      Quad1DStd quad1d;
      //shapeset->set_mode(HERMES_MODE_TRIANGLE);
      shapeset->set_mode(HERMES_MODE_QUAD);
      for (int i = 0; i < n; i++)
      {
        for (int j = i; j < n; j++)
        {
          int o = i + j + 4;
          double2* pt = quad1d.get_points(o);
          int ii = shapeset->get_edge_index(0, 0, i + nv);
          int ij = shapeset->get_edge_index(0, 0, j + nv);
          double val = 0.0;
          for (int k = 0; k < quad1d.get_num_points(o); k++)
          {
            val += pt[k][1] * shapeset->get_fn_value(ii, pt[k][0], -1.0, component)
              * shapeset->get_fn_value(ij, pt[k][0], -1.0, component);
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
      _F_
        if (e->active)
        {
          Node* en = e->en[surf_pos->surf_num];
          NodeData* nd = &ndata[en->id];
          nd->edge_bc_proj = NULL;

          if (nd->dof != H2D_UNASSIGNED_DOF && en->bnd)
            if(essential_bcs != NULL)
              if(essential_bcs->get_boundary_condition(mesh->boundary_markers_conversion.get_user_marker(en->marker)) != NULL) 
              {
                int order = get_edge_order_internal(en);
                surf_pos->marker = en->marker;
                nd->edge_bc_proj = get_bc_projection(surf_pos, order);
                extra_data.push_back(nd->edge_bc_proj);

                int i = surf_pos->surf_num, j = e->next_vert(i);
                ndata[e->vn[i]->id].vertex_bc_coef = nd->edge_bc_proj + 0;
                ndata[e->vn[j]->id].vertex_bc_coef = nd->edge_bc_proj + 1;
              }
        }
        else
        {
          int son1, son2;
          if (mesh->get_edge_sons(e, surf_pos->surf_num, son1, son2) == 2)
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
      _F_
        Element* e;
      for_all_base_elements(e, mesh)
      {
        for (unsigned int i = 0; i < e->nvert; i++)
        {
          int j = e->next_vert(i);
          if (e->vn[i]->bnd && e->vn[j]->bnd)
          {
            SurfPos surf_pos = {0, i, e, e->vn[i]->id, e->vn[j]->id, 0.0, 0.0, 1.0};
            update_edge_bc(e, &surf_pos);
          }
        }
      }
    }


    template<typename Scalar>
    void Space<Scalar>::free_extra_data()
    {
      _F_
        for (unsigned int i = 0; i < extra_data.size(); i++)
          delete [] (Scalar*) extra_data[i];
      extra_data.clear();
    }

    template HERMES_API class Space<double>;
    template HERMES_API class Space<std::complex<double> >;
  }
}