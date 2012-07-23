// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#include "neighbor.h"
#include <algorithm>

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    NeighborSearch<Scalar>::NeighborSearch(Element* el, const Mesh* mesh) :
    supported_shapes(NULL),
      mesh(mesh),
      central_transformations(LightArray<Transformations*>(4)),
      neighbor_transformations(LightArray<Transformations*>(4)),
      central_el(el),
      neighb_el(NULL),
      quad(&g_quad_2d_std)
    {
      if(central_el == NULL || central_el->active != 1)
        throw Exceptions::Exception("You must pass an active element to the NeighborSearch constructor.");
      neighbors.reserve(2);
      neighbor_edges.reserve(2);

      ignore_errors = false;
      n_neighbors = 0;
      neighborhood_type = H2D_DG_NOT_INITIALIZED;
      original_central_el_transform = 0;
    }

    template<typename Scalar>
    NeighborSearch<Scalar>::NeighborSearch(const NeighborSearch& ns) :
    supported_shapes(NULL),
      mesh(ns.mesh),
      central_transformations(LightArray<Transformations*>(4)),
      neighbor_transformations(LightArray<Transformations*>(4)),
      central_el(ns.central_el),
      neighb_el(NULL),
      neighbor_edge(ns.neighbor_edge),
      active_segment(ns.active_segment)
    {
      neighbors.reserve(2);
      neighbor_edges.reserve(2);

      for(unsigned int j = 0; j < ns.central_transformations.get_size(); j++)
        if(ns.central_transformations.present(j))
        {
          // Copy current array of transformations into new memory location.
          Transformations *tmp = new Transformations(ns.central_transformations.get(j));
          this->central_transformations.add(tmp, j);
        }
        for(unsigned int j = 0; j < ns.neighbor_transformations.get_size(); j++)
          if(ns.neighbor_transformations.present(j))
          {
            // Copy current array of transformations into new memory location.
            Transformations *tmp = new Transformations(ns.neighbor_transformations.get(j));
            this->neighbor_transformations.add(tmp, j);
          }

          if(central_el == NULL || central_el->active != 1)
            throw Exceptions::Exception("You must pass an active element to the NeighborSearch constructor.");

          for(unsigned int i = 0; i < ns.neighbors.size(); i++)
            this->neighbors.push_back(ns.neighbors[i]);
          for(unsigned int i = 0; i < ns.neighbor_edges.size(); i++)
            this->neighbor_edges.push_back(ns.neighbor_edges[i]);

          ignore_errors = ns.ignore_errors;
          n_neighbors = ns.n_neighbors;
          neighborhood_type = ns.neighborhood_type;
          original_central_el_transform = ns.original_central_el_transform;
          quad = (&g_quad_2d_std);
          active_edge = ns.active_edge;
    }

    template<typename Scalar>
    NeighborSearch<Scalar>::~NeighborSearch()
    {
      neighbor_edges.clear();
      neighbors.clear();
      clear_supported_shapes();

      for(unsigned int i = 0; i < central_transformations.get_size(); i++)
        if(this->central_transformations.present(i))
          delete this->central_transformations.get(i);
      for(unsigned int i = 0; i < neighbor_transformations.get_size(); i++)
        if(this->neighbor_transformations.present(i))
          delete this->neighbor_transformations.get(i);
    }

    template<typename Scalar>
    int NeighborSearch<Scalar>::get_num_neighbors() const
    {
      return n_neighbors;
    }

    template<typename Scalar>
    const Hermes::vector<Element*>* NeighborSearch<Scalar>::get_neighbors() const
    {
      return &neighbors;
    }

    template<typename Scalar>
    void NeighborSearch<Scalar>::clear_supported_shapes()
    {
      if(supported_shapes != NULL) delete supported_shapes; supported_shapes = NULL;
    }

    template<typename Scalar>
    void NeighborSearch<Scalar>::set_ignore_errors(bool value)
    {
      this->ignore_errors = value;
    }

    template<typename Scalar>
    void NeighborSearch<Scalar>::reset_neighb_info()
    {
      // Reset transformations.
      for(unsigned int i = 0; i < n_neighbors; i++)
      {
        if(this->central_transformations.present(i))
          this->central_transformations.get(i)->reset();
        if(this->neighbor_transformations.present(i))
          this->neighbor_transformations.get(i)->reset();
      }

      // Reset information about the neighborhood's active state.
      active_segment = 0;
      active_edge = 0;
      neighb_el = NULL;
      neighbor_edge.local_num_of_edge = 0;

      // Clear vectors with neighbor elements and their edge info for the active edge.
      neighbor_edges.clear();
      neighbors.clear();
      n_neighbors = 0;

      neighborhood_type = H2D_DG_NOT_INITIALIZED;
    }

    template<typename Scalar>
    void NeighborSearch<Scalar>::set_active_edge(int edge)
    {
      reset_neighb_info();
      active_edge = edge;

      //debug_log("central element: %d", central_el->id);
      if(central_el->en[active_edge]->bnd == 0)
      {
        neighb_el = central_el->get_neighbor(active_edge);

        // First case : The neighboring element is of the same size as the central one.
        if(neighb_el != NULL)
        {
          //debug_log("active neighbor el: %d", neighb_el->id);

          // Get local number of the edge used by the neighbor.
          for (unsigned int j = 0; j < neighb_el->get_num_surf(); j++)
            if(central_el->en[active_edge] == neighb_el->en[j])
            {
              neighbor_edge.local_num_of_edge = j;
              break;
            }

            NeighborEdgeInfo local_edge_info;
            local_edge_info.local_num_of_edge = neighbor_edge.local_num_of_edge;

            // Query the orientation of the neighbor edge relative to the central el.
            int p1 = central_el->vn[active_edge]->id;
            int p2 = central_el->vn[(active_edge + 1) % central_el->get_num_surf()]->id;
            local_edge_info.orientation = neighbor_edge_orientation(p1, p2, 0);

            neighbor_edges.push_back(local_edge_info);

            // There is only one neighbor in this case.
            n_neighbors = 1;
            neighbors.push_back(neighb_el);

            // No need for transformation, since the neighboring element is of the same size.
            neighborhood_type = H2D_DG_NO_TRANSF;
        }
        else
        {
          // Peek the vertex in the middle of the active edge (if there is none, vertex will be NULL).
          Node* vertex = mesh->peek_vertex_node(central_el->en[active_edge]->p1,  central_el->en[active_edge]->p2);

          // Endpoints of the active edge.
          int orig_vertex_id[2];
          orig_vertex_id[0] = central_el->vn[active_edge]->id;
          orig_vertex_id[1]  = central_el->vn[(active_edge + 1) % central_el->get_num_surf()]->id;

          if(vertex == NULL)
          {
            neighborhood_type = H2D_DG_GO_UP;

            Element* parent = central_el->parent;

            // Array of middle-point vertices of the intermediate parent edges that we climb up to the correct parent element.
            Node** par_mid_vertices = new Node*[Transformations::max_level];
            // Number of visited intermediate parents.
            int n_parents = 0;

            for (unsigned int j = 0; j < (unsigned) Transformations::max_level; j++)
              par_mid_vertices[j] = NULL;

            find_act_elem_up(parent, orig_vertex_id, par_mid_vertices, n_parents);

            delete [] par_mid_vertices;
          }
          else
          {
            neighborhood_type = H2D_DG_GO_DOWN;

            int sons[Transformations::max_level]; // array of virtual sons of the central el. visited on the way down to the neighbor
            int n_sons = 0; // number of used transformations

            // Start the search by going down to the first son.
            find_act_elem_down( vertex, orig_vertex_id, sons, n_sons + 1);

            //debug_log("number of neighbors on the way down: %d ", n_neighbors);
          }
        }
      }
      else
        if(!ignore_errors)
          throw Hermes::Exceptions::Exception("The given edge isn't inner");
    }

    template<typename Scalar>
    bool NeighborSearch<Scalar>::set_active_edge_multimesh(const int& edge)
    {
      Hermes::vector<unsigned int> transformations = get_transforms(original_central_el_transform);
      // Inter-element edge.
      if(is_inter_edge(edge, transformations))
      {
        set_active_edge(edge);
        update_according_to_sub_idx(transformations);
        return true;
      }
      // Intra-element edge.
      else
      {
        neighb_el = central_el;

        neighbor_transformations.add(new Transformations(transformations), 0);

        neighbor_edge.local_num_of_edge = active_edge = edge;
        NeighborEdgeInfo local_edge_info;
        local_edge_info.local_num_of_edge = neighbor_edge.local_num_of_edge;
        // The "opposite" view of the same edge has the same orientation.
        local_edge_info.orientation = 0;
        neighbor_edges.push_back(local_edge_info);

        n_neighbors = 1;
        neighbors.push_back(neighb_el);
        neighborhood_type = H2D_DG_NO_TRANSF;
        return false;
      }
    }

    template<typename Scalar>
    Hermes::vector<unsigned int> NeighborSearch<Scalar>::get_transforms(uint64_t sub_idx) const
    {
      Hermes::vector<unsigned int> transformations_backwards;
      while (sub_idx > 0)
      {
        transformations_backwards.push_back((sub_idx - 1) & 7);
        sub_idx = (sub_idx - 1) >> 3;
      }
      Hermes::vector<unsigned int> transformations;
      for(unsigned int i = 0; i < transformations_backwards.size(); i++)
        transformations.push_back(transformations_backwards[transformations_backwards.size() - 1 - i]);

      return transformations;
    }

    template<typename Scalar>
    bool NeighborSearch<Scalar>::is_inter_edge(const int& edge, const Hermes::vector<unsigned int>& transformations) const
    {
      // No subelements => of course this edge is an inter-element one.
      if(transformations.size() == 0)
        return true;

      // Triangles.
      for(unsigned int i = 0; i < transformations.size(); i++)
        if(central_el->get_mode() == HERMES_MODE_TRIANGLE)
        {
          if((edge == 0 && (transformations[i] == 2 || transformations[i] == 3)) ||
            (edge == 1 && (transformations[i] == 0 || transformations[i] == 3)) ||
            (edge == 2 && (transformations[i] == 1 || transformations[i] == 3)))
            return false;
        }
        // Quads.
        else
        {
          if((edge == 0 && (transformations[i] == 2 || transformations[i] == 3 || transformations[i] == 5)) ||
            (edge == 1 && (transformations[i] == 0 || transformations[i] == 3 || transformations[i] == 6)) ||
            (edge == 2 && (transformations[i] == 0 || transformations[i] == 1 || transformations[i] == 4)) ||
            (edge == 3 && (transformations[i] == 1 || transformations[i] == 2 || transformations[i] == 7)))
            return false;
        }
        return true;
    }

    template<typename Scalar>
    void NeighborSearch<Scalar>::update_according_to_sub_idx(const Hermes::vector<unsigned int>& transformations)
    {
      if(neighborhood_type == H2D_DG_NO_TRANSF || neighborhood_type == H2D_DG_GO_UP)
      {
        if(!neighbor_transformations.present(0)) // in case of neighborhood_type == H2D_DG_NO_TRANSF
          neighbor_transformations.add(new Transformations, 0);
        Transformations *tr = neighbor_transformations.get(0);

        for(unsigned int i = 0; i < transformations.size(); i++)
          // Triangles.
          if(central_el->get_mode() == HERMES_MODE_TRIANGLE)
            if((active_edge == 0 && transformations[i] == 0) ||
              (active_edge == 1 && transformations[i] == 1) ||
              (active_edge == 2 && transformations[i] == 2))
              tr->transf[tr->num_levels++] = (!neighbor_edge.orientation ? neighbor_edge.local_num_of_edge : (neighbor_edge.local_num_of_edge + 1) % 3);
            else
              tr->transf[tr->num_levels++] = (neighbor_edges[0].orientation ? neighbor_edge.local_num_of_edge : (neighbor_edge.local_num_of_edge + 1) % 3);
        // Quads.
          else
            if((active_edge == 0 && (transformations[i] == 0 || transformations[i] == 6)) ||
              (active_edge == 1 && (transformations[i] == 1 || transformations[i] == 4)) ||
              (active_edge == 2 && (transformations[i] == 2 || transformations[i] == 7)) ||
              (active_edge == 3 && (transformations[i] == 3 || transformations[i] == 5)))
              tr->transf[tr->num_levels++] = (!neighbor_edge.orientation ? neighbor_edge.local_num_of_edge : (neighbor_edge.local_num_of_edge + 1) % 4);
            else if((active_edge == 0 && (transformations[i] == 1 || transformations[i] == 7)) ||
              (active_edge == 1 && (transformations[i] == 2 || transformations[i] == 5)) ||
              (active_edge == 2 && (transformations[i] == 3 || transformations[i] == 6)) ||
              (active_edge == 3 && (transformations[i] == 0 || transformations[i] == 4)))
              tr->transf[tr->num_levels++] = (neighbor_edge.orientation ? neighbor_edge.local_num_of_edge : (neighbor_edge.local_num_of_edge + 1) % 4);
      }
      else handle_sub_idx_way_down(transformations);
    }

    template<typename Scalar>
    void NeighborSearch<Scalar>::handle_sub_idx_way_down(const Hermes::vector<unsigned int>& transformations)
    {
      Hermes::vector<unsigned int> neighbors_to_be_deleted;
      Hermes::vector<unsigned int> neighbors_not_to_be_deleted;

      Hermes::vector<unsigned int> updated_transformations;
      for(int i = 0; i < transformations.size(); i++)
      {
        if(! ((active_edge == 0 && transformations[i] == 4) || (active_edge == 1 && transformations[i] == 7) || (active_edge == 2 && transformations[i] == 5) || (active_edge == 3 && transformations[i] == 6)) )
        {
          if(active_edge == 0 && transformations[i] == 6)
            updated_transformations.push_back(0);
          else if(active_edge == 0 && transformations[i] == 7)
            updated_transformations.push_back(1);
          else if(active_edge == 1 && transformations[i] == 4)
            updated_transformations.push_back(1);
          else if(active_edge == 1 && transformations[i] == 5)
            updated_transformations.push_back(2);
          else if(active_edge == 2 && transformations[i] == 6)
            updated_transformations.push_back(3);
          else if(active_edge == 2 && transformations[i] == 7)
            updated_transformations.push_back(2);
          else if(active_edge == 3 && transformations[i] == 4)
            updated_transformations.push_back(0);
          else if(active_edge == 3 && transformations[i] == 5)
            updated_transformations.push_back(3);
          else
            updated_transformations.push_back(transformations[i]);
        }
      }

      // We basically identify the neighbors that are not compliant with the current sub-element mapping on the central element.
      for(unsigned int neighbor_i = 0; neighbor_i < n_neighbors; neighbor_i++)
      {
        bool deleted = false;

        Transformations* current_transforms = central_transformations.get(neighbor_i);

        for(unsigned int level = 0; level < std::min((unsigned int)updated_transformations.size(), current_transforms->num_levels); level++)
        {
          // If the found neighbor is not a neighbor of this subelement.
          if(!compatible_transformations(current_transforms->transf[level], updated_transformations[level], active_edge))
          {
            deleted = true;
            break;
          }
        }

        /*
        // If there were more sub-element transformation from the assembling than from the neighbor search.
        if(!deleted)
        {
          if((unsigned int)updated_transformations.size() > current_transforms->num_levels)
          {
            for(unsigned int level = current_transforms->num_levels; level < (unsigned int)updated_transformations.size(); level++)
            {
              // If the found neighbor is not a neighbor of this subelement.
              if(!compatible_transformations(current_transforms->transf[current_transforms->num_levels - 1], updated_transformations[level], active_edge))
              {
                deleted = true;
                break;
              }
            }
          }
        }
        */
        if(deleted)
          neighbors_to_be_deleted.push_back(neighbor_i);
        else
          neighbors_not_to_be_deleted.push_back(neighbor_i);
      }

      // Now we truly delete (in the reverse order) the neighbors.
      if(neighbors_to_be_deleted.size() > 0)
        for(unsigned int neighbors_to_be_deleted_i = neighbors_to_be_deleted.size(); neighbors_to_be_deleted_i >= 1; neighbors_to_be_deleted_i--)
          delete_neighbor(neighbors_to_be_deleted[neighbors_to_be_deleted_i - 1]);
    }

    template<typename Scalar>
    bool NeighborSearch<Scalar>::compatible_transformations(unsigned int a, unsigned int b, int edge) const
    {
      if(a == b)
        return true;
      if(edge == 0)
      {
        if((a == 0 && (b == 6 || b == 4)) ||
          (a == 1 && (b == 7 || b == 4)))
          return true;
        else
          return false;
      }
      if(edge == 1)
      {
        if((a == 1 && (b == 4 || b == 7)) ||
          (a == 2 && (b == 5 || b == 7)))
          return true;
        else
          return false;
      }
      if(edge == 2)
      {
        if((a == 2 && (b == 7 || b == 5)) ||
          (a == 3 && (b == 6 || b == 5)))
          return true;
        else
          return false;
      }
      if(edge == 3)
      {
        if((a == 3 && (b == 5 || b == 6)) ||
          (a == 0 && (b == 4 || b == 6)))
          return true;
        else
          return false;
      }
      return false;
    }

    template<typename Scalar>
    void NeighborSearch<Scalar>::clear_initial_sub_idx()
    {
      if(neighborhood_type != H2D_DG_GO_DOWN)
        return;
      // Obtain the transformations sequence.
      Hermes::vector<unsigned int> transformations = get_transforms(original_central_el_transform);

      Hermes::vector<unsigned int> updated_transformations;
      for(int i = 0; i < transformations.size(); i++)
      {
        if(! ((active_edge == 0 && transformations[i] == 4) || (active_edge == 1 && transformations[i] == 7) || (active_edge == 2 && transformations[i] == 5) || (active_edge == 3 && transformations[i] == 6)) )
        {
          if(active_edge == 0 && transformations[i] == 6)
            updated_transformations.push_back(0);
          else if(active_edge == 0 && transformations[i] == 7)
            updated_transformations.push_back(1);
          else if(active_edge == 1 && transformations[i] == 4)
            updated_transformations.push_back(1);
          else if(active_edge == 1 && transformations[i] == 5)
            updated_transformations.push_back(2);
          else if(active_edge == 2 && transformations[i] == 6)
            updated_transformations.push_back(3);
          else if(active_edge == 2 && transformations[i] == 7)
            updated_transformations.push_back(2);
          else if(active_edge == 3 && transformations[i] == 4)
            updated_transformations.push_back(0);
          else if(active_edge == 3 && transformations[i] == 5)
            updated_transformations.push_back(3);
          else
            updated_transformations.push_back(transformations[i]);
        }
      }

      // Test for active element.
      if(updated_transformations.size() == 0)
        return;

      for(unsigned int i = 0; i < n_neighbors; i++)
      {
        // Find the index where the additional subelement mapping (on top of the initial one from assembling) starts.
        unsigned int j = 0;
        // Note that we do not have to test if central_transformations is empty or how long it is, because it has to be
        // longer than transformations (and that is tested).
        // Also the function compatible_transformations() does not have to be used, as now the array central_transformations
        // has been adjusted so that it contains the array transformations.
        while(central_transformations.get(i)->transf[j] == updated_transformations[j])
          if(++j > updated_transformations.size() - 1)
            break;
        if(j > central_transformations.get(i)->num_levels)
          j = central_transformations.get(i)->num_levels;

        for(unsigned int level = central_transformations.get(i)->num_levels; level < updated_transformations.size(); level++)
        {
          if(!neighbor_transformations.present(i))
            neighbor_transformations.add(new Transformations, i);

          Transformations* neighbor_transforms = neighbor_transformations.get(i);

          // Triangles.
          if(central_el->get_mode() == HERMES_MODE_TRIANGLE)
            if((active_edge == 0 && updated_transformations[level] == 0) ||
              (active_edge == 1 && updated_transformations[level] == 1) ||
              (active_edge == 2 && updated_transformations[level] == 2))
              neighbor_transforms->transf[neighbor_transforms->num_levels++] = (!neighbor_edge.orientation ? neighbor_edge.local_num_of_edge : (neighbor_edge.local_num_of_edge + 1) % 3);
            else
              neighbor_transforms->transf[neighbor_transforms->num_levels++] = (neighbor_edge.orientation ? neighbor_edge.local_num_of_edge : (neighbor_edge.local_num_of_edge + 1) % 3);
          // Quads.
          else
            if((active_edge == 0 && (updated_transformations[level] == 0 || updated_transformations[level] == 6)) ||
              (active_edge == 1 && (updated_transformations[level] == 1 || updated_transformations[level] == 4)) ||
              (active_edge == 2 && (updated_transformations[level] == 2 || updated_transformations[level] == 7)) ||
              (active_edge == 3 && (updated_transformations[level] == 3 || updated_transformations[level] == 5)))
              neighbor_transforms->transf[neighbor_transforms->num_levels++] = (!neighbor_edge.orientation ? neighbor_edge.local_num_of_edge : (neighbor_edge.local_num_of_edge + 1) % 4);
            else if((active_edge == 0 && (updated_transformations[level] == 1 || updated_transformations[level] == 7)) ||
              (active_edge == 1 && (updated_transformations[level] == 2 || updated_transformations[level] == 5)) ||
              (active_edge == 2 && (updated_transformations[level] == 3 || updated_transformations[level] == 6)) ||
              (active_edge == 3 && (updated_transformations[level] == 0 || updated_transformations[level] == 4)))
              neighbor_transforms->transf[neighbor_transforms->num_levels++] = (neighbor_edge.orientation ? neighbor_edge.local_num_of_edge : (neighbor_edge.local_num_of_edge + 1) % 4);
        }

        central_transformations.get(i)->strip_initial_transformations(j);
      }
    }

    template<typename Scalar>
    void NeighborSearch<Scalar>::Transformations::strip_initial_transformations(unsigned int number_of_stripped)
    {
      // Create a new clear array for the remaining transformations.
      unsigned int shifted_trfs[max_level];
      memset(shifted_trfs, 0, max_level * sizeof(unsigned int));
      // Move the old one to the new one.
      for(unsigned int k = number_of_stripped; k < num_levels; k++)
        shifted_trfs[k - number_of_stripped] = transf[k];
      // Point to the new one.
      memcpy(transf, shifted_trfs, max_level*sizeof(unsigned int));
      // We also have to store the information about length of the transformation array for this neighbor.
      num_levels -= number_of_stripped;
    }

    template<typename Scalar>
    void NeighborSearch<Scalar>::delete_neighbor(unsigned int position)
    {
      for(unsigned int i = position; i < n_neighbors - 1; i++)
        central_transformations.get(i)->copy_from(central_transformations.get(i + 1));

      if(central_transformations.present(n_neighbors - 1)) // may not be true when position == n_neighbors - 1
        central_transformations.get(n_neighbors - 1)->reset();

      for(unsigned int i = position; i < n_neighbors - 1; i++)
      {
        if(neighbor_transformations.present(i + 1))
        {
          if(!neighbor_transformations.present(i))
            neighbor_transformations.add(new Transformations, i);

          neighbor_transformations.get(i)->copy_from(neighbor_transformations.get(i + 1));
        }
      }
      if(neighbor_transformations.present(n_neighbors - 1)) // may not be true when position == n_neighbors - 1
        neighbor_transformations.get(n_neighbors - 1)->reset();

      neighbor_edges.erase (neighbor_edges.begin() + position);
      neighbors.erase (neighbors.begin() + position);
      n_neighbors--;
    }

    template<typename Scalar>
    void NeighborSearch<Scalar>::find_act_elem_up( Element* elem, int* orig_vertex_id, Node** par_mid_vertices, int n_parents)
    {
      Node* edge = NULL;
      Node* vertex = NULL;

      assert(n_parents <= (int)Transformations::max_level);

      // IDs of vertices bounding the current intermediate parent edge.
      int p1 = elem->vn[active_edge]->id;
      int p2 = elem->vn[(active_edge + 1) % elem->get_num_surf()]->id;

      int id_of_par_orient_1 = p1;
      int id_of_par_orient_2 = p2;

      // Find if p1 and p2 bound a used edge (used by the neighbor element).
      edge = mesh->peek_edge_node(p1, p2);

      // Add the vertex in the middle of the parent edge to the array of intermediate parent vertices. This is for
      // consequent transformation of functions on neighbor element.
      vertex = mesh->peek_vertex_node(p1, p2);
      if(vertex != NULL)
      {
        if(n_parents == 0)
          par_mid_vertices[n_parents++] = vertex;
        else
          if(n_parents == Transformations::max_level - 1)
            throw Hermes::Exceptions::Exception("Maximum number of intermediate parents exceeded in NeighborSearch<Scalar>::finding_act_elem_up");
          else
            if(par_mid_vertices[n_parents - 1]->id != vertex->id)
              par_mid_vertices[n_parents++] = vertex;
      }

      if((edge == NULL) || (central_el->en[active_edge]->id == edge->id))
      {
        // We have not yet found the parent of the central element completely adjacent to the neighbor.
        find_act_elem_up(elem->parent, orig_vertex_id, par_mid_vertices, n_parents);
      }
      else
      {
        for (int i = 0; i < 2; i++)
        {
          // Get a pointer to the active neighbor element.
          if((edge->elem[i] != NULL) && (edge->elem[i]->active == 1))
          {
            neighb_el = edge->elem[i];  //debug_log("way up neighbor: %d", neighb_el->id);

            // Get local number of the edge used by the neighbor.
            neighbor_edge.local_num_of_edge = -1;
            for(unsigned int j = 0; j < neighb_el->get_num_surf(); j++)
              if(neighb_el->en[j] == edge)
              {
                neighbor_edge.local_num_of_edge = j;
                break;
              }
              if(neighbor_edge.local_num_of_edge == -1) throw Hermes::Exceptions::Exception("Neighbor edge wasn't found");

              Node* n = NULL;

              // Add to the array of neighbor_transformations one that transforms central el. to its parent completely
              // adjacent to the single big neighbor.
              assert(n_neighbors == 0);

              neighbor_transformations.add(new Transformations, n_neighbors);
              Transformations *neighbor_transforms = neighbor_transformations.get(n_neighbors);

              neighbor_transforms->num_levels = n_parents;

              // Go back through the intermediate inactive parents down to the central element and stack corresponding
              // neighbor_transformations into the array 'neighbor_transformations'.
              for(int j = n_parents - 1; j > 0; j-- )
              {
                n = mesh->peek_vertex_node(par_mid_vertices[j]->id, p1);
                if(n == NULL)
                {
                  neighbor_transforms->transf[n_parents - j - 1] = neighbor_edge.local_num_of_edge;
                  p1 = par_mid_vertices[j]->id;
                }
                else
                {
                  if(n->id == par_mid_vertices[j-1]->id)
                  {
                    neighbor_transforms->transf[n_parents - j - 1] = (neighbor_edge.local_num_of_edge + 1) % neighb_el->get_num_surf();
                    p2 = par_mid_vertices[j]->id;
                  }
                  else
                  {
                    neighbor_transforms->transf[n_parents - j - 1] = neighbor_edge.local_num_of_edge;
                    p1 = par_mid_vertices[j]->id;
                  }
                }
              }

              // Final transformation to the central element itself.
              if(orig_vertex_id[0] == par_mid_vertices[0]->id)
                neighbor_transforms->transf[n_parents - 1] = neighbor_edge.local_num_of_edge;
              else
                neighbor_transforms->transf[n_parents - 1] = (neighbor_edge.local_num_of_edge + 1) % neighb_el->get_num_surf();

              NeighborEdgeInfo local_edge_info;
              local_edge_info.local_num_of_edge = neighbor_edge.local_num_of_edge;

              // Query the orientation of the neighbor edge relative to the central el.
              local_edge_info.orientation = neighbor_edge_orientation(id_of_par_orient_1, id_of_par_orient_2, 0);

              neighbor_edges.push_back(local_edge_info);

              // There is only one neighbor, ...
              n_neighbors = 1;

              // ...add it to the vector of neighbors.
              neighbors.push_back(neighb_el);
          }
        }
      }
    }

    template<typename Scalar>
    void NeighborSearch<Scalar>::find_act_elem_down( Node* vertex, int* bounding_verts_id, int* sons, unsigned int n_sons)
    {
      int mid_vert = vertex->id; // ID of vertex in between vertices from par_vertex_id.
      int bnd_verts[2];
      bnd_verts[0] = bounding_verts_id[0];
      bnd_verts[1] = bounding_verts_id[1];

      assert(n_sons < (unsigned) Transformations::max_level);

      for (int i = 0; i < 2; i++)
      {
        sons[n_sons-1] = (active_edge + i) % central_el->get_num_surf();

        // Try to get a pointer to the edge between the middle vertex and one of the vertices bounding the previously
        // tested segment.
        Node* edge = mesh->peek_edge_node(mid_vert, bnd_verts[i]);

        if(edge == NULL) // The edge is not used, i.e. there is no active element on either side.
        {
          // Get the middle vertex of this edge and try again on the segments into which this vertex splits the edge.
          Node * n = mesh->peek_vertex_node(mid_vert, bnd_verts[i]);
          if(n == NULL)
            throw Hermes::Exceptions::Exception("wasn't able to find middle vertex");
          else
          {
            // Make sure the next visited segment has the same orientation as the original central element's active edge.
            if(i == 0)
              bounding_verts_id[1] = mid_vert;
            else
              bounding_verts_id[0] = mid_vert;

            find_act_elem_down( n, bounding_verts_id, sons, n_sons + 1);

            bounding_verts_id[0] = bnd_verts[0];
            bounding_verts_id[1] = bnd_verts[1];
          }
        }
        else  // We have found a used edge, the active neighbor we are looking for is on one of its sides.
        {
          for (int j = 0; j < 2; j++)
          {
            if((edge->elem[j] != NULL) && (edge->elem[j]->active == 1))
            {
              neighb_el = mesh->get_element(edge->elem[j]->id);  //debug_log("way down neighbor: %d", edge->elem[j]->id);

              // Get local number of the edge used by the neighbor.
              neighbor_edge.local_num_of_edge = -1;
              for(unsigned int k = 0; k < neighb_el->get_num_surf(); k++)
                if(neighb_el->en[k] == edge)
                {
                  neighbor_edge.local_num_of_edge = k;
                  break;
                }

                if(neighbor_edge.local_num_of_edge == -1) throw Hermes::Exceptions::Exception("Neighbor edge wasn't found");

                //assert(!central_transformations.present(n_neighbors));

                central_transformations.add(new Transformations, n_neighbors);
                Transformations *tr = central_transformations.get(n_neighbors);

                // Construct the transformation path to the current neighbor.
                for(unsigned int k = 0; k < n_sons; k++)
                  tr->transf[k] = sons[k];
                tr->num_levels = n_sons;

                NeighborEdgeInfo local_edge_info;
                local_edge_info.local_num_of_edge = neighbor_edge.local_num_of_edge;
                // Query the orientation of the neighbor edge relative to the central el.
                local_edge_info.orientation = neighbor_edge_orientation(bnd_verts[0], bnd_verts[1], i);

                neighbor_edges.push_back(local_edge_info);

                // Append the new neighbor.
                n_neighbors++;
                neighbors.push_back(neighb_el);
            }
          }
        }
      }
    }

    template<typename Scalar>
    int NeighborSearch<Scalar>::neighbor_edge_orientation(int bounding_vert1, int bounding_vert2, int segment) const
    {
      if(segment == 0)
      {
        // neighbor edge goes from parent1 to middle vertex
        if(neighb_el->vn[neighbor_edge.local_num_of_edge]->id != bounding_vert1)
          return 1; // orientation reversed
      }
      else
      {
        // neighbor edge goes from middle vertex to parent2
        if(neighb_el->vn[neighbor_edge.local_num_of_edge]->id == bounding_vert2)
          return 1; // orientation reversed
      }
      return 0;
    }

    template<typename Scalar>
    NeighborSearch<Scalar>::ExtendedShapeset::ExtendedShapeset(const ExtendedShapeset & other)
    {
      this->central_al = new AsmList<Scalar>(*other.central_al);
      this->cnt = other.cnt;
      this->dof = other.dof;
      this->neighbor_al = new AsmList<Scalar>(*other.neighbor_al);
      this->combine_assembly_lists();
    }

    template<typename Scalar>
    typename NeighborSearch<Scalar>::ExtendedShapeset* NeighborSearch<Scalar>::create_extended_asmlist(const Space<Scalar>*space, AsmList<Scalar>* al)
    {
      ExtendedShapeset* new_supp_shapes = new ExtendedShapeset(this, al, space);

      return new_supp_shapes;
    }

    template<typename Scalar>
    typename NeighborSearch<Scalar>::ExtendedShapeset* NeighborSearch<Scalar>::create_extended_asmlist_multicomponent(const Space<Scalar> *space, AsmList<Scalar>* al)
    {
      if(supported_shapes != NULL)
        delete supported_shapes;

      supported_shapes = new ExtendedShapeset(this, al, space);

      return new ExtendedShapeset(*supported_shapes);
    }

    template<typename Scalar>
    void NeighborSearch<Scalar>::set_quad_order(int order)
    {
      neighb_quad_order = quad->get_edge_points(neighbor_edge.local_num_of_edge, order, neighbors[active_segment]->get_mode());
      central_quad_order = quad->get_edge_points(active_edge, order, central_el->get_mode());
    }

    template<typename Scalar>
    int NeighborSearch<Scalar>::get_quad_eo(bool on_neighbor) const
    {
      if(on_neighbor)
        return neighb_quad_order;
      else
        return central_quad_order;
    }

    template<typename Scalar>
    int NeighborSearch<Scalar>::get_active_segment() const
    {
      return this->active_segment;
    }

    template<typename Scalar>
    void NeighborSearch<Scalar>::set_active_segment(unsigned int index)
    {
      if(index >= n_neighbors)
        throw Hermes::Exceptions::Exception("NeighborSearch<Scalar>::set_active_segment() called with an incorrect index.");

      this->active_segment = index;
      this->neighb_el = this->neighbors[index];
      this->neighbor_edge = this->neighbor_edges[index];
    }

    template<typename Scalar>
    Element* NeighborSearch<Scalar>::get_neighb_el() const
    {
      return this->neighb_el;
    }

    template<typename Scalar>
    typename NeighborSearch<Scalar>::NeighborEdgeInfo NeighborSearch<Scalar>::get_neighbor_edge() const
    {
      return this->neighbor_edge;
    }

    template<typename Scalar>
    unsigned int NeighborSearch<Scalar>::get_central_n_trans(unsigned int index) const
    {
      if(this->central_transformations.present(index))
        return this->central_transformations.get(index)->num_levels;
      else
        return 0;
    }

    template<typename Scalar>
    unsigned int NeighborSearch<Scalar>::get_central_transformations(unsigned int index_1, unsigned int index_2) const
    {
      if(!this->central_transformations.present(index_1))
        throw Hermes::Exceptions::Exception("Out of bounds of central_transformations.");
      if(index_2 >= (unsigned) Transformations::max_level)
        throw Hermes::Exceptions::Exception("Trying to access transformation deeper than allowed.");

      return this->central_transformations.get(index_1)->transf[index_2];
    }

    template<typename Scalar>
    unsigned int NeighborSearch<Scalar>::get_neighbor_n_trans(unsigned int index) const
    {
      if(this->neighbor_transformations.present(index))
        return this->neighbor_transformations.get(index)->num_levels;
      else
        return 0;
    }

    template<typename Scalar>
    unsigned int NeighborSearch<Scalar>::get_neighbor_transformations(unsigned int index_1, unsigned int index_2) const
    {
      if(!this->neighbor_transformations.present(index_1))
        throw Hermes::Exceptions::Exception("Out of bounds of neighbor_transformations.");
      if(index_2 >= (unsigned) Transformations::max_level)
        throw Hermes::Exceptions::Exception("Trying to access transformation deeper than allowed.");

      return this->neighbor_transformations.get(index_1)->transf[index_2];
    }

    template<typename Scalar>
    DiscontinuousFunc<Scalar>* NeighborSearch<Scalar>::init_ext_fn(MeshFunction<Scalar>* fu)
    {
      Func<Scalar>* fn_central = init_fn(fu, get_quad_eo(false));

      uint64_t original_transform = fu->get_transform();

      // Change the active element of the function. Note that this also resets the transformations on the function.
      fu->set_active_element(neighbors[active_segment]);

      if(neighbor_transformations.present(active_segment))
        neighbor_transformations.get(active_segment)->apply_on(fu);

      Func<Scalar>* fn_neighbor = init_fn(fu, get_quad_eo(true));

      // Restore the original function.
      fu->set_active_element(central_el);
      fu->set_transform(original_transform);

      return new DiscontinuousFunc<Scalar>(fn_central, fn_neighbor, (neighbor_edge.orientation == 1));

      //NOTE: This function is not very efficient, since it sets the active elements and possibly pushes transformations
      // for each mesh function in each cycle of the innermost assembly loop. This is neccessary because only in
      // this innermost cycle (in function DiscreteProblem::eval_form), we know the quadrature order (dependent on
      // the actual basis and test function), which is needed for creating the Func<Scalar> objects via init_fn.
      // The reason for storing the central and neighbor values of any given function in these objects is that otherwise
      // we would have to have one independent copy of the function for each of the neighboring elements. However, it
      // could unify the way PrecalcShapesets and MeshFunctions are treated in NeighborSearch and maybe these additional
      // deep memory copying, performed only after setting the active edge part (before the nested loops over basis and
      // test functions), would be actually more efficient than this. This would require implementing copy for Filters.
    }

    template<typename Scalar>
    NeighborSearch<Scalar>::ExtendedShapeset::ExtendedShapeset(NeighborSearch* neighborhood, AsmList<Scalar>* central_al, const Space<Scalar>* space) :
    central_al(central_al)
    {
      neighbor_al = new AsmList<Scalar>();
      space->get_boundary_assembly_list(neighborhood->neighb_el, neighborhood->neighbor_edge.local_num_of_edge, neighbor_al);
      combine_assembly_lists();
    }

    template<typename Scalar>
    void NeighborSearch<Scalar>::ExtendedShapeset::combine_assembly_lists()
    {
      assert(central_al != NULL && neighbor_al != NULL);
      cnt = central_al->cnt + neighbor_al->cnt;
      dof = new int[cnt];
      memcpy(dof, central_al->dof, sizeof(int)*central_al->cnt);
      memcpy(dof + central_al->cnt, neighbor_al->dof, sizeof(int)*neighbor_al->cnt);
    }

    template<typename Scalar>
    NeighborSearch<Scalar>::ExtendedShapeset::~ExtendedShapeset()
    {
      delete [] dof;
      delete neighbor_al;
    }

    template<typename Scalar>
    void NeighborSearch<Scalar>::ExtendedShapeset::free_central_al()
    {
      delete central_al;
    }

    template<typename Scalar>
    void NeighborSearch<Scalar>::ExtendedShapeset::update(NeighborSearch* neighborhood, const Space<Scalar>* space)
    {
      delete [] this->dof;
      space->get_boundary_assembly_list(neighborhood->neighb_el, neighborhood->neighbor_edge.local_num_of_edge, neighbor_al);
      combine_assembly_lists();
    }

    template<typename Scalar>
    bool NeighborSearch<Scalar>::ExtendedShapeset::has_support_on_neighbor(unsigned int index) const
    {
      return (index >= central_al->cnt);
    }

    template<typename Scalar>
    NeighborSearch<Scalar>::Transformations::Transformations() : num_levels(0)
    {
      memset(transf, 0, max_level * sizeof(int));
    }

    template<typename Scalar>
    NeighborSearch<Scalar>::Transformations::Transformations(const Transformations* t)
    {
      copy_from(t);
    }

    template<typename Scalar>
    NeighborSearch<Scalar>::Transformations::Transformations(const Hermes::vector<unsigned int>& t)
    {
      copy_from(t);
    }

    template<typename Scalar>
    void NeighborSearch<Scalar>::Transformations::copy_from(const Hermes::vector<unsigned int>& t)
    {
      num_levels = std::min<unsigned int>(t.size(), max_level);
      std::copy( t.begin(), t.begin() + num_levels, transf);
    }

    template<typename Scalar>
    void NeighborSearch<Scalar>::Transformations::copy_from(const Transformations* t)
    {
      num_levels = t->num_levels;
      memcpy(transf, t->transf, max_level * sizeof(unsigned int));
    }

    template<typename Scalar>
    void NeighborSearch<Scalar>::Transformations::copy_to(Hermes::vector<unsigned int>* t)
    {
      t->assign(transf, transf + num_levels);
    }

    template<typename Scalar>
    void NeighborSearch<Scalar>::Transformations::reset()
    {
      memset(transf, 0, num_levels * sizeof(unsigned int));
      num_levels = 0;
    }

    template<typename Scalar>
    void NeighborSearch<Scalar>::Transformations::apply_on(Transformable* tr) const
    {
      for(unsigned int i = 0; i < num_levels; i++)
        tr->push_transform(transf[i]);
    }

    template<typename Scalar>
    void NeighborSearch<Scalar>::Transformations::apply_on(const Hermes::vector<Transformable*>& tr) const
    {
      for(Hermes::vector<Transformable*>::const_iterator it = tr.begin(); it != tr.end(); ++it)
        for(unsigned int i = 0; i < num_levels; i++)
          (*it)->push_transform(transf[i]);
    }

    template class HERMES_API NeighborSearch<double>;
    template class HERMES_API NeighborSearch<std::complex<double> >;
  }
}