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

#include "discrete_problem/dg/multimesh_dg_neighbor_tree.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    void MultimeshDGNeighborTree<Scalar>::process_edge(NeighborSearch<Scalar>** neighbor_searches, int num_neighbor_searches, int& num_neighbors, bool*& processed)
    {
      MultimeshDGNeighborTreeNode root(NULL, 0);

      build_multimesh_tree(&root, neighbor_searches, num_neighbor_searches);

      // Update all NeighborSearches according to the multimesh tree.
      // After this, all NeighborSearches in neighbor_searches should have the same count
      // of neighbors and proper set of transformations
      // for the central and the neighbor element(s) alike.
      // Also check that every NeighborSearch has the same number of neighbor elements.
      num_neighbors = 0;
      for(unsigned int i = 0; i < num_neighbor_searches; i++)
      {
        NeighborSearch<Scalar>* ns = neighbor_searches[i];
        update_neighbor_search(ns, &root);
        
        if(num_neighbors == 0)
          num_neighbors = ns->n_neighbors;
        
        if(ns->n_neighbors != num_neighbors)
          throw Hermes::Exceptions::Exception("Num_neighbors of different NeighborSearches not matching in assemble_one_state().");
      }

      processed = new bool[num_neighbors];

      for(unsigned int neighbor_i = 0; neighbor_i < num_neighbors; neighbor_i++)
      {
        // If the active segment has already been processed (when the neighbor element was assembled), it is skipped.
        // We test all neighbor searches, because in the case of intra-element edge, the neighboring (the same as central) element
        // will be marked as visited, even though the edge was not calculated.
        processed[neighbor_i] = true;
        for(unsigned int i = 0; i < num_neighbor_searches; i++)
        {
          if(!neighbor_searches[i]->neighbors.at(neighbor_i)->visited)
          {
            processed[neighbor_i] = false;
            break;
          }
        }
      }
    }

    template<typename Scalar>
    void MultimeshDGNeighborTree<Scalar>::build_multimesh_tree(MultimeshDGNeighborTreeNode* root, NeighborSearch<Scalar>** neighbor_searches, int number)
    {
      for(unsigned int i = 0; i < number; i++)
      {
        NeighborSearch<Scalar>* ns = neighbor_searches[i];
        if(ns->n_neighbors == 1 && (ns->central_transformations_size == 0 || ns->central_transformations[0]->num_levels == 0))
          continue;
        for(unsigned int j = 0; j < ns->n_neighbors; j++)
          if(ns->central_transformations[j])
            insert_into_multimesh_tree(root, ns->central_transformations[j]->transf, ns->central_transformations[j]->num_levels);
      }
    }

    template<typename Scalar>
    void MultimeshDGNeighborTree<Scalar>::insert_into_multimesh_tree(MultimeshDGNeighborTreeNode* node, unsigned int* transformations, unsigned int transformation_count)
    {
      // If we are already in the leaf.
      if(transformation_count == 0)
        return;
      // Both sons are null. We have to add a new Node. Let us do it for the left sone of node.
      if(node->get_left_son() == NULL && node->get_right_son() == NULL)
      {
        node->set_left_son(new MultimeshDGNeighborTreeNode(node, transformations[0]));
        insert_into_multimesh_tree(node->get_left_son(), transformations + 1, transformation_count - 1);
      }
      // At least the left son is not null (it is impossible only for the right one to be not null, because
      // the left one always gets into the tree first, as seen above).
      else
      {
        // The existing left son is the right one to continue through.
        if(node->get_left_son()->get_transformation() == transformations[0])
          insert_into_multimesh_tree(node->get_left_son(), transformations + 1, transformation_count - 1);
        // The right one also exists, check that it is the right one, or return an error.
        else if(node->get_right_son())
        {
          if(node->get_right_son()->get_transformation() == transformations[0])
            insert_into_multimesh_tree(node->get_right_son(), transformations + 1, transformation_count - 1);
          else
            throw Hermes::Exceptions::Exception("More than two possible sons in insert_into_multimesh_tree().");
        }
        // If the right one does not exist and the left one was not correct, create a right son and continue this way.
        else
        {
          node->set_right_son(new MultimeshDGNeighborTreeNode(node, transformations[0]));
          insert_into_multimesh_tree(node->get_right_son(), transformations + 1, transformation_count - 1);
        }
      }
    }

    template<typename Scalar>
    Hermes::vector<Hermes::vector<unsigned int>*> MultimeshDGNeighborTree<Scalar>::get_multimesh_neighbors_transformations(MultimeshDGNeighborTreeNode* multimesh_tree)
    {
      // Initialize the vector.
      Hermes::vector<Hermes::vector<unsigned int>*> running_transformations;
      // Prepare the first neighbor's vector.
      running_transformations.push_back(new Hermes::vector<unsigned int>);
      // Fill the vector.
      traverse_multimesh_tree(multimesh_tree, running_transformations);
      return running_transformations;
    }

    template<typename Scalar>
    void MultimeshDGNeighborTree<Scalar>::traverse_multimesh_tree(MultimeshDGNeighborTreeNode* node,
      Hermes::vector<Hermes::vector<unsigned int>*>& running_transformations)
    {
      // If we are in the root.
      if(node->get_transformation() == 0)
      {
        if(node->get_left_son())
          traverse_multimesh_tree(node->get_left_son(), running_transformations);
        if(node->get_right_son())
          traverse_multimesh_tree(node->get_right_son(), running_transformations);
        // Delete the vector prepared by the last accessed leaf.
        delete running_transformations.back();
        running_transformations.pop_back();
        return;
      }
      // If we are in a leaf.
      if(node->get_left_son() == NULL && node->get_right_son() == NULL)
      {
        // Create a vector for the new neighbor.
        Hermes::vector<unsigned int>* new_neighbor_transformations = new Hermes::vector<unsigned int>;
        // Copy there the whole path except for this leaf.
        for(unsigned int i = 0; i < running_transformations.back()->size(); i++)
          new_neighbor_transformations->push_back((*running_transformations.back())[i]);
        // Insert this leaf into the current running transformation, thus complete it.
        running_transformations.back()->push_back(node->get_transformation());
        // Make the new_neighbor_transformations the current running transformation.
        running_transformations.push_back(new_neighbor_transformations);
        return;
      }
      else
      {
        running_transformations.back()->push_back(node->get_transformation());
        if(node->get_left_son())
          traverse_multimesh_tree(node->get_left_son(), running_transformations);
        if(node->get_right_son())
          traverse_multimesh_tree(node->get_right_son(), running_transformations);
        running_transformations.back()->pop_back();
        return;
      }
      return;
    }

    template<typename Scalar>
    void MultimeshDGNeighborTree<Scalar>::update_neighbor_search(NeighborSearch<Scalar>* ns, MultimeshDGNeighborTreeNode* multimesh_tree)
    {
      // This has to be done, because we pass ns by reference and the number of neighbors is changing.
      unsigned int num_neighbors = ns->get_num_neighbors();

      for(int i = 0; i < num_neighbors; i++)
      {
        // Find the node corresponding to this neighbor in the tree.
        MultimeshDGNeighborTreeNode* node;
        if(ns->central_transformations[i])
          node = find_node(ns->central_transformations[i]->transf, ns->central_transformations[i]->num_levels, multimesh_tree);
        else
          node = multimesh_tree;

        // Update the NeighborSearch.
        int added = update_ns_subtree(ns, node, i);
        i -= added;
        num_neighbors -= added;
      }
    }

    template<typename Scalar>
    MultimeshDGNeighborTreeNode* MultimeshDGNeighborTree<Scalar>::find_node(unsigned int* transformations,
      unsigned int transformation_count,
      MultimeshDGNeighborTreeNode* node)
    {
      // If there are no transformations left.
      if(transformation_count == 0)
        return node;
      else
      {
        if(node->get_left_son())
        {
          if(node->get_left_son()->get_transformation() == transformations[0])
            return find_node(transformations + 1, transformation_count - 1, node->get_left_son());
        }
        if(node->get_right_son())
        {
          if(node->get_right_son()->get_transformation() == transformations[0])
            return find_node(transformations + 1, transformation_count - 1, node->get_right_son());
        }
      }
      // We always should be able to empty the transformations array.
      throw
        Hermes::Exceptions::Exception("Transformation of a central element not found in the multimesh tree.");
      return NULL;
    }

    template<typename Scalar>
    int MultimeshDGNeighborTree<Scalar>::update_ns_subtree(NeighborSearch<Scalar>* ns,
      MultimeshDGNeighborTreeNode* node, unsigned int ith_neighbor)
    {
      int current_count = ns->get_num_neighbors();

      // No subtree => no work.
      // Also check the assertion that if one son is null, then the other too.
      if(node->get_left_son() == NULL)
      {
        if(node->get_right_son())
          throw Hermes::Exceptions::Exception("Only one son (right) not null in MultimeshDGNeighborTree<Scalar>::update_ns_subtree.");
        return 0;
      }

      // Key part.
      // Begin with storing the info about the current neighbor.
      Element* neighbor = ns->neighbors[ith_neighbor];
      typename NeighborSearch<Scalar>::NeighborEdgeInfo edge_info = ns->neighbor_edges[ith_neighbor];

      // Initialize the vector for central transformations->
      Hermes::vector<Hermes::vector<unsigned int>*> running_central_transformations;
      // Prepare the first new neighbor's vector. Push back the current transformations (in case of GO_DOWN neighborhood).
      running_central_transformations.push_back(new Hermes::vector<unsigned int>);
      if(ns->central_transformations[ith_neighbor])
        ns->central_transformations[ith_neighbor]->copy_to(running_central_transformations.back());

      // Initialize the vector for neighbor transformations->
      Hermes::vector<Hermes::vector<unsigned int>*> running_neighbor_transformations;
      // Prepare the first new neighbor's vector. Push back the current transformations (in case of GO_UP/NO_TRF neighborhood).
      running_neighbor_transformations.push_back(new Hermes::vector<unsigned int>);
      if(ns->neighbor_transformations[ith_neighbor])
        ns->neighbor_transformations[ith_neighbor]->copy_to(running_neighbor_transformations.back());

      // Delete the current neighbor.
      ns->delete_neighbor(ith_neighbor);

      // Move down the subtree.
      if(node->get_left_son())
        traverse_multimesh_subtree(node->get_left_son(), running_central_transformations,
        running_neighbor_transformations, edge_info, ns->active_edge,
        ns->central_el->get_mode());
      if(node->get_right_son())
        traverse_multimesh_subtree(node->get_right_son(), running_central_transformations,
        running_neighbor_transformations, edge_info, ns->active_edge,
        ns->central_el->get_mode());

      // Delete the last neighbors' info (this is a dead end, caused by the function traverse_multimesh_subtree.
      delete running_central_transformations.back();
      running_central_transformations.pop_back();
      delete running_neighbor_transformations.back();
      running_neighbor_transformations.pop_back();

      // Insert new neighbors.
      for(unsigned int i = 0; i < running_central_transformations.size(); i++)
      {
        ns->neighbors.push_back(neighbor);
        ns->neighbor_edges.push_back(edge_info);

        if(!ns->central_transformations[ns->n_neighbors])
          ns->add_central_transformations(new typename NeighborSearch<Scalar>::Transformations, ns->n_neighbors);

        if(!ns->neighbor_transformations[ns->n_neighbors])
          ns->add_neighbor_transformations(new typename NeighborSearch<Scalar>::Transformations, ns->n_neighbors);

        ns->central_transformations[ns->n_neighbors]->copy_from(*running_central_transformations[i]);
        ns->neighbor_transformations[ns->n_neighbors]->copy_from(*running_neighbor_transformations[i]);

        ns->n_neighbors++;
      }

      for(unsigned int i = 0; i < running_central_transformations.size(); i++)
        delete running_central_transformations[i];
      for(unsigned int i = 0; i < running_neighbor_transformations.size(); i++)
        delete running_neighbor_transformations[i];

      // Return the number of neighbors added/deleted.
      return ns->get_num_neighbors() - current_count;
    }

    template<typename Scalar>
    void MultimeshDGNeighborTree<Scalar>::traverse_multimesh_subtree(MultimeshDGNeighborTreeNode* node,
      Hermes::vector<Hermes::vector<unsigned int>*>& running_central_transformations,
      Hermes::vector<Hermes::vector<unsigned int>*>& running_neighbor_transformations,
      const typename NeighborSearch<Scalar>::NeighborEdgeInfo& edge_info, const int& active_edge, const int& mode)
    {
      // If we are in a leaf.
      if(node->get_left_son() == NULL && node->get_right_son() == NULL)
      {
        // Create vectors for the new neighbor.
        Hermes::vector<unsigned int>* new_neighbor_central_transformations = new Hermes::vector<unsigned int>;
        Hermes::vector<unsigned int>* new_neighbor_neighbor_transformations = new Hermes::vector<unsigned int>;

        // Copy there the whole path except for this leaf.
        for(unsigned int i = 0; i < running_central_transformations.back()->size(); i++)
          new_neighbor_central_transformations->push_back((*running_central_transformations.back())[i]);
        for(unsigned int i = 0; i < running_neighbor_transformations.back()->size(); i++)
          new_neighbor_neighbor_transformations->push_back((*running_neighbor_transformations.back())[i]);

        // Insert this leaf into the current running central transformation, thus complete it.
        running_central_transformations.back()->push_back(node->get_transformation());

        // Make the new_neighbor_central_transformations the current running central transformation.
        running_central_transformations.push_back(new_neighbor_central_transformations);

        // Take care of the neighbor transformation.
        // Insert appropriate info from this leaf into the current running neighbor transformation, thus complete it.
        if(mode == HERMES_MODE_TRIANGLE)
          if((active_edge == 0 && node->get_transformation() == 0) ||
            (active_edge == 1 && node->get_transformation() == 1) ||
            (active_edge == 2 && node->get_transformation() == 2))
            running_neighbor_transformations.back()->push_back((!edge_info.orientation ? edge_info.local_num_of_edge : (edge_info.local_num_of_edge + 1) % 3));
          else
            running_neighbor_transformations.back()->push_back((edge_info.orientation ? edge_info.local_num_of_edge : (edge_info.local_num_of_edge + 1) % 3));
        // Quads.
        else
          if((active_edge == 0 && (node->get_transformation() == 0 || node->get_transformation() == 6)) ||
            (active_edge == 1 && (node->get_transformation() == 1 || node->get_transformation() == 4)) ||
            (active_edge == 2 && (node->get_transformation() == 2 || node->get_transformation() == 7)) ||
            (active_edge == 3 && (node->get_transformation() == 3 || node->get_transformation() == 5)))
            running_neighbor_transformations.back()->push_back((!edge_info.orientation ? edge_info.local_num_of_edge : (edge_info.local_num_of_edge + 1) % H2D_MAX_NUMBER_EDGES));
          else
            running_neighbor_transformations.back()->push_back((edge_info.orientation ? edge_info.local_num_of_edge : (edge_info.local_num_of_edge + 1) % H2D_MAX_NUMBER_EDGES));

        // Make the new_neighbor_neighbor_transformations the current running neighbor transformation.
        running_neighbor_transformations.push_back(new_neighbor_neighbor_transformations);

        return;
      }
      else
      {
        // Insert this leaf into the current running central transformation, thus complete it.
        running_central_transformations.back()->push_back(node->get_transformation());

        // Insert appropriate info from this leaf into the current running neighbor transformation, thus complete it.
        // Triangles.
        if(mode == HERMES_MODE_TRIANGLE)
          if((active_edge == 0 && node->get_transformation() == 0) ||
            (active_edge == 1 && node->get_transformation() == 1) ||
            (active_edge == 2 && node->get_transformation() == 2))
            running_neighbor_transformations.back()->push_back((!edge_info.orientation ? edge_info.local_num_of_edge : (edge_info.local_num_of_edge + 1) % 3));
          else
            running_neighbor_transformations.back()->push_back((edge_info.orientation ? edge_info.local_num_of_edge : (edge_info.local_num_of_edge + 1) % 3));
        // Quads.
        else
          if((active_edge == 0 && (node->get_transformation() == 0 || node->get_transformation() == 6)) ||
            (active_edge == 1 && (node->get_transformation() == 1 || node->get_transformation() == 4)) ||
            (active_edge == 2 && (node->get_transformation() == 2 || node->get_transformation() == 7)) ||
            (active_edge == 3 && (node->get_transformation() == 3 || node->get_transformation() == 5)))
            running_neighbor_transformations.back()->push_back((!edge_info.orientation ? edge_info.local_num_of_edge : (edge_info.local_num_of_edge + 1) % H2D_MAX_NUMBER_EDGES));
          else
            running_neighbor_transformations.back()->push_back((edge_info.orientation ? edge_info.local_num_of_edge : (edge_info.local_num_of_edge + 1) % H2D_MAX_NUMBER_EDGES));

        // Move down.
        if(node->get_left_son())
          traverse_multimesh_subtree(node->get_left_son(), running_central_transformations, running_neighbor_transformations,
          edge_info, active_edge, mode);
        if(node->get_right_son())
          traverse_multimesh_subtree(node->get_right_son(), running_central_transformations, running_neighbor_transformations,
          edge_info, active_edge, mode);

        // Take this transformation out.
        running_central_transformations.back()->pop_back();
        running_neighbor_transformations.back()->pop_back();
        return;
      }
      return;
    }

    template class HERMES_API MultimeshDGNeighborTree<double>;
    template class HERMES_API MultimeshDGNeighborTree<std::complex<double> >;
  }
}