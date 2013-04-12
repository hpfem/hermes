/// This file is part of Hermes2D.
///
/// Hermes2D is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 2 of the License, or
/// (at your option) any later version.
///
/// Hermes2D is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY;without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with Hermes2D. If not, see <http:///www.gnu.org/licenses/>.

#ifndef __H2D_MULTIMESH_DG_NEIGHBOR_TREE_H
#define __H2D_MULTIMESH_DG_NEIGHBOR_TREE_H

#include "multimesh_dg_neighbor_tree_node.h"
#include "mesh/traverse.h"
#include "neighbor.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// @ingroup inner
    /// Multimesh neighbors traversal class.
    /// Internal.
    template<typename Scalar>
    class MultimeshDGNeighborTree
    {
    public:
      /// The main method, for the passed neighbor searches, it will process all multi-mesh neighbor consolidation.
      static void process_edge(NeighborSearch<Scalar>** neighbor_searches, int num_neighbor_searches, int& num_neighbors, bool*& processed);
    
    private:
      /// Initialize the tree for traversing multimesh neighbors.
      static void build_multimesh_tree(MultimeshDGNeighborTreeNode* root, NeighborSearch<Scalar>** neighbor_searches, int number);

      /// Recursive insertion function into the tree.
      static void insert_into_multimesh_tree(MultimeshDGNeighborTreeNode* node, unsigned int* transformations, unsigned int transformation_count);

      /// Return a global (unified list of central element transformations representing the neighbors on the union mesh.
      static Hermes::vector<Hermes::vector<unsigned int>*> get_multimesh_neighbors_transformations(MultimeshDGNeighborTreeNode* multimesh_tree);

      /// Traverse the multimesh tree. Used in the function get_multimesh_neighbors_transformations().
      static void traverse_multimesh_tree(MultimeshDGNeighborTreeNode* node, Hermes::vector<Hermes::vector<unsigned int>*>& running_transformations);

      /// Update the NeighborSearch according to the multimesh tree.
      static void update_neighbor_search(NeighborSearch<Scalar>* ns, MultimeshDGNeighborTreeNode* multimesh_tree);

      /// Finds a node in the multimesh tree that corresponds to the array transformations, with the length of transformation_count,
      /// starting to look for it in the MultimeshDGNeighborTreeNode node.
      static MultimeshDGNeighborTreeNode* find_node(unsigned int* transformations, unsigned int transformation_count, MultimeshDGNeighborTreeNode* node);

      /// Updates the NeighborSearch ns according to the subtree of MultimeshDGNeighborTreeNode node.
      /// Returns 0 if no neighbor was deleted, -1 otherwise.
      static int update_ns_subtree(NeighborSearch<Scalar>* ns, MultimeshDGNeighborTreeNode* node, unsigned int ith_neighbor);

      /// Traverse the multimesh subtree. Used in the function update_ns_subtree().
      static void traverse_multimesh_subtree(MultimeshDGNeighborTreeNode* node, Hermes::vector<Hermes::vector<unsigned int>*>& running_central_transformations,
        Hermes::vector<Hermes::vector<unsigned int>*>& running_neighbor_transformations, const typename NeighborSearch<Scalar>::NeighborEdgeInfo& edge_info, const int& active_edge, const int& mode);

      friend class DiscreteProblem<Scalar>;
      friend class DiscreteProblemDGAssembler<Scalar>;
      friend class KellyTypeAdapt<Scalar>;
    };
  }
}
#endif
