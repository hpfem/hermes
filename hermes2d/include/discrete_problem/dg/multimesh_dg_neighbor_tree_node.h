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

#ifndef __H2D_MULTIMESH_DG_NEIGHBOR_TREE_NODE_H
#define __H2D_MULTIMESH_DG_NEIGHBOR_TREE_NODE_H

namespace Hermes
{
  namespace Hermes2D
  {
    /// @ingroup inner
    /// Multimesh neighbors traversal inner class.
    /// Completely internal.
    class MultimeshDGNeighborTreeNode
    {
    private:
      MultimeshDGNeighborTreeNode(MultimeshDGNeighborTreeNode* parent, unsigned int transformation);
      ~MultimeshDGNeighborTreeNode();
      void set_left_son(MultimeshDGNeighborTreeNode* left_son);
      void set_right_son(MultimeshDGNeighborTreeNode* right_son);
      void set_transformation(unsigned int transformation);
      MultimeshDGNeighborTreeNode* get_left_son();
      MultimeshDGNeighborTreeNode* get_right_son();
      unsigned int get_transformation();
      MultimeshDGNeighborTreeNode* parent;
      MultimeshDGNeighborTreeNode* left_son;
      MultimeshDGNeighborTreeNode* right_son;
      unsigned int transformation;

      template<typename Scalar> friend class MultimeshDGNeighborTree;
      template<typename Scalar> friend class KellyTypeAdapt;
    };
  }
}
#endif
