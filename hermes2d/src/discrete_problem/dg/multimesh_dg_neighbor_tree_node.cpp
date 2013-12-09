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

#include "discrete_problem/dg/multimesh_dg_neighbor_tree_node.h"
#include "global.h"

namespace Hermes
{
  namespace Hermes2D
  {
    MultimeshDGNeighborTreeNode::MultimeshDGNeighborTreeNode(MultimeshDGNeighborTreeNode* parent, unsigned int transformation) : parent(parent), transformation(transformation)
    {
      left_son = right_son = NULL;
    }
    MultimeshDGNeighborTreeNode::~MultimeshDGNeighborTreeNode()
    {
      if(left_son)
      {
        delete left_son;
        left_son = NULL;
      }
      if(right_son)
      {
        delete right_son;
        right_son = NULL;
      }
    }
    void MultimeshDGNeighborTreeNode::set_left_son(MultimeshDGNeighborTreeNode* left_son)
    {
      this->left_son = left_son;
    }
    void MultimeshDGNeighborTreeNode::set_right_son(MultimeshDGNeighborTreeNode* right_son)
    {
      this->right_son = right_son;
    }
    void MultimeshDGNeighborTreeNode::set_transformation(unsigned int transformation)
    {
      this->transformation = transformation;
    }
    MultimeshDGNeighborTreeNode* MultimeshDGNeighborTreeNode::get_left_son()
    {
      return left_son;
    }
    MultimeshDGNeighborTreeNode* MultimeshDGNeighborTreeNode::get_right_son()
    {
      return right_son;
    }
    unsigned int MultimeshDGNeighborTreeNode::get_transformation()
    {
      return this->transformation;
    }
  }
}