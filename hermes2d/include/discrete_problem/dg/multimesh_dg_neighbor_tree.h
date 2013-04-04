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

namespace Hermes
{
  namespace Hermes2D
  {
    /// @ingroup inner
    /// Multimesh neighbors traversal class.
    /// Internal.
    class MultimeshDGNeighborTree
    {
    private:
      MultimeshDGNeighborTree(MultimeshDGNeighborTree* parent, unsigned int transformation);
      ~MultimeshDGNeighborTree();
      void set_left_son(MultimeshDGNeighborTree* left_son);
      void set_right_son(MultimeshDGNeighborTree* right_son);
      void set_transformation(unsigned int transformation);
      MultimeshDGNeighborTree* get_left_son();
      MultimeshDGNeighborTree* get_right_son();
      unsigned int get_transformation();
      MultimeshDGNeighborTree* parent;
      MultimeshDGNeighborTree* left_son;
      MultimeshDGNeighborTree* right_son;
      unsigned int transformation;
      template<typename Scalar> friend class DiscreteProblem;
      template<typename Scalar> friend class KellyTypeAdapt;
    };
  }
}
#endif
