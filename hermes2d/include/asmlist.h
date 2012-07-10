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

#ifndef __H2D_ASMLIST_H
#define __H2D_ASMLIST_H

#include "global.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// AsmList is a simple container for the element assembly arrays idx, dof and coef.
    /// These arrays are filled by Space::get_element_assembly_list() and used by the
    /// assembling procedure and the Solution class. The arrays are allocated and deallocated
    /// automatically by the class. The class provides a list of triples (idx, dof, coef).
    /// The triples are flattened to separate arrays of length 'cnt'.
    ///
    /// @ingroup inner
    template<typename Scalar>
    class HERMES_API AsmList
    {
    public:
      /// Constructor.
      AsmList();

      /// Destructor.
      ~AsmList();

      int* get_idx();
      int* get_dof();
      unsigned int get_cnt();
    private:
      /// Copy constructor.
      AsmList(const AsmList<Scalar> & other);

      int* idx;      ///< array of shape function indices
      int* dof;      ///< array of basis function numbers (DOFs)
      Scalar* coef;  ///< array of coefficients
      unsigned int cnt;       ///< the number of items in the arrays idx, dof and coef
      unsigned int cap;       ///< internal

      /// Adds a record for one basis function (shape functions index, basis functions index, coefficient).
      void add_triplet(int i, int d, Scalar c);

      /// Internal. Enlarges the storage capacity.
      void enlarge();
      template<typename T> friend class DiscreteProblem;
      template<typename T> friend class DiscreteProblemLinear;
      template<typename T> friend class NeighborSearch;
      template<typename T> friend class Solution;
      template<typename T> friend class Space;
      template<typename T> friend class H1Space;
      template<typename T> friend class L2Space;
      template<typename T> friend class HcurlSpace;
      template<typename T> friend class HdivSpace;
    };
  }
}
#endif