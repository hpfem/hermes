// This file is part of HermesCommon
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes2D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
/*! \file dp_interface.h
\brief Interface for DiscreteProblem required by NoxProblemInterface.
*/
#ifndef DPINTERFACE_H
#define DPINTERFACE_H

#include "../matrix.h"
#include "../tables.h"

using namespace Hermes::Algebra;

namespace Hermes
{
  namespace Solvers
  {
    /// \brief Minimalistic DiscreteProblem interface required by NoxProblemInterface.
    template<typename Scalar>
    class DiscreteProblemInterface : public Hermes::Mixins::IntegrableWithGlobalOrder, public Hermes::Mixins::SettableComputationTime
    {
    public:
      /// Get the number of unknowns.
      virtual int get_num_dofs() const = 0;

      /// Get info about presence of a matrix.
      virtual bool is_matrix_free() const = 0;

      /// Assembling.
      /// General assembling procedure for nonlinear problems. coeff_vec is the
      /// previous Newton vector. If force_diagonal_block == true, then (zero) matrix
      /// antries are created in diagonal blocks even if corresponding matrix weak
      /// forms do not exist. This is useful if the matrix is later to be merged with
      /// a matrix that has nonzeros in these blocks. The Table serves for optional
      /// weighting of matrix blocks in systems. The parameter add_dir_lift decides
      /// whether Dirichlet lift will be added while coeff_vec is converted into
      /// Solutions.
      virtual void assemble(Scalar* coeff_vec, SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs = NULL,
        bool force_diagonal_blocks = false, Table* block_weights = NULL) = 0;

      /// Assembling.
      /// Without the matrix.
      virtual void assemble(Scalar* coeff_vec, Vector<Scalar>* rhs = NULL,
        bool force_diagonal_blocks = false, Table* block_weights = NULL) = 0;

    protected:
      /// Preassembling.
      /// Precalculate matrix sparse structure.
      /// If force_diagonal_block == true, then (zero) matrix
      /// antries are created in diagonal blocks even if corresponding matrix weak
      /// forms do not exist. This is useful if the matrix is later to be merged with
      /// a matrix that has nonzeros in these blocks. The Table serves for optional
      /// weighting of matrix blocks in systems.
      virtual void create_sparse_structure() = 0;
      virtual void create_sparse_structure(SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs = NULL) = 0;

      DiscreteProblemInterface();

      template<typename T> friend class DiscreteProblemNOX;
    };
  }
}
#endif