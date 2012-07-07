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
/*! \file precond_ml.h
\brief ML (Trilinos package) preconditioners interface.
*/
#ifndef __HERMES_COMMON_PRECOND_ML_H_
#define __HERMES_COMMON_PRECOND_ML_H_
#include "../config.h"
#ifdef HAVE_ML
#include "precond.h"
#include "epetra.h"
#include <ml_MultiLevelPreconditioner.h>

namespace Hermes
{
  namespace Preconditioners
  {
    using namespace Hermes::Solvers;
    /// \brief Preconditioners built on ML.
    ///
    /// @ingroup preconds
    template <typename Scalar>
    class HERMES_API MlPrecond : public Precond<Scalar>
    {
    public:
      /// @param[in] type - type of the preconditioner[ sa | dd ]
      /// - sa = smooth aggregation
      /// - dd = domain decomposition
      MlPrecond(const char *type);
      /// Wrap ML object.
      MlPrecond(ML_Epetra::MultiLevelPreconditioner *mpc);
      virtual ~MlPrecond();

      void set_param(const char *name, const char *value);
      void set_param(const char *name, int value);
      void set_param(const char *name, double value);
    protected:
      virtual Epetra_Operator *get_obj() { return prec; }

      /// @param[in] a
      virtual void create(Matrix<Scalar> *mat);
      /// Destroy the preconditioner object.
      virtual void destroy();
      /// Compute the preconditioner.
      virtual void compute();

      void print_unused();

      // Epetra_Operator interface
      virtual int ApplyInverse(const Epetra_MultiVector &r, Epetra_MultiVector &z) const;
      virtual const Epetra_Comm &Comm() const;
      virtual const Epetra_Map &OperatorDomainMap() const;
      virtual const Epetra_Map &OperatorRangeMap() const;

      ML_Epetra::MultiLevelPreconditioner *prec;
      Teuchos::ParameterList mlist;
      EpetraMatrix<Scalar> *mat;
      unsigned owner:1;

      friend class AztecOOSolver<Scalar>;
    };
  }
}
#endif
#endif