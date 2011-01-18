// This file is part of Hermes2D
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


#ifndef __HERMES_COMMON_PRECOND_ML_H_
#define __HERMES_COMMON_PRECOND_ML_H_

#include "precond.h"
#include "epetra.h"
#ifdef HAVE_ML
  #include <ml_MultiLevelPreconditioner.h>
#endif

/// Preconditioners built on ML
///
/// @ingroup preconds
class HERMES_API MlPrecond : public Precond {
public:
  /// @param[in] type - type of the preconditioner [ sa | dd ]
  /// - sa = smooth aggregation
  /// - dd = domain decomposition
  MlPrecond(const char *type);
#ifdef HAVE_ML
  /// Wrap ML object
  MlPrecond(ML_Epetra::MultiLevelPreconditioner *mpc);
#endif
  virtual ~MlPrecond();

#ifdef HAVE_ML
  virtual Epetra_Operator *get_obj() { return prec; }
#endif

  /// @param[in] a
  virtual void create(Matrix *mat);
  /// Destroy the preconditioner object
  virtual void destroy();
  /// Compute the preconditioner
  virtual void compute();

  void set_param(const char *name, const char *value);
  void set_param(const char *name, int value);
  void set_param(const char *name, double value);

  void print_unused();

#ifdef HAVE_ML
  // Epetra_Operator interface
  virtual int ApplyInverse(const Epetra_MultiVector &r, Epetra_MultiVector &z) const;
  virtual const Epetra_Comm &Comm() const;
  virtual const Epetra_Map &OperatorDomainMap() const;
  virtual const Epetra_Map &OperatorRangeMap() const;
#endif

protected:
#ifdef HAVE_ML
  ML_Epetra::MultiLevelPreconditioner *prec;
  Teuchos::ParameterList mlist;
  EpetraMatrix *mat;
#endif
  unsigned owner:1;

  friend class AztecOOSolver;
};

#endif /* _PRECOND_ML_H_ */
