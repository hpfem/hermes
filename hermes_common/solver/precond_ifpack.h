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

#ifndef __HERMES_COMMON_PRECOND_IFPACK_H_
#define __HERMES_COMMON_PRECOND_IFPACK_H_

#include "precond.h"
#include "epetra.h"
#ifdef HAVE_IFPACK
  #include <Ifpack_Preconditioner.h>
#endif

/// Preconditioners built on IFPACK
///
/// @ingroup preconds
class HERMES_API IfpackPrecond: public Precond {
public:
  /// Constructor for relaxation methods
  /// @param[in] cls - class of the preconditioner [ point-relax | block-relax ]
  /// @param[in] name - the name of the relaxation type
  /// Possible values are: [ Jacobi | Gauss-Seidel | symmetric Gauss-Seidel ]
  IfpackPrecond(const char *cls, const char *type = "Jacobi");
  /// Constructor for domain decomposition methods
  /// @param[in] cls - class of the preconditioner [ add-schwartz ]
  /// @param[in] name: [ ic | ict | ilu | ilut ]
  /// @param[in] overlap - number
  IfpackPrecond(const char *cls, const char *type, int overlap);
#ifdef HAVE_IFPACK
  /// Wrap IFPACK object
  IfpackPrecond(Ifpack_Preconditioner *ipc);
#endif
  virtual ~IfpackPrecond();

#ifdef HAVE_IFPACK
  virtual Epetra_Operator *get_obj() { return prec; }
#endif

  virtual void create(Matrix *mat);
  virtual void destroy() { }
  virtual void compute();

  void set_param(const char *name, const char *value);
  void set_param(const char *name, int value);
  void set_param(const char *name, double value);

#ifdef HAVE_IFPACK
  // Epetra_Operator interface
  virtual int ApplyInverse(const Epetra_MultiVector &r, Epetra_MultiVector &z) const;
  virtual const Epetra_Comm &Comm() const;
  virtual const Epetra_Map &OperatorDomainMap() const;
  virtual const Epetra_Map &OperatorRangeMap() const;
#endif

protected:
  void create_point_relax(EpetraMatrix *a, const char *name);
  void create_block_relax(EpetraMatrix *a, const char *name);
  void create_add_schwartz(EpetraMatrix *a, const char *name, int overlap);
  int initialize();
  void apply_params();
#ifdef HAVE_IFPACK
  Ifpack_Preconditioner *prec;
  Teuchos::ParameterList ilist;
  EpetraMatrix *mat;
#endif
  unsigned owner:1;
  const char *cls;			// class of the preconditioner
  const char *type;
  int overlap;

  friend class AztecOOSolver;
};

#endif /* _PRECOND_IFPACK_H_ */
