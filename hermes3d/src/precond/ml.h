// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes3D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes3D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes3D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#ifndef _PRECOND_ML_H_
#define _PRECOND_ML_H_

#include "../precond.h"
#include "../solver/epetra.h"
#ifdef HAVE_ML
#include <ml_MultiLevelPreconditioner.h>
#endif

/// Preconditioners built on ML
///
/// @ingroup preconds
class MlPrecond : public Precond {
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
