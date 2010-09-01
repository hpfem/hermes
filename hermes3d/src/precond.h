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

#ifndef _PRECOND_H_
#define _PRECOND_H_

#include "matrix.h"
#include "weakform.h"

#ifdef HAVE_EPETRA
#include <Epetra_Operator.h>
#endif

/// @defgroup preconds Preconditioners
///


/// Abstract class to define interface for preconditioners
///
/// @ingroup preconds
class Precond
#ifdef HAVE_EPETRA
	: public Epetra_Operator
#endif
{
public:
	virtual void create(Matrix *mat) = 0;
	virtual void destroy() = 0;
	virtual void compute() = 0;

	double get_time() { return time; }

#ifdef HAVE_EPETRA
	virtual Epetra_Operator *get_obj() = 0;

	// Epetra_Operator interface
	virtual int SetUseTranspose(bool UseTranspose) { return 0; }
	virtual int Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const { return 0; }
	virtual int ApplyInverse(const Epetra_MultiVector &r, Epetra_MultiVector &z) const { return 0; }
	virtual double NormInf() const { return 0.0; }
	virtual const char *Label() const { return NULL; }
	virtual bool UseTranspose() const { return false; }
	virtual bool HasNormInf() const { return false; }
	virtual const Epetra_Comm &Comm() const = 0;
	virtual const Epetra_Map &OperatorDomainMap() const = 0;
	virtual const Epetra_Map &OperatorRangeMap() const = 0;
#endif

protected:
	double time;
};

#endif /* _PRECOND_H_ */
