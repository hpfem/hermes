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

#include "../h3dconfig.h"
#include "ml.h"
#include <common/callstack.h>
#include <common/error.h>
#include <common/timer.h>

#define H3D_ML_NOT_COMPILED "hermes3d was not built with ML support."

MlPrecond::MlPrecond(const char *type)
{
	_F_
#ifdef HAVE_ML
	this->time = 0.0;
	this->prec = NULL;
	this->owner = true;
	this->mat = NULL;

	if (strcmp(type, "sa") == 0) ML_Epetra::SetDefaults("SA", mlist);
	else if (strcmp(type, "dd") == 0) ML_Epetra::SetDefaults("DD", mlist);
#else
	error(H3D_ML_NOT_COMPILED);
#endif
}

#ifdef HAVE_ML
MlPrecond::MlPrecond(ML_Epetra::MultiLevelPreconditioner *mpc)
{
	_F_
	this->time = 0.0;
	this->prec = mpc;
	this->owner = false;
	this->mat = NULL;			// FIXME: get the matrix from mpc
}
#endif

MlPrecond::~MlPrecond()
{
	_F_
#ifdef HAVE_ML
	if (owner) delete prec;
#endif
}

void MlPrecond::set_param(const char *name, const char *value)
{
	_F_
#ifdef HAVE_ML
	mlist.set(name, value);
#endif
}

void MlPrecond::set_param(const char *name, int value)
{
	_F_
#ifdef HAVE_ML
	mlist.set(name, value);
#endif
}

void MlPrecond::set_param(const char *name, double value)
{
	_F_
#ifdef HAVE_ML
	mlist.set(name, value);
#endif
}

void MlPrecond::create(Matrix *m)
{
	_F_
#ifdef HAVE_ML
	EpetraMatrix *mt = dynamic_cast<EpetraMatrix *>(m);
	assert(mt != NULL);
	mat = mt;
	delete prec;
	prec = new ML_Epetra::MultiLevelPreconditioner(*mat->mat, mlist, false);
	MEM_CHECK(prec);
#endif
}

void MlPrecond::destroy()
{
	_F_
#ifdef HAVE_ML
	assert(prec != NULL);
	prec->DestroyPreconditioner();
#endif
}

void MlPrecond::compute()
{
	_F_
#ifdef HAVE_ML
	assert(prec != NULL);
	Timer tmr;
	tmr.start();
	prec->ComputePreconditioner();
	tmr.stop();
	time = tmr.get_seconds();
#endif
}

void MlPrecond::print_unused()
{
	_F_
#ifdef HAVE_ML
	assert(prec != NULL);
	prec->PrintUnused();
#endif
}

#ifdef HAVE_ML

int MlPrecond::ApplyInverse(const Epetra_MultiVector &r, Epetra_MultiVector &z) const
{
	assert(prec != NULL);
	return prec->ApplyInverse(r, z);
}

const Epetra_Comm &MlPrecond::Comm() const
{
	return mat->mat->Comm();
}

const Epetra_Map &MlPrecond::OperatorDomainMap() const
{
	return mat->mat->OperatorDomainMap();
}

const Epetra_Map &MlPrecond::OperatorRangeMap() const
{
	return mat->mat->OperatorRangeMap();
}

#endif
