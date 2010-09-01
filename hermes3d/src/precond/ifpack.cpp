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
#include "ifpack.h"
#include <common/callstack.h>
#include <common/timer.h>

#ifdef HAVE_IFPACK
#include <Ifpack_PointRelaxation.h>
#include <Ifpack_BlockRelaxation.h>
#include <Ifpack_DenseContainer.h>
#include <Ifpack_AdditiveSchwarz.h>
#include <Ifpack_ILU.h>
#include <Ifpack_ILUT.h>
#include <Ifpack_IC.h>
#include <Ifpack_ICT.h>
#include <Ifpack_Graph_Epetra_CrsGraph.h>
#endif

#define IFPACK_NOT_COMPILED "hermes3d was not built with Ifpack support."

IfpackPrecond::IfpackPrecond(const char *cls, const char *type)
{
	_F_
#ifdef HAVE_IFPACK
	this->time = 0.0;
	this->prec = NULL;
	this->owner = true;
	this->mat = NULL;

	this->cls = cls;
	this->type = type;
#else
	error(IFPACK_NOT_COMPILED);
#endif
}

IfpackPrecond::IfpackPrecond(const char *cls, const char *type, int overlap)
{
	_F_
#ifdef HAVE_IFPACK
	this->time = 0.0;
	this->prec = NULL;
	this->owner = true;
	this->mat = NULL;

	this->cls = cls;
	this->type = type;
	this->overlap = overlap;
#else
	error(IFPACK_NOT_COMPILED);
#endif
}

#ifdef HAVE_IFPACK
IfpackPrecond::IfpackPrecond(Ifpack_Preconditioner *ipc)
{
	_F_
	this->time = 0.0;
	this->prec = ipc;
	this->owner = false;
	this->mat = NULL;		// FIXME: take the matrix from ipc
}
#endif

IfpackPrecond::~IfpackPrecond()
{
	_F_
#ifdef HAVE_IFPACK
	if (owner) delete prec;
#endif
}

void IfpackPrecond::set_param(const char *name, const char *value)
{
	_F_
#ifdef HAVE_IFPACK
	ilist.set(name, value);
#endif
}

void IfpackPrecond::set_param(const char *name, int value)
{
	_F_
#ifdef HAVE_IFPACK
	ilist.set(name, value);
#endif
}

void IfpackPrecond::set_param(const char *name, double value)
{
	_F_
#ifdef HAVE_IFPACK
	ilist.set(name, value);
#endif
}

void IfpackPrecond::create(Matrix *m)
{
	_F_
#ifdef HAVE_IFPACK
	EpetraMatrix *mt = dynamic_cast<EpetraMatrix *>(m);
	assert(mt != NULL);
	mat = mt;
	if (strcmp(cls, "point-relax") == 0) {
		create_point_relax(mat, type);
		apply_params();
		initialize();
	}
	else if (strcmp(cls, "block-relax") == 0) {
		create_block_relax(mat, type);
		apply_params();
	}
	else if (strcmp(cls, "add-schwartz") == 0) {
		create_add_schwartz(mat, type, overlap);
		apply_params();
		initialize();
	}
#endif
}

void IfpackPrecond::create_point_relax(EpetraMatrix *a, const char *name)
{
	_F_
#ifdef HAVE_IFPACK
	prec = new Ifpack_PointRelaxation(a->mat);
	MEM_CHECK(prec);
	ilist.set("relaxation: type", name);
#endif
}

void IfpackPrecond::create_block_relax(EpetraMatrix *a, const char *name)
{
	_F_
#ifdef HAVE_IFPACK
	Teuchos::RCP<const Epetra_CrsGraph> rgraph = Teuchos::rcp(&(a->mat->Graph()));
	Ifpack_Graph *graph = new Ifpack_Graph_Epetra_CrsGraph(rgraph);
	MEM_CHECK(graph);

	Ifpack_Partitioner *partitioner = new Ifpack_GreedyPartitioner(graph);
	MEM_CHECK(partitioner);

	Teuchos::ParameterList list;
	list.set("partitioner: local parts", 1000);	// TODO: parametrize me
	partitioner->SetParameters(list);
	partitioner->Compute();

	prec = new Ifpack_BlockRelaxation<Ifpack_DenseContainer>(a->mat);
	MEM_CHECK(prec);
	ilist.set("relaxation: type", name);

	rgraph.release();
#endif
}

void IfpackPrecond::create_add_schwartz(EpetraMatrix *a, const char *name, int overlap)
{
	_F_
#ifdef HAVE_IFPACK
	if (strcasecmp(name, "ilu") == 0) {
		prec = new Ifpack_AdditiveSchwarz<Ifpack_ILU>(a->mat, overlap);
		MEM_CHECK(prec);
	}
	else if (strcasecmp(name, "ilut") == 0) {
		prec = new Ifpack_AdditiveSchwarz<Ifpack_ILUT>(a->mat, overlap);
		MEM_CHECK(prec);
	}
	else if (strcasecmp(name, "ic") == 0) {
		prec = new Ifpack_AdditiveSchwarz<Ifpack_IC>(a->mat, overlap);
		MEM_CHECK(prec);
	}
	else if (strcasecmp(name, "ict") == 0) {
		prec = new Ifpack_AdditiveSchwarz<Ifpack_ICT>(a->mat, overlap);
		MEM_CHECK(prec);
	}
	else
		prec = NULL;
#endif
}

int IfpackPrecond::initialize()
{
	_F_
#ifdef HAVE_IFPACK
	assert(prec != NULL);
	return prec->Initialize();
#else
	return 0;
#endif
}

void IfpackPrecond::compute()
{
	_F_
#ifdef HAVE_IFPACK
	assert(prec != NULL);
	Timer tmr;
	tmr.start();
	prec->Compute();
	tmr.stop();
	time = tmr.get_seconds();
#endif
}

void IfpackPrecond::apply_params()
{
	_F_
#ifdef HAVE_IFPACK
	prec->SetParameters(ilist);
#endif
}

#ifdef HAVE_IFPACK

int IfpackPrecond::ApplyInverse(const Epetra_MultiVector &r, Epetra_MultiVector &z) const
{
	assert(prec != NULL);
	return prec->ApplyInverse(r, z);
}


const Epetra_Comm &IfpackPrecond::Comm() const
{
	return mat->mat->Comm();
}

const Epetra_Map &IfpackPrecond::OperatorDomainMap() const
{
	return mat->mat->OperatorDomainMap();
}

const Epetra_Map &IfpackPrecond::OperatorRangeMap() const
{
	return mat->mat->OperatorRangeMap();
}
#endif
