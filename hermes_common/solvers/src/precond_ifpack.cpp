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

#include "precond_ifpack.h"

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

template<typename Scalar>
IfpackPrecond<Scalar>::IfpackPrecond(const char *cls, const char *type)
{
#ifdef HAVE_IFPACK
  this->prec = NULL;
  this->owner = true;
  this->mat = NULL;

  this->cls = cls;
  this->type = type;
#else
  error(IFPACK_NOT_COMPILED);
#endif
}

template<typename Scalar>
IfpackPrecond<Scalar>::IfpackPrecond(const char *cls, const char *type, int overlap)
{
#ifdef HAVE_IFPACK
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
template<typename Scalar>
IfpackPrecond<Scalar>::IfpackPrecond(Ifpack_Preconditioner *ipc)
{
  this->prec = ipc;
  this->owner = false;
  this->mat = NULL;		// FIXME: take the matrix from ipc
}
#endif

template<typename Scalar>
IfpackPrecond<Scalar>::~IfpackPrecond()
{

#ifdef HAVE_IFPACK
  if (owner) delete prec;
#endif
}

template<typename Scalar>
void IfpackPrecond<Scalar>::set_param(const char *name, const char *value)
{

#ifdef HAVE_IFPACK
  ilist.set(name, value);
#endif
}

template<typename Scalar>
void IfpackPrecond<Scalar>::set_param(const char *name, int value)
{

#ifdef HAVE_IFPACK
  ilist.set(name, value);
#endif
}

template<typename Scalar>
void IfpackPrecond<Scalar>::set_param(const char *name, double value)
{

#ifdef HAVE_IFPACK
  ilist.set(name, value);
#endif
}

template<typename Scalar>
void IfpackPrecond<Scalar>::create(Matrix<Scalar> *m)
{

#ifdef HAVE_IFPACK
  EpetraMatrix<Scalar> *mt = dynamic_cast<EpetraMatrix<Scalar> *>(m);
  assert(mt != NULL);
  mat = mt;
  if (strcmp(cls, "point-relax") == 0) 
  {
    create_point_relax(mat, type);
    apply_params();
    initialize();
  }
  else if (strcmp(cls, "block-relax") == 0) 
  {
    create_block_relax(mat, type);
    apply_params();
  }
  else if (strcmp(cls, "add-schwartz") == 0) 
  {
    create_add_schwartz(mat, type, overlap);
    apply_params();
    initialize();
  }
#endif
}

template<typename Scalar>
void IfpackPrecond<Scalar>::create_point_relax(EpetraMatrix<Scalar> *a, const char *name)
{

#ifdef HAVE_IFPACK
  prec = new Ifpack_PointRelaxation(a->mat);
  ilist.set("relaxation: type", name);
#endif
}

template<typename Scalar>
void IfpackPrecond<Scalar>::create_block_relax(EpetraMatrix<Scalar> *a, const char *name)
{

#ifdef HAVE_IFPACK
  Teuchos::RCP<const Epetra_CrsGraph> rgraph = Teuchos::rcp(&(a->mat->Graph()));
  Ifpack_Graph *graph = new Ifpack_Graph_Epetra_CrsGraph(rgraph);

  Ifpack_Partitioner *partitioner = new Ifpack_GreedyPartitioner(graph);

  Teuchos::ParameterList list;
  list.set("partitioner: local parts", 1000);	// TODO: parametrize me
  partitioner->SetParameters(list);
  partitioner->Compute();

  prec = new Ifpack_BlockRelaxation<Ifpack_DenseContainer>(a->mat);
  ilist.set("relaxation: type", name);

  rgraph.release();
#endif
}

template<typename Scalar>
void IfpackPrecond<Scalar>::create_add_schwartz(EpetraMatrix<Scalar> *a, const char *name, int overlap)
{

#ifdef HAVE_IFPACK
  if (strcasecmp(name, "ilu") == 0) 
  {
    prec = new Ifpack_AdditiveSchwarz<Ifpack_ILU>(a->mat, overlap);
  }
  else if (strcasecmp(name, "ilut") == 0) 
  {
    prec = new Ifpack_AdditiveSchwarz<Ifpack_ILUT>(a->mat, overlap);
  }
  else if (strcasecmp(name, "ic") == 0) 
  {
    prec = new Ifpack_AdditiveSchwarz<Ifpack_IC>(a->mat, overlap);
  }
  else if (strcasecmp(name, "ict") == 0) 
  {
    prec = new Ifpack_AdditiveSchwarz<Ifpack_ICT>(a->mat, overlap);
  }
  else
    prec = NULL;
#endif
}

template<typename Scalar>
int IfpackPrecond<Scalar>::initialize()
{

#ifdef HAVE_IFPACK
  assert(prec != NULL);
  return prec->Initialize();
#else
  return 0;
#endif
}

template<typename Scalar>
void IfpackPrecond<Scalar>::compute()
{

#ifdef HAVE_IFPACK
  assert(prec != NULL);
  prec->Compute();
#endif
}

template<typename Scalar>
void IfpackPrecond<Scalar>::apply_params()
{

#ifdef HAVE_IFPACK
  prec->SetParameters(ilist);
#endif
}

#ifdef HAVE_IFPACK
template<typename Scalar>
int IfpackPrecond<Scalar>::ApplyInverse(const Epetra_MultiVector &r, Epetra_MultiVector &z) const
{
  assert(prec != NULL);
  return prec->ApplyInverse(r, z);
}

template<typename Scalar>
const Epetra_Comm &IfpackPrecond<Scalar>::Comm() const
{
  return mat->mat->Comm();
}

template<typename Scalar>
const Epetra_Map &IfpackPrecond<Scalar>::OperatorDomainMap() const
{
  return mat->mat->OperatorDomainMap();
}

template<typename Scalar>
const Epetra_Map &IfpackPrecond<Scalar>::OperatorRangeMap() const
{
  return mat->mat->OperatorRangeMap();
}
#endif
