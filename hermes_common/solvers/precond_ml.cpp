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

#include "precond_ml.h"

template<typename Scalar>
MlPrecond<Scalar>::MlPrecond(const char *type)
{
#ifdef HAVE_ML
  this->prec = NULL;
  this->owner = true;
  this->mat = NULL;

  if (strcmp(type, "sa") == 0) ML_Epetra::SetDefaults("SA", mlist);
  else if (strcmp(type, "dd") == 0) ML_Epetra::SetDefaults("DD", mlist);
#else
  error(ML_NOT_COMPILED);
#endif
}

#ifdef HAVE_ML

template<typename Scalar>
MlPrecond<Scalar>::MlPrecond(ML_Epetra::MultiLevelPreconditioner *mpc)
{
  this->prec = mpc;
  this->owner = false;
  this->mat = NULL;			// FIXME: get the matrix from mpc
}
#endif

template<typename Scalar>
MlPrecond<Scalar>::~MlPrecond()
{
#ifdef HAVE_ML
  if (owner) delete prec;
#endif
}

template<typename Scalar>
void MlPrecond<Scalar>::set_param(const char *name, const char *value)
{
#ifdef HAVE_ML
  mlist.set(name, value);
#endif
}

template<typename Scalar>
void MlPrecond<Scalar>::set_param(const char *name, int value)
{
#ifdef HAVE_ML
  mlist.set(name, value);
#endif
}

template<typename Scalar>
void MlPrecond<Scalar>::set_param(const char *name, double value)
{
#ifdef HAVE_ML
  mlist.set(name, value);
#endif
}

template<typename Scalar>
void MlPrecond<Scalar>::create(Matrix<Scalar> *m)
{
#ifdef HAVE_ML
  EpetraMatrix<Scalar> *mt = dynamic_cast<EpetraMatrix<Scalar> *>(m);
  assert(mt != NULL);
  mat = mt;
  delete prec;
  prec = new ML_Epetra::MultiLevelPreconditioner(*mat->mat, mlist, false);
#endif
}

template<typename Scalar>
void MlPrecond<Scalar>::destroy()
{
#ifdef HAVE_ML
  assert(prec != NULL);
  prec->DestroyPreconditioner();
#endif
}

template<typename Scalar>
void MlPrecond<Scalar>::compute()
{
#ifdef HAVE_ML
  assert(prec != NULL);
  prec->ComputePreconditioner();
#endif
}

template<typename Scalar>
void MlPrecond<Scalar>::print_unused()
{
#ifdef HAVE_ML
  assert(prec != NULL);
  prec->PrintUnused();
#endif
}

#ifdef HAVE_ML

template<typename Scalar>
int MlPrecond<Scalar>::ApplyInverse(const Epetra_MultiVector &r, Epetra_MultiVector &z) const
{
  assert(prec != NULL);
  return prec->ApplyInverse(r, z);
}

template<typename Scalar>
const Epetra_Comm &MlPrecond<Scalar>::Comm() const
{
  return mat->mat->Comm();
}

template<typename Scalar>
const Epetra_Map &MlPrecond<Scalar>::OperatorDomainMap() const
{
  return mat->mat->OperatorDomainMap();
}

template<typename Scalar>
const Epetra_Map &MlPrecond<Scalar>::OperatorRangeMap() const
{
  return mat->mat->OperatorRangeMap();
}

#endif
