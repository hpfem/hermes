// This file is part of HermesCommon
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
/*! \file precond_ml.cpp
\brief ML (Trilinos package) preconditioners interface.
*/
#include "config.h"
#ifdef HAVE_ML
#include "precond_ml.h"

namespace Hermes
{
  namespace Preconditioners
  {
    template<typename Scalar>
    MlPrecond<Scalar>::MlPrecond(const char *type)
    {
      this->prec = NULL;
      this->owner = true;
      this->mat = NULL;

      if(strcmp(type, "sa") == 0) ML_Epetra::SetDefaults("SA", mlist);
      else if(strcmp(type, "dd") == 0) ML_Epetra::SetDefaults("DD", mlist);
    }

    template<typename Scalar>
    MlPrecond<Scalar>::MlPrecond(ML_Epetra::MultiLevelPreconditioner *mpc)
    {
      this->prec = mpc;
      this->owner = false;
      this->mat = NULL;      // FIXME: get the matrix from mpc
    }

    template<typename Scalar>
    MlPrecond<Scalar>::~MlPrecond()
    {
      if(owner) delete prec;
    }

    template<typename Scalar>
    void MlPrecond<Scalar>::set_param(const char *name, const char *value)
    {
      mlist.set(name, value);
    }

    template<typename Scalar>
    void MlPrecond<Scalar>::set_param(const char *name, int value)
    {
      mlist.set(name, value);
    }

    template<typename Scalar>
    void MlPrecond<Scalar>::set_param(const char *name, double value)
    {
      mlist.set(name, value);
    }

    template<typename Scalar>
    void MlPrecond<Scalar>::create(Matrix<Scalar> *m)
    {
      EpetraMatrix<Scalar> *mt = static_cast<EpetraMatrix<Scalar> *>(m);
      assert(mt != NULL);
      mat = mt;
      delete prec;
      prec = new ML_Epetra::MultiLevelPreconditioner(*mat->mat, mlist, false);
    }

    template<typename Scalar>
    void MlPrecond<Scalar>::destroy()
    {
      assert(prec != NULL);
      prec->DestroyPreconditioner();
    }

    template<typename Scalar>
    void MlPrecond<Scalar>::compute()
    {
      assert(prec != NULL);
      prec->ComputePreconditioner();
    }

    template<typename Scalar>
    void MlPrecond<Scalar>::print_unused()
    {
      assert(prec != NULL);
      prec->PrintUnused();
    }

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
    template class HERMES_API MlPrecond<double>;
    template class HERMES_API MlPrecond<std::complex<double> >;
  }
}
#endif