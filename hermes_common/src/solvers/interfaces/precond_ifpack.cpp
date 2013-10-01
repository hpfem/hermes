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
/*! \file precond_ifpack.cpp
\brief IFPACK (Trilinos package) preconditioners interface.
*/
#include "config.h"
#ifdef HAVE_IFPACK
#include "precond_ifpack.h"
#include "Ifpack_PointRelaxation.h"
#include "Ifpack_BlockRelaxation.h"
#include "Ifpack_DenseContainer.h"
#include "Ifpack_AdditiveSchwarz.h"
#include "Ifpack_ILU.h"
#include "Ifpack_ILUT.h"
#include "Ifpack_IC.h"
#include "Ifpack_ICT.h"
#include "Ifpack_Graph_Epetra_CrsGraph.h"

namespace Hermes
{
  namespace Preconditioners
{
    template<typename Scalar>
    IfpackPrecond<Scalar>::IfpackPrecond(const char *cls, const char *type)
    {
      this->prec = NULL;
      this->owner = true;
      this->mat = NULL;

      this->cls = cls;
      this->type = type;
    }

    template<typename Scalar>
    IfpackPrecond<Scalar>::IfpackPrecond(const char *cls, const char *type, int overlap)
    {
      this->prec = NULL;
      this->owner = true;
      this->mat = NULL;

      this->cls = cls;
      this->type = type;
      this->overlap = overlap;
    }

    template<typename Scalar>
    IfpackPrecond<Scalar>::IfpackPrecond(Ifpack_Preconditioner *ipc)
    {
      this->prec = ipc;
      this->owner = false;
      this->mat = NULL;    // FIXME: take the matrix from ipc
    }

    template<typename Scalar>
    IfpackPrecond<Scalar>::~IfpackPrecond()
    {
      if(owner) delete prec;
    }

    template<typename Scalar>
    void IfpackPrecond<Scalar>::set_param(const char *name, const char *value)
    {
      ilist.set(name, value);
    }

    template<typename Scalar>
    void IfpackPrecond<Scalar>::set_param(const char *name, int value)
    {
      ilist.set(name, value);
    }

    template<typename Scalar>
    void IfpackPrecond<Scalar>::set_param(const char *name, double value)
    {
      ilist.set(name, value);
    }

    template<typename Scalar>
    void IfpackPrecond<Scalar>::create(Matrix<Scalar> *m)
    {
      EpetraMatrix<Scalar> *mt = static_cast<EpetraMatrix<Scalar> *>(m);
      assert(mt != NULL);
      mat = mt;
      if(strcmp(cls, "point-relax") == 0)
      {
        create_point_relax(mat, type);
        apply_params();
        initialize();
      }
      else if(strcmp(cls, "block-relax") == 0)
      {
        create_block_relax(mat, type);
        apply_params();
      }
      else if(strcmp(cls, "add-schwartz") == 0)
      {
        create_add_schwartz(mat, type, overlap);
        apply_params();
        initialize();
      }
    }

    template<typename Scalar>
    void IfpackPrecond<Scalar>::create_point_relax(EpetraMatrix<Scalar> *a, const char *name)
    {
      prec = new Ifpack_PointRelaxation(a->mat);
      ilist.set("relaxation: type", name);
    }

    template<typename Scalar>
    void IfpackPrecond<Scalar>::create_block_relax(EpetraMatrix<Scalar> *a, const char *name)
    {
      Teuchos::RCP<const Epetra_CrsGraph> rgraph = Teuchos::rcp(&(a->mat->Graph()));
      Ifpack_Graph *graph = new Ifpack_Graph_Epetra_CrsGraph(rgraph);

      Ifpack_Partitioner *partitioner = new Ifpack_GreedyPartitioner(graph);

      Teuchos::ParameterList list;
      list.set("partitioner: local parts", 1000);  //\todo parametrize me
      partitioner->SetParameters(list);
      partitioner->Compute();

      prec = new Ifpack_BlockRelaxation<Ifpack_DenseContainer>(a->mat);
      ilist.set("relaxation: type", name);

      rgraph.release();
    }

    template<typename Scalar>
    void IfpackPrecond<Scalar>::create_add_schwartz(EpetraMatrix<Scalar> *a, const char *name, int overlap)
    {
      if(strcasecmp(name, "ilu") == 0)
      {
        prec = new Ifpack_AdditiveSchwarz<Ifpack_ILU>(a->mat, overlap);
      }
      else if(strcasecmp(name, "ilut") == 0)
      {
        prec = new Ifpack_AdditiveSchwarz<Ifpack_ILUT>(a->mat, overlap);
      }
      else if(strcasecmp(name, "ic") == 0)
      {
        prec = new Ifpack_AdditiveSchwarz<Ifpack_IC>(a->mat, overlap);
      }
      else if(strcasecmp(name, "ict") == 0)
      {
        prec = new Ifpack_AdditiveSchwarz<Ifpack_ICT>(a->mat, overlap);
      }
      else
        prec = NULL;
    }

    template<typename Scalar>
    int IfpackPrecond<Scalar>::initialize()
    {
      assert(prec != NULL);
      return prec->Initialize();
    }

    template<typename Scalar>
    void IfpackPrecond<Scalar>::compute()
    {
      assert(prec != NULL);
      prec->Compute();
    }

    template<typename Scalar>
    void IfpackPrecond<Scalar>::apply_params()
    {
      prec->SetParameters(ilist);
    }

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
    template class HERMES_API IfpackPrecond<double>;
    template class HERMES_API IfpackPrecond<std::complex<double> >;
  }
}
#endif