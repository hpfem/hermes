// This file is part of HermesCommon
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
/*! \file paralution_solver.cpp
\brief PARALUTION solver interface.
*/
#include "config.h"
#ifdef WITH_PARALUTION
#include "paralution_solver.h"

#include <limits>

namespace Hermes
{
  namespace Algebra
  {
    template<typename Scalar>
    ParalutionMatrix<Scalar>::ParalutionMatrix(ParalutionMatrixType type) 
      : CSRMatrix<Scalar>(), paralutionMatrixType(type)
    {
    }

    template<typename Scalar>
    ParalutionMatrix<Scalar>::~ParalutionMatrix()
    {
      this->paralutionMatrix.Clear();
      this->Ap = NULL;
      this->Ai = NULL;
      this->Ax = NULL;
    }

    template<typename Scalar>
    void ParalutionMatrix<Scalar>::free()
    {
      CSRMatrix<Scalar>::free();
    }

    template<typename Scalar>
    void ParalutionMatrix<Scalar>::zero()
    {
      CSRMatrix<Scalar>::zero();
      this->paralutionMatrix.Zeros();
    }

    template<typename Scalar>
    paralution::LocalMatrix<Scalar>& ParalutionMatrix<Scalar>::get_paralutionMatrix()
    {
      return this->paralutionMatrix;
    }

    template<typename Scalar>
    void ParalutionMatrix<Scalar>::alloc()
    {
      CSRMatrix<Scalar>::alloc();

      this->paralutionMatrix.SetDataPtrCSR(&this->Ap, &this->Ai, &this->Ax, "paralutionMatrix", this->nnz, this->size, this->size);
    }


    template<typename Scalar>
    ParalutionVector<Scalar>::ParalutionVector() : Vector<Scalar>(), v(NULL)
    {
    }

    template<typename Scalar>
    ParalutionVector<Scalar>::ParalutionVector(unsigned int size) : Vector<Scalar>(size), v(NULL)
    {
      this->alloc(size);
      this->paralutionVector.SetDataPtr(&this->v, "paralutionVector", this->size);
    }

    template<typename Scalar>
    void ParalutionVector<Scalar>::alloc(unsigned int n)
    {
      free();
      this->size = n;
      v = new Scalar[n];
      this->zero();
      this->paralutionVector.SetDataPtr(&this->v, "vector", this->size);
    }

    template<typename Scalar>
    ParalutionVector<Scalar>::~ParalutionVector()
    {
      this->paralutionVector.Clear();
    }

    template<typename Scalar>
    void ParalutionVector<Scalar>::free()
    {
      delete [] v;
      v = NULL;
      this->size = 0;
      this->paralutionVector.Clear();
    }

    template<typename Scalar>
    paralution::LocalVector<Scalar>& ParalutionVector<Scalar>::get_paralutionVector()
    {
      return this->paralutionVector;
    }

    template<typename Scalar>
    void ParalutionVector<Scalar>::zero()
    {
      memset(v, 0, this->size * sizeof(Scalar));
      this->paralutionVector.Zeros();
    }

    template<typename Scalar>
    void ParalutionVector<Scalar>::change_sign()
    {
      for (unsigned int i = 0; i < this->size; i++) v[i] *= -1.;
    }

    template<typename Scalar>
    void ParalutionVector<Scalar>::set(unsigned int idx, Scalar y)
    {
      v[idx] = y;
    }

    template<>
    void ParalutionVector<double>::add(unsigned int idx, double y)
    {
#pragma omp atomic
      v[idx] += y;
    }


    template<typename Scalar>
    void ParalutionVector<Scalar>::add(unsigned int n, unsigned int *idx, Scalar *y)
    {
      for (unsigned int i = 0; i < n; i++)
        v[idx[i]] += y[i];
    }

    template<typename Scalar>
    Scalar ParalutionVector<Scalar>::get(unsigned int idx)
    {
      return v[idx];
    }

    template<typename Scalar>
    void ParalutionVector<Scalar>::extract(Scalar *v) const
    {
      memcpy(v, this->v, this->size * sizeof(Scalar));
    }

    template<typename Scalar>
    void ParalutionVector<Scalar>::add_vector(Vector<Scalar>* vec)
    {
      assert(this->length() == vec->length());
      for (unsigned int i = 0; i < this->length(); i++) this->v[i] += vec->get(i);
    }

    template<typename Scalar>
    void ParalutionVector<Scalar>::add_vector(Scalar* vec)
    {
      for (unsigned int i = 0; i < this->length(); i++)
        this->v[i] += vec[i];
    }

    template<>
    bool ParalutionVector<double>::dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt, char* number_format)
    {
      switch (fmt)
      {
      case DF_MATLAB_SPARSE:
        fprintf(file, "%% Size: %dx1\n%s =[\n", this->size, var_name);
        for (unsigned int i = 0; i < this->size; i++)
        {
          Hermes::Helpers::fprint_num(file, v[i], number_format);
          fprintf(file, "\n");
        }
        fprintf(file, " ];\n");
        return true;

      case DF_HERMES_BIN:
        {
          hermes_fwrite("HERMESR\001", 1, 8, file);
          int ssize = sizeof(double);
          hermes_fwrite(&ssize, sizeof(int), 1, file);
          hermes_fwrite(&this->size, sizeof(int), 1, file);
          hermes_fwrite(v, sizeof(double), this->size, file);
          return true;
        }

      case DF_PLAIN_ASCII:
        {
          fprintf(file, "\n");
          for (unsigned int i = 0; i < size; i++)
          {
            Hermes::Helpers::fprint_num(file, v[i], number_format);
            fprintf(file, "\n");
          }

          return true;
        }

      default:
        return false;
      }
    }

    template class HERMES_API ParalutionMatrix<double>;
    template class HERMES_API ParalutionVector<double>;
  }
  namespace Solvers
  {
    template<typename Scalar>
    ParalutionLinearMatrixSolver<Scalar>::ParalutionLinearMatrixSolver(ParalutionMatrix<Scalar> *matrix, ParalutionVector<Scalar> *rhs) : IterSolver<Scalar>(), matrix(matrix), rhs(rhs), preconditioner(NULL)
    {
      paralutionSolver = new paralution::CG<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>();
      this->set_max_iters(10000);
    }

    template<typename Scalar>
    ParalutionLinearMatrixSolver<Scalar>::~ParalutionLinearMatrixSolver()
    {
      if(preconditioner)
        delete preconditioner;
    }

    template<typename Scalar>
    void ParalutionLinearMatrixSolver<Scalar>::set_solver_type(typename ParalutionLinearMatrixSolver<Scalar>::ParalutionSolverType paralutionSolverType)
    {
      if(this->paralutionSolver)
        delete this->paralutionSolver;

      switch(paralutionSolverType)
      {
        case CG:
          {
            this->paralutionSolver = new paralution::CG<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>();
          }
          break;
        case GMRES:
          {
            this->paralutionSolver = new paralution::GMRES<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>();
          }
          break;
        case BiCGStab:
          {
            this->paralutionSolver = new paralution::BiCGStab<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>();
          }
          break;
        case AMG:
          {
            this->paralutionSolver = new paralution::AMG<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>();
          }
          break;
        default:
        throw Hermes::Exceptions::Exception("A wrong Paralution solver type passed to ParalutionLinearMatrixSolver constructor.");
      }

      this->info("PARALUTION solver successfully set, you need to re-set the preconditioner should you wish to use it.");
    }

    template<typename Scalar>
    bool ParalutionLinearMatrixSolver<Scalar>::solve()
    {
      Scalar* initial_guess = new Scalar[this->get_matrix_size()];
      bool result = this->solve(initial_guess);
      delete [] initial_guess;
      return result;
    }

    template<typename Scalar>
    bool ParalutionLinearMatrixSolver<Scalar>::solve(Scalar* initial_guess)
    {
      paralution::LocalVector<Scalar> x;
      x.SetDataPtr(&initial_guess, "Initial guess", matrix->get_size());

      if(std::abs(rhs->get_paralutionVector().Norm()) < Hermes::epsilon)
      {
        x.LeaveDataPtr(&this->sln);
        x.Clear();
        return true;
      }

      rhs->get_paralutionVector().WriteFileASCII("rhs");
      assert(matrix != NULL);
      assert(rhs != NULL);
      assert(matrix->get_size() == rhs->length());

      this->paralutionSolver->SetOperator(this->matrix->get_paralutionMatrix());
      this->paralutionSolver->Build();

      this->paralutionSolver->Solve(rhs->get_paralutionVector(), &x);

      x.LeaveDataPtr(&this->sln);
      this->paralutionSolver->Clear();

      return true;
    }

    template<typename Scalar>
    int ParalutionLinearMatrixSolver<Scalar>::get_matrix_size()
    {
      return matrix->get_size();
    }

    template<typename Scalar>
    int ParalutionLinearMatrixSolver<Scalar>::get_num_iters()
    {
      return this->paralutionSolver->GetIterationCount();
    }

    template<typename Scalar>
    void ParalutionLinearMatrixSolver<Scalar>::set_max_iters(int iters)
    {
      IterSolver<Scalar>::set_max_iters(iters);
      this->paralutionSolver->InitMaxIter(iters);
    }

    template<typename Scalar>
    void ParalutionLinearMatrixSolver<Scalar>::set_tolerance(double tol, typename IterSolver<Scalar>::ToleranceType toleranceType)
    {
      IterSolver<Scalar>::set_tolerance(tol, toleranceType);
      switch(toleranceType)
      {
      case IterSolver<Scalar>::AbsoluteTolerance:
        this->paralutionSolver->InitTol(tol, 1., 1.);
        break;
      case IterSolver<Scalar>::RelativeTolerance:
        this->paralutionSolver->InitTol(std::numeric_limits<Scalar>::max(), tol, 1.);
        break;
      case IterSolver<Scalar>::DivergenceTolerance:
        this->paralutionSolver->InitTol(std::numeric_limits<Scalar>::max(), 1., tol);
        break;
      }
    }

    template<typename Scalar>
    void ParalutionLinearMatrixSolver<Scalar>::set_verbose_output(bool to_set)
    {
      Hermes::Mixins::Loggable::set_verbose_output(to_set);
      this->paralutionSolver->Verbose(10);
    }

    template<typename Scalar>
    double ParalutionLinearMatrixSolver<Scalar>::get_residual()
    {
      return 0.;
    }

    template<typename Scalar>
    void ParalutionLinearMatrixSolver<Scalar>::set_precond(Precond<Scalar> *pc)
    {
      ParalutionPrecond<Scalar>* paralutionPreconditioner = dynamic_cast<ParalutionPrecond<Scalar>*>(pc);
      if(paralutionPreconditioner)
        this->preconditioner = paralutionPreconditioner;
      else
        throw Hermes::Exceptions::Exception("A wrong preconditioner type passed to Paralution.");
    }

    template<typename Scalar>
    ParalutionPrecond<Scalar>::ParalutionPrecond(typename ParalutionPrecond<Scalar>::ParalutionPrecondType paralutionPrecondType) : Precond<Scalar>()
    {
      switch(paralutionPrecondType)
      {
        case Jacobi:
          {
            this->paralutionPreconditioner = new paralution::Jacobi<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>();
          }
          break;
        case ILU:
          {
            this->paralutionPreconditioner = new paralution::ILU<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>();
          }
          break;
        case MultiColoredILU:
          {
            this->paralutionPreconditioner = new paralution::MultiColoredILU<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>();
          }
          break;
        case MultiColoredSGS:
          {
            this->paralutionPreconditioner = new paralution::MultiColoredSGS<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>();
          }
          break;
        case IC:
          {
            this->paralutionPreconditioner = new paralution::IC<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>();
          }
          break;
        case AIChebyshev:
          {
            this->paralutionPreconditioner = new paralution::AIChebyshev<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>();
          }
          break;
        default:
        throw Hermes::Exceptions::Exception("A wrong Paralution preconditioner type passed to ParalutionPrecond constructor.");
      }
    }

    template<typename Scalar>
    paralution::Preconditioner<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>& ParalutionPrecond<Scalar>::get_paralutionPreconditioner()
    {
      return (*this->paralutionPreconditioner);
    }

    template class HERMES_API ParalutionLinearMatrixSolver<double>;
    template class HERMES_API ParalutionPrecond<double>;
  }
}
#endif
