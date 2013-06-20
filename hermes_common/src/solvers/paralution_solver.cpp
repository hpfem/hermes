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
      this->paralutionMatrix.Clear();
      this->Ap = NULL;
      this->Ai = NULL;
      this->Ax = NULL;
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
      this->paralutionVector.Clear();
      v = NULL;
      this->size = 0;
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
    Scalar ParalutionVector<Scalar>::get(unsigned int idx) const
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
    ParalutionLinearMatrixSolver<Scalar>::ParalutionLinearMatrixSolver(ParalutionMatrix<Scalar> *matrix, ParalutionVector<Scalar> *rhs) : IterSolver<Scalar>(), matrix(matrix), rhs(rhs), preconditioner(NULL), paralutionSolverType(CG)
    {
      this->set_max_iters(1000);
      this->set_precond(new Preconditioners::ParalutionPrecond<Scalar>(Hermes::Preconditioners::ParalutionPrecond<Scalar>::ILU));
      this->set_tolerance(1e-8, IterSolver<double>::AbsoluteTolerance);
    }

    template<typename Scalar>
    ParalutionLinearMatrixSolver<Scalar>::~ParalutionLinearMatrixSolver()
    {
      if(preconditioner)
        delete preconditioner;
      this->sln = NULL;
    }

    template<typename Scalar>
    void ParalutionLinearMatrixSolver<Scalar>::set_solver_type(typename ParalutionLinearMatrixSolver<Scalar>::ParalutionSolverType paralutionSolverType)
    {
      this->paralutionSolverType = paralutionSolverType;
    }

    template<typename Scalar>
    paralution::IterativeLinearSolver<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>* ParalutionLinearMatrixSolver<Scalar>::create_paralutionSolver()
    {
      switch(this->paralutionSolverType)
      {
      case CG:
        {
          return new paralution::CG<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>();
        }
        break;
      case GMRES:
        {
          return new paralution::GMRES<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>();
        }
        break;
      case BiCGStab:
        {
          return new paralution::BiCGStab<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>();
        }
        break;
      case AMG:
        {
          throw Hermes::Exceptions::MethodNotImplementedException("Algebraic multigrid interface.");
          return NULL;
        }
        break;
      default:
        throw Hermes::Exceptions::Exception("A wrong Paralution solver type detected.");
      }
    }


    template<typename Scalar>
    bool ParalutionLinearMatrixSolver<Scalar>::solve()
    {
      if(this->sln)
        delete [] this->sln;
      this->sln = new Scalar[this->get_matrix_size()];
      memset(this->sln, Scalar(0), this->get_matrix_size() * sizeof(Scalar));
      return this->solve(this->sln);
    }

    template<typename Scalar>
    bool ParalutionLinearMatrixSolver<Scalar>::solve(Scalar* initial_guess)
    {
      // Handle sln.
      if(this->sln)
          delete [] this->sln;
      this->sln = new Scalar[this->get_matrix_size()];

      // Create initial guess.
      if(initial_guess)
        memcpy(this->sln, initial_guess, this->get_matrix_size() * sizeof(Scalar));
      else
        memset(this->sln, Scalar(0), this->get_matrix_size() * sizeof(Scalar));
      
      paralution::LocalVector<Scalar> x;
      x.SetDataPtr(&this->sln, "Initial guess", matrix->get_size());

      // Handle the situation when rhs == 0(vector).
      if(std::abs(rhs->get_paralutionVector().Norm()) < Hermes::epsilon)
      {
        x.LeaveDataPtr(&this->sln);
        x.Clear();
        return true;
      }

      // Create a solver according to the current type settings.
      paralution::IterativeLinearSolver<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>* paralutionSolver = this->create_paralutionSolver();

      // Set tolerances.
      switch(this->toleranceType)
      {
      case IterSolver<Scalar>::AbsoluteTolerance:
        paralutionSolver->InitTol(this->tolerance, 0., std::numeric_limits<Scalar>::max());
        break;
      case IterSolver<Scalar>::RelativeTolerance:
        paralutionSolver->InitTol(std::numeric_limits<Scalar>::max(), this->tolerance, std::numeric_limits<Scalar>::max());
        break;
      case IterSolver<Scalar>::DivergenceTolerance:
        paralutionSolver->InitTol(std::numeric_limits<Scalar>::max(), 0., this->tolerance);
        break;
      }

      // Set max iters.
      paralutionSolver->InitMaxIter(this->max_iters);

      // Set verbose_level.
      if(this->get_verbose_output())
        paralutionSolver->Verbose(10);

      // Set operator, preconditioner, build, solve.
      paralutionSolver->SetOperator(this->matrix->get_paralutionMatrix());
      if(this->preconditioner)
        paralutionSolver->SetPreconditioner(this->preconditioner->get_paralutionPreconditioner());
      paralutionSolver->Build();
      paralutionSolver->Solve(rhs->get_paralutionVector(), &x);

      // Store num_iters.
      num_iters = paralutionSolver->GetIterationCount();

      // Store final_residual
      final_residual = paralutionSolver->GetCurrentResidual();

      // Destroy the paralution vector, keeping the data in sln.
      x.LeaveDataPtr(&this->sln);

      // Delete the paralution solver.
      delete paralutionSolver;

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
      return this->num_iters;
    }

    template<typename Scalar>
    void ParalutionLinearMatrixSolver<Scalar>::set_max_iters(int iters)
    {
      IterSolver<Scalar>::set_max_iters(iters);
    }

    template<typename Scalar>
    void ParalutionLinearMatrixSolver<Scalar>::set_tolerance(double tol, typename IterSolver<Scalar>::ToleranceType toleranceType)
    {
      IterSolver<Scalar>::set_tolerance(tol, toleranceType);
      
    }

    template<typename Scalar>
    void ParalutionLinearMatrixSolver<Scalar>::set_verbose_output(bool to_set)
    {
      Hermes::Mixins::Loggable::set_verbose_output(to_set);
    }

    template<typename Scalar>
    double ParalutionLinearMatrixSolver<Scalar>::get_residual()
    {
      return final_residual;
    }

    template<typename Scalar>
    void ParalutionLinearMatrixSolver<Scalar>::set_precond(Precond<Scalar> *pc)
    {
      if(this->preconditioner)
        delete this->preconditioner;

      Preconditioners::ParalutionPrecond<Scalar>* paralutionPreconditioner = dynamic_cast<Preconditioners::ParalutionPrecond<Scalar>*>(pc);
      if(paralutionPreconditioner)
        this->preconditioner = paralutionPreconditioner;
      else
        throw Hermes::Exceptions::Exception("A wrong preconditioner type passed to Paralution.");
    }

    template class HERMES_API ParalutionLinearMatrixSolver<double>;
  }

  namespace Preconditioners
  {
    template<typename Scalar>
    ParalutionPrecond<Scalar>::ParalutionPrecond(typename ParalutionPrecond<Scalar>::ParalutionPreconditionerType paralutionPrecondType) : Precond<Scalar>()
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

    template class HERMES_API ParalutionPrecond<double>;
  }
}
#endif
