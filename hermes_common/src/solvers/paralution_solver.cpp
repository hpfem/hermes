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
    ParalutionVector<Scalar>::ParalutionVector() : SimpleVector<Scalar>()
    {
    }

    template<typename Scalar>
    ParalutionVector<Scalar>::ParalutionVector(unsigned int size) : SimpleVector<Scalar>(size)
    {
      this->alloc(size);
      this->paralutionVector.SetDataPtr(&this->v, "paralutionVector", this->size);
    }

    template<typename Scalar>
    void ParalutionVector<Scalar>::alloc(unsigned int n)
    {
      free();
      this->size = n;
      this->v = new Scalar[n];
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
      this->v = NULL;
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
      memset(this->v, 0, this->size * sizeof(Scalar));
      this->paralutionVector.Zeros();
    }

    template class HERMES_API ParalutionMatrix<double>;
    template class HERMES_API ParalutionVector<double>;
  }
  namespace Solvers
  {
    template<typename Scalar>
    IterativeParalutionLinearMatrixSolver<Scalar>::IterativeParalutionLinearMatrixSolver() : IterSolver<Scalar>(), matrix(NULL), rhs(NULL), preconditioner(NULL), paralutionSolverType(CG), paralutionSolver(NULL)
    {
      this->set_max_iters(1000);
      this->set_precond(new Preconditioners::ParalutionPrecond<Scalar>(Hermes::Preconditioners::ParalutionPrecond<Scalar>::ILU));
      this->set_tolerance(1e-8, IterSolver<double>::AbsoluteTolerance);
    }

    template<typename Scalar>
    IterativeParalutionLinearMatrixSolver<Scalar>::IterativeParalutionLinearMatrixSolver(ParalutionMatrix<Scalar> *matrix, ParalutionVector<Scalar> *rhs) : IterSolver<Scalar>(), matrix(matrix), rhs(rhs), preconditioner(NULL), paralutionSolverType(CG), paralutionSolver(NULL)
    {
      this->set_max_iters(1000);
      this->set_precond(new Preconditioners::ParalutionPrecond<Scalar>(Hermes::Preconditioners::ParalutionPrecond<Scalar>::ILU));
      this->set_tolerance(1e-8, IterSolver<double>::AbsoluteTolerance);
    }

    template<typename Scalar>
    void IterativeParalutionLinearMatrixSolver<Scalar>::reset_internal_solver()
    {
      if(this->paralutionSolver)
        delete this->paralutionSolver;
      this->paralutionSolver = NULL;
    }

    template<typename Scalar>
    IterativeParalutionLinearMatrixSolver<Scalar>::~IterativeParalutionLinearMatrixSolver()
    {
      if(this->paralutionSolver)
        delete this->paralutionSolver;
      if(preconditioner)
        delete preconditioner;
      this->sln = NULL;
    }

    template<typename Scalar>
    void IterativeParalutionLinearMatrixSolver<Scalar>::set_solver_type(typename IterativeParalutionLinearMatrixSolver<Scalar>::ParalutionSolverType paralutionSolverType)
    {
      this->paralutionSolverType = paralutionSolverType;
      this->reset_internal_solver();
    }

    template<typename Scalar>
    paralution::IterativeLinearSolver<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>*
      IterativeParalutionLinearMatrixSolver<Scalar>::return_paralutionSolver(typename IterativeParalutionLinearMatrixSolver<Scalar>::ParalutionSolverType type)
    {
      switch(type)
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
      default:
        throw Hermes::Exceptions::Exception("A wrong Paralution solver type detected.");
        return NULL;
      }
    }

    template<typename Scalar>
    void IterativeParalutionLinearMatrixSolver<Scalar>::solve()
    {
      if(this->sln)
        delete [] this->sln;
      this->sln = new Scalar[this->get_matrix_size()];
      memset(this->sln, Scalar(0), this->get_matrix_size() * sizeof(Scalar));
      this->solve(this->sln);
    }

    template<typename Scalar>
    void IterativeParalutionLinearMatrixSolver<Scalar>::init_internal_solver()
    {
      // Reset if the matrix has changed.
      if(this->reuse_scheme != HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY)
        this->reset_internal_solver();

      // Create a solver according to the current type settings.
      if(!this->paralutionSolver)
      {
        this->paralutionSolver = this->return_paralutionSolver(this->paralutionSolverType);

        // Set operator, preconditioner, build.
        if(this->preconditioner)
          paralutionSolver->SetPreconditioner(this->preconditioner->get_paralutionPreconditioner());
        paralutionSolver->SetOperator(this->matrix->get_paralutionMatrix());
        paralutionSolver->Build();
      }

      // Set verbose_level.
      if(this->get_verbose_output())
        this->paralutionSolver->Verbose(10);

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
    }

    template<typename Scalar>
    void IterativeParalutionLinearMatrixSolver<Scalar>::solve(Scalar* initial_guess)
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
      }

      // Init.
      this->init_internal_solver();

      // Solve.
      paralutionSolver->Solve(rhs->get_paralutionVector(), &x);

      // Store num_iters.
      num_iters = paralutionSolver->GetIterationCount();

      // Store final_residual
      final_residual = paralutionSolver->GetCurrentResidual();

      // Destroy the paralution vector, keeping the data in sln.
      x.LeaveDataPtr(&this->sln);
    }

    template<typename Scalar>
    int IterativeParalutionLinearMatrixSolver<Scalar>::get_matrix_size()
    {
      return matrix->get_size();
    }

    template<typename Scalar>
    int IterativeParalutionLinearMatrixSolver<Scalar>::get_num_iters()
    {
      return this->num_iters;
    }

    template<typename Scalar>
    void IterativeParalutionLinearMatrixSolver<Scalar>::set_verbose_output(bool to_set)
    {
      Hermes::Mixins::Loggable::set_verbose_output(to_set);
    }

    template<typename Scalar>
    double IterativeParalutionLinearMatrixSolver<Scalar>::get_residual()
    {
      return final_residual;
    }

    template<typename Scalar>
    void IterativeParalutionLinearMatrixSolver<Scalar>::set_precond(Precond<Scalar> *pc)
    {
      if(this->preconditioner)
        delete this->preconditioner;

      Preconditioners::ParalutionPrecond<Scalar>* paralutionPreconditioner = dynamic_cast<Preconditioners::ParalutionPrecond<Scalar>*>(pc);
      if(paralutionPreconditioner)
        this->preconditioner = paralutionPreconditioner;
      else
        throw Hermes::Exceptions::Exception("A wrong preconditioner type passed to Paralution.");
    }

    template<typename Scalar>
    AMGParalutionLinearMatrixSolver<Scalar>::AMGParalutionLinearMatrixSolver(ParalutionMatrix<Scalar> *matrix, ParalutionVector<Scalar> *rhs) : AMGSolver<Scalar>(), matrix(matrix), rhs(rhs), paralutionSolver(NULL)
    {
      this->set_max_iters(1000);
      this->set_tolerance(1e-8, LoopSolver<Scalar>::AbsoluteTolerance);
      this->smootherSolverType = IterativeParalutionLinearMatrixSolver<Scalar>::CG;
      this->smootherPreconditionerType = ParalutionPrecond<Scalar>::MultiColoredSGS;
    }

    template<typename Scalar>
    void AMGParalutionLinearMatrixSolver<Scalar>::reset_internal_solver()
    {
      if(this->paralutionSolver)
        delete this->paralutionSolver;
      this->paralutionSolver = NULL;
    }

    template<typename Scalar>
    AMGParalutionLinearMatrixSolver<Scalar>::~AMGParalutionLinearMatrixSolver()
    {
      delete this->paralutionSolver;
      this->sln = NULL;
    }

    template<typename Scalar>
    void AMGParalutionLinearMatrixSolver<Scalar>::solve()
    {
      if(this->sln)
        delete [] this->sln;
      this->sln = new Scalar[this->get_matrix_size()];
      memset(this->sln, Scalar(0), this->get_matrix_size() * sizeof(Scalar));
      this->solve(this->sln);
    }

    template<typename Scalar>
    void AMGParalutionLinearMatrixSolver<Scalar>::init_internal_solver()
    {
      // Reset if the matrix has changed.
      if(this->reuse_scheme != HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY)
        this->reset_internal_solver();

      // Create a solver according to the current type settings.
      if(!this->paralutionSolver)
      {
        this->paralutionSolver = new paralution::AMG<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>();
        this->paralutionSolver->SetManualSmoothers(true);
        this->paralutionSolver->SetOperator(this->matrix->get_paralutionMatrix());
        this->paralutionSolver->BuildHierarchy();

        // Set operator, smoother, build.
        int levels = this->paralutionSolver->GetNumLevels();
        paralution::IterativeLinearSolver<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar >** smoothers = new paralution::IterativeLinearSolver<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar >*[levels-1];
        paralution::Preconditioner<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar >** preconditioners = new paralution::Preconditioner<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar >*[levels-1];

        for (int i = 0; i < levels - 1; ++i)
        {
          smoothers[i] = IterativeParalutionLinearMatrixSolver<Scalar>::return_paralutionSolver(this->smootherSolverType);
          preconditioners[i] = ParalutionPrecond<Scalar>::return_paralutionPreconditioner(this->smootherPreconditionerType);

          smoothers[i]->SetPreconditioner(*preconditioners[i]);
          smoothers[i]->Verbose(0);
        }

        paralutionSolver->SetSmoother(smoothers);
        paralutionSolver->SetSmootherPreIter(1);
        paralutionSolver->SetSmootherPostIter(2);

        paralutionSolver->Build();
      }

      // Set verbose_level.
      if(this->get_verbose_output())
        this->paralutionSolver->Verbose(10);

      // Set tolerances.
      switch(this->toleranceType)
      {
      case AMGSolver<Scalar>::AbsoluteTolerance:
        paralutionSolver->InitTol(this->tolerance, 0., std::numeric_limits<Scalar>::max());
        break;
      case AMGSolver<Scalar>::RelativeTolerance:
        paralutionSolver->InitTol(std::numeric_limits<Scalar>::max(), this->tolerance, std::numeric_limits<Scalar>::max());
        break;
      case AMGSolver<Scalar>::DivergenceTolerance:
        paralutionSolver->InitTol(std::numeric_limits<Scalar>::max(), 0., this->tolerance);
        break;
      }

      // Set max iters.
      paralutionSolver->InitMaxIter(this->max_iters);

      if(HermesCommonApi.get_integral_param_value(useAccelerators))
      {
        this->paralutionSolver->MoveToAccelerator();
        this->matrix->get_paralutionMatrix().MoveToAccelerator();
        this->rhs->get_paralutionVector().MoveToAccelerator();
      }
    }

    template<typename Scalar>
    void AMGParalutionLinearMatrixSolver<Scalar>::solve(Scalar* initial_guess)
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
      }

      // Init.
      this->init_internal_solver();

      // Solve.
      paralutionSolver->Solve(rhs->get_paralutionVector(), &x);

      // Store num_iters.
      num_iters = paralutionSolver->GetIterationCount();

      // Store final_residual
      final_residual = paralutionSolver->GetCurrentResidual();

      // Destroy the paralution vector, keeping the data in sln.
      x.LeaveDataPtr(&this->sln);
    }

    template<typename Scalar>
    int AMGParalutionLinearMatrixSolver<Scalar>::get_matrix_size()
    {
      return matrix->get_size();
    }

    template<typename Scalar>
    int AMGParalutionLinearMatrixSolver<Scalar>::get_num_iters()
    {
      return this->num_iters;
    }

    template<typename Scalar>
    void AMGParalutionLinearMatrixSolver<Scalar>::set_verbose_output(bool to_set)
    {
      Hermes::Mixins::Loggable::set_verbose_output(to_set);
    }

    template<typename Scalar>
    double AMGParalutionLinearMatrixSolver<Scalar>::get_residual()
    {
      return final_residual;
    }

    template<typename Scalar>
    void AMGParalutionLinearMatrixSolver<Scalar>::set_smoother(typename IterativeParalutionLinearMatrixSolver<Scalar>::ParalutionSolverType solverType_, typename ParalutionPrecond<Scalar>::ParalutionPreconditionerType preconditionerType_)
    {
      this->smootherPreconditionerType = preconditionerType_;
      this->smootherSolverType = solverType_;
    }

    template class HERMES_API IterativeParalutionLinearMatrixSolver<double>;
    template class HERMES_API AMGParalutionLinearMatrixSolver<double>;
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

    template<typename Scalar>
    paralution::Preconditioner<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>* ParalutionPrecond<Scalar>::return_paralutionPreconditioner(typename ParalutionPrecond<Scalar>::ParalutionPreconditionerType type)
    {
      switch(type)
      {
      case Jacobi:
        {
          return new paralution::Jacobi<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>();
        }
        break;
      case ILU:
        {
          return new paralution::ILU<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>();
        }
        break;
      case MultiColoredILU:
        {
          return new paralution::MultiColoredILU<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>();
        }
        break;
      case MultiColoredSGS:
        {
          return new paralution::MultiColoredSGS<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>();
        }
        break;
      case IC:
        {
          return new paralution::IC<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>();
        }
        break;
      case AIChebyshev:
        {
          return new paralution::AIChebyshev<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>();
        }
        break;
      default:
        throw Hermes::Exceptions::Exception("A wrong Paralution preconditioner type passed to ParalutionPrecond constructor.");
        return NULL;
      }
    }

    template class HERMES_API ParalutionPrecond<double>;
  }
}
#endif
