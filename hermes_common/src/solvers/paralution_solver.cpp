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
    ParalutionVector<Scalar>::ParalutionVector() : SimpleVector<Scalar>(), paralutionVector(new paralution::LocalVector<Scalar>)
    {
    }

    template<typename Scalar>
    ParalutionVector<Scalar>::ParalutionVector(unsigned int size) : SimpleVector<Scalar>(size), paralutionVector(new paralution::LocalVector<Scalar>)
    {
      this->alloc(size);
      this->paralutionVector->SetDataPtr(&this->v, "paralutionVector", this->size);
    }

    template<typename Scalar>
    void ParalutionVector<Scalar>::alloc(unsigned int n)
    {
      SimpleVector<Scalar>::alloc(n);
      this->paralutionVector->Clear();
      this->paralutionVector->SetDataPtr(&this->v, "vector", this->size);
    }

    template<typename Scalar>
    ParalutionVector<Scalar>::~ParalutionVector()
    {
      free();
      /*
         Temporarily suspended - introduced a memory leak
         Reason: heap corruption upon/after execution in Agros.
         Temporarily suspended - so far attempts to replicate in Hermes failed.
      */
      // delete this->paralutionVector;
    }

    template<typename Scalar>
    void ParalutionVector<Scalar>::free()
    {
      this->paralutionVector->Clear();
      this->v = NULL;
      SimpleVector<Scalar>::free();
    }

    template<typename Scalar>
    paralution::LocalVector<Scalar>* ParalutionVector<Scalar>::get_paralutionVector()
    {
      return this->paralutionVector;
    }

    template<typename Scalar>
    void ParalutionVector<Scalar>::zero()
    {
      memset(this->v, 0, this->size * sizeof(Scalar));
    }

    template class HERMES_API ParalutionMatrix<double>;
    template class HERMES_API ParalutionVector<double>;
  }
  namespace Solvers
  {
    template<typename Scalar>
    AbstractParalutionLinearMatrixSolver<Scalar>::AbstractParalutionLinearMatrixSolver() : LoopSolver<Scalar>(), matrix(NULL), rhs(NULL), paralutionSolver(NULL)
    {
      this->set_max_iters(1000);
      this->set_tolerance(1e-8, AbsoluteTolerance);
    }

    template<typename Scalar>
    AbstractParalutionLinearMatrixSolver<Scalar>::AbstractParalutionLinearMatrixSolver(ParalutionMatrix<Scalar> *matrix, ParalutionVector<Scalar> *rhs) : LoopSolver<Scalar>(), matrix(matrix), rhs(rhs), paralutionSolver(NULL)
    {
      this->set_max_iters(1000);
      this->set_tolerance(1e-8, AbsoluteTolerance);
    }

    template<typename Scalar>
    AbstractParalutionLinearMatrixSolver<Scalar>::~AbstractParalutionLinearMatrixSolver()
    {
      delete this->paralutionSolver;
      this->sln = NULL;
    }

    template<typename Scalar>
    void AbstractParalutionLinearMatrixSolver<Scalar>::reset_internal_solver()
    {
      if(this->paralutionSolver)
        delete this->paralutionSolver;
      this->paralutionSolver = NULL;
    }

    template<typename Scalar>
    void AbstractParalutionLinearMatrixSolver<Scalar>::solve()
    {
      this->solve(NULL);
    }

    
    template<typename Scalar>
    void AbstractParalutionLinearMatrixSolver<Scalar>::presolve_init()
    {
      // Reinit iff matrix has changed.
      if(this->reuse_scheme != HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY)
        this->reset_internal_solver();
      this->init_internal_solver();
      
      // Set verbose_level.
      if(this->get_verbose_output())
        this->paralutionSolver->Verbose(10);
      else
        this->paralutionSolver->Verbose(0);

      // Set tolerances.
      switch(this->toleranceType)
      {
      case AbsoluteTolerance:
        paralutionSolver->InitTol(this->tolerance, 0., std::numeric_limits<Scalar>::max());
        break;
      case RelativeTolerance:
        paralutionSolver->InitTol(0., this->tolerance, std::numeric_limits<Scalar>::max());
        break;
      case DivergenceTolerance:
        paralutionSolver->InitTol(0., 0., this->tolerance);
        break;
      }

      // Set max iters.
      paralutionSolver->InitMaxIter(this->max_iters);
    }

    template<typename Scalar>
    void AbstractParalutionLinearMatrixSolver<Scalar>::solve(Scalar* initial_guess)
    {
      // Handle sln.
      if(this->sln && this->sln != initial_guess)
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
      if(std::abs(rhs->get_paralutionVector()->Norm()) < Hermes::epsilon)
        x.LeaveDataPtr(&this->sln);

      // (Re-)init.
      this->presolve_init();

      // Move to Accelerators.
      if(HermesCommonApi.get_integral_param_value(useAccelerators))
      {
        this->paralutionSolver->MoveToAccelerator();
        this->matrix->get_paralutionMatrix().MoveToAccelerator();
        this->rhs->get_paralutionVector()->MoveToAccelerator();
        x.MoveToAccelerator();
      }

      // Solve.
      paralutionSolver->Solve(*rhs->get_paralutionVector(), &x);

      // Store num_iters.
      num_iters = paralutionSolver->GetIterationCount();

      // Store final_residual
      final_residual = paralutionSolver->GetCurrentResidual();

      // Destroy the paralution vector, keeping the data in sln.
      x.LeaveDataPtr(&this->sln);
    }

    template<typename Scalar>
    int AbstractParalutionLinearMatrixSolver<Scalar>::get_matrix_size()
    {
      return matrix->get_size();
    }

    template<typename Scalar>
    int AbstractParalutionLinearMatrixSolver<Scalar>::get_num_iters()
    {
      return this->num_iters;
    }

    template<typename Scalar>
    void AbstractParalutionLinearMatrixSolver<Scalar>::set_verbose_output(bool to_set)
    {
      Hermes::Mixins::Loggable::set_verbose_output(to_set);
    }

    template<typename Scalar>
    double AbstractParalutionLinearMatrixSolver<Scalar>::get_residual_norm()
    {
      return final_residual;
    }

    template<typename Scalar>
    IterativeParalutionLinearMatrixSolver<Scalar>::IterativeParalutionLinearMatrixSolver() : AbstractParalutionLinearMatrixSolver<Scalar>(), IterSolver<Scalar>(), preconditioner(NULL)
    {
      this->set_precond(new Preconditioners::ParalutionPrecond<Scalar>(ILU));
    }

    template<typename Scalar>
    IterativeParalutionLinearMatrixSolver<Scalar>::IterativeParalutionLinearMatrixSolver(ParalutionMatrix<Scalar> *matrix, ParalutionVector<Scalar> *rhs) : AbstractParalutionLinearMatrixSolver<Scalar>(matrix, rhs), IterSolver<Scalar>(), preconditioner(NULL)
    {
      this->set_precond(new Preconditioners::ParalutionPrecond<Scalar>(ILU));
    }

    template<typename Scalar>
    IterativeParalutionLinearMatrixSolver<Scalar>::~IterativeParalutionLinearMatrixSolver()
    {
      if(preconditioner)
        delete preconditioner;
    }

    template<typename Scalar>
    void IterativeParalutionLinearMatrixSolver<Scalar>::set_solver_type(IterSolverType iterSolverType)
    {
      IterSolver<double>::set_solver_type(iterSolverType);
      this->reset_internal_solver();
    }

    template<typename Scalar>
    paralution::IterativeLinearSolver<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>*
      IterativeParalutionLinearMatrixSolver<Scalar>::return_paralutionSolver(IterSolverType type)
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
        throw Hermes::Exceptions::Exception("A wrong solver type detected in PARALUTION.");
        return NULL;
      }
    }

    template<typename Scalar>
    void IterativeParalutionLinearMatrixSolver<Scalar>::init_internal_solver()
    {
      // Create a solver according to the current type settings.
      if(!this->paralutionSolver)
      {
        this->paralutionSolver = this->return_paralutionSolver(this->iterSolverType);

        // Set operator, preconditioner, build.
        if(this->preconditioner)
          this->paralutionSolver->SetPreconditioner(this->preconditioner->get_paralutionPreconditioner());
        this->paralutionSolver->SetOperator(this->matrix->get_paralutionMatrix());
        this->paralutionSolver->Build();
      }
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
    AMGParalutionLinearMatrixSolver<Scalar>::AMGParalutionLinearMatrixSolver(ParalutionMatrix<Scalar> *matrix, ParalutionVector<Scalar> *rhs) : AbstractParalutionLinearMatrixSolver<Scalar>(matrix, rhs), AMGSolver<Scalar>()
    {
    }

    template<typename Scalar>
    void AMGParalutionLinearMatrixSolver<Scalar>::set_smoother(IterSolverType solverType_, PreconditionerType preconditionerType_)
    {
      AMGSolver<double>::set_smoother(solverType_, preconditionerType_);
    }

    template<typename Scalar>
    AMGParalutionLinearMatrixSolver<Scalar>::~AMGParalutionLinearMatrixSolver()
    {
    }

    template<typename Scalar>
    void AMGParalutionLinearMatrixSolver<Scalar>::init_internal_solver()
    {
      // Create a solver according to the current type settings.
      if(!this->paralutionSolver)
      {
        paralution::AMG<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>* AMG_solver;
        this->paralutionSolver = AMG_solver = new paralution::AMG<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>();
        AMG_solver->SetManualSmoothers(true);
        AMG_solver->SetOperator(this->matrix->get_paralutionMatrix());
        AMG_solver->BuildHierarchy();

        // Set operator, smoother, build.
        int levels = AMG_solver->GetNumLevels();
        paralution::IterativeLinearSolver<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar >** smoothers = new paralution::IterativeLinearSolver<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar >*[levels-1];
        paralution::Preconditioner<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar >** preconditioners = new paralution::Preconditioner<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar >*[levels-1];

        for (int i = 0; i < levels - 1; ++i)
        {
          smoothers[i] = IterativeParalutionLinearMatrixSolver<Scalar>::return_paralutionSolver(this->smootherSolverType);
          preconditioners[i] = ParalutionPrecond<Scalar>::return_paralutionPreconditioner(this->smootherPreconditionerType);

          smoothers[i]->SetPreconditioner(*preconditioners[i]);
          if(this->get_verbose_output())
            smoothers[i]->Verbose(10);
          else
            smoothers[i]->Verbose(0);
        }

        AMG_solver->SetSmoother(smoothers);
        AMG_solver->SetSmootherPreIter(3);
        AMG_solver->SetSmootherPostIter(3);

        AMG_solver->Build();
      }
    }

    template class HERMES_API IterativeParalutionLinearMatrixSolver<double>;
    template class HERMES_API AMGParalutionLinearMatrixSolver<double>;
  }

  namespace Preconditioners
  {
    template<typename Scalar>
    ParalutionPrecond<Scalar>::ParalutionPrecond(PreconditionerType preconditionerType) : Precond<Scalar>()
    {
      switch(preconditionerType)
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
        throw Hermes::Exceptions::Exception("A wrong paralution preconditioner type passed to ParalutionPrecond constructor.");
      }
    }

    template<typename Scalar>
    paralution::Preconditioner<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>& ParalutionPrecond<Scalar>::get_paralutionPreconditioner()
    {
      return (*this->paralutionPreconditioner);
    }

    template<typename Scalar>
    paralution::Preconditioner<paralution::LocalMatrix<Scalar>, paralution::LocalVector<Scalar>, Scalar>* ParalutionPrecond<Scalar>::return_paralutionPreconditioner(PreconditionerType type)
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
        throw Hermes::Exceptions::Exception("A wrong preconditioner type passed to ParalutionPrecond constructor.");
        return NULL;
      }
    }

    template class HERMES_API ParalutionPrecond<double>;
  }
}
#endif
