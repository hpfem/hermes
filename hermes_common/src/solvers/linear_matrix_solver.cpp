// This file is part of HermesCommon
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes2D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
/*! \file linear_solver.cpp
\brief General linear solver functionality.
*/
#include "linear_matrix_solver.h"
#include "umfpack_solver.h"
#include "superlu_solver.h"
#include "amesos_solver.h"
#include "petsc_solver.h"
#include "mumps_solver.h"
#include "aztecoo_solver.h"
#include "paralution_solver.h"
#include "api.h"
#include "exceptions.h"

using namespace Hermes::Algebra;

namespace Hermes
{
  namespace Solvers
  {
    template<typename Scalar>
    LinearMatrixSolver<Scalar>::LinearMatrixSolver(MatrixStructureReuseScheme reuse_scheme) : reuse_scheme(reuse_scheme)
    {
      sln = NULL;
      time = -1.0;
      n_eq = 1;
      node_wise_ordering = false;
    }

    template<typename Scalar>
    LinearMatrixSolver<Scalar>::~LinearMatrixSolver()
    {
      if(sln != NULL)
        delete [] sln;
    }

    template<typename Scalar>
    Scalar *LinearMatrixSolver<Scalar>::get_sln_vector()
    {
      return sln;
    }

    template<typename Scalar>
    double LinearMatrixSolver<Scalar>::get_time()
    {
      return time;
    }

    template<typename Scalar>
    void LinearMatrixSolver<Scalar>::set_reuse_scheme()
    {
      set_reuse_scheme(HERMES_CREATE_STRUCTURE_FROM_SCRATCH);
    }

    template<typename Scalar>
    MatrixStructureReuseScheme LinearMatrixSolver<Scalar>::get_used_reuse_scheme() const 
    { 
      return reuse_scheme; 
    }

    template<typename Scalar>
    void LinearMatrixSolver<Scalar>::set_reuse_scheme(MatrixStructureReuseScheme reuse_scheme) 
    {
      this->reuse_scheme = reuse_scheme; 
    }

    template<typename Scalar>
    void LinearMatrixSolver<Scalar>::use_node_wise_ordering(unsigned int num_pdes)
    {
      this->reuse_scheme = HERMES_CREATE_STRUCTURE_FROM_SCRATCH;
      this->n_eq = num_pdes;
      this->node_wise_ordering = true;
    }

    template<typename Scalar>
    void LinearMatrixSolver<Scalar>::use_equations_wise_ordering()
    {
      this->reuse_scheme = HERMES_CREATE_STRUCTURE_FROM_SCRATCH; 
      this->node_wise_ordering = false;
    }

    template<typename Scalar>
    double LinearMatrixSolver<Scalar>::get_residual_norm()
    {
      throw Exceptions::MethodNotOverridenException("LinearMatrixSolver<Scalar>::get_residual_norm");
      return 0.;
    }

    template<typename Scalar>
    DirectSolver<Scalar>* LinearMatrixSolver<Scalar>::as_DirectSolver() const
    {
      DirectSolver<Scalar>* solver = dynamic_cast<DirectSolver<Scalar>*>(const_cast<LinearMatrixSolver<Scalar>*>(this));
      if(solver)
        return solver;
      else
      {
        throw Hermes::Exceptions::LinearMatrixSolverException("Can not cast to DirectSolver.");
      }
    }

    template<typename Scalar>
    LoopSolver<Scalar>* LinearMatrixSolver<Scalar>::as_LoopSolver() const
    {
      LoopSolver<Scalar>* solver = dynamic_cast<LoopSolver<Scalar>*>(const_cast<LinearMatrixSolver<Scalar>*>(this));
      if(solver)
        return solver;
      else
      {
        throw Hermes::Exceptions::LinearMatrixSolverException("Can not cast to LoopSolver.");
      }
    }

    template<typename Scalar>
    IterSolver<Scalar>* LinearMatrixSolver<Scalar>::as_IterSolver() const
    {
      IterSolver<Scalar>* solver = dynamic_cast<IterSolver<Scalar>*>(const_cast<LinearMatrixSolver<Scalar>*>(this));
      if(solver)
        return solver;
      else
      {
        throw Hermes::Exceptions::LinearMatrixSolverException("Can not cast to IterSolver.");
      }
    }

    template<typename Scalar>
    AMGSolver<Scalar>* LinearMatrixSolver<Scalar>::as_AMGSolver() const
    {
      AMGSolver<Scalar>* solver = dynamic_cast<AMGSolver<Scalar>*>(const_cast<LinearMatrixSolver<Scalar>*>(this));
      if(solver)
        return solver;
      else
      {
        throw Hermes::Exceptions::LinearMatrixSolverException("Can not cast to AMGSolver.");
      }
    }


    template<typename Scalar>
    ExternalSolver<Scalar>* static_create_external_solver(CSCMatrix<Scalar> *m, SimpleVector<Scalar> *rhs)
    {
      throw Exceptions::MethodNotOverridenException("ExternalSolver<Scalar>::create_external_solver");
      return NULL;
    }

    template<>
    ExternalSolver<double>::creation ExternalSolver<double>::create_external_solver = static_create_external_solver<double>;
    template<>
    ExternalSolver<std::complex<double> >::creation ExternalSolver<std::complex<double> >::create_external_solver = static_create_external_solver<std::complex<double> >;

    template<>
    HERMES_API LinearMatrixSolver<double>* create_linear_solver(Matrix<double>* matrix, Vector<double>* rhs, bool use_direct_solver)
    {
      Vector<double>* rhs_dummy = NULL;
      switch (use_direct_solver ? Hermes::HermesCommonApi.get_integral_param_value(Hermes::directMatrixSolverType) : Hermes::HermesCommonApi.get_integral_param_value(Hermes::matrixSolverType))
      {
      case Hermes::SOLVER_EXTERNAL:
        {
          if(rhs != NULL)
            return ExternalSolver<double>::create_external_solver(static_cast<CSCMatrix<double>*>(matrix), static_cast<SimpleVector<double>*>(rhs));
          else
            return ExternalSolver<double>::create_external_solver(static_cast<CSCMatrix<double>*>(matrix), static_cast<SimpleVector<double>*>(rhs_dummy));
        }
      case Hermes::SOLVER_AZTECOO:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver AztecOO selected as a direct solver.");
#if defined HAVE_AZTECOO && defined HAVE_EPETRA
          if(rhs != NULL) return new AztecOOSolver<double>(static_cast<EpetraMatrix<double>*>(matrix), static_cast<EpetraVector<double>*>(rhs));
          else return new AztecOOSolver<double>(static_cast<EpetraMatrix<double>*>(matrix), static_cast<EpetraVector<double>*>(rhs_dummy));
#else
          throw Hermes::Exceptions::Exception("AztecOO not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_AMESOS:
        {
#if defined HAVE_AMESOS && defined HAVE_EPETRA
          if(rhs != NULL) return new AmesosSolver<double>("Amesos_Klu", static_cast<EpetraMatrix<double>*>(matrix), static_cast<EpetraVector<double>*>(rhs));
          else return new AmesosSolver<double>("Amesos_Klu", static_cast<EpetraMatrix<double>*>(matrix), static_cast<EpetraVector<double>*>(rhs_dummy));
#else
          throw Hermes::Exceptions::Exception("Amesos not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_MUMPS:
        {
#ifdef WITH_MUMPS
          if(rhs != NULL) return new MumpsSolver<double>(static_cast<MumpsMatrix<double>*>(matrix), static_cast<SimpleVector<double>*>(rhs));
          else return new MumpsSolver<double>(static_cast<MumpsMatrix<double>*>(matrix), static_cast<SimpleVector<double>*>(rhs_dummy));
#else
          throw Hermes::Exceptions::Exception("MUMPS was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_PETSC:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver PETSc selected as a direct solver.");
#ifdef WITH_PETSC
          if(rhs != NULL) return new PetscLinearMatrixSolver<double>(static_cast<PetscMatrix<double>*>(matrix), static_cast<PetscVector<double>*>(rhs));
          else return new PetscLinearMatrixSolver<double>(static_cast<PetscMatrix<double>*>(matrix), static_cast<PetscVector<double>*>(rhs_dummy));
#else
          throw Hermes::Exceptions::Exception("PETSc not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_UMFPACK:
        {
#ifdef WITH_UMFPACK
          if(rhs != NULL) return new UMFPackLinearMatrixSolver<double>(static_cast<CSCMatrix<double>*>(matrix), static_cast<SimpleVector<double>*>(rhs));
          else return new UMFPackLinearMatrixSolver<double>(static_cast<CSCMatrix<double>*>(matrix), static_cast<SimpleVector<double>*>(rhs_dummy));
#else
          throw Hermes::Exceptions::Exception("UMFPACK was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_PARALUTION_ITERATIVE:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver PARALUTION selected as a direct solver.");
#ifdef WITH_PARALUTION
          return new IterativeParalutionLinearMatrixSolver<double>(static_cast<ParalutionMatrix<double>*>(matrix), static_cast<ParalutionVector<double>*>(rhs));
#else
          throw Hermes::Exceptions::Exception("PARALUTION was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_PARALUTION_AMG:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The AMG solver PARALUTION selected as a direct solver.");
#ifdef WITH_PARALUTION
          return new AMGParalutionLinearMatrixSolver<double>(static_cast<ParalutionMatrix<double>*>(matrix), static_cast<ParalutionVector<double>*>(rhs));
#else
          throw Hermes::Exceptions::Exception("PARALUTION was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_SUPERLU:
        {
#ifdef WITH_SUPERLU
          if(rhs != NULL) return new SuperLUSolver<double>(static_cast<CSCMatrix<double>*>(matrix), static_cast<SimpleVector<double>*>(rhs));
          else return new SuperLUSolver<double>(static_cast<CSCMatrix<double>*>(matrix), static_cast<SimpleVector<double>*>(rhs_dummy));
#else
          throw Hermes::Exceptions::Exception("SuperLU was not installed.");
#endif
          break;
        }
      default:
        throw Hermes::Exceptions::Exception("Unknown matrix solver requested in create_linear_solver().");
      }
      return NULL;
    }

    template<>
    HERMES_API LinearMatrixSolver<std::complex<double> >* create_linear_solver(Matrix<std::complex<double> >* matrix, Vector<std::complex<double> >* rhs, bool use_direct_solver)
    {
      Vector<std::complex<double> >* rhs_dummy = NULL;
      switch (use_direct_solver ? Hermes::HermesCommonApi.get_integral_param_value(Hermes::directMatrixSolverType) : Hermes::HermesCommonApi.get_integral_param_value(Hermes::matrixSolverType))
      {
      case Hermes::SOLVER_EXTERNAL:
        {
          if(rhs != NULL)
            return ExternalSolver<std::complex<double> >::create_external_solver(static_cast<CSCMatrix<std::complex<double> >*>(matrix), static_cast<SimpleVector<std::complex<double> >*>(rhs));
          else
            return ExternalSolver<std::complex<double> >::create_external_solver(static_cast<CSCMatrix<std::complex<double> >*>(matrix), static_cast<SimpleVector<std::complex<double> >*>(rhs_dummy));
        }
      case Hermes::SOLVER_AZTECOO:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver AztecOO selected as a direct solver.");
#if defined HAVE_AZTECOO && defined HAVE_EPETRA
          if(rhs != NULL) return new AztecOOSolver<std::complex<double> >(static_cast<EpetraMatrix<std::complex<double> >*>(matrix), static_cast<EpetraVector<std::complex<double> >*>(rhs));
          else return new AztecOOSolver<std::complex<double> >(static_cast<EpetraMatrix<std::complex<double> >*>(matrix), static_cast<EpetraVector<std::complex<double> >*>(rhs_dummy));
#else
          throw Hermes::Exceptions::Exception("AztecOO not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_AMESOS:
        {
#if defined HAVE_AMESOS && defined HAVE_EPETRA
          if(rhs != NULL) return new AmesosSolver<std::complex<double> >("Amesos_Klu", static_cast<EpetraMatrix<std::complex<double> >*>(matrix), static_cast<EpetraVector<std::complex<double> >*>(rhs));
          else return new AmesosSolver<std::complex<double> >("Amesos_Klu", static_cast<EpetraMatrix<std::complex<double> >*>(matrix), static_cast<EpetraVector<std::complex<double> >*>(rhs_dummy));
#else
          throw Hermes::Exceptions::Exception("Amesos not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_MUMPS:
        {
#ifdef WITH_MUMPS
          if(rhs != NULL) return new MumpsSolver<std::complex<double> >(static_cast<MumpsMatrix<std::complex<double> >*>(matrix), static_cast<SimpleVector<std::complex<double> >*>(rhs));
          else return new MumpsSolver<std::complex<double> >(static_cast<MumpsMatrix<std::complex<double> >*>(matrix), static_cast<SimpleVector<std::complex<double> >*>(rhs_dummy));
#else
          throw Hermes::Exceptions::Exception("MUMPS was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_PETSC:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver PETSc selected as a direct solver.");
#ifdef WITH_PETSC
          if(rhs != NULL) return new PetscLinearMatrixSolver<std::complex<double> >(static_cast<PetscMatrix<std::complex<double> >*>(matrix), static_cast<PetscVector<std::complex<double> >*>(rhs));
          else return new PetscLinearMatrixSolver<std::complex<double> >(static_cast<PetscMatrix<std::complex<double> >*>(matrix), static_cast<PetscVector<std::complex<double> >*>(rhs_dummy));
#else
          throw Hermes::Exceptions::Exception("PETSc not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_UMFPACK:
        {
#ifdef WITH_UMFPACK
          if(rhs != NULL) return new UMFPackLinearMatrixSolver<std::complex<double> >(static_cast<CSCMatrix<std::complex<double> >*>(matrix), static_cast<SimpleVector<std::complex<double> >*>(rhs));
          else return new UMFPackLinearMatrixSolver<std::complex<double> >(static_cast<CSCMatrix<std::complex<double> >*>(matrix), static_cast<SimpleVector<std::complex<double> >*>(rhs_dummy));
#else
          throw Hermes::Exceptions::Exception("UMFPACK was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_PARALUTION_ITERATIVE:
      case Hermes::SOLVER_PARALUTION_AMG:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver PARALUTION selected as a direct solver.");
#ifdef WITH_PARALUTION
          throw Hermes::Exceptions::Exception("PARALUTION only works for real problems.");
#else
          throw Hermes::Exceptions::Exception("PARALUTION was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_SUPERLU:
        {
#ifdef WITH_SUPERLU
          if(rhs != NULL) return new SuperLUSolver<std::complex<double> >(static_cast<CSCMatrix<std::complex<double> >*>(matrix), static_cast<SimpleVector<std::complex<double> >*>(rhs));
          else return new SuperLUSolver<std::complex<double> >(static_cast<CSCMatrix<std::complex<double> >*>(matrix), static_cast<SimpleVector<std::complex<double> >*>(rhs_dummy));
#else
          throw Hermes::Exceptions::Exception("SuperLU was not installed.");
#endif
          break;
        }
      default:
        throw Hermes::Exceptions::Exception("Unknown matrix solver requested in create_linear_solver().");
      }
      return NULL;
    }

    template <typename Scalar>
    ExternalSolver<Scalar>::ExternalSolver(CSCMatrix<Scalar> *m, SimpleVector<Scalar> *rhs) : LinearMatrixSolver<Scalar>(HERMES_CREATE_STRUCTURE_FROM_SCRATCH), m(m), rhs(rhs)
    {
    }

    template <typename Scalar>
    SimpleExternalSolver<Scalar>::SimpleExternalSolver(CSCMatrix<Scalar> *m, SimpleVector<Scalar> *rhs) : ExternalSolver<Scalar>(m, rhs)
    {
    }

    template <typename Scalar>
    void SimpleExternalSolver<Scalar>::solve()
    {
      solve(NULL);
    }

    template <typename Scalar>
    void SimpleExternalSolver<Scalar>::solve(Scalar* initial_guess)
    {
      // Output.
      this->process_matrix_output(this->m);
      this->process_vector_output(this->rhs);

      // External process.
      std::string resultFileName = this->command();

      // Handling of the result file.
      this->sln = new Scalar[this->m->get_size()];
      SimpleVector<Scalar> temp;
      temp.alloc(this->m->get_size());
      temp.import_from_file((char*)resultFileName.c_str(), "x", EXPORT_FORMAT_PLAIN_ASCII );
      memcpy(this->sln, temp.v, this->m->get_size() * sizeof(Scalar));
    }

    template <typename Scalar>
    DirectSolver<Scalar>::DirectSolver(MatrixStructureReuseScheme reuse_scheme) : LinearMatrixSolver<Scalar>(reuse_scheme)
    {
    }

    template <typename Scalar>
    void DirectSolver<Scalar>::solve(Scalar* initial_guess)
    {
      this->solve();
    }

    template <typename Scalar>
    LoopSolver<Scalar>::LoopSolver(MatrixStructureReuseScheme reuse_scheme) : LinearMatrixSolver<Scalar>(reuse_scheme), max_iters(10000), tolerance(1e-8)
    {
    }

    template<typename Scalar>
    void LoopSolver<Scalar>::set_tolerance(double tol)
    {
      this->tolerance = tol;
      this->toleranceType = AbsoluteTolerance;
    }

    template<typename Scalar>
    void LoopSolver<Scalar>::set_tolerance(double tol, LoopSolverToleranceType toleranceType)
    {
      this->tolerance = tol;
      this->toleranceType = toleranceType;
    }

    template<typename Scalar>
    void LoopSolver<Scalar>::set_max_iters(int iters)
    {
      this->max_iters = iters;
    }

    template <typename Scalar>
    IterSolver<Scalar>::IterSolver(MatrixStructureReuseScheme reuse_scheme) : LoopSolver<Scalar>(reuse_scheme), precond_yes(false), iterSolverType(CG)
    {
    }

    template<typename Scalar>
    void IterSolver<Scalar>::set_solver_type(IterSolverType iterSolverType)
    {
      this->iterSolverType = iterSolverType;
    }

    template <typename Scalar>
    AMGSolver<Scalar>::AMGSolver(MatrixStructureReuseScheme reuse_scheme) : LoopSolver<Scalar>(reuse_scheme), smootherSolverType(CG), smootherPreconditionerType(MultiColoredSGS)
    {
    }

    template<typename Scalar>
    void AMGSolver<Scalar>::set_smoother(IterSolverType solverType_, PreconditionerType preconditionerType_)
    {
      this->smootherPreconditionerType = preconditionerType_;
      this->smootherSolverType = solverType_;
    }

    template class HERMES_API LinearMatrixSolver<double>;
    template class HERMES_API LinearMatrixSolver<std::complex<double> >;
    template class HERMES_API SimpleExternalSolver<double>;
    template class HERMES_API SimpleExternalSolver<std::complex<double> >;
    template class HERMES_API ExternalSolver<double>;
    template class HERMES_API ExternalSolver<std::complex<double> >;
    template class HERMES_API DirectSolver<double>;
    template class HERMES_API DirectSolver<std::complex<double> >;
    template class HERMES_API LoopSolver<double>;
    template class HERMES_API LoopSolver<std::complex<double> >;
    template class HERMES_API IterSolver<double>;
    template class HERMES_API IterSolver<std::complex<double> >;
    template class HERMES_API AMGSolver<double>;
    template class HERMES_API AMGSolver<std::complex<double> >;
  }
}