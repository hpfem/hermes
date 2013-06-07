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

using namespace Hermes::Algebra;

namespace Hermes
{
  namespace Solvers
  {
    template<typename Scalar>
    LinearMatrixSolver<Scalar>::LinearMatrixSolver()
    {
      sln = NULL;
      time = -1.0;
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
    void LinearMatrixSolver<Scalar>::set_factorization_scheme()
    {
      set_factorization_scheme(HERMES_FACTORIZE_FROM_SCRATCH);
    }

    template<>
    HERMES_API LinearMatrixSolver<double>* create_linear_solver(Matrix<double>* matrix, Vector<double>* rhs, bool use_direct_solver)
    {
      Vector<double>* rhs_dummy = NULL;
      switch (use_direct_solver ? Hermes::HermesCommonApi.get_integral_param_value(Hermes::directMatrixSolverType) : Hermes::HermesCommonApi.get_integral_param_value(Hermes::matrixSolverType))
      {
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
          if(rhs != NULL) return new MumpsSolver<double>(static_cast<MumpsMatrix<double>*>(matrix), static_cast<MumpsVector<double>*>(rhs));
          else return new MumpsSolver<double>(static_cast<MumpsMatrix<double>*>(matrix), static_cast<MumpsVector<double>*>(rhs_dummy));
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
          if(rhs != NULL) return new UMFPackLinearMatrixSolver<double>(static_cast<UMFPackMatrix<double>*>(matrix), static_cast<UMFPackVector<double>*>(rhs));
          else return new UMFPackLinearMatrixSolver<double>(static_cast<UMFPackMatrix<double>*>(matrix), static_cast<UMFPackVector<double>*>(rhs_dummy));
#else
          throw Hermes::Exceptions::Exception("UMFPACK was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_PARALUTION:
        {
          if(use_direct_solver)
            throw Hermes::Exceptions::Exception("The iterative solver PARALUTION selected as a direct solver.");
#ifdef WITH_PARALUTION
          return new ParalutionLinearMatrixSolver<double>(static_cast<ParalutionMatrix<double>*>(matrix), static_cast<ParalutionVector<double>*>(rhs));
#else
          throw Hermes::Exceptions::Exception("PARALUTION was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_SUPERLU:
        {
#ifdef WITH_SUPERLU
          if(rhs != NULL) return new SuperLUSolver<double>(static_cast<SuperLUMatrix<double>*>(matrix), static_cast<SuperLUVector<double>*>(rhs));
          else return new SuperLUSolver<double>(static_cast<SuperLUMatrix<double>*>(matrix), static_cast<SuperLUVector<double>*>(rhs_dummy));
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
          if(rhs != NULL) return new MumpsSolver<std::complex<double> >(static_cast<MumpsMatrix<std::complex<double> >*>(matrix), static_cast<MumpsVector<std::complex<double> >*>(rhs));
          else return new MumpsSolver<std::complex<double> >(static_cast<MumpsMatrix<std::complex<double> >*>(matrix), static_cast<MumpsVector<std::complex<double> >*>(rhs_dummy));
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
          if(rhs != NULL) return new UMFPackLinearMatrixSolver<std::complex<double> >(static_cast<UMFPackMatrix<std::complex<double> >*>(matrix), static_cast<UMFPackVector<std::complex<double> >*>(rhs));
          else return new UMFPackLinearMatrixSolver<std::complex<double> >(static_cast<UMFPackMatrix<std::complex<double> >*>(matrix), static_cast<UMFPackVector<std::complex<double> >*>(rhs_dummy));
#else
          throw Hermes::Exceptions::Exception("UMFPACK was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_PARALUTION:
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
          if(rhs != NULL) return new SuperLUSolver<std::complex<double> >(static_cast<SuperLUMatrix<std::complex<double> >*>(matrix), static_cast<SuperLUVector<std::complex<double> >*>(rhs));
          else return new SuperLUSolver<std::complex<double> >(static_cast<SuperLUMatrix<std::complex<double> >*>(matrix), static_cast<SuperLUVector<std::complex<double> >*>(rhs_dummy));
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

    template<typename Scalar>
    void DirectSolver<Scalar>::set_factorization_scheme(FactorizationScheme reuse_scheme)
    {
      factorization_scheme = reuse_scheme;
    }

    template<typename Scalar>
    void IterSolver<Scalar>::set_tolerance(double tol)
    {
      this->tolerance = tol;
      this->toleranceType = AbsoluteTolerance;
    }

    template<typename Scalar>
    void IterSolver<Scalar>::set_tolerance(double tol, typename IterSolver<Scalar>::ToleranceType toleranceType)
    {
      this->tolerance = tol;
      this->toleranceType = toleranceType;
    }

    template<typename Scalar>
    void IterSolver<Scalar>::set_max_iters(int iters)
    {
      this->max_iters = iters;
    }

    template class HERMES_API LinearMatrixSolver<double>;
    template class HERMES_API LinearMatrixSolver<std::complex<double> >;
    template class HERMES_API DirectSolver<double>;
    template class HERMES_API DirectSolver<std::complex<double> >;
    template class HERMES_API IterSolver<double>;
    template class HERMES_API IterSolver<std::complex<double> >;
  }
}