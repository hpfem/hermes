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
#include "linear_solver.h"
#include "umfpack_solver.h"
#include "superlu_solver.h"
#include "amesos_solver.h"
#include "petsc_solver.h"
#include "mumps_solver.h"
#include "nox_solver.h"
#include "aztecoo_solver.h"

using namespace Hermes::Algebra;

namespace Hermes
{
  namespace Solvers 
  {
    template<typename Scalar>
    LinearSolver<Scalar>* create_linear_solver(Hermes::MatrixSolverType matrix_solver_type, Matrix<Scalar>* matrix, Vector<Scalar>* rhs)
    {
      _F_;
      Vector<Scalar>* rhs_dummy = NULL;
      switch (matrix_solver_type) 
      {
      case Hermes::SOLVER_AZTECOO:
        {
#if defined HAVE_AZTECOO && defined HAVE_EPETRA
          info("Using AztecOO."); 
          if (rhs != NULL) return new AztecOOSolver<Scalar>(static_cast<EpetraMatrix<Scalar>*>(matrix), static_cast<EpetraVector<Scalar>*>(rhs));
          else return new AztecOOSolver<Scalar>(static_cast<EpetraMatrix<Scalar>*>(matrix), static_cast<EpetraVector<Scalar>*>(rhs_dummy));
#else
          error("AztecOO not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_AMESOS:
        {
#if defined HAVE_AMESOS && defined HAVE_EPETRA
          info("Using Amesos.");         
          if (rhs != NULL) return new AmesosSolver<Scalar>("Amesos_Klu", static_cast<EpetraMatrix<Scalar>*>(matrix), static_cast<EpetraVector<Scalar>*>(rhs));
          else return new AmesosSolver<Scalar>("Amesos_Klu", static_cast<EpetraMatrix<Scalar>*>(matrix), static_cast<EpetraVector<Scalar>*>(rhs_dummy));
#else
          error("Amesos not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_MUMPS: 
        {
#ifdef WITH_MUMPS
          info("Using Mumps.");         
          if (rhs != NULL) return new MumpsSolver<Scalar>(static_cast<MumpsMatrix<Scalar>*>(matrix), static_cast<MumpsVector<Scalar>*>(rhs)); 
          else return new MumpsSolver<Scalar>(static_cast<MumpsMatrix<Scalar>*>(matrix), static_cast<MumpsVector<Scalar>*>(rhs_dummy)); 
#else
          error("MUMPS was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_PETSC: 
        {
#ifdef WITH_PETSC
          info("Using PETSc.");        
          if (rhs != NULL) return new PetscLinearSolver<Scalar>(static_cast<PetscMatrix<Scalar>*>(matrix), static_cast<PetscVector<Scalar>*>(rhs)); 
          else return new PetscLinearSolver<Scalar>(static_cast<PetscMatrix<Scalar>*>(matrix), static_cast<PetscVector<Scalar>*>(rhs_dummy)); 
#else
          error("PETSc not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_UMFPACK: 
        {
#ifdef WITH_UMFPACK
          info("Using UMFPack.");
          if (rhs != NULL) return new UMFPackLinearSolver<Scalar>(static_cast<UMFPackMatrix<Scalar>*>(matrix), static_cast<UMFPackVector<Scalar>*>(rhs)); 
          else return new UMFPackLinearSolver<Scalar>(static_cast<UMFPackMatrix<Scalar>*>(matrix), static_cast<UMFPackVector<Scalar>*>(rhs_dummy));  
#else
          error("UMFPACK was not installed.");
#endif
          break;
        }
      case Hermes::SOLVER_SUPERLU: 
        {
#ifdef WITH_SUPERLU
          info("Using SuperLU.");       
          if (rhs != NULL) return new SuperLUSolver<Scalar>(static_cast<SuperLUMatrix<Scalar>*>(matrix), static_cast<SuperLUVector<Scalar>*>(rhs)); 
          else return new SuperLUSolver<Scalar>(static_cast<SuperLUMatrix<Scalar>*>(matrix), static_cast<SuperLUVector<Scalar>*>(rhs_dummy)); 
#else
          error("SuperLU was not installed.");
#endif
          break;
        }
      default: 
        error("Unknown matrix solver requested.");
      }
      return NULL;
    }

    template HERMES_API LinearSolver<double>*  create_linear_solver(Hermes::MatrixSolverType matrix_solver, 
      Matrix<double>* matrix, Vector<double>* rhs);

    template HERMES_API LinearSolver<std::complex<double> >*  create_linear_solver(Hermes::MatrixSolverType matrix_solver, 
      Matrix<std::complex<double> >* matrix, Vector<std::complex<double> >* rhs);

    template class HERMES_API LinearSolver<double>;
    template class HERMES_API LinearSolver<std::complex<double> >;
    template class HERMES_API DirectSolver<double>;
    template class HERMES_API DirectSolver<std::complex<double> >;
    template class HERMES_API IterSolver<double>;
    template class HERMES_API IterSolver<std::complex<double> >;
  }
}
