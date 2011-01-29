// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes3D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes3D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes3D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#include "amesos.h"
#include "../callstack.h"

#ifdef HAVE_AMESOS
  #include <Amesos_ConfigDefs.h>
#endif

#ifdef HAVE_AMESOS
  Amesos AmesosSolver::factory;
#endif

// Amesos solver ///////////////////////////////////////////////////////////////////////////////////

AmesosSolver::AmesosSolver(const char *solver_type, EpetraMatrix *m, EpetraVector *rhs)
  : LinearSolver(HERMES_FACTORIZE_FROM_SCRATCH), m(m), rhs(rhs)
{
  _F_
#ifdef HAVE_AMESOS
  solver = factory.Create(solver_type, problem);
  assert(solver != NULL);
  // WARNING: Amesos does not use RCP to allocate the Amesos_BaseSolver, 
  //          so don't forget to delete it!
  //          ( Amesos.cpp, line 88, called from factory.Create(): 
  //            return new Amesos_Klu(LinearProblem); )
#else
  error(AMESOS_NOT_COMPILED);
#endif
}

AmesosSolver::~AmesosSolver()
{
  _F_
#ifdef HAVE_AMESOS
  delete solver;
#endif
}

bool AmesosSolver::is_available(const char *name)
{
  _F_
#ifdef HAVE_AMESOS
  return factory.Query(name);
#else
  return false;
#endif
}

void AmesosSolver::set_use_transpose(bool use_transpose)
{
  _F_
#ifdef HAVE_AMESOS
  solver->SetUseTranspose(use_transpose);
#endif
}

bool AmesosSolver::use_transpose()
{
  _F_
#ifdef HAVE_AMESOS
  return solver->UseTranspose();
#else
  return false;
#endif
}

bool AmesosSolver::solve()
{
  _F_
#ifdef HAVE_AMESOS
  assert(m != NULL);
  assert(rhs != NULL);
  
  assert(m->size == rhs->size);
  
  TimePeriod tmr;  

#if defined(H1D_COMPLEX) || defined(H2D_COMPLEX) || defined (H3D_COMPLEX)
  error("AmesosSolver::solve() not yet implemented for complex problems");
#else
  problem.SetOperator(m->mat);
  problem.SetRHS(rhs->vec);
  Epetra_Vector x(*rhs->std_map);
  problem.SetLHS(&x);
#endif

  if (!setup_factorization())
  {
    warning("AmesosSolver: LU factorization could not be completed");
    return false;
  }

  int status = solver->Solve();
  if (status != 0) 
  {
    error("AmesosSolver: Solution failed.");
    return false;
  }
  
  tmr.tick();
  time = tmr.accumulated();

  delete [] sln;
  sln = new scalar[m->size]; MEM_CHECK(sln);
  // copy the solution into sln vector
  memset(sln, 0, m->size * sizeof(scalar));
  
#if defined(H1D_COMPLEX) || defined(H2D_COMPLEX) || defined (H3D_COMPLEX)
#else 
  for (unsigned int i = 0; i < m->size; i++) sln[i] = x[i];
#endif

  return true;
#else
  return false;
#endif
}

bool AmesosSolver::setup_factorization()
{
  _F_
#ifdef HAVE_AMESOS
  // Perform both factorization phases for the first time.
  int eff_fact_scheme;
  if (factorization_scheme != HERMES_FACTORIZE_FROM_SCRATCH && 
      solver->NumSymbolicFact() == 0 && solver->NumNumericFact() == 0)
    eff_fact_scheme = HERMES_FACTORIZE_FROM_SCRATCH;
  else
    eff_fact_scheme = factorization_scheme;
  
  int status;
  switch(eff_fact_scheme)
  {
    case HERMES_FACTORIZE_FROM_SCRATCH:
      //debug_log("Factorizing symbolically.");
      status = solver->SymbolicFactorization();
      if (status != 0)
      {
        warning("Symbolic factorization failed.");
        return false;
      }
      
    case HERMES_REUSE_MATRIX_REORDERING:
    case HERMES_REUSE_MATRIX_REORDERING_AND_SCALING:
      status = solver->NumericFactorization();
      if (status != 0) 
      {
        warning("Numeric factorization failed.");
        return false;
      }
  }
  
  return true;
#else
  return false;
#endif
}