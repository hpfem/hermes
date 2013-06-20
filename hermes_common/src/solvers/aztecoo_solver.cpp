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
/*! \file aztecoo_solver.cpp
\brief AztecOOSolver class as an interface to AztecOO.
*/
#include "config.h"
#ifdef HAVE_AZTECOO
#include "aztecoo_solver.h"
#include "callstack.h"
#include "Komplex_LinearProblem.h"

namespace Hermes
{
  namespace Solvers
  {
    template<typename Scalar>
    AztecOOSolver<Scalar>::AztecOOSolver(EpetraMatrix<Scalar> *m, EpetraVector<Scalar> *rhs)
      : IterSolver<Scalar>(), m(m), rhs(rhs)
    {
      pc = NULL;
    }

    template<typename Scalar>
    AztecOOSolver<Scalar>::~AztecOOSolver()
    {
      }

    template<typename Scalar>
    void AztecOOSolver<Scalar>::set_solver(const char *name)
    {
      int az_solver;
      if(name)
      {
        if(strcasecmp(name, "gmres") == 0) az_solver = AZ_gmres;
        else if(strcasecmp(name, "cg") == 0) az_solver = AZ_cg;
        else if(strcasecmp(name, "cgs") == 0) az_solver = AZ_cgs;
        else if(strcasecmp(name, "tfqmr") == 0) az_solver = AZ_tfqmr;
        else if(strcasecmp(name, "bicgstab") == 0) az_solver = AZ_bicgstab;
        else az_solver = AZ_gmres;
      }
      else
      {
        az_solver = AZ_gmres;
      }

      aztec.SetAztecOption(AZ_solver, az_solver);
    }

    template<typename Scalar>
    void AztecOOSolver<Scalar>::set_tolerance(double tol)
    {
      this->tolerance = tol;
    }

    template<typename Scalar>
    void AztecOOSolver<Scalar>::set_max_iters(int iters)
    {
      this->max_iters = iters;
    }

    template<typename Scalar>
    void AztecOOSolver<Scalar>::set_precond(Precond<Scalar> *pc)
    {
      this->pc = dynamic_cast<EpetraPrecond<Scalar>*>(pc);
      if (this->pc)
        this->precond_yes = true;
      // TODO: else warn that a wrong type of preconditioner has been used.
    }

    template<typename Scalar>
    void AztecOOSolver<Scalar>::set_precond(const char *name)
    {
      int az_precond;
      if(name)
      {
        if(strcasecmp(name, "none") == 0)
          az_precond = AZ_none;
        else if(strcasecmp(name, "jacobi") == 0)
          az_precond = AZ_Jacobi;
        else if(strcasecmp(name, "neumann") == 0)
          az_precond = AZ_Neumann;
        else if(strcasecmp(name, "least-squares") == 0)
          az_precond = AZ_ls;
        else
          az_precond = AZ_none;
      }
      else
        az_precond = AZ_none; //asi by to melo byt nastaveno

      this->precond_yes = (az_precond != AZ_none);
      aztec.SetAztecOption(AZ_precond, az_precond);
    }

    template<typename Scalar>
    void AztecOOSolver<Scalar>::set_option(int option, int value)
    {
      aztec.SetAztecOption(option, value);
    }

    template<typename Scalar>
    void AztecOOSolver<Scalar>::set_param(int param, double value)
    {
      aztec.SetAztecParam(param, value);
    }

    template<typename Scalar>
    int AztecOOSolver<Scalar>::get_matrix_size()
    {
      return m->size;
    }
    
    template<>
    bool AztecOOSolver<double>::solve()
    {
      assert(m != NULL);
      assert(rhs != NULL);
      assert(m->size == rhs->size);

      // no output
      aztec.SetAztecOption(AZ_output, AZ_none);  // AZ_all | AZ_warnings | AZ_last | AZ_summary

      // setup the problem
      aztec.SetUserMatrix(m->mat);
      aztec.SetRHS(rhs->vec);
      Epetra_Vector x(*rhs->std_map);
     
      aztec.SetLHS(&x);

      if(pc != NULL)
      {
        Epetra_Operator *op = pc->get_obj();
        assert(op != NULL);    // can work only with Epetra_Operators
        aztec.SetPrecOperator(op);
      }

      // solve it
      aztec.Iterate(this->max_iters, this->tolerance);

      this->tick();
      this->time = this->accumulated();

      delete [] this->sln;
      this->sln = new double[m->size];
      memset(this->sln, 0, m->size * sizeof(double));

      // copy the solution into sln vector
      for (unsigned int i = 0; i < m->size; i++) this->sln[i] = x[i];
      return true;
    }


    template<>
    bool AztecOOSolver<double>::solve(double *initial_guess)
    {
      assert(m != NULL);
      assert(rhs != NULL);
      assert(m->size == rhs->size);

      // no output
      aztec.SetAztecOption(AZ_output, AZ_none);  // AZ_all | AZ_warnings | AZ_last | AZ_summary

      // setup the problem
      aztec.SetUserMatrix(m->mat);
      aztec.SetRHS(rhs->vec);
      Epetra_Vector x(*rhs->std_map);
      
      if (initial_guess)
        for (unsigned int i = 0; i < m->size; i++)
          x[i] = initial_guess[i];
      
      aztec.SetLHS(&x);

      if(pc != NULL)
      {
        Epetra_Operator *op = pc->get_obj();
        assert(op != NULL);    // can work only with Epetra_Operators
        aztec.SetPrecOperator(op);
      }

      // solve it
      aztec.Iterate(this->max_iters, this->tolerance);

      this->tick();
      this->time = this->accumulated();

      delete [] this->sln;
      this->sln = new double[m->size];
      memset(this->sln, 0, m->size * sizeof(double));

      // copy the solution into sln vector
      for (unsigned int i = 0; i < m->size; i++) this->sln[i] = x[i];
      return true;
    }

    template<>
    bool AztecOOSolver<std::complex<double> >::solve()
    {
      assert(m != NULL);
      assert(rhs != NULL);
      assert(m->size == rhs->size);

      // no output
      aztec.SetAztecOption(AZ_output, AZ_none);  // AZ_all | AZ_warnings | AZ_last | AZ_summary

      double c0r = 1.0, c0i = 0.0;
      double c1r = 0.0, c1i = 1.0;

      Epetra_Vector xr(*rhs->std_map);
      Epetra_Vector xi(*rhs->std_map);

      Komplex_LinearProblem kp(c0r, c0i, *m->mat, c1r, c1i, *m->mat_im, xr, xi, *rhs->vec, *rhs->vec_im);
      Epetra_LinearProblem *lp = kp.KomplexProblem();
      aztec.SetProblem(*lp,true);

      // solve it
      aztec.Iterate(this->max_iters, this->tolerance);

      kp.ExtractSolution(xr, xi);

      delete [] this->sln;
      this->sln = new std::complex<double>[m->size];
      memset(this->sln, 0, m->size * sizeof(std::complex<double>));

      // copy the solution into sln vector
      for (unsigned int i = 0; i < m->size; i++) this->sln[i] = std::complex<double>(xr[i], xi[i]);
      return true;
    }
    
    template<>
    bool AztecOOSolver<std::complex<double> >::solve(std::complex<double>* initial_guess)
    {
      assert(m != NULL);
      assert(rhs != NULL);
      assert(m->size == rhs->size);

      // no output
      aztec.SetAztecOption(AZ_output, AZ_none);  // AZ_all | AZ_warnings | AZ_last | AZ_summary

      double c0r = 1.0, c0i = 0.0;
      double c1r = 0.0, c1i = 1.0;

      Epetra_Vector xr(*rhs->std_map);
      Epetra_Vector xi(*rhs->std_map);
      
      if (initial_guess)
      {
        for (unsigned int i = 0; i < m->size; i++)
        {
          xr[i] = initial_guess[i].real();
          xi[i] = initial_guess[i].imag();
        }
      }

      Komplex_LinearProblem kp(c0r, c0i, *m->mat, c1r, c1i, *m->mat_im, xr, xi, *rhs->vec, *rhs->vec_im);
      Epetra_LinearProblem *lp = kp.KomplexProblem();
      aztec.SetProblem(*lp,true);

      // solve it
      aztec.Iterate(this->max_iters, this->tolerance);

      kp.ExtractSolution(xr, xi);

      delete [] this->sln;
      this->sln = new std::complex<double>[m->size];
      memset(this->sln, 0, m->size * sizeof(std::complex<double>));

      // copy the solution into sln vector
      for (unsigned int i = 0; i < m->size; i++) this->sln[i] = std::complex<double>(xr[i], xi[i]);
      return true;
    }

    template<typename Scalar>
    int AztecOOSolver<Scalar>::get_num_iters()
    {
      return aztec.NumIters();
    }

    template<typename Scalar>
    double AztecOOSolver<Scalar>::get_residual()
    {
      return aztec.TrueResidual();
    }

    template class HERMES_API AztecOOSolver<double>;
    template class HERMES_API AztecOOSolver<std::complex<double> >;
  }
}
#endif