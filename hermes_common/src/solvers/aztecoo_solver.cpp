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
      : IterSolver<Scalar>(), m(m), rhs(rhs), final_matrix(NULL), P(NULL), Q(NULL), row_perm(NULL), col_perm(NULL)
    {
      pc = NULL;
    }

    template<typename Scalar>
    AztecOOSolver<Scalar>::~AztecOOSolver()
    {
      free_permutation_data();
    }
    
    template<typename Scalar>
    void AztecOOSolver<Scalar>::free_permutation_data()
    {
      if (row_perm)
      {
        assert(col_perm);
        assert(Q);
        assert(final_matrix);
      
        delete row_perm;
        delete col_perm;
        delete P;
        delete Q;
        delete final_matrix;
      }
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
    void AztecOOSolver<Scalar>::set_reuse_scheme(MatrixStructureReuseScheme reuse_scheme)
    {
      LinearMatrixSolver<Scalar>::set_reuse_scheme(reuse_scheme);
      /*
      if (reuse_scheme == HERMES_REUSE_MATRIX_STRUCTURE_COMPLETELY || reuse_scheme == HERMES_REUSE_MATRIX_REORDERING_AND_SCALING)
        aztec.SetAztecOption(AZ_pre_calc, AZ_reuse);
      else if (reuse_scheme == HERMES_REUSE_MATRIX_REORDERING)
        aztec.SetAztecOption(AZ_pre_calc, AZ_recalc);
      else 
        aztec.SetAztecOption(AZ_pre_calc, AZ_calc);
      */
    }

    template<typename Scalar>
    int AztecOOSolver<Scalar>::get_matrix_size()
    {
      return m->size;
    }
    
    template<typename Scalar>
    void AztecOOSolver<Scalar>::use_node_wise_ordering(unsigned int num_pdes) 
    { 
      LinearMatrixSolver<Scalar>::use_node_wise_ordering(num_pdes);
      this->free_permutation_data();
    }
    
    template<typename Scalar>
    void AztecOOSolver<Scalar>::use_equations_wise_ordering() 
    { 
      LinearMatrixSolver<Scalar>::use_equations_wise_ordering();
      this->free_permutation_data();
    }
    
    template<typename Scalar>
    void AztecOOSolver<Scalar>::create_permutation_vectors()
    {
      int ndof = m->size;
      int ndof_per_pde = ndof/this->n_eq;
      int j = 0, jc = 0;
      int k = 0, kc = 0;
      
      this->row_perm = new int [ndof];
      this->col_perm = new int [ndof];
      
      for (int i = 0; i < ndof; i++)
      {
        this->row_perm[i] = j%ndof + jc;
        this->col_perm[i] = k%ndof + kc;
        
        j += ndof_per_pde;
        if (j == ndof)
        {
          jc++;
          j = 0;
        }
        
        k += this->n_eq;
        if (k == ndof)
        {
          kc++;
          k = 0;
        }
      }
    }
    
    template<>
    void AztecOOSolver<double>::solve()
    {
      assert(m != NULL);
      assert(rhs != NULL);
      assert(m->size == rhs->size);

      // no output
      aztec.SetAztecOption(AZ_output, AZ_summary);  // AZ_all | AZ_warnings | AZ_last | AZ_summary

      // setup the problem
      if (reuse_scheme == HERMES_CREATE_STRUCTURE_FROM_SCRATCH)
      //if (aztec.GetAztecOption(AZ_pre_calc) == AZ_calc)
      {
        if (node_wise_ordering)
        {
          if (row_perm)
            free_permutation_data();
          
          create_permutation_vectors();
          
          // NOTE: RowMap() == RangeMap() == ColMap() == DomainMap()
          P = new EpetraExt::Permutation<Epetra_CrsMatrix>(Copy, m->mat->RowMap(), row_perm);
          Q = new EpetraExt::Permutation<Epetra_CrsMatrix>(Copy, m->mat->ColMap(), col_perm);          
          
          final_matrix = new EpetraMatrix<double>((*P)((*Q)(*m->mat,true)));
          
          // NOTE: According to Trilinos docs, final_matrix should be fill_completed by now. However, when row permutation is performed, 
          // the input matrix is not finalized (unlike the case of column permutation) - this may possibly be a bug in Trilinos.
          // Also, when doing a row permutation, all maps of the source matrix are set to a temporary permutation
          // map in order to export the values from the source matrix to the target matrix in a permuted way. Then,
          // according to code documentation in EpetraExt_Permutation_impl.h, the row indexing for the new permuted
          // matrix (specified by the target matrix' RowMap_) is set to the original indexing (the original matrix' 
          // RowMap_). However, the other maps of the new matrix are not reset (and still point to the temporary
          // permutation map), which causes another set of problems. Hence, replaceMap() implementation in EpetraExt_Permutation_impl.h
          // has to be changed to call FillComplete with correct maps as input (calling FillComplete here doesn't suffice, since it
          // does not prevent creating a non-trivial Exporter, which is the root cause of the problems (and which cannot be deleted
          // from outside.)
        }
        else
          final_matrix = m;
        
        aztec.SetUserMatrix(final_matrix->mat);
      }
            
      EpetraExt::Permutation<Epetra_MultiVector> Pv(Copy, rhs->vec->Map(), row_perm);
      Epetra_Vector *final_rhs; 
      if (row_perm)
        final_rhs = Pv(*rhs->vec)(0);
      else
        final_rhs = rhs->vec;

      aztec.SetRHS(final_rhs);

      Epetra_Vector x(final_matrix->mat->DomainMap());
      aztec.SetLHS(&x);

      if(pc != NULL)
      {
        if(reuse_scheme == HERMES_CREATE_STRUCTURE_FROM_SCRATCH)
        //if(aztec.GetAztecOption(AZ_pre_calc) == AZ_calc)
        {
          pc->create(final_matrix); 
          pc->compute();
          aztec.SetPrecOperator(pc->get_obj());
        }
        else if (reuse_scheme == HERMES_REUSE_MATRIX_REORDERING || reuse_scheme == HERMES_REUSE_MATRIX_REORDERING_AND_SCALING)
        //else if (aztec.GetAztecOption(AZ_pre_calc) == AZ_recalc)
        {
          pc->recompute();
        }
        else
        {
          assert(pc->get_obj());
#ifdef HAVE_ML
          ML_Epetra::MultiLevelPreconditioner* op_ml = dynamic_cast<ML_Epetra::MultiLevelPreconditioner*>(pc->get_obj());
          assert(op_ml->IsPreconditionerComputed());
#endif            
        }
      }
      
      // solve it
      aztec.Iterate(this->max_iters, this->tolerance);

      this->tick();
      this->time = this->accumulated();

      delete [] this->sln;
      this->sln = new double[final_matrix->size];
      memset(this->sln, 0, final_matrix->size * sizeof(double));

      // copy the solution into sln vector
      if (col_perm)
        for (unsigned int i = 0; i < final_matrix->size; i++) this->sln[i] = x[col_perm[i]];
      else
        for (unsigned int i = 0; i < final_matrix->size; i++) this->sln[i] = x[i];                
    }


    template<>
    void AztecOOSolver<double>::solve(double *initial_guess)
    {
      assert(m != NULL);
      assert(rhs != NULL);
      assert(m->size == rhs->size);

      // no output
      aztec.SetAztecOption(AZ_output, AZ_summary);  // AZ_all | AZ_warnings | AZ_last | AZ_summary
      
      if (this->get_verbose_output())
        aztec.SetAztecOption(AZ_output, AZ_all);

      // setup the problem
      if (reuse_scheme == HERMES_CREATE_STRUCTURE_FROM_SCRATCH)
      //if (aztec.GetAztecOption(AZ_pre_calc) == AZ_calc)
      {
        if (node_wise_ordering)
        {
          if (row_perm)
            free_permutation_data();
          
          create_permutation_vectors();
          
          // NOTE: RowMap() == RangeMap() == ColMap() == DomainMap()
          P = new EpetraExt::Permutation<Epetra_CrsMatrix>(Copy, m->mat->RowMap(), row_perm);
          Q = new EpetraExt::Permutation<Epetra_CrsMatrix>(Copy, m->mat->ColMap(), col_perm);          
          
          final_matrix = new EpetraMatrix<double>((*P)((*Q)(*m->mat,true)));
          
          // NOTE: According to Trilinos docs, final_matrix should be fill_completed by now. However, when row permutation is performed, 
          // the input matrix is not finalized (unlike the case of column permutation) - this may possibly be a bug in Trilinos.
          // Also, when doing a row permutation, all maps of the source matrix are set to a temporary permutation
          // map in order to export the values from the source matrix to the target matrix in a permuted way. Then,
          // according to code documentation in EpetraExt_Permutation_impl.h, the row indexing for the new permuted
          // matrix (specified by the target matrix' RowMap_) is set to the original indexing (the original matrix' 
          // RowMap_). However, the other maps of the new matrix are not reset (and still point to the temporary
          // permutation map), which causes another set of problems. Hence, replaceMap() implementation in EpetraExt_Permutation_impl.h
          // has to be changed to call FillComplete with correct maps as input (calling FillComplete here doesn't suffice, since it
          // does not prevent creating a non-trivial Exporter, which is the root cause of the problems (and which cannot be deleted
          // from outside.)
        }
        else
          final_matrix = m;
        
        aztec.SetUserMatrix(final_matrix->mat);
      }
            
      EpetraExt::Permutation<Epetra_MultiVector> Pv(Copy, rhs->vec->Map(), row_perm);
      Epetra_Vector *final_rhs; 
      if (row_perm)
        final_rhs = Pv(*rhs->vec)(0);
      else
        final_rhs = rhs->vec;

      aztec.SetRHS(final_rhs);

      Epetra_Vector x(final_matrix->mat->DomainMap());
      
      if (initial_guess)
      {
        if (row_perm)
          for (unsigned int i = 0; i < m->size; i++)
            x[i] = initial_guess[row_perm[i]];
        else
          for (unsigned int i = 0; i < m->size; i++)
            x[i] = initial_guess[i];
      }
      
      aztec.SetLHS(&x);
      
      if(pc != NULL)
      {
        if(reuse_scheme == HERMES_CREATE_STRUCTURE_FROM_SCRATCH)
        //if(aztec.GetAztecOption(AZ_pre_calc) == AZ_calc)
        {
          pc->create(final_matrix); 
          pc->compute();
          aztec.SetPrecOperator(pc->get_obj());
        }
        else if (reuse_scheme == HERMES_REUSE_MATRIX_REORDERING || reuse_scheme == HERMES_REUSE_MATRIX_REORDERING_AND_SCALING)
        //else if (aztec.GetAztecOption(AZ_pre_calc) == AZ_recalc)
        {
          pc->recompute();
        }
        else
        {
          assert(pc->get_obj());
#ifdef HAVE_ML
          ML_Epetra::MultiLevelPreconditioner* op_ml = dynamic_cast<ML_Epetra::MultiLevelPreconditioner*>(pc->get_obj());
          assert(op_ml->IsPreconditionerComputed());
#endif            
        }
      }
      
      // solve it
      aztec.Iterate(this->max_iters, this->tolerance);

      this->tick();
      this->time = this->accumulated();

      delete [] this->sln;
      this->sln = new double[final_matrix->size];
      memset(this->sln, 0, final_matrix->size * sizeof(double));
      
      // copy the solution into sln vector
      if (col_perm)
        for (unsigned int i = 0; i < final_matrix->size; i++) this->sln[i] = x[col_perm[i]];
      else
        for (unsigned int i = 0; i < final_matrix->size; i++) this->sln[i] = x[i];  
    }

    template<>
    void AztecOOSolver<std::complex<double> >::solve()
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
      for (unsigned int i = 0; i < m->size; i++)
        this->sln[i] = std::complex<double>(xr[i], xi[i]);
    }
    
    template<>
    void AztecOOSolver<std::complex<double> >::solve(std::complex<double>* initial_guess)
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
      for (unsigned int i = 0; i < m->size; i++)
        this->sln[i] = std::complex<double>(xr[i], xi[i]);
    }

    template<typename Scalar>
    int AztecOOSolver<Scalar>::get_num_iters()
    {
      return aztec.NumIters();
    }

    template<typename Scalar>
    double AztecOOSolver<Scalar>::get_residual_norm()
    {
      return aztec.TrueResidual();
    }

    template class HERMES_API AztecOOSolver<double>;
    template class HERMES_API AztecOOSolver<std::complex<double> >;
  }
}
#endif