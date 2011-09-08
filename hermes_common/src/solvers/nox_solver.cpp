// This file is part of HermesCommon
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.
/*! \file nox_solver.cpp
\brief NOX (nonliner) solver interface.
*/
#include "nox_solver.h"
#if (defined HAVE_NOX && defined HAVE_EPETRA && defined HAVE_TEUCHOS)

namespace Hermes 
{
  namespace Solvers 
  {
    static Epetra_SerialComm seq_comm;

    template<typename Scalar>
    NoxDiscreteProblem<Scalar>::NoxDiscreteProblem(DiscreteProblemInterface<Scalar>* problem) : dp(problem)
    {
      this->precond = Teuchos::null;
      if(!this->dp->is_matrix_free()) 
        this->dp->create_sparse_structure(&jacobian);
    }

    template<typename Scalar>
    bool NoxDiscreteProblem<Scalar>::computeF(const Epetra_Vector &x, Epetra_Vector &f, FillType flag)
    {
      EpetraVector<Scalar> xx(x);  // wrap our structures around core Epetra objects
      EpetraVector<Scalar> rhs(f); 

      rhs.zero(); 

      Scalar* coeff_vec = new Scalar[xx.length()];
      xx.extract(coeff_vec);
      this->dp->assemble(coeff_vec, NULL, &rhs); // NULL is for the global matrix.
      delete [] coeff_vec;

      return true;
    }

    template<typename Scalar>
    bool NoxDiscreteProblem<Scalar>::computeJacobian(const Epetra_Vector &x, Epetra_Operator &op)
    {
      Epetra_RowMatrix *jac = dynamic_cast<Epetra_RowMatrix *>(&op);
      assert(jac != NULL);

      EpetraVector<Scalar> xx(x);			// wrap our structures around core Epetra objects
      EpetraMatrix<Scalar> jacob(*jac);

      jacob.zero();

      Scalar* coeff_vec = new Scalar[xx.length()];
      xx.extract(coeff_vec);  
      this->dp->assemble(coeff_vec, &jacob, NULL); // NULL is for the right-hand side.
      delete [] coeff_vec;
      //jacob.finish();

      return true;
    }

    template<typename Scalar>
    bool NoxDiscreteProblem<Scalar>::computePreconditioner(const Epetra_Vector &x, Epetra_Operator &m,
      Teuchos::ParameterList *precParams)
    {
      assert(precond != Teuchos::null);
      EpetraVector<Scalar> xx(x);			// wrap our structures around core Epetra objects

      jacobian.zero();

      Scalar* coeff_vec = new Scalar[xx.length()];
      xx.extract(coeff_vec);
      this->dp->assemble(coeff_vec, &jacobian, NULL);  // NULL is for the right-hand side.
      delete [] coeff_vec;
      //jacobian.finish();

      precond->create(&jacobian);
      precond->compute();
      m = *precond->get_obj();

      return true;
    }

    template<typename Scalar>
    NoxSolver<Scalar>::NoxSolver(DiscreteProblemInterface<Scalar>* problem) : NonlinearSolver<Scalar>(problem),ndp(problem)
    {
      // default values
      // convergence test
      conv.max_iters = 10;
      conv.abs_resid = 1.0e-6;
      conv.rel_resid = 1.0e-2;
      conv.norm_type = NOX::Abstract::Vector::TwoNorm;
      conv.stype = NOX::StatusTest::NormF::Scaled;
      conv.update = 1.0e-5;
      conv.wrms_rtol = 1.0e-2;
      conv.wrms_atol = 1.0e-8;

      conv_flag.absresid = 1;
      conv_flag.relresid = 0;
      conv_flag.update = 0;
      conv_flag.wrms = 0;

      nl_pars = Teuchos::rcp(new Teuchos::ParameterList);
      // Set the nonlinear solver method
      nl_pars->set("Nonlinear Solver", "Line Search Based");
      nl_pars->sublist("Printing").set("Output Information", NOX::Utils::Error);

      // Sublist for line search
      Teuchos::ParameterList &search_pars = nl_pars->sublist("Line Search");
      search_pars.set("Method", "Full Step");

      // Sublist for direction
      Teuchos::ParameterList &dir_pars = nl_pars->sublist("Direction");
      dir_pars.set("Method", "Newton");
      Teuchos::ParameterList &newton_pars = dir_pars.sublist("Newton");

      newton_pars.set("Forcing Term Method", "Constant");

      // Sublist for linear solver for the Newton method
      Teuchos::ParameterList &ls_pars = newton_pars.sublist("Linear Solver");
      ls_pars.set("Aztec Solver", "GMRES");
      ls_pars.set("Max Iterations", 800);
      ls_pars.set("Tolerance", 1e-8);
      ls_pars.set("Size of Krylov Subspace", 50);
      ls_pars.set("Preconditioner Reuse Policy", "Recompute");
      ls_pars.set("Output Frequency", AZ_all);
      ls_pars.set("Max Age Of Prec", 999);
    }

    template<typename Scalar>
    NoxSolver<Scalar>::~NoxSolver()
    {
      // FIXME: this does not destroy the "interface_", and Trilinos 
      // complains at closing main.cpp.
      this->dp->invalidate_matrix();
    }

    template<typename Scalar>
    void NoxDiscreteProblem<Scalar>::set_precond(Teuchos::RCP<Precond<Scalar> > &pc)
    {
      precond = pc;
      this->dp->create_sparse_structure(&jacobian);
    }

    template<typename Scalar>
    void NoxSolver<Scalar>::set_precond(Precond<Scalar> &pc)
    {
      Teuchos::RCP<Precond<Scalar> > tpc = Teuchos::rcpFromRef(pc);
      ndp.set_precond(tpc);
      nl_pars->sublist("Direction").sublist("Newton").sublist("Linear Solver").set("Preconditioner", "User Defined");
    }

    template<typename Scalar>
    void NoxSolver<Scalar>::set_precond(const char *pc)
    {
      nl_pars->sublist("Direction").sublist("Newton").sublist("Linear Solver").set("Preconditioner",pc);
    }
    
    template<typename Scalar>
    void NoxSolver<Scalar>::set_output_flags(int flags)
    {
      nl_pars->sublist("Printing").set("Output Information", flags);
    }

    template<typename Scalar>
    void NoxSolver<Scalar>::set_ls_type(const char *type){
      nl_pars->sublist("Direction").sublist("Newton").sublist("Linear Solver").set("Aztec Solver",type);
    }

    template<typename Scalar>
    void NoxSolver<Scalar>::set_ls_max_iters(int iters)
    { 
      nl_pars->sublist("Direction").sublist("Newton").sublist("Linear Solver").set("Max Iterations",iters);
    }

    template<typename Scalar>
    void NoxSolver<Scalar>::set_ls_tolerance(double tolerance) 
    { 
      nl_pars->sublist("Direction").sublist("Newton").sublist("Linear Solver").set("Tolerance",tolerance);
    }

    template<typename Scalar>
    void NoxSolver<Scalar>::set_ls_sizeof_krylov_subspace(int size) 
    { 
      nl_pars->sublist("Direction").sublist("Newton").sublist("Linear Solver").set("Size of Krylov Subspace",size);
    }

    template<typename Scalar>
    void NoxSolver<Scalar>::set_precond_reuse(const char * pc_reuse)
    {
      nl_pars->sublist("Direction").sublist("Newton").sublist("Linear Solver").set("Preconditioner Reuse Policy", pc_reuse);
    }
    template<typename Scalar>
    void NoxSolver<Scalar>::set_precond_max_age(int max_age)
    {
      nl_pars->sublist("Direction").sublist("Newton").sublist("Linear Solver").set("Max Age Of Prec", max_age);
    }

    template<typename Scalar>
    bool NoxSolver<Scalar>::solve(Scalar* coeff_vec)
    {
      // Put the initial coeff_vec into the inner structure for the initial guess.
      Hermes::Algebra::EpetraVector<Scalar> temp_init_sln;
      temp_init_sln.alloc(this->dp->get_num_dofs());
      for (int i = 0; i < this->dp->get_num_dofs(); i++)
        temp_init_sln.set(i, coeff_vec[i]);
      NOX::Epetra::Vector init_sln(*temp_init_sln.vec);

      // prepare variables
      // print settings
      Teuchos::ParameterList &print_pars = nl_pars->sublist("Printing");
      // linear system settings
      Teuchos::ParameterList &ls_pars = nl_pars->sublist("Direction").sublist("Newton").sublist("Linear Solver");
      // preconditioner
      Teuchos::RCP<Precond<Scalar> > precond = ndp.get_precond();
      //linear system     
      Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> lin_sys;
      // problem
      Teuchos::RCP<NOX::Epetra::Interface::Required> i_req = Teuchos::rcpFromRef(ndp);
      // jacobian matrix
      Teuchos::RCP<Epetra_RowMatrix> jac_mat;

      // Create linear system 
      if(this->dp->is_matrix_free()) 
      {
        if(precond == Teuchos::null) 
        { //Matrix free without preconditioner
          lin_sys = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(print_pars, ls_pars, i_req, init_sln));
        }
        else 
        { //Matrix free with preconditioner
          Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> i_prec = Teuchos::rcpFromRef(ndp);
          lin_sys = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(print_pars, ls_pars, i_req,
            i_prec, precond, init_sln));
        }
      }
      else
      {  // not matrix Free (with jacobian)
        Teuchos::RCP<NOX::Epetra::Interface::Jacobian> i_jac = Teuchos::rcpFromRef(ndp);
        jac_mat = Teuchos::rcp(ndp.get_jacobian()->mat);
        if(precond == Teuchos::null)
        { //Matrix  without preconditioner
          lin_sys = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(print_pars, ls_pars, i_req,
            i_jac, jac_mat, init_sln));
        }
        else
        { //Matrix  with preconditioner
          Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> i_prec = Teuchos::rcpFromRef(ndp);
          lin_sys = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(print_pars, ls_pars, 
            i_jac, jac_mat, i_prec, precond, init_sln));
        }
      }

      // Create the Group
      Teuchos::RCP<NOX::Epetra::Group> grp = 
        Teuchos::rcp(new NOX::Epetra::Group(print_pars, i_req, init_sln, lin_sys));

      // Create convergence tests
      Teuchos::RCP<NOX::StatusTest::Combo> converged =
        Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));

      if(conv_flag.absresid) 
      {
        Teuchos::RCP<NOX::StatusTest::NormF> absresid =
          Teuchos::rcp(new NOX::StatusTest::NormF(conv.abs_resid, conv.norm_type, conv.stype));
        converged->addStatusTest(absresid);
      }

      if(conv_flag.relresid) 
      {
        Teuchos::RCP<NOX::StatusTest::NormF> relresid = 
          Teuchos::rcp(new NOX::StatusTest::NormF(*grp.get(), conv.rel_resid));
        converged->addStatusTest(relresid);
      }

      if(conv_flag.update) 
      {
        Teuchos::RCP<NOX::StatusTest::NormUpdate> update = 
          Teuchos::rcp(new NOX::StatusTest::NormUpdate(conv.update));
        converged->addStatusTest(update);
      }

      if(conv_flag.wrms) 
      {
        Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms =
          Teuchos::rcp(new NOX::StatusTest::NormWRMS(conv.wrms_rtol, conv.wrms_atol));
        converged->addStatusTest(wrms);
      }

      Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters = 
        Teuchos::rcp(new NOX::StatusTest::MaxIters(conv.max_iters));

      Teuchos::RCP<NOX::StatusTest::FiniteValue> fv = 
        Teuchos::rcp(new NOX::StatusTest::FiniteValue);

      Teuchos::RCP<NOX::StatusTest::Combo> cmb = 
        Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));

      cmb->addStatusTest(fv);
      cmb->addStatusTest(converged);
      cmb->addStatusTest(maxiters);

      // Output parameters
      Teuchos::RCP<Teuchos::ParameterList> final_pars = nl_pars;

      // Create the solver
      Teuchos::RCP<NOX::Solver::Generic> solver = NOX::Solver::buildSolver(grp, cmb, final_pars);

      /// Solve.
      NOX::StatusTest::StatusType status = solver->solve();


      if(!this->dp->is_matrix_free())
        jac_mat.release();	// release the ownership (we take care of jac_mat by ourselves)

      bool success;
      if(status == NOX::StatusTest::Converged) 
      {
        // get result informations
        num_iters = solver->getNumIterations();
        residual = solver->getSolutionGroup().getNormF();
        num_lin_iters = final_pars->sublist("Direction").sublist("Newton").sublist("Linear Solver").sublist("Output").get("Total Number of Linear Iterations", -1);
        achieved_tol = final_pars->sublist("Direction").sublist("Newton").sublist("Linear Solver").sublist("Output").get("Achieved Tolerance", 0.0);

        // Get the Epetra_Vector with the final solution from the solver
        const NOX::Epetra::Group &f_grp =
          dynamic_cast<const NOX::Epetra::Group &>(solver->getSolutionGroup());
        const Epetra_Vector &f_sln =
          (dynamic_cast<const NOX::Epetra::Vector &>(f_grp.getX())).getEpetraVector();
        // extract solution
        int n = this->dp->get_num_dofs();
        delete [] this->sln_vector;
        this->sln_vector = new double[n];
        memset(this->sln_vector, 0, n * sizeof(double));
        f_sln.ExtractCopy(this->sln_vector);

        success = true;
      }
      else { // not converged
        num_iters = -1;
        success = false;
      }
      return success;
    }
    template class HERMES_API NoxSolver<double>;
    //template class HERMES_API NoxSolver<std::complex<double> >; //complex version of nox solver is not implemented 
  }
}
#endif
