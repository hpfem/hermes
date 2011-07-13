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
    void NoxSolver<Scalar>::prealloc_jacobian()
    {
      this->dp->create_sparse_structure(&jacobian);
    }

    template<typename Scalar>
    bool NoxSolver<Scalar>::computeF(const Epetra_Vector &x, Epetra_Vector &f, FillType flag)
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
    bool NoxSolver<Scalar>::computeJacobian(const Epetra_Vector &x, Epetra_Operator &op)
    {
      Epetra_RowMatrix *jac = dynamic_cast<Epetra_RowMatrix *>(&op);
      assert(jac != NULL);

      EpetraVector<Scalar> xx(x);			// wrap our structures around core Epetra objects
      EpetraMatrix<Scalar> jacobian(*jac);

      jacobian.zero();

      Scalar* coeff_vec = new Scalar[xx.length()];
      xx.extract(coeff_vec);
      this->dp->assemble(coeff_vec, &jacobian, NULL); // NULL is for the right-hand side.
      delete [] coeff_vec;
      //jacobian.finish();

      return true;
    }

    template<typename Scalar>
    bool NoxSolver<Scalar>::computePreconditioner(const Epetra_Vector &x, Epetra_Operator &m,
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
    NoxSolver<Scalar>::NoxSolver(DiscreteProblemInterface<Scalar>* problem) : NonlinearSolver<Scalar>(problem)
    {
      if(!this->dp->is_matrix_free()) 
        prealloc_jacobian();

      this->precond = Teuchos::null;

      // default values
      nl_dir = "Newton";
      output_flags = NOX::Utils::Error;
      // linear solver settings
      ls_type = "GMRES";
      ls_max_iters = 800;
      ls_tolerance = 1e-8;
      ls_sizeof_krylov_subspace = 50;
      precond_type = "None";
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
    }

    template<typename Scalar>
    NoxSolver<Scalar>::NoxSolver(DiscreteProblemInterface<Scalar> *problem, unsigned message_type, const char* ls_type, const char* nl_dir, 
      double ls_tolerance,
      const char* precond_type,
      unsigned flag_absresid,
      double abs_resid,
      unsigned flag_relresid,
      double rel_resid,
      int max_iters,
      double update,
      int ls_max_iters,
      int ls_sizeof_krylov_subspace,
      NOX::Abstract::Vector::NormType norm_type,
      NOX::StatusTest::NormF::ScaleType stype,
      double wrms_rtol,
      double wrms_atol,
      unsigned flag_update,
      unsigned flag_wrms
      ) : NonlinearSolver<Scalar>(problem)
    {
      if(!this->dp->is_matrix_free()) 
        prealloc_jacobian();

      this->precond = Teuchos::null;

      // default values
      this->nl_dir = nl_dir;
      output_flags = message_type;
      // linear solver settings
      this->ls_type = ls_type;
      this->ls_max_iters = ls_max_iters;
      this->ls_tolerance = ls_tolerance;
      this->ls_sizeof_krylov_subspace = ls_sizeof_krylov_subspace;
      this->precond_type = precond_type;
      // convergence test
      this->conv.max_iters = max_iters;
      this->conv.abs_resid = abs_resid;
      this->conv.rel_resid = rel_resid;
      this->conv.norm_type = norm_type;
      this->conv.stype = stype;
      this->conv.update = update;
      this->conv.wrms_rtol = wrms_rtol;
      this->conv.wrms_atol = wrms_atol;

      this->conv_flag.absresid = flag_absresid;
      this->conv_flag.relresid = flag_relresid;
      this->conv_flag.update = flag_update;
      this->conv_flag.wrms = flag_wrms;
    }

    template<typename Scalar>
    NoxSolver<Scalar>::~NoxSolver()
    {
      // FIXME: this does not destroy the "interface_", and Trilinos 
      // complains at closing main.cpp.
      this->dp->invalidate_matrix();
    }

    template<typename Scalar>
    void NoxSolver<Scalar>::set_precond(Teuchos::RCP<Precond<Scalar> > &pc)
    {
      this->precond_yes = true;
      precond = pc;
      prealloc_jacobian();
    }

    template<typename Scalar>
    void NoxSolver<Scalar>::set_precond(const char *pc)
    {
      this->precond_yes = true;
      precond_type = pc;
    }

    template<>
    void NoxSolver<double>::get_final_solution(Teuchos::RCP<NOX::Solver::Generic> & solver)
    {
      const NOX::Epetra::Group &f_grp =
        dynamic_cast<const NOX::Epetra::Group &>(solver->getSolutionGroup());
      const Epetra_Vector &f_sln =
        (dynamic_cast<const NOX::Epetra::Vector &>(f_grp.getX())).getEpetraVector();
      // extract solution
      int n = this->dp->get_num_dofs();
      delete [] sln_vector;
      sln_vector = new double[n];
      memset(sln_vector, 0, n * sizeof(double));
      f_sln.ExtractCopy(sln_vector);
    }

    template<>
    void NoxSolver<std::complex<double> >::get_final_solution(Teuchos::RCP<NOX::Solver::Generic> & solver)
    {
      // extract solution
      int n = this->dp->get_num_dofs();
      delete [] sln_vector;
      sln_vector = new std::complex<double>[n];
      memset(sln_vector, 0, n * sizeof(std::complex<double>));
    }

    template<typename Scalar>
    bool NoxSolver<Scalar>::solve(Scalar* coeff_vec)
    {
      // Put the initial coeff_vec into the inner structure for the initial guess.
      /// \todo Put this into a separate method.
      Hermes::Algebra::EpetraVector<Scalar> temp_init_sln;
      temp_init_sln.alloc(this->dp->get_num_dofs());
      for (int i = 0; i < this->dp->get_num_dofs(); i++)
        temp_init_sln.set(i, coeff_vec[i]);
      NOX::Epetra::Vector init_sln(*temp_init_sln.vec);

      // Create the top level parameter list
      Teuchos::RCP<Teuchos::ParameterList> nl_pars_ptr = Teuchos::rcp(new Teuchos::ParameterList);
      Teuchos::ParameterList &nl_pars = *nl_pars_ptr.get();

      // Set the nonlinear solver method
      nl_pars.set("Nonlinear Solver", "Line Search Based");

      // Set the printing parameters in the "Printing" sublist
      Teuchos::ParameterList &print_pars = nl_pars.sublist("Printing");
      print_pars.set("Output Information", output_flags);

      // Sublist for line search
      Teuchos::ParameterList &search_pars = nl_pars.sublist("Line Search");
      search_pars.set("Method", "Full Step");

      // Sublist for direction
      Teuchos::ParameterList &dir_pars = nl_pars.sublist("Direction");
      dir_pars.set("Method", nl_dir);
      Teuchos::ParameterList &newton_pars = dir_pars.sublist(nl_dir);

      if(strcmp(nl_dir, "Newton") == 0)
        newton_pars.set("Forcing Term Method", "Constant");

      // Sublist for linear solver for the Newton method
      Teuchos::ParameterList &ls_pars = newton_pars.sublist("Linear Solver");
      ls_pars.set("Aztec Solver", ls_type);
      ls_pars.set("Max Iterations", ls_max_iters);
      ls_pars.set("Tolerance", ls_tolerance);
      ls_pars.set("Size of Krylov Subspace", ls_sizeof_krylov_subspace);
      /// \todo parametrize me.
      ls_pars.set("Preconditioner Reuse Policy", "Recompute");
      ls_pars.set("Output Frequency", AZ_all);

      // precond stuff
      Teuchos::RCP<Precond<Scalar> > precond = this->get_precond();
      if(this->precond_yes == false)
        ls_pars.set("Preconditioner", "None");
      else 
        if(this->dp->is_matrix_free()) 
          ls_pars.set("Preconditioner", "User Defined");
        else 
        {
          if(strcasecmp(precond_type, "ML") == 0)
            ls_pars.set("Preconditioner", "ML");
          else 
            if(strcasecmp(precond_type, "Ifpack") == 0)
              ls_pars.set("Preconditioner", "Ifpack");
            else 
            {
              warn("Unsupported type of preconditioner.");
              ls_pars.set("Preconditioner", "None");
            }
        }

        /// \todo Parametrize me.
        ls_pars.set("Max Age Of Prec", 999);

        Teuchos::RCP<NOX::Epetra::Interface::Required> i_req = Teuchos::rcp(this);
        Teuchos::RCP<NOX::Epetra::Interface::Jacobian> i_jac = Teuchos::rcp(this);
        Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> i_prec = Teuchos::rcp(this);
        Teuchos::RCP<Epetra_RowMatrix> jac_mat;
        Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> lin_sys;


        if(this->dp->is_matrix_free()) 
        {
          // Matrix<Scalar>-Free (Epetra_Operator)
          if(precond == Teuchos::null) 
          {
            Teuchos::RCP<NOX::Epetra::MatrixFree> mf = 
              Teuchos::rcp(new NOX::Epetra::MatrixFree(print_pars, Teuchos::rcp(this), init_sln));
            i_jac = mf;
            lin_sys = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(print_pars, ls_pars, i_req,
              i_jac, mf, init_sln));
          }
          else 
          {
            const Teuchos::RCP<Epetra_Operator> pc = precond;
            lin_sys = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(print_pars, ls_pars, i_req,
              i_prec, pc, init_sln));
          }
        }
        else {  // not Matrix<Scalar> Free
          // Create the Epetra_RowMatrix.
          jac_mat = Teuchos::rcp(this->get_jacobian()->mat);
          i_jac = Teuchos::rcp(this);
          lin_sys = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(print_pars, ls_pars, i_req,
            i_jac, jac_mat, init_sln));
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

        Teuchos::RCP<Teuchos::ParameterList> final_pars = nl_pars_ptr;

        // Create the method
        Teuchos::RCP<NOX::Solver::Generic> solver = NOX::Solver::buildSolver(grp, cmb, final_pars);


        /// Solve.
        NOX::StatusTest::StatusType status = solver->solve();


        if(!this->dp->is_matrix_free())
          jac_mat.release();	// release the ownership (we take care of jac_mat by ourselves)

        bool success;
        if(status == NOX::StatusTest::Converged) 
        {
          num_iters = solver->getNumIterations();
          residual = solver->getSolutionGroup().getNormF();
          num_lin_iters = final_pars->sublist("Direction").sublist(nl_dir).sublist("Linear Solver").sublist("Output").get("Total Number of Linear Iterations", -1);
          achieved_tol = final_pars->sublist("Direction").sublist(nl_dir).sublist("Linear Solver").sublist("Output").get("Achieved Tolerance", 0.0);

          // Get the Epetra_Vector with the final solution from the solver
          get_final_solution(solver);
          success = true;
        }
        else { // not converged
          num_iters = -1;
          success = false;
        }
        return success;
    }
    template class HERMES_API NoxSolver<double>;
    template class HERMES_API NoxSolver<std::complex<double> >;
  }
}
#endif
