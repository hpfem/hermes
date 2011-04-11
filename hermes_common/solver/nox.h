// This file is part of Hermes2D
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

#ifndef __HERMES_COMMON_NOX_SOLVER_H_
#define __HERMES_COMMON_NOX_SOLVER_H_

#include "solver.h"
#include "epetra.h"

#ifdef HAVE_NOX
  #include <NOX.H>
    #ifdef _POSIX_C_SOURCE
        # undef _POSIX_C_SOURCE	// pyconfig.h included by NOX_Epetra defines it
    #endif
    #ifdef _XOPEN_SOURCE
        # undef _XOPEN_SOURCE	// pyconfig.h included by NOX_Epetra defines it
    #endif
  #include <NOX_Epetra.H>
#endif

class NoxProblemInterface;

/// Encapsulation of NOX nonlinear solver
///
/// @ingroup solvers
class HERMES_API NoxSolver : public IterSolver
{
public:
  // Basic constructor.
  NoxSolver(DiscreteProblemInterface *problem);
#ifdef HAVE_NOX
  // Enhanced constructor.
  // For details of the parameter message_type, please see NOX_Utils.H, enum MsgType.
  // http://trilinos.sandia.gov/packages/docs/r4.0/packages/nox/doc/html/parameters.html
  // http://trilinos.sandia.gov/packages/docs/r7.0/packages/nox/doc/html/classNOX_1_1Epetra_1_1LinearSystemAztecOO.html
  NoxSolver(DiscreteProblemInterface *problem, unsigned message_type, const char* ls_type = "GMRES", const char* nl_dir = "Newton", 
    double ls_tolerance = 1e-8,
    const char* precond_type = "None",
    unsigned flag_absresid = 1,
    double abs_resid = 1.0e-6,
    unsigned flag_relresid = 0,
    double rel_resid = 1.0e-2,
    int max_iters = 10,
    double update = 1.0e-5,
    int ls_max_iters = 800,
    int ls_sizeof_krylov_subspace = 50,
    NOX::Abstract::Vector::NormType norm_type = NOX::Abstract::Vector::TwoNorm,
    NOX::StatusTest::NormF::ScaleType stype = NOX::StatusTest::NormF::Scaled,
    double wrms_rtol = 1.0e-2,
    double wrms_atol = 1.0e-8,
    unsigned flag_update = 0,
    unsigned flag_wrms = 0
   );
#endif


  virtual ~NoxSolver();

  bool set_init_sln(double *ic);
  bool set_init_sln(EpetraVector *ic);
  virtual bool solve();

  virtual int get_num_iters() { return num_iters; }
  virtual double get_residual()  { return residual; }
  int get_num_lin_iters() { return num_lin_iters; }
  double get_achieved_tol()  { return achieved_tol; }

  // settings for the solver
  void set_nl_method(const char *par);
  void set_output_flags(int flags) { output_flags = flags; }

  // linear solver setters
  void set_ls_type(const char *type) { ls_type = type; }
  void set_ls_max_iters(int iters) { ls_max_iters = iters; }
  void set_ls_tolerance(double tolerance) { ls_tolerance = tolerance; }
  void set_ls_sizeof_krylov_subspace(int size) { ls_sizeof_krylov_subspace = size; }

  // convergence params
#ifdef HAVE_NOX
  void set_norm_type(NOX::Abstract::Vector::NormType type)  { conv.norm_type = type; }
  void set_scale_type(NOX::StatusTest::NormF::ScaleType type)  { conv.stype = type; }
#endif
  void set_conv_iters(int iters)        { conv.max_iters = iters; }
  void set_conv_abs_resid(double resid) { conv_flag.absresid = 1; conv.abs_resid = resid; }
  void set_conv_rel_resid(double resid) { conv_flag.relresid = 1; conv.rel_resid = resid; }
  void disable_abs_resid() { conv_flag.absresid = 0; }
  void disable_rel_resid() { conv_flag.relresid = 0; }
  void set_conv_update(double update)   { conv_flag.update = 1; conv.update = update; }
  void set_conv_wrms(double rtol, double atol) 
  {
    conv_flag.wrms = 1;
    conv.wrms_rtol = rtol;
    conv.wrms_atol = atol;
  }
  
#ifdef HAVE_TEUCHOS
  virtual void set_precond(Teuchos::RCP<Precond> &pc);
  virtual void unset_precond();
#else
  virtual void set_precond(Precond* pc) 
  { 
    warning("Teuchos is currently required to use IFPACK/ML preconditioners. No preconditioning will be performed.");  
  }
#endif
  virtual void set_precond(const char *pc);

protected:
#ifdef HAVE_NOX
  Teuchos::RCP<NoxProblemInterface> interface_;
#endif
  int num_iters;
  double residual;
  int num_lin_iters;
  double achieved_tol;  
  const char *nl_dir;

  int output_flags;
  const char *ls_type;
  int ls_max_iters;
  double ls_tolerance;
  int ls_sizeof_krylov_subspace;
  
  const char* precond_type;
  
  // convergence params
  struct conv_t {
    int max_iters;
    double abs_resid;
    double rel_resid;
#ifdef HAVE_NOX
    NOX::Abstract::Vector::NormType norm_type;
    NOX::StatusTest::NormF::ScaleType stype;
#endif
    double update;
    double wrms_rtol;
    double wrms_atol;
  } conv;

  struct conv_flag_t {
    unsigned absresid:1;
    unsigned relresid:1;
    unsigned wrms:1;
    unsigned update:1;
  } conv_flag;
};

#endif
