#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

class WeakFormHeatTransferNewtonTimedepSDIRKStage1 : public WeakForm
{
public:
  WeakFormHeatTransferNewtonTimedepSDIRKStage1(double alpha, double tau, 
      Solution* sln_prev_time, double butcher_a_11, double gamma, double butcher_c_1) : WeakForm(1) {
    add_matrix_form(new MatrixFormVolHeatTransfer(0, 0, alpha, tau, butcher_a_11));

    VectorFormVolHeatTransfer* vector_form = new VectorFormVolHeatTransfer(0, alpha, tau, gamma, butcher_c_1);
    vector_form->ext.push_back(sln_prev_time);
    add_vector_form(vector_form);
  };

private:
  class MatrixFormVolHeatTransfer : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolHeatTransfer(int i, int j, double alpha, double tau, double butcher_a_11) : 
        WeakForm::MatrixFormVol(i, j), alpha(alpha), tau(tau), butcher_a_11(butcher_a_11) {
      sym = HERMES_NONSYM; 
    }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      Scalar result = 0;
      Func<Scalar>* u_prev_newton = u_ext[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->val[i] * v->val[i] / tau
                           + butcher_a_11 * (dlam_du<Real>(u_prev_newton->val[i]) * u->val[i] 
                                      * (u_prev_newton->dx[i] * v->dx[i] + u_prev_newton->dy[i] * v->dy[i])
                                      + lam<Real>(u_prev_newton->val[i]) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i])));
      return result;
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    // Thermal conductivity (temperature-dependent)
    // Note: for any u, this function has to be positive.
    template<typename Real>
    Real lam(Real u) { 
      return 1 + pow(u, alpha); 
    }

    // Derivative of the thermal conductivity with respect to 'u'.
    template<typename Real>
    Real dlam_du(Real u) { 
      return alpha*pow(u, alpha - 1); 
    }
    
    // Members.
    double alpha;
    double tau;
    double butcher_a_11;
  };

  class VectorFormVolHeatTransfer : public WeakForm::VectorFormVol
  {
  public:
  VectorFormVolHeatTransfer(int i, double alpha, double tau, double gamma, double butcher_c_1) :
      WeakForm::VectorFormVol(i), alpha(alpha), tau(tau), gamma(gamma), butcher_c_1(butcher_c_1) { }

  template<typename Real, typename Scalar>
  Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
    Scalar result = 0;
    Func<Scalar>* Y1_prev_newton = u_ext[0];
    Func<Scalar>* u_prev_time = ext->fn[0];
    for (int i = 0; i < n; i++) {
      result += wt[i] * (Y1_prev_newton->val[i] - u_prev_time->val[i]) * v->val[i] / tau;
    }
    result += gamma * res_ss(n, wt, u_ext, v, e, ext, wf->get_current_time() + butcher_c_1 * tau);
    return result;
  }

  scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
    return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
  }

  Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
    return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
  }

  template<typename Real, typename Scalar>
  Scalar res_ss(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext, double t) {
    Scalar result = 0;
    Func<Scalar>* Y_prev_newton = u_ext[0];
    for (int i = 0; i < n; i++) {
      result += wt[i] * (lam<Real>(Y_prev_newton->val[i]) * (Y_prev_newton->dx[i] * v->dx[i] + Y_prev_newton->dy[i] * v->dy[i])
                          - heat_src<Real>(e->x[i], e->y[i]) * v->val[i]);
    }
    return result;
  }

  // Heat sources (can be a general function of 'x' and 'y').
  template<typename Real>
  Real heat_src(Real x, Real y) {
    return 1.0;
  }

  // Thermal conductivity (temperature-dependent)
  // Note: for any u, this function has to be positive.
  template<typename Real>
  Real lam(Real u)  { 
    return 1 + pow(u, alpha); 
  }

  // Members.
  double alpha;
  double tau;
  double gamma;
  double butcher_c_1;
  };
};



class WeakFormHeatTransferNewtonTimedepSDIRKStage2 : public WeakForm
{
public:
  WeakFormHeatTransferNewtonTimedepSDIRKStage2(double alpha, double tau, Solution* sln_prev_time, Solution* sln_stage_1,
    double butcher_a_11, double gamma, double butcher_b_1, double butcher_b_2, double butcher_c_1, double butcher_c_2) : WeakForm(1) 
  {
    add_matrix_form(new MatrixFormVolHeatTransfer(0, 0, alpha, tau, butcher_a_11));

    VectorFormVolHeatTransfer* vector_form = new VectorFormVolHeatTransfer(0, alpha, tau, gamma, butcher_b_1,
                                             butcher_b_2, butcher_c_1, butcher_c_2);
    vector_form->ext.push_back(sln_prev_time);
    vector_form->ext.push_back(sln_stage_1);
    add_vector_form(vector_form);
  };

private:
  class MatrixFormVolHeatTransfer : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolHeatTransfer(int i, int j, double alpha, double tau, double butcher_a_11) : 
        WeakForm::MatrixFormVol(i, j), alpha(alpha), tau(tau), butcher_a_11(butcher_a_11) {
      sym = HERMES_NONSYM; 
    }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      Scalar result = 0;
      Func<Scalar>* u_prev_newton = u_ext[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->val[i] * v->val[i] / tau
                           + butcher_a_11 * (dlam_du<Real>(u_prev_newton->val[i]) * u->val[i] 
                                      * (u_prev_newton->dx[i] * v->dx[i] + u_prev_newton->dy[i] * v->dy[i])
                                      + lam<Real>(u_prev_newton->val[i]) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i])));
      return result;
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    // Thermal conductivity (temperature-dependent)
    // Note: for any u, this function has to be positive.
    template<typename Real>
    Real lam(Real u) { 
      return 1 + pow(u, alpha); 
    }

    // Derivative of the thermal conductivity with respect to 'u'.
    template<typename Real>
    Real dlam_du(Real u) { 
      return alpha*pow(u, alpha - 1); 
    }
    
    // Members.
    double alpha;
    double tau;
    double butcher_a_11;
  };

  class VectorFormVolHeatTransfer : public WeakForm::VectorFormVol
  {
  public:
    VectorFormVolHeatTransfer(int i, double alpha, double tau, double gamma, double butcher_b_1, 
        double butcher_b_2, double butcher_c_1, double butcher_c_2) :
        WeakForm::VectorFormVol(i), alpha(alpha), tau(tau), gamma(gamma), butcher_b_1(butcher_b_1),
        butcher_b_2(butcher_b_2), butcher_c_1(butcher_c_1), butcher_c_2(butcher_c_2) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      Scalar result = 0;
      Func<Scalar>* Y2_prev_newton = u_ext[0];
      Func<Scalar>* u_prev_time = ext->fn[0];
      Func<Scalar>* Y1[] = {ext->fn[1]};
      for (int i = 0; i < n; i++) {
        result += wt[i] * (Y2_prev_newton->val[i] - u_prev_time->val[i]) * v->val[i] / tau;
      }
      result += butcher_b_1 * res_ss(n, wt, Y1, v, e, ext, wf->get_current_time() + butcher_c_1 * tau);
      result += butcher_b_2 * res_ss(n, wt, u_ext, v, e, ext, wf->get_current_time() + butcher_c_2 * tau);
      return result;
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    template<typename Real, typename Scalar>
    Scalar res_ss(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext, double t) {
      Scalar result = 0;
      Func<Scalar>* Y_prev_newton = u_ext[0];
      for (int i = 0; i < n; i++) {
        result += wt[i] * (lam<Real>(Y_prev_newton->val[i]) * (Y_prev_newton->dx[i] * v->dx[i] + Y_prev_newton->dy[i] * v->dy[i])
                            - heat_src<Real>(e->x[i], e->y[i]) * v->val[i]);
      }
      return result;
    }

    // Heat sources (can be a general function of 'x' and 'y').
    template<typename Real>
    Real heat_src(Real x, Real y) {
      return 1.0;
    }

    // Thermal conductivity (temperature-dependent)
    // Note: for any u, this function has to be positive.
    template<typename Real>
    Real lam(Real u)  { 
      return 1 + pow(u, alpha); 
    }

    // Members.
    double alpha;
    double tau;
    double gamma;
    double butcher_b_1;
    double butcher_b_2;
    double butcher_c_1;
    double butcher_c_2;
  };
};



class EssentialBCNonConst : public EssentialBC {
public:
  EssentialBCNonConst(std::string marker) : EssentialBC(Hermes::vector<std::string>())
  {
    markers.push_back(marker);
  }

  ~EssentialBCNonConst() {};

  inline EssentialBCValueType get_value_type() const { return EssentialBC::BC_FUNCTION; }

  scalar function(double x, double y) const
  {
    return (x+10)*(y+10)/100.;
  }
};