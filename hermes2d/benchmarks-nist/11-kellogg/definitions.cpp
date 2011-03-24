#include "weakform/weakform.h"
#include "weakform_library/laplace.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh, double sigma, double tau, double rho) 
          : ExactSolutionScalar(mesh), sigma(sigma), tau(tau), rho(rho) {
  };

  virtual double value(double x, double y) const {
    double theta = atan2(y,x);
    if (theta < 0) theta = theta + 2.*M_PI;
    double r = sqrt(x*x + y*y);
    double mu;

    if (theta <= M_PI/2.)
      mu = cos((M_PI/2. - sigma)*tau) * cos((theta - M_PI/2. + rho)*tau);
    else if (theta <= M_PI)
      mu = cos(rho*tau) * cos((theta - M_PI + sigma)*tau);
    else if (theta <= 3.*M_PI/2.)
      mu = cos(sigma*tau) * cos((theta - M_PI - rho)*tau);
    else
      mu = cos((M_PI/2. - rho)*tau) * cos((theta - 3.*M_PI/2. - sigma)*tau);

    return pow(r, tau) * mu;
  };

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    double theta = atan2(y,x);
    if (theta < 0) theta = theta + 2*M_PI;
    double r = sqrt(x*x + y*y);
    
    // x-derivative
    if (theta <= M_PI/2.) 
      dx = tau*x*pow(r, (2.*(-1 + tau/2.))) * cos((M_PI/2. - sigma)*tau) * cos(tau*(-M_PI/2. + rho + theta)) 
          + (tau*y*pow(r, tau)*cos((M_PI/2. - sigma)*tau) * sin(tau*(-M_PI/2. + rho + theta))/(r*r));
    else if (theta <= M_PI)
      dx = tau*x * pow(r, (2.*(-1 + tau/2.))) * cos(rho*tau) * cos(tau*(-M_PI + sigma + theta)) 
          + (tau*y * pow(r, tau) * cos(rho*tau) * sin(tau*(-M_PI + sigma + theta))/(r*r));
    else if (theta <= 3.*M_PI/2.)
      dx = tau*x * pow(r, (2.*(-1 + tau/2.))) * cos(sigma*tau) * cos(tau*(-M_PI - rho + theta)) 
          + (tau*y * pow(r, tau) * cos(sigma*tau) * sin(tau*(-M_PI - rho + theta))/(r*r));
    else
      dx = tau*x* pow(r, (2*(-1 + tau/2.))) * cos((M_PI/2. - rho)*tau) * cos(tau*(-3.*M_PI/2. - sigma + theta)) 
          + (tau*y*pow(r, tau) * cos((M_PI/2. - rho)*tau) * sin(tau*(-3.*M_PI/2. - sigma + theta))/(r*r));
    
    // y-derivative
    if (theta <= M_PI/2.)
      dy = tau*y * pow(r, (2*(-1 + tau/2.))) * cos((M_PI/2. - sigma)*tau) * cos(tau*(-M_PI/2. + rho + theta)) 
          - (tau * pow(r, tau) * cos((M_PI/2. - sigma)*tau) *sin(tau*(-M_PI/2. + rho + theta))*x/(r*r));
    else if (theta <= M_PI)
      dy = tau*y* pow(r, (2*(-1 + tau/2.))) * cos(rho*tau) * cos(tau*(-M_PI + sigma + theta)) 
          - (tau * pow(r, tau) * cos(rho*tau) * sin(tau*(-M_PI + sigma + theta))*x/(r*r));
    else if (theta <= 3.*M_PI/2.)
      dy = tau*y * pow(r, (2*(-1 + tau/2.))) * cos(sigma*tau) * cos(tau*(-M_PI - rho + theta)) 
          - (tau * pow(r, tau) * cos(sigma*tau) * sin(tau*(-M_PI - rho + theta))*x/(r*r));
    else 
      dy = tau*y * pow(r, (2*(-1 + tau/2.))) * cos((M_PI/2. - rho)*tau) * cos(tau*(-3.*M_PI/2. - sigma + theta)) 
           - (tau * pow(r, tau) * cos((M_PI/2. - rho)*tau) * sin(tau*((-3.*M_PI)/2. - sigma + theta))*x/(r*r));

  };

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(6);
  }

  double sigma, tau, rho;
};

/* Weak forms */

class CustomWeakFormPoisson : public WeakForm
{
public:
  CustomWeakFormPoisson(double r) : WeakForm(1)
  {
    add_matrix_form(new CustomMatrixFormVol_I_III(0, 0, r));
    add_matrix_form(new CustomMatrixFormVol_II_IV(0, 0));
  };

private:
  class CustomMatrixFormVol_I_III : public WeakForm::MatrixFormVol
  {
  public:
    CustomMatrixFormVol_I_III(int i, int j, double r) 
      : WeakForm::MatrixFormVol(i, j, HERMES_SYM, "0"), r(r) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      return r * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    double r;
  };

  class CustomMatrixFormVol_II_IV : public WeakForm::MatrixFormVol
  {
  public:
    CustomMatrixFormVol_II_IV(int i, int j) 
          : WeakForm::MatrixFormVol(i, j, HERMES_SYM, "1") { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }
  };
};
