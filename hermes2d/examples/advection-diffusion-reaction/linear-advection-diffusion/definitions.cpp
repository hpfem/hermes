#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

class WeakFormLinearAdvectionDiffusion : public WeakForm
{
public:
  // Problem parameters.
  double const_f;

  WeakFormLinearAdvectionDiffusion(bool stabilization_on, bool shock_capturing_on, double b1, double b2, double epsilon) 
      : WeakForm(1), stabilization_on(stabilization_on), shock_capturing_on(shock_capturing_on) {
    add_matrix_form(new MatrixFormVolAdvectionDiffusion(0, 0, b1, b2, epsilon));
  };

private:
  class MatrixFormVolAdvectionDiffusion : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolAdvectionDiffusion(int i, int j, double b1, double b2, double epsilon) 
          : WeakForm::MatrixFormVol(i, j), b1(b1), b2(b2), epsilon(epsilon) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      Scalar result = 0;
      Real h_e = e->diam;
      Real s_c = 0.9;
      for (int i=0; i < n; i++) {
        result += wt[i] * (epsilon * (u->dx[i]*v->dx[i] + u->dy[i]*v->dy[i])
                                   - (b1 * u->val[i] * v->dx[i] + b2 * u->val[i] * v->dy[i]));
      
        if(static_cast<WeakFormLinearAdvectionDiffusion*>(wf)->shock_capturing_on) {
          Real R_squared = pow(b1 * u->dx[i] + b2 * u->dy[i], 2.);
          Real R = sqrt(R_squared); //This just does fabs(b1 * u->dx[i] + b2 * u->dy[i]); but it can be parsed
          result += wt[i] * s_c * 0.5 * h_e * R *
                    (u->dx[i]*v->dx[i] + u->dy[i]*v->dy[i]) /
                    (sqrt(pow(u->dx[i], 2) + pow(u->dy[i], 2)) + 1.e-8);
        }
        
        if(static_cast<WeakFormLinearAdvectionDiffusion*>(wf)->stabilization_on) {
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
          double b_norm = sqrt(b1*b1 + b2*b2);
          Real tau = 1. / sqrt(9*pow(4*epsilon/pow(h_e, 2), 2) + pow(2*b_norm/h_e, 2));
          result += wt[i] * tau * (-b1 * v->dx[i] - b2 * v->dy[i] + epsilon * v->laplace[i])
                                * (-b1 * u->dx[i] - b2 * u->dy[i] + epsilon * u->laplace[i]);
#else
          error("Define H2D_SECOND_DERIVATIVES_ENABLED in h2d_common.h if you want to use second derivatives of shape functions in weak forms.");
#endif
        }
      }
      return result;
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }
    // Members.
    double b1, b2, epsilon;
  };

  // Members.
  bool stabilization_on;
  bool shock_capturing_on;
};

class EssentialBCNonConst : public EssentialBC {
public:
  EssentialBCNonConst(std::string marker) : EssentialBC(Hermes::vector<std::string>()) {
    markers.push_back(marker);
  }

  ~EssentialBCNonConst() { }

  inline EssentialBCValueType get_value_type() const { return EssentialBC::BC_FUNCTION; }

  scalar function(double x, double y) const {
    return 2 - pow(x, 0.1) - pow(y, 0.1);
  }
};