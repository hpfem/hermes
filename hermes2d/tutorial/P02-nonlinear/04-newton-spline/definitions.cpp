#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

/* Weak forms */

class CustomWeakFormHeatTransferNewton : public WeakForm
{
public:
  CustomWeakFormHeatTransferNewton(CubicSpline* cspline) : WeakForm(1) {
    // Jacobian.
    add_matrix_form(new MatrixFormVolHeatTransfer(0, 0, cspline));
    // Residual.
    add_vector_form(new VectorFormVolHeatTransfer(0, cspline));
  };

private:
  class MatrixFormVolHeatTransfer : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolHeatTransfer(int i, int j, CubicSpline* cspline) 
      : WeakForm::MatrixFormVol(i, j, HERMES_NONSYM), cspline(cspline) { }


    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                 Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      Scalar result = 0;
      Func<Scalar>* u_prev = u_ext[0];
      bool success;
      for (int i = 0; i < n; i++) {
        result += wt[i] * (cspline->get_derivative(u_prev->val[i], success) * u->val[i] * 
                           (u_prev->dx[i] * v->dx[i] + u_prev->dy[i] * v->dy[i])
                           + cspline->get_value(u_prev->val[i], success) * (u->dx[i] * v->dx[i] 
                           + u->dy[i] * v->dy[i]));
      }
      return result;
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) {
      return Ord(10); //matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    // Spline representing temperature-dependent thermal conductivity.
    CubicSpline* cspline;
  };

  class VectorFormVolHeatTransfer : public WeakForm::VectorFormVol
  {
  public:
    VectorFormVolHeatTransfer(int i, CubicSpline* cspline) 
          : WeakForm::VectorFormVol(i), cspline(cspline) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) {
      Scalar result = 0;
      Func<Scalar>* u_prev = u_ext[0];
      bool success;
      for (int i = 0; i < n; i++)
        result += wt[i] * (cspline->get_value(u_prev->val[i], success) * (u_prev->dx[i] 
                           * v->dx[i] + u_prev->dy[i] * v->dy[i])
		           - heat_src<Real>(e->x[i], e->y[i]) * v->val[i]);
      return result;
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                 Geom<double> *e, ExtData<scalar> *ext) {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) {
      return Ord(10); //vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    // Heat sources (can be a general function of 'x' and 'y').
    template<typename Real>
    Real heat_src(Real x, Real y) {
      return 1.0;
    }

    // Spline representing temperature-dependent thermal conductivity.
    CubicSpline* cspline;
  };
};

/* Initial consition for the Newton's method */

class CustomInitialSolutionHeatTransfer : public ExactSolutionScalar
{
public:
  CustomInitialSolutionHeatTransfer(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    dx = (y+10)/100.;
    dy = (x+10)/100.;
  };

  virtual scalar value (double x, double y) const {
    return (x+10)*(y+10)/100. + 2;
  };

  virtual Ord ord(Ord x, Ord y) const {
    return x*y;
  }
};

/* Essential boundary conditions */

class CustomEssentialBCNonConst : public EssentialBoundaryCondition {
public:
  CustomEssentialBCNonConst(std::string marker) 
           : EssentialBoundaryCondition(Hermes::vector<std::string>())
  {
    markers.push_back(marker);
  }

  ~CustomEssentialBCNonConst() {};

  inline EssentialBCValueType get_value_type() const { 
    return EssentialBoundaryCondition::BC_FUNCTION; 
  }

  scalar value(double x, double y) const
  {
    return (x+10)*(y+10)/100.;
  }
};


