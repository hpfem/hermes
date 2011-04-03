#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

/* Weak forms */

class CustomWeakFormHeatTransferNewton : public WeakForm
{
public:
  CustomWeakFormHeatTransferNewton() : WeakForm(1) {
    // Jacobian.
    add_matrix_form(new MatrixFormVolHeatTransfer(0, 0));
    // Residual.
    add_vector_form(new VectorFormVolHeatTransfer(0));
  };

private:
  class MatrixFormVolHeatTransfer : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolHeatTransfer(int i, int j) 
          : WeakForm::MatrixFormVol(i, j, HERMES_NONSYM) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      Scalar result = 0;
      Func<Scalar>* u_prev = u_ext[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * (dlam_du<Real>(u_prev->val[i]) * u->val[i] * 
                           (u_prev->dx[i] * v->dx[i] + u_prev->dy[i] * v->dy[i])
                           + lam<Real>(u_prev->val[i]) * (u->dx[i] * v->dx[i] 
                           + u->dy[i] * v->dy[i]));
      return result;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    // Thermal conductivity (temperature-dependent)
    // Note: for any u, this function has to be positive.
    template<typename Real>
    Real lam(Real u) const { 
      return 1 + pow(u, 4); 
    }

    // Derivative of the thermal conductivity with respect to 'u'.
    template<typename Real>
    Real dlam_du(Real u) const { 
      return 4*pow(u, 3); 
    }
  };

  class VectorFormVolHeatTransfer : public WeakForm::VectorFormVol
  {
  public:
    VectorFormVolHeatTransfer(int i) : WeakForm::VectorFormVol(i) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) const {
      Scalar result = 0;
      Func<Scalar>* u_prev = u_ext[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * (lam<Real>(u_prev->val[i]) * (u_prev->dx[i] 
                           * v->dx[i] + u_prev->dy[i] * v->dy[i])
		           - heat_src<Real>(e->x[i], e->y[i]) * v->val[i]);
      return result;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                 Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    // Heat sources (can be a general function of 'x' and 'y').
    template<typename Real>
    Real heat_src(Real x, Real y) const {
      return 1.0;
    }

    // Thermal conductivity (temperature-dependent)
    // For any u, this function has to be positive.
    template<typename Real>
    Real lam(Real u) const { 
      return 1 + pow(u, 4); 
    }
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

  virtual scalar value(double x, double y, double n_x, double n_y, double t_x, double t_y) const
  {
    return (x+10)*(y+10)/100.;
  }
};


