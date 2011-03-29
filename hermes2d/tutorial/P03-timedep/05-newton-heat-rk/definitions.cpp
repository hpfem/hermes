#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

class WeakFormHeatTransferRungeKuttaTimedep : public WeakForm
{
public:
  WeakFormHeatTransferRungeKuttaTimedep(double alpha, Solution* sln_prev_time) : WeakForm(1) {
    MatrixFormVolHeatTransfer* matrix_form = new MatrixFormVolHeatTransfer(0, 0, alpha);
    matrix_form->ext.push_back(sln_prev_time);
    add_matrix_form(matrix_form);

    VectorFormVolHeatTransfer* vector_form = new VectorFormVolHeatTransfer(0, alpha);
    vector_form->ext.push_back(sln_prev_time);
    add_vector_form(vector_form);
  };

private:
  class MatrixFormVolHeatTransfer : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolHeatTransfer(int i, int j, double alpha) : WeakForm::MatrixFormVol(i, j), alpha(alpha) {
      sym = HERMES_NONSYM; 
    }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      Func<Scalar>* K_sln = u_ext[0];
      Func<Scalar>* sln_prev_time = ext->fn[0];

      // Stationary part of the Jacobian matrix (time derivative term left out).
      Scalar result1 = 0, result2 = 0;
      for (int i = 0; i < n; i++) {
        Scalar sln_val_i = sln_prev_time->val[i] + K_sln->val[i];
        Scalar sln_dx_i = sln_prev_time->dx[i] + K_sln->dx[i];
        Scalar sln_dy_i = sln_prev_time->dy[i] + K_sln->dy[i];
        result1 += -wt[i] * dlam_du<Real>(sln_val_i) * u->val[i] * (sln_dx_i * v->dx[i] + sln_dy_i * v->dy[i]);
        result2 += -wt[i] * lam<Real>(sln_val_i) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
      }

      return result1 + result2;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    virtual WeakForm::MatrixFormVol* clone() {
      return new MatrixFormVolHeatTransfer(*this);
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
    
    // Member.
    double alpha;
  };

  class VectorFormVolHeatTransfer : public WeakForm::VectorFormVol
  {
  public:
    VectorFormVolHeatTransfer(int i, double alpha) : WeakForm::VectorFormVol(i), alpha(alpha) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      Func<Scalar>* K_sln = u_ext[0];
      Func<Scalar>* sln_prev_time = ext->fn[0];

      // Stationary part of the residual (time derivative term left out).
      Scalar result1 = 0, result2 = 0;
      for (int i = 0; i < n; i++) {
        Scalar sln_val_i = sln_prev_time->val[i] + K_sln->val[i];
        Scalar sln_dx_i = sln_prev_time->dx[i] + K_sln->dx[i];
        Scalar sln_dy_i = sln_prev_time->dy[i] + K_sln->dy[i];
        result1 = result1 - wt[i] * lam<Real>(sln_val_i) * (sln_dx_i * v->dx[i] + sln_dy_i * v->dy[i]);
        result2 = result2 + wt[i] * heat_src<Real>(e->x[i], e->y[i]) * v->val[i];
      }

      return result1 + result2;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    virtual WeakForm::VectorFormVol* clone() {
      return new VectorFormVolHeatTransfer(*this);
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

    // Member.
    double alpha;
  };
};

class EssentialBCNonConst : public EssentialBoundaryCondition {
public:
  EssentialBCNonConst(std::string marker) : EssentialBoundaryCondition(Hermes::vector<std::string>())
  {
    markers.push_back(marker);
  }

  ~EssentialBCNonConst() {};

  inline EssentialBCValueType get_value_type() const { return EssentialBoundaryCondition::BC_FUNCTION; }

  virtual scalar value(double x, double y) const
  {
    return (x+10)*(y+10)/100.;
  }
};

class InitialSolutionHeatTransfer : public ExactSolutionScalar
{
public:
  InitialSolutionHeatTransfer(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  // Function representing an exact one-dimension valued solution.
  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    dx = (y+10)/10.;
    dy = (x+10)/10.;
  };

  virtual scalar value (double x, double y) const {
    return (x+10)*(y+10)/100.;
  };

  virtual Ord ord(Ord x, Ord y) const {
    return (x+10)*(y+10)/100.;
  }
};
