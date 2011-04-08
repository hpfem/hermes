#include "weakform/weakform.h"
#include "weakform_library/h1.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

/* Initial condition */

class CustomInitialConditionWave : public ExactSolutionScalar
{
public:
  CustomInitialConditionWave(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  virtual scalar value (double x, double y) const {
    return exp(-x*x - y*y);
  }

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    dx = exp(-x*x - y*y) * (-2*x);
    dy = exp(-x*x - y*y) * (-2*y);
  };

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(10);
  }
};

/* Weak forms */

class CustomWeakFormWave : public WeakForm
{
public:

  CustomWeakFormWave(double tau, double c_squared, Solution* u_prev_sln, Solution* v_prev_sln) : WeakForm(2) {
    add_matrix_form(new MatrixFormVolWave_0_1);
    add_matrix_form(new MatrixFormVolWave_1_0(c_squared));

    add_vector_form(new VectorFormVolWave_0());
    
    add_vector_form(new VectorFormVolWave_1(c_squared));
  };

private:
  // This form is custom because of the clone() method.
  class MatrixFormVolWave_0_1 : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolWave_0_1() : WeakForm::MatrixFormVol(0, 1) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      return int_u_v<Real, Scalar>(n, wt, u, v);
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    virtual WeakForm::MatrixFormVol* clone() {
      return new MatrixFormVolWave_0_1(*this);
    }
  };

  // This form is custom because of the clone() method.
  class MatrixFormVolWave_1_0 : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolWave_1_0(double c_squared) 
          : WeakForm::MatrixFormVol(1, 0), c_squared(c_squared) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      return - c_squared * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    virtual WeakForm::MatrixFormVol* clone() {
      return new MatrixFormVolWave_1_0(*this);
    }

    double c_squared;
  };

  class VectorFormVolWave_0 : public WeakForm::VectorFormVol
  {
  public:
    VectorFormVolWave_0() : WeakForm::VectorFormVol(0) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) const {
      Scalar result = 0;

      for (int i = 0; i < n; i++)
        result += wt[i] * u_ext[1]->val[i] * v->val[i];

      return result;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                 Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
            ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    virtual WeakForm::VectorFormVol* clone() {
      return new VectorFormVolWave_0(*this);
    }
  };

  class VectorFormVolWave_1 : public WeakForm::VectorFormVol
  {
  public:
    VectorFormVolWave_1(double c_squared) 
          : WeakForm::VectorFormVol(1), c_squared(c_squared) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) const {
      Scalar result = 0;
      
      for (int i = 0; i < n; i++)
        result += wt[i] * (u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i]);
      
      return - c_squared * result;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                 Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
            ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    virtual WeakForm::VectorFormVol* clone() {
      return new VectorFormVolWave_1(*this);
    }

    double c_squared;
  };
};
