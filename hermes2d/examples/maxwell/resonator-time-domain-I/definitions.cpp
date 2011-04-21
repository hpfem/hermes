#include "weakform/weakform.h"
#include "weakform_library/hcurl.h"
#include "integrals/h1.h"
#include "boundaryconditions/essential_bcs.h"

using namespace WeakFormsHcurl;
using namespace WeakFormsHcurl::VolumetricMatrixForms;

/* Initial condition for E */

class CustomInitialConditionWave : public ExactSolutionVector
{
public:
  CustomInitialConditionWave(Mesh* mesh) : ExactSolutionVector(mesh) {};

  virtual scalar2 value (double x, double y) const {
    return scalar2(sin(x) * cos(y), -cos(x) * sin(y));
  }

  virtual void derivatives (double x, double y, scalar2& dx, scalar2& dy) const {
    dx[0] = cos(x) * cos(y);
    dx[1] = sin(x) * sin(y);
    dy[0] = -sin(x) * sin(y);
    dy[1] = -cos(x) * cos(y);
  }

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(10);
  }
};

/* Weak forms */

class CustomWeakFormWave : public WeakForm
{
public:

  CustomWeakFormWave(double c_squared) : WeakForm(2) {
    // Jacobian.
    add_matrix_form(new MatrixFormVolWave_0_1(c_squared));
    add_matrix_form(new MatrixFormVolWave_1_0);

    // Residual.
    add_vector_form(new VectorFormVolWave_0(c_squared));
    add_vector_form(new VectorFormVolWave_1());
  };

private:
  class MatrixFormVolWave_0_1 : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolWave_0_1(double c_squared) 
      : WeakForm::MatrixFormVol(0, 1, HERMES_ANY, HERMES_NONSYM), c_squared(c_squared) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      Scalar result = 0;
      for (int i=0; i < n; i++) {
        result += wt[i] * (u->dy[i] * v->val0[i] - u->dx[i] * v->val1[i]);
      }
      return c_squared * result;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    virtual WeakForm::MatrixFormVol* clone() {
      return new MatrixFormVolWave_0_1(*this);
    }

    double c_squared;
  };

  class MatrixFormVolWave_1_0 : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolWave_1_0() 
      : WeakForm::MatrixFormVol(1, 0, HERMES_ANY, HERMES_NONSYM) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      Scalar result = 0;
      for (int i=0; i < n; i++) {
        result -= wt[i] * u->curl[i] * v->val[i];
      }

      return result;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    virtual WeakForm::MatrixFormVol* clone() {
      return new MatrixFormVolWave_1_0(*this);
    }
  };

  class VectorFormVolWave_0 : public WeakForm::VectorFormVol
  {
  public:
    VectorFormVolWave_0(double c_squared) 
          : WeakForm::VectorFormVol(0), c_squared(c_squared) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) const {
      Scalar result = 0;
      Func<Scalar>* sln_prev = u_ext[1];

      for (int i = 0; i < n; i++) {
        result += wt[i] * (sln_prev->dy[i] * v->val0[i] - sln_prev->dx[i] * v->val1[i]);
      }
      return c_squared * result;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                 Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
            ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    virtual WeakForm::VectorFormVol* clone() {
      return new VectorFormVolWave_0(*this);
    }

    double c_squared;
  };

  class VectorFormVolWave_1 : public WeakForm::VectorFormVol
  {
  public:
    VectorFormVolWave_1() 
          : WeakForm::VectorFormVol(1) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) const {
      Scalar result = 0;
      Func<Scalar>* sln_prev = u_ext[0];
      
      for (int i = 0; i < n; i++) {
        result -= wt[i] * sln_prev->curl[i] * v->val[i];
      }
      return result;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                 Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    virtual WeakForm::VectorFormVol* clone() {
      return new VectorFormVolWave_1(*this);
    }
  };
};
