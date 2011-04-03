#include "weakform/weakform.h"
#include "weakform_library/hcurl.h"
#include "integrals/integrals_h1.h"
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

  CustomWeakFormWave(double tau, double c_squared, Solution* E_prev_sln, Solution* B_prev_sln) : WeakForm(2) {
    add_matrix_form(new MatrixFormVolWave_0_1(c_squared));
    add_matrix_form(new MatrixFormVolWave_1_0);

    VectorFormVolWave_0* vector_form_0 = new VectorFormVolWave_0(c_squared);
    vector_form_0->ext.push_back(B_prev_sln);
    add_vector_form(vector_form_0);

    VectorFormVolWave_1* vector_form_1 = new VectorFormVolWave_1();
    vector_form_1->ext.push_back(E_prev_sln);
    add_vector_form(vector_form_1);
  };

private:
  class MatrixFormVolWave_0_1 : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolWave_0_1(double c_squared) 
          : WeakForm::MatrixFormVol(0, 1, HERMES_NONSYM), c_squared(c_squared) { }

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
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    virtual WeakForm::MatrixFormVol* clone() {
      return new MatrixFormVolWave_0_1(*this);
    }

    // Member.
    double c_squared;
  };

  class MatrixFormVolWave_1_0 : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolWave_1_0() 
      : WeakForm::MatrixFormVol(1, 0, HERMES_NONSYM) { }

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
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
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
      Func<Scalar>* K_sln = u_ext[1];
      Func<Scalar>* sln_prev_time = ext->fn[0];

      for (int i = 0; i < n; i++) {
        Scalar sln_dx_i = sln_prev_time->dx[i] + K_sln->dx[i];
        Scalar sln_dy_i = sln_prev_time->dy[i] + K_sln->dy[i];
        result += wt[i] * (sln_dy_i * v->val0[i] - sln_dx_i * v->val1[i]);
      }
      return c_squared * result;
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
      Func<Scalar>* K_sln = u_ext[0];
      Func<Scalar>* sln_prev_time = ext->fn[0];
      
      for (int i = 0; i < n; i++) {
        Scalar sln_curl_i = sln_prev_time->curl[i] + K_sln->curl[i];
        result -= wt[i] * sln_curl_i * v->val[i];
      }
      return result;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                 Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    virtual WeakForm::VectorFormVol* clone() {
      return new VectorFormVolWave_1(*this);
    }
  };
};
