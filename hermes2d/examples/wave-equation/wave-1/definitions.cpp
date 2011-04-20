#include "weakform/weakform.h"
#include "weakform_library/h1.h"
#include "integrals/h1.h"
#include "boundaryconditions/essential_bcs.h"

using namespace WeakFormsH1::VolumetricMatrixForms;

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

  CustomWeakFormWave(double tau, double c_squared, Solution* u_prev_sln, 
                     Solution* v_prev_sln) : WeakForm(2) {
    // Volumetric matrix forms.
    add_matrix_form(new DefaultLinearMass(0, 1, 1.0, HERMES_NONSYM));
    add_matrix_form(new DefaultLinearDiffusion(1, 0, -c_squared, HERMES_NONSYM));

    // Volumetric surface forms.
    add_vector_form(new VectorFormVolWave_0());
    add_vector_form(new VectorFormVolWave_1(c_squared));
  };

private:
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
      return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
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
      return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
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
