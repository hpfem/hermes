#include "weakform/weakform.h"
#include "weakform_library/h1.h"
#include "integrals/h1.h"
#include "boundaryconditions/essential_bcs.h"

using namespace WeakFormsH1;

/* Weak forms */

// NOTE: The linear problem in the Picard's method is 
//       solved using the Newton's method, but this is 
//       not a Newton's method.

// Thermal conductivity (temperature-dependent)
// For any u, this function has to be positive.
template<typename Real>
Real lam(Real u) {
  return 1 + pow(u, 4);
}

class CustomWeakFormPicard : public WeakForm
{
public:
  CustomWeakFormPicard(Solution* prev_iter_sln, double heat_src) : WeakForm(1) 
  {
    // Jacobian.
    CustomJacobian* matrix_form = new CustomJacobian(0, 0);
    matrix_form->ext.push_back(prev_iter_sln);
    add_matrix_form(matrix_form);

    // Residual.
    CustomResidual* vector_form = new CustomResidual(0, heat_src);
    vector_form->ext.push_back(prev_iter_sln);
    add_vector_form(vector_form);
  };

private:
  class CustomJacobian : public WeakForm::MatrixFormVol
  {
  public:
    CustomJacobian(int i, int j) : WeakForm::MatrixFormVol(i, j) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      Scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * lam<Real>(ext->fn[0]->val[i]) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);

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
  };

  class CustomResidual : public WeakForm::VectorFormVol
  {
  public:
    CustomResidual(int i, double heat_src) : WeakForm::VectorFormVol(i), heat_src(heat_src) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      Scalar result = 0;
      for (int i = 0; i < n; i++) {
        result += wt[i] * lam<Real>(ext->fn[0]->val[i]) * (u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i]);
        result -= wt[i] * heat_src * v->val[i];
      }

      return result;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[],
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

  private:
    double heat_src;
  };
};
