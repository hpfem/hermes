#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/boundaryconditions.h"

class WeakFormNeumann : public WeakForm
{
public:
  // Problem parameters.
  double const_f;

  WeakFormNeumann(double const_f,
                  double const_gamma_bottom,
                  double const_gamma_outer,
                  double const_gamma_left) : WeakForm(1), const_f(const_f)
  {
    add_matrix_form(new MatrixFormVolNeumann(0, 0));
    add_vector_form(new VectorFormVolNeumann(0));

    add_vector_form_surf(new VectorFormSurfNeumann(0, BDY_BOTTOM, const_gamma_bottom));
    add_vector_form_surf(new VectorFormSurfNeumann(0, BDY_OUTER, const_gamma_outer));
    add_vector_form_surf(new VectorFormSurfNeumann(0, BDY_LEFT, const_gamma_left));
  };

private:
  class MatrixFormVolNeumann : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolNeumann(int i, int j) : WeakForm::MatrixFormVol(i, j) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }
  };

  class VectorFormVolNeumann : public WeakForm::VectorFormVol
  {
  public:
    VectorFormVolNeumann(int i) : WeakForm::VectorFormVol(i) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      return static_cast<WeakFormNeumann *>(wf)->const_f * int_v<Real, Scalar>(n, wt, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }
  };

  class VectorFormSurfNeumann : public WeakForm::VectorFormSurf
  {
  private:
      double constant;
  public:
    VectorFormSurfNeumann(int i, std::string area, double constant) : WeakForm::VectorFormSurf(i, area), constant(constant) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      return constant * int_v<Real, Scalar>(n, wt, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
        return matrix_form_surf<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
        return matrix_form_surf<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }
  };
};
