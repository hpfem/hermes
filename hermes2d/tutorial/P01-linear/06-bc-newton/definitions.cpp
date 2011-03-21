#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

class WeakFormNewton : public WeakForm
{
public:
  // Problem parameters.
  WeakFormNewton(double h, double T0, std::string natural_bc_bnd_part) : WeakForm(1)
  {
    add_matrix_form(new MatrixFormVolNewton(0, 0));

    add_matrix_form_surf(new MatrixFormSurfNewton(0, 0, natural_bc_bnd_part, h));
    add_vector_form_surf(new VectorFormSurfNewton(0, natural_bc_bnd_part, h, T0));
  };

private:
  class MatrixFormVolNewton : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolNewton(int i, int j) : WeakForm::MatrixFormVol(i, j) { }

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

  class MatrixFormSurfNewton : public WeakForm::MatrixFormSurf
  {
  private:
      double h;
  public:
    MatrixFormSurfNewton(int i, int j, std::string area, double h) : WeakForm::MatrixFormSurf(i, j, area), h(h) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      return h * int_u_v<Real, Scalar>(n, wt, u, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
        return matrix_form_surf<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
        return matrix_form_surf<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }
  };

  class VectorFormSurfNewton : public WeakForm::VectorFormSurf
  {
  private:
      double h;
  public:
    VectorFormSurfNewton(int i, std::string area, double h, double T0) : WeakForm::VectorFormSurf(i, area), h(h) { }

    template<typename Real, typename Scalar>
    Scalar vector_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      return h * T0 * int_v<Real, Scalar>(n, wt, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
        return vector_form_surf<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
        return vector_form_surf<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }
  };
};

