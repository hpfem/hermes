#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/boundaryconditions.h"

class WeakFormLameEquations : public WeakForm
{
public:
  // Problem parameters.
  double const_f;

  WeakFormLameEquations(double mu, double lambda, std::string non_zero_neumann_bnd, double f_0, double f_1) : WeakForm(2) {
    add_matrix_form(new MatrixFormVolLameEquations_0_0(mu, lambda));
    add_matrix_form(new MatrixFormVolLameEquations_0_1(mu, lambda));
    add_matrix_form(new MatrixFormVolLameEquations_1_1(mu, lambda));

    add_vector_form_surf(new VectorFormSurfLameEquations_0(non_zero_neumann_bnd, f_0));
    add_vector_form_surf(new VectorFormSurfLameEquations_1(non_zero_neumann_bnd, f_1));
  };

private:
  class MatrixFormVolLameEquations_0_0 : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolLameEquations_0_0(double mu, double lambda) : WeakForm::MatrixFormVol(0, 0), mu(mu), lambda(lambda) { 
      sym = HERMES_SYM;
    }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      return (lambda + 2*mu) * int_dudx_dvdx<Real, Scalar>(n, wt, u, v) +
                      mu * int_dudy_dvdy<Real, Scalar>(n, wt, u, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    // Members.
    double mu;
    double lambda;
  };

  class MatrixFormVolLameEquations_0_1 : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolLameEquations_0_1(double mu, double lambda) : WeakForm::MatrixFormVol(0, 1), mu(mu), lambda(lambda) { 
      sym = HERMES_SYM;
    }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      return lambda * int_dudy_dvdx<Real, Scalar>(n, wt, u, v) +
             mu * int_dudx_dvdy<Real, Scalar>(n, wt, u, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    // Members.
    double mu;
    double lambda;
  };

  class MatrixFormVolLameEquations_1_1 : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolLameEquations_1_1(double mu, double lambda) : WeakForm::MatrixFormVol(1, 1), mu(mu), lambda(lambda) { 
      sym = HERMES_SYM;
    }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      return mu * int_dudx_dvdx<Real, Scalar>(n, wt, u, v) +
         (lambda + 2*mu) * int_dudy_dvdy<Real, Scalar>(n, wt, u, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    // Members.
    double mu;
    double lambda;
  };

  class VectorFormSurfLameEquations_0 : public WeakForm::VectorFormSurf
  {
  public:
    VectorFormSurfLameEquations_0(std::string marker, double f_0) : WeakForm::VectorFormSurf(0, marker), f_0(f_0) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      return f_0 * int_v<Real, Scalar>(n, wt, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    // Member.
    double f_0;
  };

  class VectorFormSurfLameEquations_1 : public WeakForm::VectorFormSurf
  {
  public:
    VectorFormSurfLameEquations_1(std::string marker, double f_1) : WeakForm::VectorFormSurf(1, marker), f_1(f_1) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      return f_1 * int_v<Real, Scalar>(n, wt, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    // Member.
    double f_1;
  };
};