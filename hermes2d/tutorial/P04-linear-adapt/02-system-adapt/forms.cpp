#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

// Exact solution.
#include "exact_solution.cpp"

class WeakFormFitzHughNagumo : public WeakForm
{
public:
  WeakFormFitzHughNagumo(ExactSolutionFitzHughNagumo1* exact_solution1, ExactSolutionFitzHughNagumo2* exact_solution2)
          : WeakForm(2) {
    add_matrix_form(new MatrixFormVolFitzHughNagumo_0_0(exact_solution1->d_u));
    add_matrix_form(new MatrixFormVolFitzHughNagumo_0_1(exact_solution1->sigma));
    add_matrix_form(new MatrixFormVolFitzHughNagumo_1_0);
    add_matrix_form(new MatrixFormVolFitzHughNagumo_1_1(exact_solution2->d_v));

    add_vector_form(new VectorFormVolFitzHughNagumo_1(exact_solution1, exact_solution2));
    add_vector_form(new VectorFormVolFitzHughNagumo_2(exact_solution1, exact_solution2));
  }

private:
  class MatrixFormVolFitzHughNagumo_0_0 : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolFitzHughNagumo_0_0(double D_u) : WeakForm::MatrixFormVol(0, 0), d_u(d_u) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      return D_u * D_u * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v) - int_u_v<Real, Scalar>(n, wt, u, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext)
    {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
    {
       return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    // Member.
    double d_u;
  };

  class MatrixFormVolFitzHughNagumo_0_1 : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolFitzHughNagumo_0_1(double sigma) : WeakForm::MatrixFormVol(0, 1), sigma(sigma) { 
    }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      return SIGMA * int_u_v<Real, Scalar>(n, wt, u, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
       return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    // Member.
    double sigma;
  };

  class MatrixFormVolFitzHughNagumo_1_0 : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolFitzHughNagumo_1_0() : WeakForm::MatrixFormVol(1, 0) { 
    }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      return -int_u_v<Real, Scalar>(n, wt, u, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
       return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }
  };

  class MatrixFormVolFitzHughNagumo_1_1 : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolFitzHughNagumo_1_1(double d_v) : WeakForm::MatrixFormVol(1, 1), d_v(d_v) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      return D_v * D_v * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v) + int_u_v<Real, Scalar>(n, wt, u, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext)
    {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
    {
       return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    // Member.
    double d_v;
  };

  class VectorFormVolFitzHughNagumo_1 : public WeakForm::VectorFormVol
  {
  public:
    VectorFormVolFitzHughNagumo_1(ExactSolutionFitzHughNagumo1* exact_solution1, ExactSolutionFitzHughNagumo2* exact_solution2)
      : WeakForm::VectorFormVol(0), exact_solution1(exact_solution1), exact_solution2(exact_solution2), d_u(exact_solution1->d_u) {}
    
    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      Scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * (g_1(e->x[i], e->y[i]) * v->val[i]);
      return result;
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return Ord(30);
    }

    double g_1(double x, double y) {
      double Laplace_u = exact_solution1->ddUdtt(x) * exact_solution1->U(y) + exact_solution1->U(x) * exact_solution1->ddUdtt(y);
      double u = exact_solution1->U(x) * exact_solution1->U(y);
      double v = exact_solution2->V(x) * exact_solution2->V(y);
      return -d_u * d_u * Laplace_u - u + SIGMA * v;
    }

    // Member.
    ExactSolutionFitzHughNagumo1 * exact_solution1;
    ExactSolutionFitzHughNagumo2 * exact_solution2;
    double d_u;
  };

  class VectorFormVolFitzHughNagumo_2 : public WeakForm::VectorFormVol
  {
  public:
    VectorFormVolFitzHughNagumo_2(ExactSolutionFitzHughNagumo1* exact_solution1, ExactSolutionFitzHughNagumo2* exact_solution2)
          : WeakForm::VectorFormVol(1), exact_solution1(exact_solution1), exact_solution2(exact_solution2), d_v(exact_solution2->d_v) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      Scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * (g_2(e->x[i], e->y[i]) * v->val[i]);
      return result;
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return Ord(30);
    }

    double g_2(double x, double y) {
      double Laplace_v = exact_solution2->ddVdtt(x) * exact_solution2->V(y) + exact_solution2->V(x) * exact_solution2->ddVdtt(y);
      double u = exact_solution1->U(x) * exact_solution1->U(y);
      double v = exact_solution2->V(x) * exact_solution2->V(y);
      return -d_v*d_v * Laplace_v - u + v;
    }

    // Member.
    ExactSolutionFitzHughNagumo1 * exact_solution1;
    ExactSolutionFitzHughNagumo2 * exact_solution2;
    double d_v;
  };
};