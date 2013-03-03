#include "definitions.h"

const std::string MatrixSolverNames[6] = {
  "UMFPACK",
  "PETSc",
  "MUMPS",
  "SuperLU",
  "Trilinos/Amesos",
  "Trilinos/AztecOO"
};

double CustomRightHandSide::value(double x, double y) const
{
  return - kx(x, y) * dudx(x, y) - ky(x, y) * dudy(x, y) - k(x, y) * (dudxx(x, y) + dudyy(x, y));
}

Ord CustomRightHandSide::ord(Ord x, Ord y) const
{
  return - kx(x, y) * dudx(x, y) - ky(x, y) * dudy(x, y) - k(x, y) * (dudxx(x, y) + dudyy(x, y));
}

template<typename Real>
Real CustomRightHandSide::u(Real x, Real y) const
{
  return (x - x*x) * (y - y*y);
}

template<typename Real>
Real CustomRightHandSide::dudx(Real x, Real y) const
{
  return (1- 2*x) * y * (1 - y);
}

template<typename Real>
Real CustomRightHandSide::dudy(Real x, Real y) const
{
  return (1- 2*y) * x * (1 - x);
}

template<typename Real>
Real CustomRightHandSide::dudxx(Real x, Real y) const
{
  return -2.0 * (y-y*y);
}

template<typename Real>
Real CustomRightHandSide::dudyy(Real x, Real y) const
{
  return -2.0 * (x-x*x);
}

template<typename Real>
Real CustomRightHandSide::dudxy(Real x, Real y) const
{
  return (1- 2*y) * (1 - 2*x);
}

template<typename Real>
Real CustomRightHandSide::k(Real x, Real y) const
{
  return 1.0 / Hermes::sqrt(1.0 + sqr(dudx(x, y)) + sqr(dudy(x, y)));
}

template<typename Real>
Real CustomRightHandSide::kx(Real x, Real y) const
{
  return -0.5 * Hermes::pow(1.0 + sqr(dudx(x, y)) + sqr(dudy(x, y)), -1.5) *
                   (2.0 * dudx(x, y) * dudxx(x, y) + 2.0 * dudy(x, y) * dudxy(x, y));
}

template<typename Real>
Real CustomRightHandSide::ky(Real x, Real y) const
{
  return -0.5 * Hermes::pow(1.0 + sqr(dudx(x, y)) + sqr(dudy(x, y)), -1.5) *
                   (2.0 * dudx(x, y) * dudxy(x, y) + 2.0 * dudy(x, y) * dudyy(x, y));
}

double CustomExactSolution::value(double x, double y) const
{
  return  x * y * (1-x) * (1-y);
}

void CustomExactSolution::derivatives (double x, double y, double& dx, double& dy) const
{
  dx = (1- 2*x) * y * (1 - y);
  dy = (1- 2*y) * x * (1 - x);
}

Ord CustomExactSolution::ord(Ord x, Ord y) const
{
  return (1- 2*x) * y * (1 - y);
}

CustomWeakForm::CustomWeakForm(bool JFNK, bool precondition_jacobian, bool precondition_jacobian_approx) : WeakForm<double>(1, JFNK)
{
  // Jacobian forms - volumetric.
  if(!JFNK || precondition_jacobian)
    add_matrix_form(new JacobianFormVol(0, 0));

  // Residual forms - volumetric.
  add_vector_form(new ResidualFormVol(0, new CustomRightHandSide()));

  // Preconditioning form.
  if(JFNK && precondition_jacobian_approx)
    add_matrix_form(new PrecondFormVol(0, 0));
}

double CustomWeakForm::JacobianFormVol::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
                                              Func<double> *v, Geom<double> *e, Func<double> **ext) const
{
  double result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * ( -2.75 * Hermes::pow(1.0 + sqr(u_ext[0]->dx[i]) + sqr(u_ext[0]->dy[i]), -1.5) *
                       (2.0 * u_ext[0]->dx[i] * u->dx[i] + 2.0 * u_ext[0]->dx[i] * u->dx[i])
                       * (u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i]) +
                       (Hermes::pow(1.0 + sqr(u_ext[0]->dx[i]) + sqr(u_ext[0]->dy[i]), -0.5))
                       * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]) );
  return result;
}

Ord CustomWeakForm::JacobianFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                         Geom<Ord> *e, Func<Ord> **ext) const
{
  // Returning the sum of the degrees of the basis and test function plus two.
  return Ord(10);
}

double CustomWeakForm::ResidualFormVol::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                                              Geom<double> *e, Func<double> **ext) const
{
  Func<double>* u = u_ext[0];
  double result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * ((Hermes::pow(1.0 + sqr(u->dx[i]) + sqr(u->dy[i]), -0.5)) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i])
                       - rhs->value(e->x[i], e->y[i]) * v->val[i] );
  return result;
}

Ord CustomWeakForm::ResidualFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                                         Geom<Ord> *e, Func<Ord> **ext) const
{
  // Returning the sum of the degrees of the test function and solution plus two.
  return Ord(10);
}

double CustomWeakForm::PrecondFormVol::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
                                             Func<double> *v, Geom<double> *e, Func<double> **ext) const
{
  double result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * ( u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
  return result;
}

Ord CustomWeakForm::PrecondFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                        Geom<Ord> *e, Func<Ord> **ext) const
{
  // Returning the sum of the degrees of the basis and test function plus two.
  return Ord(10);
}