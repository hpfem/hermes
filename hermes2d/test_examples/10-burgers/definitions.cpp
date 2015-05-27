// ext[0] is the previous solution.
// u is the solution function, v is the test function currently in hand.
// u, v, and ext[0] have values (->val[i]), derivatives (->dx[i], dy[i]), [i] is the index of the integration point (there are n integration points, in this case just one)
// wt are integration weights
// ext[0]->fn_central, ext[0]->fn_neighbor: on internal edges, this represents the value on the currently assembled element (central), and the neighbor.
// 
#include "definitions.h"

// u = -x - y + t - 1
void CustomExactSolution::derivatives(double x, double y, double& dx, double& dy) const
{
  dx = -1.;
  dy = -1.;
};

double CustomExactSolution::value(double x, double y) const
{
  return -x - y + t - 1.;
}

Ord CustomExactSolution::ord(double x, double y) const
{
  return Hermes::Ord(0);
}

MeshFunction<double>* CustomExactSolution::clone() const
{
  return new CustomExactSolution(this->mesh);
}

CustomExactSolution::~CustomExactSolution()
{
}

void CustomExactSolution::setTime(double time)
{
  this->t = time;
}

CustomWeakForm::CustomWeakForm(MeshFunctionSharedPtr<double> previousSln, MeshSharedPtr mesh) : WeakForm<double>(1), previousSln(previousSln), mesh(mesh)
{
  this->set_ext(previousSln);
  add_matrix_form(new CustomMatrixFormVol(0, 0));
  add_vector_form(new CustomVectorFormVol(0));

  add_vector_form_DG(new CustomVectorFormInterface(0));
  add_vector_form_surf(new CustomVectorFormSurface(0));
}

WeakForm<double>* CustomWeakForm::clone() const
{
  return new CustomWeakForm(*this);
}

template<typename Real, typename Scalar>
Scalar CustomWeakForm::CustomMatrixFormVol::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v,
  GeomVol<Real> *e, Func<Scalar> **ext) const
{
  Scalar result = Scalar(0);
  for (int i = 0; i < n; i++)
    result += wt[i] * u->val[i] * v->val[i];
  return result / this->wf->get_current_time_step();
}

double CustomWeakForm::CustomMatrixFormVol::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
  GeomVol<double> *e, Func<double> **ext) const
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord CustomWeakForm::CustomMatrixFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
  GeomVol<Ord> *e, Func<Ord> **ext) const
{
  return Ord(0);
}

MatrixFormVol<double>* CustomWeakForm::CustomMatrixFormVol::clone() const
{
  return new CustomWeakForm::CustomMatrixFormVol(*this);
}

template<typename Real, typename Scalar>
Scalar CustomWeakForm::CustomVectorFormVol::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
  GeomVol<Real> *e, Func<Scalar> **ext) const
{
  Scalar result = Scalar(0);
  for (int i = 0; i < n; i++)
    result += wt[i] * ext[0]->val[i] * v->val[i];

  result /= this->wf->get_current_time_step();

  for (int i = 0; i < n; i++)
    result += wt[i] * 0.5 * ext[0]->val[i] * ext[0]->val[i] * (v->dx[i] + v->dy[i]);

  return result;
}

double CustomWeakForm::CustomVectorFormVol::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
  GeomVol<double> *e, Func<double> **ext) const
{
  return vector_form<double, double>(n, wt, u_ext, v, e, ext);
}

Ord CustomWeakForm::CustomVectorFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
  GeomVol<Ord> *e, Func<Ord> **ext) const
{
  return Ord(0);
}

VectorFormVol<double>* CustomWeakForm::CustomVectorFormVol::clone() const
{
  return new CustomWeakForm::CustomVectorFormVol(*this);
}

template<typename Real, typename Scalar>
Scalar CustomWeakForm::CustomVectorFormInterface::vector_form(int n, double *wt, DiscontinuousFunc<Scalar>** u_ext, Func<Real> *v,
  InterfaceGeom<Real> *e, DiscontinuousFunc<Scalar> **ext) const
{
  Scalar result = Scalar(0);

  for (int i = 0; i < n; i++)
  {
    result += wt[i] * 0.25 * (ext[0]->fn_central->val[i] * ext[0]->fn_central->val[i] + ext[0]->fn_neighbor->val[i] * ext[0]->fn_neighbor->val[i]) * v->val[i] * (e->nx[i] + e->ny[i]);
    
    result -= wt[i] * 0.5 * std::max(ext[0]->fn_central->val[i], ext[0]->fn_neighbor->val[i]) * (ext[0]->fn_neighbor->val[i] - ext[0]->fn_central->val[i]) * v->val[i];
  }
  return result;
}

double CustomWeakForm::CustomVectorFormInterface::value(int n, double *wt, DiscontinuousFunc<double> **u_ext, Func<double> *v,
  InterfaceGeom<double> *e, DiscontinuousFunc<double> **ext) const
{
  return vector_form<double, double>(n, wt, u_ext, v, e, ext);
}

Ord CustomWeakForm::CustomVectorFormInterface::ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, Func<Ord> *v,
  InterfaceGeom<Ord> *e, DiscontinuousFunc<Ord> **ext) const
{
  return Ord(0);
}

VectorFormDG<double>* CustomWeakForm::CustomVectorFormInterface::clone() const
{
  return new CustomWeakForm::CustomVectorFormInterface(*this);
}

double CustomWeakForm::CustomVectorFormSurface::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
  GeomSurf<double> *e, Func<double> **ext) const
{
  double result = 0;

  for (int i = 0; i < n; i++)
  {
    double exact_sln_value = -e->x[i] - e->y[i] + this->wf->get_current_time() - 1.;
    result += wt[i] * 0.25 * (ext[0]->val[i] * ext[0]->val[i] + exact_sln_value * exact_sln_value) * v->val[i] * (e->nx[i] + e->ny[i]);
    result -= wt[i] * 0.5 * std::max(ext[0]->val[i], exact_sln_value) * (exact_sln_value - ext[0]->val[i]) * v->val[i];
  }
  return result;
}

Ord CustomWeakForm::CustomVectorFormSurface::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, GeomSurf<Ord> *e, Func<Ord> **ext) const
{
  return Ord(0);
}

VectorFormSurf<double>* CustomWeakForm::CustomVectorFormSurface::clone() const
{
  return new CustomWeakForm::CustomVectorFormSurface(*this);
}