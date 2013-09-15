#include "definitions.h"

double CustomExactFunction1::val(double x) 
{
  return Hermes::cos(M_PI*x/2);
}
  
double CustomExactFunction1::dx(double x) 
{
  return -Hermes::sin(M_PI*x/2)*(M_PI/2.);
}
  
double CustomExactFunction1::ddxx(double x) 
{
  return -Hermes::cos(M_PI*x/2)*(M_PI/2.)*(M_PI/2.);
}


double CustomExactFunction2::val(double x) 
{
  return 1. - (exp(K*x) + exp(-K*x))/(exp(K) + exp(-K));
}
  
double CustomExactFunction2::dx(double x) 
{
  return -K*(exp(K*x) - exp(-K*x))/(exp(K) + exp(-K));
}
  
double CustomExactFunction2::ddxx(double x) 
{
  return -K*K*(exp(K*x) + exp(-K*x))/(exp(K) + exp(-K));
}


CustomRightHandSide1::CustomRightHandSide1(double K, double d_u, double sigma)
  : Hermes2DFunction<double>(), d_u(d_u), sigma(sigma) 
{
  cef1 = new CustomExactFunction1();
  cef2 = new CustomExactFunction2(K);
}

double CustomRightHandSide1::value(double x, double y) const 
{
  double Laplace_u = cef1->ddxx(x) * cef1->val(y)
                     + cef1->val(x) * cef1->ddxx(y);
  double u = cef1->val(x) * cef1->val(y);
  double v = cef2->val(x) * cef2->val(y);
  return -d_u * d_u * Laplace_u - u + sigma * v;
}

Ord CustomRightHandSide1::value(Ord x, Ord y) const 
{
  return Ord(10);
}

CustomRightHandSide1::~CustomRightHandSide1() 
{ 
  delete cef1; 
  delete cef2;
}

CustomRightHandSide2::CustomRightHandSide2(double K, double d_v)
      : Hermes2DFunction<double>(), d_v(d_v) 
{
  cef1 = new CustomExactFunction1();
  cef2 = new CustomExactFunction2(K);
}

double CustomRightHandSide2::value(double x, double y) const 
{
  double Laplace_v = cef2->ddxx(x) * cef2->val(y)
                     + cef2->val(x) * cef2->ddxx(y);
  double u = cef1->val(x) * cef1->val(y);
  double v = cef2->val(x) * cef2->val(y);
  return -d_v*d_v * Laplace_v - u + v;
}

Ord CustomRightHandSide2::value(Ord x, Ord y) const 
{
  return Ord(10);
}

CustomRightHandSide2::~CustomRightHandSide2() 
{ 
  delete cef1; 
  delete cef2;
}


ExactSolutionFitzHughNagumo1::ExactSolutionFitzHughNagumo1(MeshSharedPtr mesh)
     : ExactSolutionScalar<double>(mesh) 
{
  cef1 = new CustomExactFunction1();
}

double ExactSolutionFitzHughNagumo1::value(double x, double y) const 
{
  return cef1->val(x)*cef1->val(y);
}

void ExactSolutionFitzHughNagumo1::derivatives(double x, double y, double& dx, double& dy) const 
{
  dx = cef1->dx(x)*cef1->val(y);
  dy = cef1->val(x)*cef1->dx(y);
}

Ord ExactSolutionFitzHughNagumo1::ord(double x, double y) const 
{
  return Ord(10);
}

ExactSolutionFitzHughNagumo1::~ExactSolutionFitzHughNagumo1() 
{
  delete cef1;
}

MeshFunction<double>* ExactSolutionFitzHughNagumo1::clone() const
{
	return new ExactSolutionFitzHughNagumo1(this->mesh);
}


ExactSolutionFitzHughNagumo2::ExactSolutionFitzHughNagumo2(MeshSharedPtr mesh, double K)
     : ExactSolutionScalar<double>(mesh), K(K)
{
  cef2 = new CustomExactFunction2(K);
}

MeshFunction<double>* ExactSolutionFitzHughNagumo2::clone() const
{
	return new ExactSolutionFitzHughNagumo2(this->mesh, this->K);
}

double ExactSolutionFitzHughNagumo2::value(double x, double y) const 
{
  return cef2->val(x)*cef2->val(y);
}

void ExactSolutionFitzHughNagumo2::derivatives(double x, double y, double& dx, double& dy) const 
{
  dx = cef2->dx(x)*cef2->val(y);
  dy = cef2->val(x)*cef2->dx(y);
}

Ord ExactSolutionFitzHughNagumo2::ord(double x, double y) const 
{
  return Ord(10);
}

ExactSolutionFitzHughNagumo2::~ExactSolutionFitzHughNagumo2() 
{
  delete cef2;
}

double CustomResidual1::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                              Geom<double> *e, Func<double> **ext) const
{
   double result = 0;
   for (int i = 0; i < n; i++) 
   {
     result += wt[i] * (    d_u*d_u * (u_ext[0]->dx[i]*v->dx[i] + u_ext[0]->dy[i]*v->dy[i]) 
                          - u_ext[0]->val[i]*v->val[i] 
                          + sigma*u_ext[1]->val[i]*v->val[i]
                          - g1->value(e->x[i], e->y[i])*v->val[i]
                       );
   }
 
   return result;
}

Ord CustomResidual1::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                         Geom<Ord> *e, Func<Ord> **ext) const 
{
  return Ord(20);
}

VectorFormVol<double>* CustomResidual1::clone() const 
{
  return new CustomResidual1(this->d_u, this->sigma, new CustomRightHandSide1(g1->cef2->K, this->d_u, this->sigma));
}

double CustomResidual2::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                              Geom<double> *e, Func<double> **ext) const
{
   double result = 0;
   for (int i = 0; i < n; i++) 
   {
     result += wt[i] * (    d_v*d_v * (u_ext[1]->dx[i]*v->dx[i] + u_ext[1]->dy[i]*v->dy[i]) 
                          - u_ext[0]->val[i]*v->val[i] 
                          + u_ext[1]->val[i]*v->val[i]
                          - g2->value(e->x[i], e->y[i])*v->val[i]
                       );
   }
 
   return result;
  }

Ord CustomResidual2::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                         Geom<Ord> *e, Func<Ord> **ext) const 
{
  return Ord(20);
} 

VectorFormVol<double>* CustomResidual2::clone() const 
{
  return new CustomResidual2(this->d_v, new CustomRightHandSide2(g2->cef2->K, g2->d_v));
}

CustomWeakForm::CustomWeakForm(CustomRightHandSide1* g1, CustomRightHandSide2* g2) : WeakForm<double>(2) 
{
  // Jacobian.
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<double>(0, 0, HERMES_ANY, new Hermes1DFunction<double>(g1->d_u * g1->d_u)));
  add_matrix_form(new WeakFormsH1::DefaultMatrixFormVol<double>(0, 0, HERMES_ANY, new Hermes2DFunction<double>(-1.0)));
  add_matrix_form(new WeakFormsH1::DefaultMatrixFormVol<double>(0, 1, HERMES_ANY, new Hermes2DFunction<double>(g1->sigma), HERMES_NONSYM));
  add_matrix_form(new WeakFormsH1::DefaultMatrixFormVol<double>(1, 0, HERMES_ANY, new Hermes2DFunction<double>(-1.0), HERMES_NONSYM));
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<double>(1, 1, HERMES_ANY, new Hermes1DFunction<double>(g2->d_v * g2->d_v)));
  add_matrix_form(new WeakFormsH1::DefaultMatrixFormVol<double>(1, 1, HERMES_ANY, new Hermes2DFunction<double>(1.0)));

  // Residual.
  add_vector_form(new CustomResidual1(g1->d_u, g1->sigma, g1));
  add_vector_form(new CustomResidual2(g2->d_v, g2));
}
