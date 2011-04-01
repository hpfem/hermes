#include "weakform/weakform.h"
#include "weakform_library/h1.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

//using namespace Laplace;
using namespace WeakFormsH1::VolumetricMatrixForms;
using namespace WeakFormsH1::VolumetricVectorForms;
//using namespace WeakFormsH1::SurfaceMatrixForms;
//using namespace WeakFormsH1::SurfaceVectorForms;
using namespace WeakFormsH1::RightHandSides;

/* Custom function that is used in the exact solution and in right-hand side */

class CustomExactFunction1
{
public:
  CustomExactFunction1() { };

  double val(double x) {
    return cos(M_PI*x/2);
  }
  double dx(double x) {
    return -sin(M_PI*x/2)*(M_PI/2.);
  }
  double ddxx(double x) {
    return -cos(M_PI*x/2)*(M_PI/2.)*(M_PI/2.);
  }
};

class CustomExactFunction2
{
public:
  CustomExactFunction2(double K) : K(K) {};

  double val(double x) {
    return 1. - (exp(K*x) + exp(-K*x))/(exp(K) + exp(-K));
  }
  double dx(double x) {
    return -K*(exp(K*x) - exp(-K*x))/(exp(K) + exp(-K));
  }
  double ddxx(double x) {
    return -K*K*(exp(K*x) + exp(-K*x))/(exp(K) + exp(-K));
  }
  
  double K;
};

/* Right-hand side */

class CustomRightHandSide1: public DefaultNonConstRightHandSide
{
public:
  CustomRightHandSide1(double K, double d_u, double sigma) 
    : DefaultNonConstRightHandSide(), d_u(d_u), sigma(sigma) {
    cef1 = new CustomExactFunction1();
    cef2 = new CustomExactFunction2(K);
  };

  virtual scalar value(double x, double y) const {
    double Laplace_u = cef1->ddxx(x) * cef1->val(y) 
                       + cef1->val(x) * cef1->ddxx(y);
    double u = cef1->val(x) * cef1->val(y);
    double v = cef2->val(x) * cef2->val(y);
    return -d_u * d_u * Laplace_u - u + sigma * v;
  }

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(20);
  }

  ~CustomRightHandSide1() { delete cef1; delete cef2;}

  CustomExactFunction1* cef1;
  CustomExactFunction2* cef2;
  double d_u, sigma;
};

class CustomRightHandSide2: public DefaultNonConstRightHandSide
{
public:
  CustomRightHandSide2(double K, double d_v) 
    : DefaultNonConstRightHandSide(), d_v(d_v) {
    cef1 = new CustomExactFunction1();
    cef2 = new CustomExactFunction2(K);
  };

  virtual scalar value(double x, double y) const {
    double Laplace_v = cef2->ddxx(x) * cef2->val(y) 
                       + cef2->val(x) * cef2->ddxx(y);
    double u = cef1->val(x) * cef1->val(y);
    double v = cef2->val(x) * cef2->val(y);
    return -d_v*d_v * Laplace_v - u + v;
  }

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(20);
  }

  ~CustomRightHandSide2() { delete cef1; delete cef2;}

  CustomExactFunction1* cef1;
  CustomExactFunction2* cef2;
  double d_v;
};

/* Exact solution */

class ExactSolutionFitzHughNagumo1 : public ExactSolutionScalar
{
public:
  ExactSolutionFitzHughNagumo1(Mesh* mesh) 
       : ExactSolutionScalar(mesh) {
    cef1 = new CustomExactFunction1();
  }

  virtual scalar value (double x, double y) const {
    return cef1->val(x)*cef1->val(y);
  }

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    dx = cef1->dx(x)*cef1->val(y);
    dy = cef1->val(x)*cef1->ddxx(y);
  }

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(10);
  }

  ~ExactSolutionFitzHughNagumo1() {
    delete cef1;
  }

  CustomExactFunction1* cef1;
};

class ExactSolutionFitzHughNagumo2 : public ExactSolutionScalar
{
public:
  ExactSolutionFitzHughNagumo2(Mesh* mesh, double K) 
       : ExactSolutionScalar(mesh) {
    cef2 = new CustomExactFunction2(K);
  }
  virtual scalar value (double x, double y) const {
    return cef2->val(x)*cef2->val(y);
  }

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    dx = cef2->dx(x)*cef2->val(y);
    dy = cef2->val(x)*cef2->dx(y);
  }

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(20);
  }

  ~ExactSolutionFitzHughNagumo2() {
    delete cef2;
  }

  CustomExactFunction2* cef2;
};

/* Weak forms */

class WeakFormFitzHughNagumo : public WeakForm
{
public:
  WeakFormFitzHughNagumo(CustomRightHandSide1* rhs_1, CustomRightHandSide2* rhs_2)
          : WeakForm(2) {
    add_matrix_form(new DefaultLinearDiffusion(0, 0, D_u * D_u));
    add_matrix_form(new DefaultLinearMass(0, 0, -1.0));
    add_matrix_form(new DefaultLinearMass(0, 1, rhs_1->sigma, HERMES_NONSYM));
    add_matrix_form(new DefaultLinearMass(1, 0, -1.0, HERMES_NONSYM));     
    add_matrix_form(new DefaultLinearDiffusion(1, 1, D_v * D_v));
    add_matrix_form(new DefaultLinearMass(1, 1, 1.0));
    

    add_vector_form(new DefaultVectorFormNonConst(0, rhs_1));
    add_vector_form(new DefaultVectorFormNonConst(1, rhs_2));
  }
};
