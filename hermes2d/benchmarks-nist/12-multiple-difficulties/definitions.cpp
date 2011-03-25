#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform_library/laplace.h"

/* Right-hand side */

class CustomRightHandSide: public DefaultNonConstRightHandSide
{
public:
  CustomRightHandSide(double alpha_p, double x_p, double y_p, double alpha_w, double x_w, 
                      double y_w, double omega_c, double r_0, double epsilon)
    : DefaultNonConstRightHandSide(),
      alpha_p(alpha_p), x_p(x_p), y_p(y_p), alpha_w(alpha_w), x_w(x_w), y_w(y_w), 
      omega_c(omega_c), r_0(r_0), epsilon(epsilon) { }

  double value(double x, double y) const {
    //For more elegant form please execute file "generate_rhs.py" 

    double a_P = (-alpha_p * pow((x - x_p), 2) - alpha_p * pow((y - y_p), 2));

    double a_W = pow(x - x_w, 2);
    double b_W = pow(y - y_w, 2);
    double c_W = sqrt(a_W + b_W);
    double d_W = ((alpha_w * x - (alpha_w * x_w)) * (2 * x - (2 * x_w)));
    double e_W = ((alpha_w * y - (alpha_w * y_w)) * (2 * y - (2 * y_w)));
    double f_W = (pow(alpha_w * c_W - (alpha_w * r_0), 2) + 1.0);
    double g_W = (alpha_w * c_W - (alpha_w * r_0));

    return -(4 * exp(a_P) * alpha_p * (alpha_p * (x - x_p) * (x - x_p) + alpha_p * (y - y_p) * (y - y_p) - 1)
           + ((alpha_w/(c_W * f_W)) - (d_W/(2 * pow(a_W + b_W, 1.5) * f_W)) - ((alpha_w * d_W * g_W)/((a_W + b_W) * pow(f_W, 2))) 
           + (alpha_w/(c_W * f_W)) - (e_W/(2 * pow(a_W + b_W, 1.5) * f_W)) - ((alpha_w * e_W * g_W)/((a_W + b_W) * pow(f_W, 2))))
	     + (1.0 / epsilon) * (1.0 / epsilon) * exp(-(1 + y) / epsilon));  
  }

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(10);
  }

  double alpha_p, x_p, y_p, alpha_w, x_w, y_w, omega_c, r_0, epsilon;
};

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh, double alpha_p, double x_p, double y_p, double alpha_w, 
                  double x_w, double y_w, double omega_c, double r_0, double epsilon)
    : ExactSolutionScalar(mesh), alpha_p(alpha_p), x_p(x_p), y_p(y_p), alpha_w(alpha_w), 
      x_w(x_w), y_w(y_w), omega_c(omega_c), r_0(r_0), epsilon(epsilon) { }

  double get_angle(double y, double x) const {
    double theta = atan2(y, x);
    if (theta < 0)
      theta += 2 * M_PI;
    return theta;
  }

  virtual double value(double x, double y) const {
    double ALPHA_C = (M_PI/ omega_c);

    return exp(-alpha_p * (pow((x - x_p), 2) + pow((y - y_p), 2)))
           + (pow(sqrt(x*x + y*y), ALPHA_C) * sin(ALPHA_C * get_angle(y, x)))
           + atan(alpha_w * (sqrt(pow(x - x_w, 2) + pow(y - y_w, 2)) - r_0))
           + exp(-(1 + y) / epsilon);
  };

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    double a_P = -alpha_p * ( (x - x_p) * (x - x_p) + (y - y_p) * (y - y_p));

    double ALPHA_C = (M_PI/ omega_c);
    double a_C = sqrt(x*x + y*y);
    double b_C = pow(a_C, (ALPHA_C - 1.0));
    double c_C = pow(a_C, ALPHA_C);
    double d_C = ((y*y)/(x*x) + 1.0 );

    double a_W = pow(x - x_w, 2);
    double b_W = pow(y - y_w, 2);
    double c_W = sqrt(a_W + b_W);
    double d_W = (alpha_w * x - (alpha_w * x_w));
    double e_W = (alpha_w * y - (alpha_w * y_w));
    double f_W = (pow(alpha_w * c_W - (alpha_w * r_0), 2) + 1.0);

    dx = -exp(a_P) * (2 * alpha_p * (x - x_p))
         + (((ALPHA_C* x* sin(ALPHA_C * get_angle(y,x)) *b_C)/a_C) 
         - ((ALPHA_C *y *cos(ALPHA_C * get_angle(y, x)) * c_C)/(pow(x, 2.0) *d_C)))
         + (d_W / (c_W * f_W));
    dy = -exp(a_P) * (2 * alpha_p * (y - y_p))
         + (((ALPHA_C* cos(ALPHA_C* get_angle(y, x)) *c_C)/(x * d_C)) 
         + ((ALPHA_C* y* sin(ALPHA_C* get_angle(y, x)) *b_C)/a_C))
         + (e_W / (c_W * f_W))
         + (-1) * (1.0 / epsilon) * exp(-(1 + y) / epsilon); 
  };

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(10);
  }

  double alpha_p, x_p, y_p, alpha_w, x_w, y_w, omega_c, r_0;
  double epsilon;
};

/* Weak forms */

class CustomWeakFormPoisson : public WeakForm
{
public:
  CustomWeakFormPoisson(DefaultNonConstRightHandSide* rhs) : WeakForm(1) {
    add_matrix_form(new DefaultMatrixFormStiffness(0, 0));
    add_vector_form(new DefaultVectorFormVolNonConst(0, rhs));
  };
};

