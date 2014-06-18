#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace WeakFormsH1;
using namespace Hermes::Hermes2D::Views;
using Hermes::Ord;

/* Right-hand side */

class CustomRightHandSide : public Hermes::Hermes2DFunction<double>
{
public:
  CustomRightHandSide(double alpha_w, double alpha_p, double x_w, double y_w,
    double r_0, double omega_c, double epsilon, double x_p, double y_p)
    : Hermes::Hermes2DFunction<double>(), alpha_w(alpha_w), alpha_p(alpha_p),
    x_w(x_w), y_w(y_w), r_0(r_0), omega_c(omega_c), epsilon(epsilon), x_p(x_p), y_p(y_p) {};

  virtual double value(double x, double y) const;

  virtual Ord value(Ord x, Ord y) const;

  double alpha_w;
  double alpha_p;
  double x_w;
  double y_w;
  double r_0;
  double omega_c;
  double epsilon;
  double x_p;
  double y_p;
};

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar<double>
{
public:
  CustomExactSolution(MeshSharedPtr mesh, double alpha_w, double alpha_p, double x_w, double y_w,
    double r_0, double omega_c, double epsilon, double x_p, double y_p)
    : ExactSolutionScalar<double>(mesh), alpha_w(alpha_w), alpha_p(alpha_p),
    x_w(x_w), y_w(y_w), r_0(r_0), omega_c(omega_c), epsilon(epsilon), x_p(x_p), y_p(y_p) {};

  virtual double value(double x, double y) const;

  virtual void derivatives(double x, double y, double& dx, double& dy) const;

  virtual Ord ord(double x, double y) const;

  double get_angle(double y, double x) const;

  MeshFunction<double>* clone() const;

  double alpha_w;
  double alpha_p;
  double x_w;
  double y_w;
  double r_0;
  double omega_c;
  double epsilon;
  double x_p;
  double y_p;
};
