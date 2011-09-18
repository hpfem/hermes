#include "hermes2d.h"

using namespace Hermes::Hermes2D;
using namespace WeakFormsH1;
using Hermes::Ord;

/*  Exact solution */

class CustomExactSolution : public ExactSolutionScalar<double>
{
public:
  CustomExactSolution(Mesh* mesh) : ExactSolutionScalar<double>(mesh)
  {
  }

  virtual void derivatives(double x, double y, double& dx, double& dy) const;

  virtual double value(double x, double y) const;

  virtual Ord ord(Ord x, Ord y) const;
};

/* Custom function f */

class CustomFunction : public Hermes::Hermes2DFunction<double>
{
public:
  CustomFunction() : Hermes::Hermes2DFunction<double>()
  {
  }

  virtual double value(double x, double y) const;

  virtual Ord value(Ord x, Ord y) const;
};
