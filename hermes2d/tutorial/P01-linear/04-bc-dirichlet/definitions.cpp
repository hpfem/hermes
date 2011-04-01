#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform_library/h1.h"

using namespace WeakFormsH1::VolumetricMatrixForms;
using namespace WeakFormsH1::VolumetricVectorForms;

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh, double const_f) 
            : ExactSolutionScalar(mesh), const_f(const_f) {};

  double value(double x, double y) const {
    return (-const_f/4.0) * (x*x + y*y);
  }

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    dx = -const_f*x/2.0;
    dy = -const_f*y/2.0;
  };

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(2);
  }

  double const_f;
};

/* Weak forms */

class CustomWeakFormPoisson : public WeakForm
{
public:
  CustomWeakFormPoisson(double const_f) : WeakForm(1) {
    add_matrix_form(new DefaultLinearDiffusion(0, 0));
    add_vector_form(new DefaultVectorFormConst(0, const_f));
  };
};
