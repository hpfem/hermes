#include "weakform/weakform.h"
#include "weakform_library/h1.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

using namespace WeakFormsH1::VolumetricMatrixForms;
using namespace WeakFormsH1::VolumetricVectorForms;

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  virtual scalar value (double x, double y) const {
    return x*x +y*y;
  }

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    dx = 2*x;
    dy = 2*y;
  };

  virtual Ord ord(Ord x, Ord y) const {
    return x*x + y*y;
  }
};

/* Weak forms */

class WeakFormPoisson : public WeakForm
{
public:
  WeakFormPoisson() : WeakForm(1) {
    add_matrix_form(new DefaultLinearDiffusion(0, 0));
    add_vector_form(new DefaultVectorFormConst(0, -4.0));
  };
};

/* Weak forms for NOX */

class WeakFormPoissonNox : public WeakForm
{
public:
  WeakFormPoissonNox(bool is_matfree) : WeakForm(1) {
    this->is_matfree = is_matfree;
    add_matrix_form(new DefaultLinearDiffusion(0, 0));
    add_vector_form(new DefaultResidualLinearDiffusion(0));
    add_vector_form(new DefaultVectorFormConst(0, 4.0));
  };
};
