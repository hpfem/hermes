#include "weakform/weakform.h"
#include "weakform_library/h1.h"
#include "integrals/h1.h"
#include "boundaryconditions/essential_bcs.h"

using namespace WeakFormsH1::VolumetricMatrixForms;
using namespace WeakFormsH1::VolumetricVectorForms;
using namespace WeakFormsH1::SurfaceVectorForms;

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
    CustomExactSolution(Mesh* mesh) : ExactSolutionScalar(mesh) { };

    virtual void derivatives(double x, double y, scalar& dx, scalar& dy) const {
        dx = cos(x);
        dy = 0;
    };

    virtual double value(double x, double y) const {
        return sin(x);
    }

    virtual Ord ord(Ord x, Ord y) const {
        return Ord(7);
    }
};

/* Right-hand side */

class CustomRightHandSide: public DefaultFunction
{
public:
    CustomRightHandSide() : DefaultFunction() { }

    virtual scalar value(double x, double y) const {
        return sin(x);
    };

    virtual Ord ord(Ord x, Ord y) const {
        return Ord(7);
    }
};

/* Weak forms */

class CustomWeakFormPoisson : public WeakForm
{
public:
    CustomWeakFormPoisson() : WeakForm(1) {
        add_matrix_form(new DefaultJacobianDiffusion(0, 0));
        add_vector_form(new DefaultVectorFormVol(0, HERMES_ANY, 1.0, new CustomRightHandSide));
        add_vector_form_surf(new DefaultVectorFormSurf(0, BDY_RIGHT, -1.0));
    }
};
