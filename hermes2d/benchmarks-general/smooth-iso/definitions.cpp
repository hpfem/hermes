#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform_library/h1.h"

using namespace WeakFormsH1::VolumetricMatrixForms;
using namespace WeakFormsH1::VolumetricVectorForms;
using namespace WeakFormsH1::RightHandSides;

/*  Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
    CustomExactSolution(Mesh* mesh) : ExactSolutionScalar(mesh) { }

    virtual void derivatives(double x, double y, scalar& dx, scalar& dy) const {
        dx = cos(x)*sin(y);
        dy = sin(x)*cos(y);
    };

    virtual double value(double x, double y) const {
        return sin(x)*sin(y);
    }

    virtual Ord ord(Ord x, Ord y) const {
        return Ord(7);
    }
};

/* Right-hand side */

class CustomRightHandSide: public DefaultNonConstRightHandSide
{
public:
    CustomRightHandSide() : DefaultNonConstRightHandSide() { }

    virtual scalar value(double x, double y) const {
        return 2*sin(x)*sin(y);
    };

    virtual Ord ord(Ord x, Ord y) const {
        return Ord(7);
    }
};

/* Weak forms */

class CustomWeakFormPoisson : public WeakForm
{
public:
    CustomWeakFormPoisson(DefaultNonConstRightHandSide* rhs) : WeakForm(1) {
        add_matrix_form(new DefaultLinearDiffusion(0, 0));
        add_vector_form(new DefaultVectorFormNonConst(0, rhs));
    }
};
