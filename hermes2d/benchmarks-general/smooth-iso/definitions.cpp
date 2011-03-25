#include "weakform/weakform.h"
#include "weakform_library/laplace.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"


// Exact solution to the 2D equation.
class CustomExactSolution : public ExactSolutionScalar
{
public:
    CustomExactSolution(Mesh* mesh) : ExactSolutionScalar(mesh)
    {
    };

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



// Right-hand side for the 2D equation -Laplace u  = f
class CustomRightHandSide: public DefaultNonConstRightHandSide
{
public:
    CustomRightHandSide() : DefaultNonConstRightHandSide(){}

    virtual scalar value(double x, double y) const {
        return 2*sin(x)*sin(y);
    };

    virtual Ord ord(Ord x, Ord y) const {
        return Ord(7);
    }
};


// Weak forms for the 2D equation with Dirichlet boundary conditions.
class CustomWeakFormPoisson : public WeakForm
{
public:
    CustomWeakFormPoisson(DefaultNonConstRightHandSide* rhs) : WeakForm(1) {
        add_matrix_form(new DefaultMatrixFormStiffness(0, 0));
        add_vector_form(new DefaultVectorFormVolNonConst(0, rhs));
    }
};
