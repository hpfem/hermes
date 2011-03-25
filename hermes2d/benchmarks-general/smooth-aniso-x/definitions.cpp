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



// Right-hand side for the 2D equation -Laplace u  = f
class CustomRightHandSide: public DefaultNonConstRightHandSide
{
public:
    CustomRightHandSide() : DefaultNonConstRightHandSide(){}

    virtual scalar value(double x, double y) const {
        return sin(x);
    };

    virtual Ord ord(Ord x, Ord y) const {
        return Ord(7);
    }
};

class VectorFormSurfAnisoX : public WeakForm::VectorFormSurf
{
public:
    VectorFormSurfAnisoX(int i, std::string area) : VectorFormSurf(i, area){}

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<scalar> *v, Geom<scalar> *e, ExtData<scalar> *ext)
    {
      return -int_v<scalar, scalar>(n, wt, v);
    }


    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
    {
      return Ord(7);  // returning the polynomial degree of the test function plus two
    }
};


// Weak forms for the 2D equation with Dirichlet boundary conditions.
class CustomWeakFormPoisson : public WeakForm
{
public:
    CustomWeakFormPoisson(DefaultNonConstRightHandSide* rhs) : WeakForm(1) {
        add_matrix_form(new DefaultMatrixFormStiffness(0, 0));
        add_vector_form(new DefaultVectorFormVolNonConst(0, rhs));
        add_vector_form_surf(new VectorFormSurfAnisoX(0, BDY_RIGHT));
    }
};
