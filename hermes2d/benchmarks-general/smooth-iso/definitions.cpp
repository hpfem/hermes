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

    virtual void derivatives(double x, double y, scalar& dx, scalar& dy) {
        dx = cos(x)*sin(y);
        dy = sin(x)*cos(y);
    };

    virtual double value(double x, double y) {
        return sin(x)*sin(y);
    }

    virtual Ord ord(Ord x, Ord y) {
        return Ord(15);
    }
};


// Right-hand side for the 2D equation -Laplace u  = f
class CustomRightHandSide: public WeakForm::VectorFormVol
{
public:
    CustomRightHandSide(unsigned i) : WeakForm::VectorFormVol(i) {}

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                         Geom<double> *e, ExtData<scalar> *ext) {

        scalar result = 0;
        for (int i = 0; i < n; i++)
            result += wt[i] * (2*sin(e->x[i])*sin(e->y[i])*v->val[i]);
        return result;
    }


    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e,
                    ExtData<Ord> *ext) {
        return Ord(15);
    }
};


// Weak forms for the 2D equation with Dirichlet boundary conditions.
class CustomWeakFormPoisson : public WeakForm
{
public:
    CustomWeakFormPoisson() : WeakForm(1) {
        add_matrix_form(new DefaultMatrixFormVolConst(0, 0));
        add_vector_form(new CustomRightHandSide(0));
    }
};
