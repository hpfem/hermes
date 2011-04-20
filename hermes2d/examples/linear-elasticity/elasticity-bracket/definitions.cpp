#include "weakform_library/elasticity.h"
#include "integrals/h1.h"
#include "boundaryconditions/essential_bcs.h"

class CustomWeakForm : public DefaultWeakFormLinearElasticity
{
public:
  CustomWeakForm(double E, double nu, double rho_g, std::string non_zero_neumann_bnd, double f0, double f1) 
            : DefaultWeakFormLinearElasticity(E, nu, rho_g) {
    double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));  // First Lame constant.
    double mu = E / (2*(1 + nu));                        // Second Lame constant.

    add_vector_form_surf(new VectorFormSurfForce_0(non_zero_neumann_bnd, f0));
    add_vector_form_surf(new VectorFormSurfForce_1(non_zero_neumann_bnd, f1));
  };

private:
  class VectorFormSurfForce_0 : public WeakForm::VectorFormSurf
  {
  public:
    VectorFormSurfForce_0(std::string marker, double f0) : WeakForm::VectorFormSurf(0, marker), f0(f0) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      return f0 * int_v<Real, Scalar>(n, wt, v);
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    // Member.
    double f0;
  };

  class VectorFormSurfForce_1 : public WeakForm::VectorFormSurf
  {
  public:
    VectorFormSurfForce_1(std::string marker, double f1) : WeakForm::VectorFormSurf(1, marker), f1(f1) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      return f1 * int_v<Real, Scalar>(n, wt, v);
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    // Member.
    double f1;
  };
};
