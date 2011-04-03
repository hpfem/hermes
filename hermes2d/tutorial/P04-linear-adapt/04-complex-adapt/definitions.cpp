#include "weakform/weakform.h"
#include "weakform_library/h1.h"

using namespace WeakFormsH1::VolumetricMatrixForms;
using namespace WeakFormsH1::VolumetricVectorForms;
  
/* Weak forms */
class CustomWeakFormMagnetics : public WeakForm
{ 
public:
  CustomWeakFormMagnetics(std::string mat_air, double mu_air,
                          std::string mat_iron, double mu_iron, double gamma_iron,
                          std::string mat_wire, double mu_wire, double j_ext, double omega)
  : WeakForm(1) {
    scalar ii =  cplx(0.0, 1.0);

    add_matrix_form(new DefaultLinearDiffusion(0, 0, mat_air, ii * 1.0/mu_air));
    add_matrix_form(new DefaultLinearDiffusion(0, 0, mat_iron, 1.0/mu_iron));
    add_matrix_form(new DefaultLinearDiffusion(0, 0, mat_wire, 1.0/mu_wire));
    add_matrix_form(new DefaultLinearMass(0, 0, mat_iron, omega * gamma_iron));

    add_vector_form(new DefaulVectorFormConst(0, mat_wire, -j_ext));
  };
};

/*
template<typename Real, typename Scalar>
Scalar bilinear_form_iron(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  scalar ii = cplx(0.0, 1.0);
  return 1./mu_iron * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v) + ii*omega*gamma_iron*int_u_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_wire(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return 1./mu_0 * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_air(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return 1./mu_0 * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v); // conductivity gamma is zero
}

template<typename Real, typename Scalar>
Scalar linear_form_wire(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return J_wire * int_v<Real, Scalar>(n, wt, v);
}
*/
