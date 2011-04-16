/*
// Bilinear form (material 1)
template<typename Real, typename Scalar>
Scalar bilinear_form_1(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return D_1 * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v)
         + SIGMA_A_1 * int_u_v<Real, Scalar>(n, wt, u, v);
}

// Bilinear form (material 2)
template<typename Real, typename Scalar>
Scalar bilinear_form_2(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return D_2 * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v)
         + SIGMA_A_2 * int_u_v<Real, Scalar>(n, wt, u, v);
}

// Bilinear form (material 3)
template<typename Real, typename Scalar>
Scalar bilinear_form_3(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return D_3 * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v)
         + SIGMA_A_3 * int_u_v<Real, Scalar>(n, wt, u, v);
}

// Bilinear form (material 4)
template<typename Real, typename Scalar>
Scalar bilinear_form_4(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return D_4 * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v)
         + SIGMA_A_4 * int_u_v<Real, Scalar>(n, wt, u, v);
}

// Bilinear form (material 5)
template<typename Real, typename Scalar>
Scalar bilinear_form_5(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return D_5 * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v)
         + SIGMA_A_5 * int_u_v<Real, Scalar>(n, wt, u, v);
}

// Integration order for the bilinear forms
Ord bilinear_form_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return u->val[0] * v->val[0]; // returning the sum of the degrees of the basis
                                // and test function
}

// Linear form (material 1)
template<typename Real, typename Scalar>
Scalar linear_form_1(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return Q_EXT_1 * int_v<Real, Scalar>(n, wt, v);
}

// Linear form (material 3)
template<typename Real, typename Scalar>
Scalar linear_form_3(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return Q_EXT_3 * int_v<Real, Scalar>(n, wt, v);
}

// Integration order for the linear forms
Ord linear_form_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return v->val[0];  // q_ext is piecewise constant, thus
                     // returning the polynomial degree of the test function;
}
*/
#include "weakform/weakform.h"
#include "weakform_library/h1.h"

using namespace WeakFormsH1::VolumetricMatrixForms;
using namespace WeakFormsH1::VolumetricVectorForms;

class CustomWeakFormNeutronics : public WeakForm
{ 
public:
  CustomWeakFormNeutronics(Hermes::vector<std::string> layers, 
                                    Hermes::vector<double> D_map, 
                                    Hermes::vector<double> Sigma_a_map, 
                                    Hermes::vector<double> Sources_map)
  : WeakForm(1) {
    
    for (unsigned int i = 0; i < layers.size(); i++)
    {
      // Diffusion.
      add_matrix_form(new DefaultLinearDiffusion(0, 0, layers[i], D_map[i], HERMES_SYM));
    
      // Absorption.
      add_matrix_form(new DefaultLinearMass(0, 0, layers[i], Sigma_a_map[i], HERMES_SYM));
      
      // Sources.
      add_vector_form(new DefaultVectorFormConst(0, layers[i], Sources_map[i]));
    }
  }
};