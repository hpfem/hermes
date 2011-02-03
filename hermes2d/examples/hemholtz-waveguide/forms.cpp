#include <integrals/integrals_h1.h>

template<typename Real, typename Scalar>
Scalar matrix_form_real_real(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v) - sqr(omega) * mur * mu0 * epsr * eps0 * int_u_v<Real, Scalar>(n, wt, u, v);    
}

template<typename Real, typename Scalar>
Scalar matrix_form_real_imag(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return -omega * mur * mu0 * sigma * int_u_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar matrix_form_imag_real(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return  omega * mur * mu0 * sigma * int_u_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar matrix_form_imag_imag(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v) - sqr(omega) * mur * mu0 * epsr * eps0 * int_u_v<Real, Scalar>(n, wt, u, v);   
}

template<typename Real, typename Scalar>
Scalar matrix_form_surface_imag_real(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return - beta*int_u_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar matrix_form_surface_real_imag(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return beta * int_u_v<Real, Scalar>(n, wt, u, v);
}
