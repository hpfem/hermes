template<typename real>
real rhs(real x, real y, real z)
{
  real ddxx = 0;
  real ddyy = 0;
  real ddzz = 0;

  switch (ANISO_TYPE) {
    case ANISO_X: ddxx = -sin(x); break;
    case ANISO_Y: ddyy = -sin(y); break;
    case ANISO_Z: ddzz = -sin(z); break;

    case ANISO_X | ANISO_Y:
      ddxx = - sin(x) * sin(y);
      ddyy = - sin(x) * sin(y);
      break;

    case ANISO_X | ANISO_Z:
      ddxx = - sin(x) * sin(z);
      ddzz = - sin(x) * sin(z);
      break;

    case ANISO_Y | ANISO_Z:
      ddyy = - sin(y) * sin(z);
      ddzz = - sin(y) * sin(z);
      break;

    case ANISO_X | ANISO_Y | ANISO_Z:
      ddxx = - sin(x) * sin(y) * sin(z);
      ddyy = - sin(x) * sin(y) * sin(z);
      ddzz = - sin(x) * sin(y) * sin(z);
      break;
  }

  return -ddxx - ddyy - ddzz;
}

template<typename real, typename scalar>
scalar bilinear_form(int n, double *wt, Func<scalar> *u_ext[], Func<real> *u, Func<real> *v, Geom<real> *e, ExtData<scalar> *data)
{
  return int_grad_u_grad_v<real, scalar>(n, wt, u, v, e);
}

template<typename real, typename scalar>
scalar linear_form(int n, double *wt, Func<scalar> *u_ext[], Func<real> *u, Geom<real> *e, ExtData<scalar> *data)
{
  return int_F_v<real, scalar>(n, wt, rhs, u, e);
}
