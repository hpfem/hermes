/* Exact solution */

double fn(double x, double y, double z)
{
  switch (ANISO_TYPE) {
    case ANISO_X: return sin(x);
    case ANISO_Y: return sin(y);
    case ANISO_Z: return sin(z);

    case ANISO_X | ANISO_Y: return sin(x) * sin(y);
    case ANISO_X | ANISO_Z: return sin(x) * sin(z);
    case ANISO_Y | ANISO_Z: return sin(y) * sin(z);

    case ANISO_X | ANISO_Y | ANISO_Z: return sin(x) * sin(y) * sin(z);

    default:
    return 0;
  }
}

// Needed for calculation of norms and used by visualizer.
double fndd(double x, double y, double z, 
                      double &dx, double &dy, double &dz)
{
  switch (ANISO_TYPE) {
    case ANISO_X:
      dx = cos(x);
      dy = 0;
      dz = 0;
      break;

    case ANISO_Y:
      dx = 0;
      dy = cos(y);
      dz = 0;
      break;

    case ANISO_Z:
      dx = 0;
      dy = 0;
      dz = cos(z);
      break;

    case ANISO_X | ANISO_Y:
      dx = cos(x) * sin(y);
      dy = sin(x) * cos(y);
      dz = 0;
      break;

    case ANISO_X | ANISO_Z:
      dx = cos(x) * sin(z);
      dy = 0;
      dz = sin(x) * cos(z);
      break;

    case ANISO_Y | ANISO_Z:
      dx = 0;
      dy = cos(y) * sin(z);
      dz = sin(y) * cos(z);
      break;

    case ANISO_X | ANISO_Y | ANISO_Z:
      dx = cos(x) * sin(y) * sin(z);
      dy = sin(x) * cos(y) * sin(z);
      dz = sin(x) * sin(y) * cos(z);
      break;
    }

  return fn(x, y, z);
}

/* Weak forms */

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

