//////  Bilinear and linear forms - axisymmetric arrangement  ////////////////////////////////////////////////

// NOTE: The global variable 'k_eff' from main.cpp is used in the linear forms.

template<typename Real, typename Scalar>
Scalar int_x_u_v(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (e->x[i] * u->val[i] * v->val[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar int_x_grad_u_grad_v(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * e->x[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar projection_biform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_x_u_v<Real, Scalar>(n, wt, u, v, e) + int_x_grad_u_grad_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar projection_liform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_x_u_v<Real, Scalar>(n, wt, ext->fn[0], v, e) + int_x_grad_u_grad_v<Real, Scalar>(n, wt, ext->fn[0], v, e);
}

//////////   Eq 1   /////////////////////////////////////////////////////////////////////////////////////////

template<typename Real, typename Scalar>
Scalar biform_0_0(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (D[e->elem_marker - 1][0]) * int_x_grad_u_grad_v<Real, Scalar>(n, wt, u, v, e) +
         (Sr[e->elem_marker - 1][0]) * int_x_u_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar biform_surf_0_0(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (0.5) * int_x_u_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar liform_0(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (chi[e->elem_marker - 1][0] / k_eff) * (nu[e->elem_marker - 1][0] * Sf[e->elem_marker - 1][0] * ext->fn[0]->val[i] +
                                                   nu[e->elem_marker - 1][1] * Sf[e->elem_marker - 1][1] * ext->fn[1]->val[i] +
                                                   nu[e->elem_marker - 1][2] * Sf[e->elem_marker - 1][2] * ext->fn[2]->val[i] +
                                                   nu[e->elem_marker - 1][3] * Sf[e->elem_marker - 1][3] * ext->fn[3]->val[i])
                                      * e->x[i] * v->val[i];
  return result;
}

//////////   Eq 2   /////////////////////////////////////////////////////////////////////////////////////////

template<typename Real, typename Scalar>
Scalar biform_1_1(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (D[e->elem_marker - 1][1]) * int_x_grad_u_grad_v<Real, Scalar>(n, wt, u, v, e) +
         (Sr[e->elem_marker - 1][1]) * int_x_u_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar biform_surf_1_1(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (0.5) * int_x_u_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar biform_1_0(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (- Ss[e->elem_marker - 1][1][0]) * int_x_u_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar liform_1(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (chi[e->elem_marker - 1][1] / k_eff) * (nu[e->elem_marker - 1][0] * Sf[e->elem_marker - 1][0] * ext->fn[0]->val[i] +
                                                   nu[e->elem_marker - 1][1] * Sf[e->elem_marker - 1][1] * ext->fn[1]->val[i] +
                                                   nu[e->elem_marker - 1][2] * Sf[e->elem_marker - 1][2] * ext->fn[2]->val[i] +
                                                   nu[e->elem_marker - 1][3] * Sf[e->elem_marker - 1][3] * ext->fn[3]->val[i])
                                      * e->x[i] * v->val[i];
  return result;
}

//////////   Eq 3   /////////////////////////////////////////////////////////////////////////////////////////

template<typename Real, typename Scalar>
Scalar biform_2_2(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (D[e->elem_marker - 1][2]) * int_x_grad_u_grad_v<Real, Scalar>(n, wt, u, v, e) +
         (Sr[e->elem_marker - 1][2]) * int_x_u_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar biform_surf_2_2(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (0.5) * int_x_u_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar biform_2_1(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (- Ss[e->elem_marker - 1][2][1]) * int_x_u_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar liform_2(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (chi[e->elem_marker - 1][2] / k_eff) * (nu[e->elem_marker - 1][0] * Sf[e->elem_marker - 1][0] * ext->fn[0]->val[i] +
                                                   nu[e->elem_marker - 1][1] * Sf[e->elem_marker - 1][1] * ext->fn[1]->val[i] +
                                                   nu[e->elem_marker - 1][2] * Sf[e->elem_marker - 1][2] * ext->fn[2]->val[i] +
                                                   nu[e->elem_marker - 1][3] * Sf[e->elem_marker - 1][3] * ext->fn[3]->val[i])
                                      * e->x[i] * v->val[i];
  return result;
}

//////////   Eq 4   /////////////////////////////////////////////////////////////////////////////////////////

template<typename Real, typename Scalar>
Scalar biform_3_3(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (D[e->elem_marker - 1][3]) * int_x_grad_u_grad_v<Real, Scalar>(n, wt, u, v, e) +
         (Sr[e->elem_marker - 1][3]) * int_x_u_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar biform_surf_3_3(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (0.5) * int_x_u_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar biform_3_2(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return (- Ss[e->elem_marker - 1][3][2]) * int_x_u_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar liform_3(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (chi[e->elem_marker - 1][3] / k_eff) * (nu[e->elem_marker - 1][0] * Sf[e->elem_marker - 1][0] * ext->fn[0]->val[i] +
                                                   nu[e->elem_marker - 1][1] * Sf[e->elem_marker - 1][1] * ext->fn[1]->val[i] +
                                                   nu[e->elem_marker - 1][2] * Sf[e->elem_marker - 1][2] * ext->fn[2]->val[i] +
                                                   nu[e->elem_marker - 1][3] * Sf[e->elem_marker - 1][3] * ext->fn[3]->val[i])
                                      * e->x[i] * v->val[i];
  return result;
}

//////  Determining the quadrature order used for integrating the respective forms.  ////////////////////////////////////////////////

#define DIAG_BIFORM_VOL_ORD(i)\
    template<>\
    Ord biform_##i##_##i(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)\
    {\
      return int_grad_u_grad_v<Ord, Ord>(n, wt, u, v) * int_u_v<Ord, Ord>(n, wt, u, v);\
    }
#define DIAG_BIFORM_SURF_ORD(i)\
    template<>\
    Ord biform_surf_##i##_##i(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)\
    {\
      return int_u_v<Ord, Ord>(n, wt, u, v);\
    }
#define OFFDIAG_BIFORM_VOL_ORD(i,j)\
    template<>\
    Ord biform_##i##_##j(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)\
    {\
      return int_u_v<Ord, Ord>(n, wt, u, v);\
    }
#define LIFORM_VOL_ORD(i)\
    template<>\
    Ord liform_##i(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)\
    {\
      return (ext->fn[0]->val[0] + ext->fn[1]->val[0] + ext->fn[2]->val[0] + ext->fn[3]->val[0])*e->x[0]*v->val[0];\
    }
    
DIAG_BIFORM_VOL_ORD(0)
DIAG_BIFORM_VOL_ORD(1)
DIAG_BIFORM_VOL_ORD(2)
DIAG_BIFORM_VOL_ORD(3)

DIAG_BIFORM_SURF_ORD(0)
DIAG_BIFORM_SURF_ORD(1)
DIAG_BIFORM_SURF_ORD(2)
DIAG_BIFORM_SURF_ORD(3)

OFFDIAG_BIFORM_VOL_ORD(1,0)
OFFDIAG_BIFORM_VOL_ORD(2,1)
OFFDIAG_BIFORM_VOL_ORD(3,2)

LIFORM_VOL_ORD(0)
LIFORM_VOL_ORD(1)
LIFORM_VOL_ORD(2)
LIFORM_VOL_ORD(3)
