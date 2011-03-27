// Bessel function of the first kind, order n, defined in bessel.cpp
double jv(double n, double x);

static void exact_sol_val(double x, double y, double z, scalar &e0, scalar &e1)
{
  double t1 = x*x;
  double t2 = y*y;
  double t4 = sqrt(t1+t2);
  double t5 = jv(-1.0/3.0,t4);
  double t6 = 1/t4;
  double t7 = jv(2.0/3.0,t4);
  double t11 = (t5-2.0/3.0*t6*t7)*t6;
  double t12 = atan2(y,x);
  if (t12 < 0) t12 += 2.0*M_PI;
  double t13 = 2.0/3.0*t12;
  double t14 = cos(t13);
  double t17 = sin(t13);
  double t18 = t7*t17;
  double t20 = 1/t1;
  double t23 = 1/(1.0+t2*t20);
  e0 = t11*y*t14-2.0/3.0*t18/x*t23;
  e1 = -t11*x*t14-2.0/3.0*t18*y*t20*t23;
}

static void exact_sol(double x, double y, double z, scalar &e0, scalar &e1, scalar &e1dx, scalar &e0dy)
{
  exact_sol_val(x, y, z, e0, e1);

  double t1 = x*x;
  double t2 = y*y;
  double t3 = t1+t2;
  double t4 = sqrt(t3);
  double t5 = jv(2.0/3.0,t4);
  double t6 = 1/t4;
  double t7 = jv(-1.0/3.0,t4);
  double t11 = (-t5-t6*t7/3.0)*t6;
  double t14 = 1/t4/t3;
  double t15 = t14*t5;
  double t21 = t7-2.0/3.0*t6*t5;
  double t22 = 1/t3*t21;
  double t27 = atan2(y,x);
  if (t27 < 0) t27 += 2.0*M_PI;
  double t28 = 2.0/3.0*t27;
  double t29 = cos(t28);
  double t32 = t21*t14;
  double t35 = t21*t6;
  double t36 = t35*t29;
  double t39 = sin(t28);
  double t41 = 1/t1;
  double t43 = 1.0+t2*t41;
  double t44 = 1/t43;
  double t47 = 4.0/3.0*t35/x*t39*y*t44;
  double t48 = t5*t29;
  double t49 = t1*t1;
  double t52 = t43*t43;
  double t53 = 1/t52;
  double t57 = t5*t39;
  double t59 = 1/t1/x;
  e1dx =-(t11*x+2.0/3.0*t15*x-2.0/3.0*t22*x)
         *t6*x*t29+t32*t1*t29-t36-t47+4.0/9.0*t48*t2/t49*t53+4.0/3.0*t57*y*t59*t44-4.0/3.0*t57*t2*y/t49/x*t53;
  e0dy = (t11*y+2.0/3.0*t15*y-2.0/3.0*t22*y)*t6*y*t29-t32*t2*t29+t36-t47-4.0/9.0*t48*t41*t53+4.0/3.0*t57*t59*t53*y;
}

// Exact solution.
scalar3 exact(double x, double y, double z, scalar3 &dx, scalar3 &dy, scalar3 &dz)
{
  scalar3 ex(0.0, 0.0, 0.0);

  exact_sol(x, y, z, ex[0], ex[1], dx[1], dy[0]);
  return ex;
}

// Weak forms.
template<typename real, typename scalar>
scalar biform(int n, double *wt, Func<scalar> *u_ext[], Func<real> *u, Func<real> *v, Geom<real> *e, ExtData<scalar> *ext)
{
  return 1.0/mu_r * hcurl_int_curl_u_curl_v<real, scalar>(n, wt, u, v, e)
    - sqr(kappa) * hcurl_int_u_v<real, scalar>(n, wt, u, v, e);
}

Ord biform_surf_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(v->val[0].get_max_order());
}

scalar biform_surf(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext)
{
  // j * kappa * E_T * F_T
  // E_T = nu x E x nu  (nu is outer normal)
  cplx ii = cplx(0.0, 1.0);
  scalar result = 0;
  for (int i = 0; i < n; i++) {
    scalar uu[3] = { u->val0[i], u->val1[i], u->val2[i] };
    scalar tpu[3];
    calc_tan_proj(e->nx[i], e->ny[i], e->nz[i], uu, tpu);

    scalar vv[3] = { v->val0[i], v->val1[i], v->val2[i] };
    scalar tpv[3];
    calc_tan_proj(e->nx[i], e->ny[i], e->nz[i], vv, tpv);

    result += wt[i] * (uu[0] * vv[0] + uu[1] * vv[1] + uu[2] * vv[2]);
  }

  return ii * (-kappa) * result;
}

scalar liform_surf(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext)
{
  cplx ii = cplx(0.0, 1.0);
  scalar result = 0;
  for (int i = 0; i < n; i++) {
    scalar3 dx0(0, 0, 0), dy0(0, 0, 0), dz0(0, 0, 0);
    scalar3 ev0(0.0, 0.0, 0.0);
    ev0 = exact(e->x[i], e->y[i], e->z[i], dx0, dy0, dz0);

    scalar curl_e[3];
    scalar dx[3], dy[3], dz[3];
    dx[0] = dx0[0]; dx[1] = dx0[1]; dx[2] = dx0[2]; 
    dy[0] = dy0[0]; dy[1] = dy0[1]; dy[2] = dy0[2]; 
    dz[0] = dz0[0]; dz[1] = dz0[1]; dz[2] = dz0[2]; 
    calc_curl(dx, dy, dz, curl_e);
    scalar tpe[3];
    scalar ev[3];
    ev[0] = ev0[0]; ev[1] = ev0[1]; ev[2] = ev0[2];
    calc_tan_proj(e->nx[i], e->ny[i], e->nz[i], ev, tpe);

    scalar g[3] = {
      (e->nz[i] * curl_e[1] - e->ny[i] * curl_e[2]) - ii * kappa * tpe[0],
      (e->nx[i] * curl_e[2] - e->nz[i] * curl_e[0]) - ii * kappa * tpe[1],
      (e->ny[i] * curl_e[0] - e->nx[i] * curl_e[1]) - ii * kappa * tpe[2],
    };

    // tpv is tangencial projection of v (test function)
    scalar vv[3] = { v->val0[i], v->val1[i], v->val2[i] };
    scalar tpv[3];
    calc_tan_proj(e->nx[i], e->ny[i], e->nz[i], vv, tpv);

    result += wt[i] * (g[0] * tpv[0] + g[1] * tpv[1] + g[2] * tpv[2]);
  }

  return result;
}

// Maximal polynomial order to integrate surface linear form.
Ord liform_surf_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(v->val[0].get_max_order());
}
