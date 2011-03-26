// Bilinear form, block 0, 0
template<typename Real, typename Scalar>
Scalar bilinear_form_0_0(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                         Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i=0; i < n; i++) {
    result += wt[i] * (u->val0[i] * v->val0[i] + u->val1[i] * v->val1[i]);
  }

  return result / time_step / (SPEED_OF_LIGHT * SPEED_OF_LIGHT);
}

// Bilinear form, block 0, 1
template<typename Real, typename Scalar>
Scalar bilinear_form_0_1(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                         Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i=0; i < n; i++) {
    result += wt[i] * (-u->dy[i] * v->val0[i] + u->dx[i] * v->val1[i]);
  }

  return result;
}

// Bilinear form, block 1, 0
template<typename Real, typename Scalar>
Scalar bilinear_form_1_0(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                         Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i=0; i < n; i++) {
    result += wt[i] * u->curl[i] * v->val[i];
  }

  return result;
}

// Bilinear form, block 1, 1
template<typename Real, typename Scalar>
Scalar bilinear_form_1_1(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                         Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i=0; i < n; i++) {
    result += wt[i] * u->val[i] * v->val[i];
  }

  return result / time_step;
}

// Linear form, block 0
template<typename Real, typename Scalar>
Scalar linear_form_0(int n, double *wt, Func<Scalar> *u_ext[], 
                     Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Func<Scalar>* u_prev = ext->fn[0];  // E on previous time level.

  Scalar result = 0;
  for (int i=0; i < n; i++) {
    result += wt[i] * (u_prev->val0[i] * v->val0[i] + u_prev->val1[i] * v->val1[i]);
  }

  return result / time_step / (SPEED_OF_LIGHT * SPEED_OF_LIGHT);
}

// Linear form, block 1
template<typename Real, typename Scalar>
Scalar linear_form_1(int n, double *wt, Func<Scalar> *u_ext[], 
                     Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Func<Scalar>* u_prev = ext->fn[0];  // B on previous time level.

  Scalar result = 0;
  for (int i=0; i < n; i++) {
    result += wt[i] * u_prev->val[i] * v->val[i];
  }

  return result / time_step;
}





