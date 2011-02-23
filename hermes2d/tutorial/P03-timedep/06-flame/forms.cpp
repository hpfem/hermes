// definition of reaction rate omega

void omega_fn(int n, Hermes::vector<scalar*> values, Hermes::vector<scalar*> dx, Hermes::vector<scalar*> dy,
                      scalar* out, scalar* outdx, scalar* outdy)
{
  for (int i = 0; i < n; i++)
  {
    scalar t1 = values.at(0)[i] - 1.0;
    scalar t2 = t1 * beta;
    scalar t3 = 1.0 + t1 * alpha;
    scalar t4 = sqr(beta) / (2.0*Le) * exp(t2 / t3);
    scalar t5 = (beta / (t3 * t3)) * values.at(1)[i];
    out[i] = t4 * values.at(1)[i];
    outdx[i] = t4 * (dx.at(1)[i] + dx.at(0)[i] * t5);
    outdy[i] = t4 * (dy.at(1)[i] + dy.at(0)[i] * t5);
  }
}

void omega_dt_fn(int n, Hermes::vector<scalar*> values, Hermes::vector<scalar*> dx, Hermes::vector<scalar*> dy,
                        scalar* out, scalar* outdx, scalar* outdy)
{
  for (int i = 0; i < n; i++)
  {
    scalar t1 = values.at(0)[i] - 1.0;
    scalar t2 = t1 * beta;
    scalar t3 = 1.0 + t1 * alpha;
    scalar t4 = sqr(beta) / (2.0*Le) * exp(t2 / t3);
    scalar t5 = (beta / (t3 * t3));
    out[i] = t4 * t5 * values.at(1)[i];
    outdx[i] = 0.0;
    outdy[i] = 0.0; // not important
  }
}

void omega_dc_fn(int n, Hermes::vector<scalar*> values, Hermes::vector<scalar*> dx, Hermes::vector<scalar*> dy,
                        scalar* out, scalar* outdx, scalar* outdy)
{
  for (int i = 0; i < n; i++)
  {
    scalar t1 = values.at(0)[i] - 1.0;
    scalar t2 = t1 * beta;
    scalar t3 = 1.0 + t1 * alpha;
    scalar t4 = sqr(beta) / (2.0*Le) * exp(t2 / t3);
    out[i] = t4;
    outdx[i] = 0.0;
    outdy[i] = 0.0; // not important
  }
}

// weak forms for the Newton's method

template<typename Real, typename Scalar>
Scalar newton_bilinear_form_0_0(int n, double *wt, Func<Real> *u_ext[], 
                                Func<Real> *vj, Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Real>* domegadt = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (  1.5 * vj->val[i] * vi->val[i] / time_step
                      +  vj->dx[i] * vi->dx[i] + vj->dy[i] * vi->dy[i]
                      - domegadt->val[i] * vj->val[i] * vi->val[i] );
  return result;
}

template<typename Real, typename Scalar>
Scalar newton_bilinear_form_0_0_surf(int n, double *wt, Func<Real> *u_ext[], 
                                     Func<Real> *vj, Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (kappa * vj->val[i] * vi->val[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar newton_bilinear_form_0_1(int n, double *wt, Func<Real> *u_ext[], 
                                Func<Real> *vj, Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Real>* domegady = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (- domegady->val[i] * vj->val[i] * vi->val[i] );
  return result;
}

template<typename Real, typename Scalar>
Scalar newton_bilinear_form_1_0(int n, double *wt, Func<Real> *u_ext[], 
                                Func<Real> *vj, Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Real>* domegadt = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * ( domegadt->val[i] * vj->val[i] * vi->val[i] );
  return result;
}

template<typename Real, typename Scalar>
Scalar newton_bilinear_form_1_1(int n, double *wt, Func<Real> *u_ext[], 
                                Func<Real> *vj, Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Real>* domegady = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (  1.5 * vj->val[i] * vi->val[i] / time_step
                      +  (vj->dx[i] * vi->dx[i] + vj->dy[i] * vi->dy[i]) / Le
                      + domegady->val[i] * vj->val[i] * vi->val[i] );
  return result;
}


template<typename Real, typename Scalar>
Scalar newton_linear_form_0(int n, double *wt, Func<Real> *u_ext[], 
                            Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Real>* titer = u_ext[0];
  Func<Real>* t_prev_time_1 = ext->fn[0];
  Func<Real>* t_prev_time_2 = ext->fn[1];
  Func<Real>* omega = ext->fn[2];
  for (int i = 0; i < n; i++)
    result += wt[i] * ( (3.0 * titer->val[i] - 4.0 * t_prev_time_1->val[i] 
                         + t_prev_time_2->val[i]) * vi->val[i] / (2.0 * time_step) +
                        (titer->dx[i] * vi->dx[i] + titer->dy[i] * vi->dy[i]) -
                        omega->val[i] * vi->val[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar newton_linear_form_0_surf(int n, double *wt, Func<Real> *u_ext[], 
                                 Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Real>* t_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (kappa * t_prev_newton->val[i] * vi->val[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar newton_linear_form_1(int n, double *wt, Func<Real> *u_ext[], 
                            Func<Real> *vi, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Real>* c_prev_newton = u_ext[1];
  Func<Real>* c_prev_time_1 = ext->fn[0];
  Func<Real>* c_prev_time_2 = ext->fn[1];
  Func<Real>* omega = ext->fn[2];
  for (int i = 0; i < n; i++)
		result += wt[i] * ( (3.0 * c_prev_newton->val[i] - 4.0 * c_prev_time_1->val[i] + c_prev_time_2->val[i])
                         * vi->val[i] / (2.0 * time_step) +
                        (c_prev_newton->dx[i] * vi->dx[i] + c_prev_newton->dy[i] * vi->dy[i]) / Le +
                        omega->val[i] * vi->val[i]);
  return result;
}
