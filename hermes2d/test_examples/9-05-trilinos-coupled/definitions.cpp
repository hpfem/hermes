using namespace Hermes;

// definition of reaction rate omega
void omega_fn(int n, Hermes::vector<double*> values, Hermes::vector<double*> dx, Hermes::vector<double*> dy,
                      double* out, double* outdx, double* outdy)
{
  for (int i = 0; i < n; i++)
  {
    double t1 = std::max(values.at(0)[i],0.0) - 1.0;
    double t2 = t1 * beta;
    double t3 = 1.0 + t1 * alpha;
    double t4 = sqr(beta) / (2.0*Le) * exp(t2 / t3);
    double t5 = (beta / (t3 * t3)) * values.at(1)[i];
    out[i] = t4 * values.at(1)[i];
    outdx[i] = t4 * (dx.at(1)[i] + dx.at(0)[i] * t5);
    outdy[i] = t4 * (dy.at(1)[i] + dy.at(0)[i] * t5);
  }
}

void omega_dt_fn(int n, Hermes::vector<double*> values, Hermes::vector<double*> dx, Hermes::vector<double*> dy,
                        double* out, double* outdx, double* outdy)
{
  for (int i = 0; i < n; i++)
  {
    double t1 = std::max(values.at(0)[i],0.0) - 1.0;
    double t2 = t1 * beta;
    double t3 = 1.0 + t1 * alpha;
    double t4 = sqr(beta) / (2.0*Le) * exp(t2 / t3);
    double t5 = (beta / (t3 * t3));
    out[i] = t4 * t5 * values.at(1)[i];
    outdx[i] = 0.0;
    outdy[i] = 0.0; // not important
  }
}

void omega_dc_fn(int n, Hermes::vector<double*> values, Hermes::vector<double*> dx, Hermes::vector<double*> dy,
                        double* out, double* outdx, double* outdy)
{
  for (int i = 0; i < n; i++)
  {
    double t1 = std::max(values.at(0)[i],0.0) - 1.0;
    double t2 = t1 * beta;
    double t3 = 1.0 + t1 * alpha;
    double t4 = sqr(beta) / (2.0*Le) * exp(t2 / t3);
    out[i] = t4;
    outdx[i] = 0.0;
    outdy[i] = 0.0; // not important
  }
}

// weak forms for the Newton's method

template<typename Real, typename double>
double newton_bilinear_form_0_0(int n, double *wt, Func<Real> *u_ext[], 
                                Func<Real> *vj, Func<Real> *vi, Geom<Real> *e, ExtData<double> *ext)
{
  double result = 0;
  Func<Real>* domegadt = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (  1.5 * vj->val[i] * vi->val[i] / TAU
                      +  vj->dx[i] * vi->dx[i] + vj->dy[i] * vi->dy[i]
                      - domegadt->val[i] * vj->val[i] * vi->val[i] );
  return result;
}

template<typename Real, typename double>
double newton_bilinear_form_0_0_surf(int n, double *wt, Func<Real> *u_ext[], 
                                     Func<Real> *vj, Func<Real> *vi, Geom<Real> *e, ExtData<double> *ext)
{
  double result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (kappa * vj->val[i] * vi->val[i]);
  return result;
}

template<typename Real, typename double>
double newton_bilinear_form_0_1(int n, double *wt, Func<Real> *u_ext[], 
                                Func<Real> *vj, Func<Real> *vi, Geom<Real> *e, ExtData<double> *ext)
{
  double result = 0;
  Func<Real>* domegady = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (- domegady->val[i] * vj->val[i] * vi->val[i] );
  return result;
}

template<typename Real, typename double>
double newton_bilinear_form_1_0(int n, double *wt, Func<Real> *u_ext[], 
                                Func<Real> *vj, Func<Real> *vi, Geom<Real> *e, ExtData<double> *ext)
{
  double result = 0;
  Func<Real>* domegadt = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * ( domegadt->val[i] * vj->val[i] * vi->val[i] );
  return result;
}

template<typename Real, typename double>
double newton_bilinear_form_1_1(int n, double *wt, Func<Real> *u_ext[], 
                                Func<Real> *vj, Func<Real> *vi, Geom<Real> *e, ExtData<double> *ext)
{
  double result = 0;
  Func<Real>* domegady = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (  1.5 * vj->val[i] * vi->val[i] / TAU
                      +  (vj->dx[i] * vi->dx[i] + vj->dy[i] * vi->dy[i]) / Le
                      + domegady->val[i] * vj->val[i] * vi->val[i] );
  return result;
}


template<typename Real, typename double>
double newton_linear_form_0(int n, double *wt, Func<Real> *u_ext[], 
                            Func<Real> *vi, Geom<Real> *e, ExtData<double> *ext)
{
  double result = 0;
  Func<Real>* titer = u_ext[0];
  Func<Real>* t_prev_time_1 = ext->fn[0];
  Func<Real>* t_prev_time_2 = ext->fn[1];
  Func<Real>* omega = ext->fn[2];
  for (int i = 0; i < n; i++)
    result += wt[i] * ( (3.0 * titer->val[i] - 4.0 * t_prev_time_1->val[i] 
                         + t_prev_time_2->val[i]) * vi->val[i] / (2.0 * TAU) +
                        (titer->dx[i] * vi->dx[i] + titer->dy[i] * vi->dy[i]) -
                        omega->val[i] * vi->val[i]);
  return result;
}

template<typename Real, typename double>
double newton_linear_form_0_surf(int n, double *wt, Func<Real> *u_ext[], 
                                 Func<Real> *vi, Geom<Real> *e, ExtData<double> *ext)
{
  double result = 0;
  Func<Real>* t_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (kappa * t_prev_newton->val[i] * vi->val[i]);
  return result;
}

template<typename Real, typename double>
double newton_linear_form_1(int n, double *wt, Func<Real> *u_ext[], 
                            Func<Real> *vi, Geom<Real> *e, ExtData<double> *ext)
{
  double result = 0;
  Func<Real>* c_prev_newton = u_ext[1];
  Func<Real>* c_prev_time_1 = ext->fn[0];
  Func<Real>* c_prev_time_2 = ext->fn[1];
  Func<Real>* omega = ext->fn[2];
  for (int i = 0; i < n; i++)
		result += wt[i] * ( (3.0 * c_prev_newton->val[i] - 4.0 * c_prev_time_1->val[i] + c_prev_time_2->val[i])
                         * vi->val[i] / (2.0 * TAU) +
                        (c_prev_newton->dx[i] * vi->dx[i] + c_prev_newton->dy[i] * vi->dy[i]) / Le +
                        omega->val[i] * vi->val[i]);
  return result;
}

// Preconditioner weak forms.
template<typename Real, typename double>
double precond_0_0(int n, double *wt, Func<double>* u_ext[], Func<Real> *vj, 
                   Func<Real> *vi, Geom<Real> *e, ExtData<double> *ext)
{
  double result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (  1.5 * vj->val[i] * vi->val[i] / TAU
                      +  vj->dx[i] * vi->dx[i] + vj->dy[i] * vi->dy[i]);
  return result;
}

template<typename Real, typename double>
double precond_1_1(int n, double *wt,  Func<double>* u_ext[], Func<Real> *vj, 
                   Func<Real> *vi, Geom<Real> *e, ExtData<double> *ext)
{
  double result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (  1.5 * vj->val[i] * vi->val[i] / TAU
                      +  (vj->dx[i] * vi->dx[i] + vj->dy[i] * vi->dy[i]) / Le );
  return result;
}
