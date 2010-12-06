// Heat conduction equation
template<typename Real, typename Scalar>
Scalar jac_TT(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *uj, Func<Real> *ui, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0.0;
  Func<Scalar>* T_prev_newton = u_ext[0];
  
  for (int i = 0; i < n; i++)
    result += wt[i] * (rho * cp * uj->val[i] * ui->val[i] / TAU
            + dk_dT(T_prev_newton->val[i]) * uj->val[i] * (T_prev_newton->dx[i] * ui->dx[i]
            + T_prev_newton->dy[i] * ui->dy[i])
            + k(T_prev_newton->val[i]) * (uj->dx[i] * ui->dx[i]
            + uj->dy[i] * ui->dy[i]));
  return result;
}

Ord jac_TT_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}


template<typename Real, typename Scalar>
Scalar jac_Tphi(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *uj, Func<Real> *ui, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0.0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (-kappa * xsfiss * uj->val[i] * ui->val[i]);
  return result;
}

Ord jac_Tphi_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}

template<typename Real, typename Scalar>
Scalar res_T(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *ui, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0.0;
  Func<Scalar>* T_prev_newton = u_ext[0];
  Func<Scalar>* phi_prev_newton = u_ext[1];
  Func<Scalar>* T_prev_time = ext->fn[0];
  
  for (int i = 0; i < n; i++)
    result += wt[i] * (rho * cp * (T_prev_newton->val[i] - T_prev_time->val[i]) / TAU * ui->val[i]
            + k(T_prev_newton->val[i]) * (T_prev_newton->dx[i] * ui->dx[i]
            + T_prev_newton->dy[i] * ui->dy[i])
            - kappa * xsfiss * phi_prev_newton->val[i] * ui->val[i]
            - qT(e->x[i], e->y[i]) * ui->val[i]);
  return result;
}

Ord res_T_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}

// Neutronics equation
template<typename Real, typename Scalar>
Scalar jac_phiphi(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *uj, Func<Real> *ui, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0.0;
  Func<Scalar>* T_prev_newton = u_ext[0];
  
  for (int i = 0; i < n; i++)
    result += wt[i] * (invvel * uj->val[i] * ui->val[i] / TAU
            + xsdiff * (uj->dx[i] * ui->dx[i]
            + uj->dy[i] * ui->dy[i])
            + xsrem(T_prev_newton->val[i]) * uj->val[i] * ui->val[i]
            - nu * xsfiss * uj->val[i] * ui->val[i]);
  return result;
}

Ord jac_phiphi_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}

template<typename Real, typename Scalar>
Scalar jac_phiT(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *uj, Func<Real> *ui, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0.0;
  Func<Scalar>* T_prev_newton = u_ext[0];
  Func<Scalar>* phi_prev_newton = u_ext[1];
  
  for (int i = 0; i < n; i++)
    result += wt[i] * (dxsrem_dT(T_prev_newton->val[i]) * uj->val[i] * phi_prev_newton->val[i] * ui->val[i]);
  return result;
}

Ord jac_phiT_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}

template<typename Real, typename Scalar>
Scalar res_phi(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *ui, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0.0;
  Func<Scalar>* T_prev_newton = u_ext[0];
  Func<Scalar>* phi_prev_newton = u_ext[1];
  Func<Scalar>* phi_prev_time = ext->fn[0];
  
  for (int i = 0; i < n; i++)
    result += wt[i] * (invvel * (phi_prev_newton->val[i] - phi_prev_time->val[i]) / TAU * ui->val[i]
            + xsdiff * (phi_prev_newton->dx[i] * ui->dx[i]
            + phi_prev_newton->dy[i] * ui->dy[i])
            + xsrem(T_prev_newton->val[i]) * phi_prev_newton->val[i] * ui->val[i] 
            - nu * xsfiss * phi_prev_newton->val[i] * ui->val[i]
            - q(e->x[i], e->y[i]) * ui->val[i]);
  return result;
}

Ord res_phi_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}

