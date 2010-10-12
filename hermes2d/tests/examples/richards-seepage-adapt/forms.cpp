// Jacobian matrix - volumetric part
double jac_form_vol_euler(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                          Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
  double x = e->x[0];
  double y = e->x[1];
  if (is_in_mat_1(x,y)) {
    K_S = K_S_1;
    ALPHA = ALPHA_1;
    THETA_R = THETA_R_1;
    THETA_S = THETA_R_1;
    N = N_1;
    M = M_1;
  }
  if (is_in_mat_2(x,y)) {
    K_S = K_S_2;
    ALPHA = ALPHA_2;
    THETA_R = THETA_R_2;
    THETA_S = THETA_R_2;
    N = N_2;
    M = M_2;
  }
  if (is_in_mat_3(x,y)) {
    K_S = K_S_3;
    ALPHA = ALPHA_3;
    THETA_R = THETA_R_3;
    THETA_S = THETA_R_3;
    N = N_3;
    M = M_3;
  }
  if (is_in_mat_4(x,y)) {
    K_S = K_S_4;
    ALPHA = ALPHA_4;
    THETA_R = THETA_R_4;
    THETA_S = THETA_R_4;
    N = N_4;
    M = M_4;
  }

  double result = 0;
  Func<double>* h_prev_newton = u_ext[0];
  Func<double>* h_prev_time = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (
		         C(h_prev_newton->val[i]) * u->val[i] * v->val[i] / TAU
		         + dCdh(h_prev_newton->val[i]) * u->val[i] * h_prev_newton->val[i] * v->val[i] / TAU
		         - dCdh(h_prev_newton->val[i]) * u->val[i] * h_prev_time->val[i] * v->val[i] / TAU
			 + K(h_prev_newton->val[i]) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i])
                         + dKdh(h_prev_newton->val[i]) * u->val[i] * 
                           (h_prev_newton->dx[i]*v->dx[i] + h_prev_newton->dy[i]*v->dy[i])
                         - dKdh(h_prev_newton->val[i]) * u->dy[i] * v->val[i]
                         - ddKdhh(h_prev_newton->val[i]) * u->val[i] * h_prev_newton->dy[i] * v->val[i]
                      );
  return result;
}

// Jacobian matrix - volumetric part
double jac_form_vol_cranic(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                           Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
  double x = e->x[0];
  double y = e->x[1];
  if (is_in_mat_1(x,y)) {
    K_S = K_S_1;
    ALPHA = ALPHA_1;
    THETA_R = THETA_R_1;
    THETA_S = THETA_R_1;
    N = N_1;
    M = M_1;
  }
  if (is_in_mat_2(x,y)) {
    K_S = K_S_2;
    ALPHA = ALPHA_2;
    THETA_R = THETA_R_2;
    THETA_S = THETA_R_2;
    N = N_2;
    M = M_2;
  }
  if (is_in_mat_3(x,y)) {
    K_S = K_S_3;
    ALPHA = ALPHA_3;
    THETA_R = THETA_R_3;
    THETA_S = THETA_R_3;
    N = N_3;
    M = M_3;
  }
  if (is_in_mat_4(x,y)) {
    K_S = K_S_4;
    ALPHA = ALPHA_4;
    THETA_R = THETA_R_4;
    THETA_S = THETA_R_4;
    N = N_4;
    M = M_4;
  }

  double result = 0;
  Func<double>* h_prev_newton = u_ext[0];
  Func<double>* h_prev_time = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * 0.5 * ( // implicit Euler part:
		         C(h_prev_newton->val[i]) * u->val[i] * v->val[i] / TAU
		         + dCdh(h_prev_newton->val[i]) * u->val[i] * h_prev_newton->val[i] * v->val[i] / TAU
		         - dCdh(h_prev_newton->val[i]) * u->val[i] * h_prev_time->val[i] * v->val[i] / TAU
			 + K(h_prev_newton->val[i]) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i])
                         + dKdh(h_prev_newton->val[i]) * u->val[i] * 
                           (h_prev_newton->dx[i]*v->dx[i] + h_prev_newton->dy[i]*v->dy[i])
                         - dKdh(h_prev_newton->val[i]) * u->dy[i] * v->val[i]
                         - ddKdhh(h_prev_newton->val[i]) * u->val[i] * h_prev_newton->dy[i] * v->val[i]
                       )
            + wt[i] * 0.5 * ( // explicit Euler part, 
		         C(h_prev_time->val[i]) * u->val[i] * v->val[i] / TAU
                       );
  return result;
}

// Integration order for Jacobian matrix - volumetric part
Ord jac_form_vol_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                     Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}

// Jacobian matrix - surface part on bdy 1
double jac_form_surf_1_euler(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                       Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
  double result = 0;
  Func<double>* h_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++) {
    result += wt[i] * dKdh(h_prev_newton->val[i]) * u->val[i] * v->val[i];
  }

  return result;
}

// Jacobian matrix - surface part on bdy 1
double jac_form_surf_1_cranic(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                       Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
  double result = 0;
  Func<double>* h_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++) {
    // Just the implicit Euler contributes:
    result += wt[i] * 0.5 * dKdh(h_prev_newton->val[i]) * u->val[i] * v->val[i];
  }

  return result;
}

// Integration order for Jacobian matrix - surface part on bdy 1
Ord jac_form_surf_1_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
                        Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}

// Jacobian matrix - surface part on bdy 4
double jac_form_surf_4_euler(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                       Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
  double result = 0;
  Func<double>* h_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++) {
    result -= wt[i] * dKdh(h_prev_newton->val[i]) * u->val[i] * v->val[i];
  }

  return result;
}

// Jacobian matrix - surface part on bdy 4
double jac_form_surf_4_cranic(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                       Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
  double result = 0;
  Func<double>* h_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++) {
    // Just the implicit Euler contributes:
    result -= wt[i] * 0.5 * dKdh(h_prev_newton->val[i]) * u->val[i] * v->val[i];
  }

  return result;
}

// Integration order for Jacobian matrix - surface part on bdy 4
Ord jac_form_surf_4_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
                        Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}

// Jacobian matrix - surface part on bdy 6
double jac_form_surf_6_euler(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                       Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
  double result = 0;
  Func<double>* h_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++) {
    result += wt[i] * dKdh(h_prev_newton->val[i]) * u->val[i] * v->val[i];
  }

  return result;
}

// Jacobian matrix - surface part on bdy 6
double jac_form_surf_6_cranic(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                       Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
  double result = 0;
  Func<double>* h_prev_newton = u_ext[0];
  // Just the implicit Euler contributes:
  for (int i = 0; i < n; i++) {
    result += wt[i] * 0.5 * dKdh(h_prev_newton->val[i]) * u->val[i] * v->val[i];
  }

  return result;
}

// Integration order for Jacobian matrix - surface part on bdy 6
Ord jac_form_surf_6_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
                        Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}

// Residual vector - volumetric part
double res_form_vol_euler(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                    Geom<double> *e, ExtData<double> *ext)
{
  double result = 0;
  Func<double>* h_prev_newton = u_ext[0];
  Func<double>* h_prev_time = ext->fn[0];
  for (int i = 0; i < n; i++) {
    result += wt[i] * (
		       C(h_prev_newton->val[i]) * (h_prev_newton->val[i] - h_prev_time->val[i]) * v->val[i] / TAU
                       + K(h_prev_newton->val[i]) * (h_prev_newton->dx[i] * v->dx[i] + h_prev_newton->dy[i] * v->dy[i])
                       - dKdh(h_prev_newton->val[i]) * h_prev_newton->dy[i] * v->val[i]
                      );
  }
  return result;
}

// Residual vector - volumetric part
double res_form_vol_cranic(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                    Geom<double> *e, ExtData<double> *ext)
{
  double result = 0;
  Func<double>* h_prev_newton = u_ext[0];
  Func<double>* h_prev_time = ext->fn[0];
  for (int i = 0; i < n; i++) {
    result += wt[i] * 0.5 * ( // implicit Euler part
		       C(h_prev_newton->val[i]) * (h_prev_newton->val[i] - h_prev_time->val[i]) * v->val[i] / TAU
                       + K(h_prev_newton->val[i]) * (h_prev_newton->dx[i] * v->dx[i] + h_prev_newton->dy[i] * v->dy[i])
                       - dKdh(h_prev_newton->val[i]) * h_prev_newton->dy[i] * v->val[i]
                      )
            + wt[i] * 0.5 * ( // explicit Euler part
		       C(h_prev_time->val[i]) * (h_prev_newton->val[i] - h_prev_time->val[i]) * v->val[i] / TAU
                       + K(h_prev_time->val[i]) * (h_prev_time->dx[i] * v->dx[i] + h_prev_time->dy[i] * v->dy[i])
                       - dKdh(h_prev_time->val[i]) * h_prev_time->dy[i] * v->val[i]
		       );
  }
  return result;
}

// Integration order for residual vector - volumetric part
Ord res_form_vol_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                     Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}

// Residual vector - surface part on bdy 1
double res_form_surf_1_euler(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                       Geom<double> *e, ExtData<double> *ext)
{
  double result = 0;
  Func<double>* h_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++) {
    result += wt[i] * K(h_prev_newton->val[i]) * v->val[i];
  }
  return result;
}

// Residual vector - surface part on bdy 1
double res_form_surf_1_cranic(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                       Geom<double> *e, ExtData<double> *ext)
{
  double result = 0;
  Func<double>* h_prev_newton = u_ext[0];
  Func<double>* h_prev_time = ext->fn[0];
  for (int i = 0; i < n; i++) {
    result += wt[i] * 0.5 * (K(h_prev_newton->val[i]) + K(h_prev_time->val[i])) * v->val[i];
  }
  return result;
}

// Integration order for residual vector - surface part on bdy 1
Ord res_form_surf_1_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                        Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}

// Residual vector - surface part on bdy 4
double res_form_surf_4_euler(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                       Geom<double> *e, ExtData<double> *ext)
{
  double result = 0;
  Func<double>* h_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++) {
    result -= wt[i] * K(h_prev_newton->val[i]) * v->val[i];
  }
  return result;
}

// Residual vector - surface part on bdy 4
double res_form_surf_4_cranic(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                       Geom<double> *e, ExtData<double> *ext)
{
  double result = 0;
  Func<double>* h_prev_newton = u_ext[0];
  Func<double>* h_prev_time = ext->fn[0];
  for (int i = 0; i < n; i++) {
    result -= wt[i] * 0.5 * (K(h_prev_newton->val[i]) + K(h_prev_time->val[i]))* v->val[i];
  }
  return result;
}

// Integration order for residual vector - surface part on bdy 4
Ord res_form_surf_4_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                        Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}

// Residual vector - surface part on bdy 6
double res_form_surf_6_euler(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                       Geom<double> *e, ExtData<double> *ext)
{
  double result = 0;
  Func<double>* h_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++) {
    result += wt[i] * (q_function() + K(h_prev_newton->val[i])) * v->val[i];
  }
  return result;
}

// Residual vector - surface part on bdy 6
double res_form_surf_6_cranic(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                       Geom<double> *e, ExtData<double> *ext)
{
  double result = 0;
  Func<double>* h_prev_newton = u_ext[0];
  Func<double>* h_prev_time = ext->fn[0];
  for (int i = 0; i < n; i++) {
    result += wt[i] * (q_function() + 0.5 * (K(h_prev_newton->val[i]) + K(h_prev_time->val[i]))) * v->val[i];
  }
  return result;
}

// Integration order for residual vector - surface part on bdy 6
Ord res_form_surf_6_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                        Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}

