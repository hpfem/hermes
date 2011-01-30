// The first part of the file dontains forms for the Newton's
// method. Forms for Picard are in the second part.

/*** NEWTON ***/

// Jacobian matrix - volumetric part

double jac_form_vol_euler(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                          Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
  int elem_marker = e->elem_marker;
  double result = 0;
  Func<double>* h_prev_newton = u_ext[0];
  Func<double>* h_prev_time = ext->fn[0];

  for (int i = 0; i < n; i++)
    result += wt[i] * (
		         C(h_prev_newton->val[i], elem_marker) * u->val[i] * v->val[i] / time_step
		         + dCdh(h_prev_newton->val[i], elem_marker) * u->val[i] * h_prev_newton->val[i] * v->val[i] / time_step
		         - dCdh(h_prev_newton->val[i], elem_marker) * u->val[i] * h_prev_time->val[i] * v->val[i] / time_step
			 + K(h_prev_newton->val[i], elem_marker) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i])
                         + dKdh(h_prev_newton->val[i], elem_marker) * u->val[i] * 
                           (h_prev_newton->dx[i]*v->dx[i] + h_prev_newton->dy[i]*v->dy[i])
                         - dKdh(h_prev_newton->val[i], elem_marker) * u->dy[i] * v->val[i]
                         - ddKdhh(h_prev_newton->val[i], elem_marker) * u->val[i] * h_prev_newton->dy[i] * v->val[i]
                      );
  return result;
}

// Jacobian matrix - volumetric part
double jac_form_vol_cranic(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                           Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
  double x = e->x[0];
  double y = e->x[1];


  int elem_marker = e->elem_marker;
  double result = 0;
  Func<double>* h_prev_newton = u_ext[0];
  Func<double>* h_prev_time = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * 0.5 * ( // implicit Euler part:
		         C(h_prev_newton->val[i], elem_marker) * u->val[i] * v->val[i] / time_step
		         + dCdh(h_prev_newton->val[i], elem_marker) * u->val[i] * h_prev_newton->val[i] * v->val[i] / time_step
		         - dCdh(h_prev_newton->val[i], elem_marker) * u->val[i] * h_prev_time->val[i] * v->val[i] / time_step
			 + K(h_prev_newton->val[i], elem_marker) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i])
                         + dKdh(h_prev_newton->val[i], elem_marker) * u->val[i] * 
                           (h_prev_newton->dx[i]*v->dx[i] + h_prev_newton->dy[i]*v->dy[i])
                         - dKdh(h_prev_newton->val[i], elem_marker) * u->dy[i] * v->val[i]
                         - ddKdhh(h_prev_newton->val[i], elem_marker) * u->val[i] * h_prev_newton->dy[i] * v->val[i]
                       )
            + wt[i] * 0.5 * ( // explicit Euler part, 
		         C(h_prev_time->val[i], elem_marker) * u->val[i] * v->val[i] / time_step
                       );
  return result;
}

// Integration order for Jacobian matrix - volumetric part
Ord jac_form_vol_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                     Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}

// Residual vector - volumetric part
double res_form_vol_euler(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                    Geom<double> *e, ExtData<double> *ext)
{
  int elem_marker = e->elem_marker;
  double result = 0;
  Func<double>* h_prev_newton = u_ext[0];
  Func<double>* h_prev_time = ext->fn[0];
  for (int i = 0; i < n; i++) {
    result += wt[i] * (
		       C(h_prev_newton->val[i], elem_marker) * (h_prev_newton->val[i] - h_prev_time->val[i]) * v->val[i] / time_step
                       + K(h_prev_newton->val[i], elem_marker) * (h_prev_newton->dx[i] * v->dx[i] + h_prev_newton->dy[i] * v->dy[i])
                       - dKdh(h_prev_newton->val[i], elem_marker) * h_prev_newton->dy[i] * v->val[i]
                      );
  }
  return result;
}

// Residual vector - volumetric part
double res_form_vol_cranic(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                    Geom<double> *e, ExtData<double> *ext)
{
  double result = 0;
  int elem_marker = e->elem_marker;
  Func<double>* h_prev_newton = u_ext[0];
  Func<double>* h_prev_time = ext->fn[0];
  for (int i = 0; i < n; i++) {
    result += wt[i] * 0.5 * ( // implicit Euler part
		       C(h_prev_newton->val[i], elem_marker) * (h_prev_newton->val[i] - h_prev_time->val[i]) * v->val[i] / time_step
                       + K(h_prev_newton->val[i], elem_marker) * (h_prev_newton->dx[i] * v->dx[i] + h_prev_newton->dy[i] * v->dy[i])
                       - dKdh(h_prev_newton->val[i], elem_marker) * h_prev_newton->dy[i] * v->val[i]
                      )
            + wt[i] * 0.5 * ( // explicit Euler part
		       C(h_prev_time->val[i], elem_marker) * (h_prev_newton->val[i] - h_prev_time->val[i]) * v->val[i] / time_step
                       + K(h_prev_time->val[i], elem_marker) * (h_prev_time->dx[i] * v->dx[i] + h_prev_time->dy[i] * v->dy[i])
                       - dKdh(h_prev_time->val[i], elem_marker) * h_prev_time->dy[i] * v->val[i]
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

/*** PICARD ***/

// Bilinear form for implicit Euler.
double bilinear_form_picard_euler(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                                  Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
  int elem_marker = e->elem_marker;
  double result = 0;
  Func<double>* h_prev_picard = ext->fn[0];

  for (int i = 0; i < n; i++) {
    result += wt[i] * (
                         C(h_prev_picard->val[i], elem_marker) * u->val[i] * v->val[i] / time_step
                         + K(h_prev_picard->val[i], elem_marker) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i])
                         - dKdh(h_prev_picard->val[i], elem_marker) * u->dy[i] * v->val[i]
		       );
  }

  return result;
}

Ord bilinear_form_picard_euler_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
                             Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}

// Right-hand side for implicit Euler.
double linear_form_picard_euler(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                                Geom<double> *e, ExtData<double> *ext)
{
  int elem_marker = e->elem_marker;
  double result = 0;
  Func<double>* h_prev_picard = ext->fn[0];
  Func<double>* h_prev_time = ext->fn[1];
  for (int i = 0; i < n; i++) result += wt[i] * C(h_prev_picard->val[i], elem_marker) * h_prev_time->val[i] * v->val[i] / time_step;
  return result;
}

Ord linear_form_picard_euler_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                                 Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}

// TODO:

// Bilinear form for Crank-Nicolson.
// Right-hand side for Crank-Nicolson.
