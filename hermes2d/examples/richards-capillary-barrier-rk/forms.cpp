// Weak forms for the Newton's method. 

// Jacobian matrix - volumetric part

double jac_form_vol(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                    Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
  int elem_marker = e->elem_marker;
  double result = 0;
  Func<double>* h_prev_newton = u_ext[0];

  for (int i = 0; i < n; i++)
    result += wt[i] * (
		       //C(h_prev_newton->val[i], elem_marker) * u->val[i] * v->val[i] / TAU
		       //+ dCdh(h_prev_newton->val[i], elem_marker) * u->val[i] * h_prev_newton->val[i] * v->val[i] / TAU
		       //- dCdh(h_prev_newton->val[i], elem_marker) * u->val[i] * h_prev_time->val[i] * v->val[i] / TAU
			 + K(h_prev_newton->val[i], elem_marker) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i])
                         + dKdh(h_prev_newton->val[i], elem_marker) * u->val[i] * 
                           (h_prev_newton->dx[i]*v->dx[i] + h_prev_newton->dy[i]*v->dy[i])
                         - dKdh(h_prev_newton->val[i], elem_marker) * u->dy[i] * v->val[i]
                         - ddKdhh(h_prev_newton->val[i], elem_marker) * u->val[i] * h_prev_newton->dy[i] * v->val[i]
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
double res_form_vol(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                    Geom<double> *e, ExtData<double> *ext)
{
  int elem_marker = e->elem_marker;
  double result = 0;
  Func<double>* h_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++) {
    result += -wt[i] * (
                       K(h_prev_newton->val[i], elem_marker) * (h_prev_newton->dx[i] * v->dx[i] + h_prev_newton->dy[i] * v->dy[i])
                       - dKdh(h_prev_newton->val[i], elem_marker) * h_prev_newton->dy[i] * v->val[i]
		       ) / C(h_prev_newton->val[i], elem_marker);
  }
  return result;
}

// Integration order for residual vector - volumetric part
Ord res_form_vol_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                     Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}
