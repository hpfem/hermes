/*** Definition of residiual vectors ***/

template<class Real, class Scalar>
Scalar Fc_euler(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
	Scalar result = 0;
	Func<Scalar>* C_prev_time = ext->fn[0];
	Func<Scalar>* C_prev_newton = u_ext[0];
	Func<Scalar>* phi_prev_newton = u_ext[1];
	for (int i = 0; i < n; i++) {
		result += wt[i] * ((C_prev_newton->val[i] - C_prev_time->val[i]) * v->val[i] / (*TAU) +
				D * (C_prev_newton->dx[i] * v->dx[i] + C_prev_newton->dy[i] * v->dy[i]) +
				K * C_prev_newton->val[i] * (phi_prev_newton->dx[i] * v->dx[i] + phi_prev_newton->dy[i] * v->dy[i]));
	}
	return result;
}

template<class Real, class Scalar>
Scalar Fphi_euler(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
	Scalar result = 0;
	Func<Scalar>* C_prev_newton = u_ext[0];
	Func<Scalar>* phi_prev_newton = u_ext[1];

	for (int i = 0; i < n; i++) {

	  result += wt[i] * ((phi_prev_newton->dx[i] * v->dx[i] + phi_prev_newton->dy[i] * v->dy[i]) +
					L * v->val[i] * (C0 - C_prev_newton->val[i]));
	}

	return result;
}

// matrix 0_0
template<class Real, class Scalar>
Scalar J_euler_DFcDYc(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
	Scalar result = 0;
	Func<Scalar>* phi_prev_newton = u_ext[1];
	for (int i = 0; i < n; i++) {
		result += wt[i] * (u->val[i] * v->val[i] / (*TAU) +
				D * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]) +
				K * u->val[i] * (phi_prev_newton->dx[i] * v->dx[i] + phi_prev_newton->dy[i] * v->dy[i]));
	}
	return result;
}

//matrix 0_1
template<class Real, class Scalar>
Scalar J_euler_DFcDYphi(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
	Scalar result = 0;
	Func<Scalar>* C_prev_newton = u_ext[0];
	for (int i = 0; i < n; i++) {
		result += wt[i] * K * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]) * C_prev_newton->val[i];
	}
	return result;
}

//matrix 1_0
template<class Real, class Scalar>
Scalar J_euler_DFphiDYc(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
	Scalar result = 0;
	for (int i = 0; i < n; i++) {
		result += wt[i] * ( -L * u->val[i] * v->val[i]);
	}
	return result;
}

//matrix 1_1
template<class Real, class Scalar>
Scalar J_euler_DFphiDYphi(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
	Scalar result = 0;
	for (int i = 0; i < n; i++) {
		result += wt[i] * ( u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
	}
	return result;
}


// Cranck-Nicholson forms

template<class Real, class Scalar>
Scalar Fc_cranic(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
  Scalar result = 0;
  Func<Scalar>* C_prev_time = ext->fn[0];
  Func<Scalar>* phi_prev_time = ext->fn[1];
  Func<Scalar>* C_prev_newton = u_ext[0];
  Func<Scalar>* phi_prev_newton = u_ext[1];
  for (int i = 0; i < n; i++) {
    result += wt[i] * ((C_prev_newton->val[i] - C_prev_time->val[i]) * v->val[i] / (*TAU) +
        0.5 * D * (C_prev_newton->dx[i] * v->dx[i] + C_prev_newton->dy[i] * v->dy[i]) +
        0.5 * D * (C_prev_time->dx[i] * v->dx[i] + C_prev_time->dy[i] * v->dy[i]) +
        0.5 * K * C_prev_newton->val[i] * (phi_prev_newton->dx[i] * v->dx[i] + phi_prev_newton->dy[i] * v->dy[i]) +
        0.5 * K * C_prev_time->val[i] * (phi_prev_time->dx[i] * v->dx[i] + phi_prev_time->dy[i] * v->dy[i]));
  }
  return result;
}

template<class Real, class Scalar>
Scalar Fphi_cranic(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
    Scalar result = 0;
    Func<Scalar>* C_prev_newton = u_ext[0];
    Func<Scalar>* phi_prev_newton = u_ext[1];
    for (int i = 0; i < n; i++) {
      result += wt[i] * ((phi_prev_newton->dx[i] * v->dx[i] + phi_prev_newton->dy[i] * v->dy[i]) +
            L * v->val[i] * (C0 - C_prev_newton->val[i]));
    }
    return result;
}


template<class Real, class Scalar>
Scalar J_cranic_DFcDYc(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
  Scalar result = 0;
  Func<Scalar>* phi_prev_newton = u_ext[1];
  for (int i = 0; i < n; i++) {
    result += wt[i] * (u->val[i] * v->val[i] / (*TAU) +
        0.5 * D * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]) +
        0.5 * K * u->val[i] * (phi_prev_newton->dx[i] * v->dx[i] + phi_prev_newton->dy[i] * v->dy[i]));
  }
  return result;
}

template<class Real, class Scalar>
Scalar J_cranic_DFcDYphi(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
  Scalar result = 0;
  Func<Scalar>* C_prev_newton = u_ext[0];
	for (int i = 0; i < n; i++) {
	  result += wt[i] * (0.5 * K * C_prev_newton->val[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
	}
	return result;
}

template<class Real, class Scalar>
Scalar J_cranic_DFphiDYc(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
  Scalar result = 0;
  for (int i = 0; i < n; i++) {
    result += wt[i] * ( -L * u->val[i] * v->val[i]);
  }
  return result;
}

template<class Real, class Scalar>
Scalar J_cranic_DFphiDYphi(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
  Scalar result = 0;
  for (int i = 0; i < n; i++) {
    result += wt[i] * ( u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
  }
  return result;
}
