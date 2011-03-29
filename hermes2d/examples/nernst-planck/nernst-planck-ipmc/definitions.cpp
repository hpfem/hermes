class CustomWeakFormNernstPlanckEuler : public WeakForm
{
public:
  CustomWeakFormNernstPlanckEuler(double* tau, scalar C0, double lin_force_coup, double mech_lambda, double mech_mu, double K, double L, double D, Solution* C_prev_time) : WeakForm(4) {
    for(unsigned int i = 0; i < 4; i++) {
      CustomWeakFormNernstPlanckEuler::VectorFormVol* vector_form = new CustomWeakFormNernstPlanckEuler::VectorFormVol(i, tau, C0, lin_force_coup, mech_lambda, mech_mu, K, L, D);
      if(i == 0)
        vector_form->ext.push_back(C_prev_time);
      add_vector_form(vector_form);
      for(unsigned int j = 0; j < 4; j++)
        add_matrix_form(new CustomWeakFormNernstPlanckEuler::MatrixFormVol(i, j, tau, C0, lin_force_coup, mech_lambda, mech_mu, K, L, D));
    }
  };

private:
  class MatrixFormVol : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVol(int i, int j, double* tau, double C0, double lin_force_coup, double mech_lambda, double mech_mu, double K, double L, double D) 
      : WeakForm::MatrixFormVol(i, j), i(i), j(j), tau(tau), C0(C0), lin_force_coup(lin_force_coup), mech_lambda(mech_lambda), mech_mu(mech_mu), K(K), L(L), D(D) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      Scalar result = 0;
	    Func<Scalar>* prev_newton;
      switch(i * 10 + j) {
        case 0:
	        prev_newton = u_ext[1];
	        for (int i = 0; i < n; i++) {
		        result += wt[i] * (u->val[i] * v->val[i] / *(this->tau) +
				        this->D * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]) +
				        this->K * u->val[i] * (prev_newton->dx[i] * v->dx[i] + prev_newton->dy[i] * v->dy[i]));
	        }
	        return result;
          break;
        case 1:
	        prev_newton = u_ext[0];
	        for (int i = 0; i < n; i++) {
		        result += wt[i] * this->K * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]) * prev_newton->val[i];
	        }
	        return result;
          break;
        case 2:
          return result;
          break;
        case 3:
          return result;
        break;
        case 10:
	        for (int i = 0; i < n; i++) {
		        result += wt[i] * ( -this->L * u->val[i] * v->val[i]);
	        }
	        return result;
          break;
        case 11:
	        for (int i = 0; i < n; i++) {
		        result += wt[i] * ( u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
	        }
	        return result;
          break;
        case 12:
          return result;
          break;
        case 13:
          return result;
        break;
        case 20:
          for (int i = 0; i < n; i++) {
            result += wt[i] * this->lin_force_coup * u->val[i] * v->val[i];
          }
          return result;
          break;
        case 21:
          return result;
        break;
        case 22:
          for (int i = 0; i < n; i++) {
            result += wt[i] * ((2 * this->mech_mu + this->mech_lambda) * u->dx[i] * v->dx[i] + this->mech_mu * u->dy[i] * v->dy[i]);
          }
          return result;
        break;
        case 23:
          for (int i = 0; i < n; i++) {
            result += wt[i] * (this->mech_mu * u->dx[i] * v->dy[i] + this->mech_lambda * u->dy[i] * v->dx[i]);
          }
          return result;
        break;
        case 30:
          return result;
        break;
        case 31:
          return result;
        break;
        case 32:
          for (int i = 0; i < n; i++) {
            result += wt[i] * (this->mech_mu * u->dy[i] * v->dx[i] + this->mech_lambda * u->dx[i] * v->dy[i]);
          }
          return result;
        break;
        case 33:
          
          for (int i = 0; i < n; i++) {
            result += wt[i] * ((2 * this->mech_mu + this->mech_lambda) * u->dy[i] * v->dy[i] + this->mech_mu * u->dx[i] * v->dx[i]);
          }
          return result;
          break;
        default:
          
          return result;
      }
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    // Members.
    int i, j;
    double* tau;
    double C0;
    double lin_force_coup;
    double mech_lambda;
    double mech_mu;
    double K;
    double L;
    double D;
  };

  class VectorFormVol : public WeakForm::VectorFormVol
  {
  public:
    VectorFormVol(int i, double* tau, double C0, double lin_force_coup, double mech_lambda, double mech_mu, double K, double L, double D) 
      : WeakForm::VectorFormVol(i), i(i), tau(tau), C0(C0), lin_force_coup(lin_force_coup), mech_lambda(mech_lambda), mech_mu(mech_mu), K(K), L(L), D(D) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], 
                        Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      Scalar result = 0;
      Func<Scalar>* C_prev_time;
	    Func<Scalar>* C_prev_newton;
	    Func<Scalar>* phi_prev_newton;
      Func<Scalar>* u1_prev_newton;
      Func<Scalar>* u2_prev_newton;
      switch(i) {
        case 0:
	        C_prev_time = ext->fn[0];
	        C_prev_newton = u_ext[0];
	        phi_prev_newton = u_ext[1];
	        for (int i = 0; i < n; i++) {
		        result += wt[i] * ((C_prev_newton->val[i] - C_prev_time->val[i]) * v->val[i] / *(this->tau) +
				        this->D * (C_prev_newton->dx[i] * v->dx[i] + C_prev_newton->dy[i] * v->dy[i]) +
				        this->K * C_prev_newton->val[i] * (phi_prev_newton->dx[i] * v->dx[i] + phi_prev_newton->dy[i] * v->dy[i]));
	        }
	        return result;
          break;
        case 1:
	        C_prev_newton = u_ext[0];
	        phi_prev_newton = u_ext[1];
	        for (int i = 0; i < n; i++) {
	          result += wt[i] * ((phi_prev_newton->dx[i] * v->dx[i] + phi_prev_newton->dy[i] * v->dy[i]) +
					        this->L * v->val[i] * (this->C0 - C_prev_newton->val[i]));
	        }
	        return result;
          break;
        case 2:
          C_prev_newton = u_ext[0];
          u1_prev_newton = u_ext[2];
          u2_prev_newton = u_ext[3];
          for (int i = 0; i < n; i++) {
            result += wt[i] * ((2*this->mech_mu + this->mech_lambda) * u1_prev_newton->dx[i] * v->dx[i] + this->mech_mu * u1_prev_newton->dy[i] * v->dy[i] +
                this->mech_mu * u2_prev_newton->dx[i] * v->dy[i] + this->mech_lambda * u2_prev_newton->dy[i] * v->dx[i] -
                this->lin_force_coup * (C_prev_newton->val[i] - this->C0));
          }
          return result;
          break;
        case 3:
          u1_prev_newton = u_ext[2];
          u2_prev_newton = u_ext[3];
          for (int i = 0; i < n; i++) {
            result += wt[i] * ((2*this->mech_mu + this->mech_lambda) * u2_prev_newton->dy[i] * v->dy[i] + this->mech_mu * u2_prev_newton->dx[i] * v->dx[i] +
                this->mech_mu * u1_prev_newton->dy[i] * v->dx[i] + this->mech_lambda * u1_prev_newton->dx[i] * v->dy[i]);
          }
          return result;
          break;
        default:
          return result;
      }
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    // Members.
    int i;
    double* tau;
    double C0;
    double lin_force_coup;
    double mech_lambda;
    double mech_mu;
    double K;
    double L;
    double D;
  };
};

// Cranck-Nicholson forms
// TODO: finish this & use it.
template<class Real, class Scalar>
Scalar Fc_cranic(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
  Scalar result = 0;
  Func<Scalar>* C_prev_time = ext->fn[0];
  Func<Scalar>* phi_prev_time = ext->fn[1];
  Func<Scalar>* C_prev_newton = u_ext[0];
  Func<Scalar>* phi_prev_newton = u_ext[1];
  for (int i = 0; i < n; i++) {
    result += wt[i] * ((C_prev_newton->val[i] - C_prev_time->val[i]) * v->val[i] / *(this->tau) +
        0.5 * this->D * (C_prev_newton->dx[i] * v->dx[i] + C_prev_newton->dy[i] * v->dy[i]) +
        0.5 * this->D * (C_prev_time->dx[i] * v->dx[i] + C_prev_time->dy[i] * v->dy[i]) +
        0.5 * this->K * C_prev_newton->val[i] * (phi_prev_newton->dx[i] * v->dx[i] + phi_prev_newton->dy[i] * v->dy[i]) +
        0.5 * this->K * C_prev_time->val[i] * (phi_prev_time->dx[i] * v->dx[i] + phi_prev_time->dy[i] * v->dy[i]));
  }
  return result;
}

template<class Real, class Scalar>
Scalar Fphi_cranic(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
    Scalar result = 0;
    Func<Scalar>* C_prev_newton = u_ext[0];
    Func<Scalar>* phi_prev_newton = u_ext[1];
    for (int i = 0; i < n; i++) {
      result += wt[i] * ((phi_prev_newton->dx[i] * v->dx[i] + phi_prev_newton->dy[i] * v->dy[i]) +
            this->L * v->val[i] * (this->C0 - C_prev_newton->val[i]));
    }
    return result;
}


template<class Real, class Scalar>
Scalar J_cranic_DFcDYc(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
  Scalar result = 0;
  Func<Scalar>* phi_prev_newton = u_ext[1];
  for (int i = 0; i < n; i++) {
    result += wt[i] * (u->val[i] * v->val[i] / *(this->tau) +
        0.5 * this->D * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]) +
        0.5 * this->K * u->val[i] * (phi_prev_newton->dx[i] * v->dx[i] + phi_prev_newton->dy[i] * v->dy[i]));
  }
  return result;
}

template<class Real, class Scalar>
Scalar J_cranic_DFcDYphi(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
  Scalar result = 0;
  Func<Scalar>* C_prev_newton = u_ext[0];
	for (int i = 0; i < n; i++) {
	  result += wt[i] * (0.5 * this->K * C_prev_newton->val[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
	}
	return result;
}

template<class Real, class Scalar>
Scalar J_cranic_DFphiDYc(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
  Scalar result = 0;
  for (int i = 0; i < n; i++) {
    result += wt[i] * ( -this->L * u->val[i] * v->val[i]);
  }
  return result;
}

template<class Real, class Scalar>
Scalar J_cranic_DFphiDYphi(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
  Scalar result = 0;
  for (int i = 0; i < n; i++) {
    result += wt[i] * ( u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
  }
  return result;
}
