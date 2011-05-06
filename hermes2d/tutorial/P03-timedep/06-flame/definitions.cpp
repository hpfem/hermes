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

class CustomWeakForm : public WeakForm
{
public:
  CustomWeakForm(Hermes::vector<std::string> newton_boundaries, double time_step, 
    double Le, double kappa, Solution* sln_prev_time, bool JFNK = false) : WeakForm(2, JFNK)
  {
    // Jacobian forms - volumetric.
    add_matrix_form(new JacobianFormVol(0, 0, heatcap, rho, lambda, tau));

    // Jacobian forms - surface.
    add_matrix_form_surf(new JacobianFormSurf(0, 0, newton_boundaries, alpha, lambda));


    // Residual forms - volumetric.
    ResidualFormVol* res_form = new ResidualFormVol(0, heatcap, rho, lambda, tau);
    res_form->ext.push_back(sln_prev_time);
    add_vector_form(res_form);

    // Residual forms - surface.
    add_vector_form_surf(new ResidualFormSurf(0, newton_boundaries, alpha, lambda, temp_ext));
  }

  ~CustomWeakForm() {}

private:
  class newton_bilinear_form_0_0 : public WeakForm::MatrixFormVol
  {
  public:
    newton_bilinear_form_0_0(int i, int j, double time_step) : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_SYM),
    time_step(time_step) { }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      Func<scalar>* domegadt = ext->fn[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * (  1.5 * u->val[i] * v->val[i] / time_step
                        +  u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]
                        - domegadt->val[i] * u->val[i] * v->val[i] );
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      // Returning the sum of the degrees of the basis and test function plus two.
      return Ord(10);
    }

    double time_step;
  };

  class newton_bilinear_form_0_1 : public WeakForm::MatrixFormVol
  {
  public:
    newton_bilinear_form_0_1(int i, int j) : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_SYM) { }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      Func<Real>* domegady = ext->fn[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * (- domegady->val[i] * u->val[i] * v->val[i] );
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      // Returning the sum of the degrees of the basis and test function plus two.
      return Ord(10);
    }

  };

  class newton_bilinear_form_1_0 : public WeakForm::MatrixFormVol
  {
  public:
    newton_bilinear_form_1_0(int i, int j) : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_SYM) { }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      Func<Real>* domegadt = ext->fn[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * ( domegadt->val[i] * u->val[i] * v->val[i] );
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      // Returning the sum of the degrees of the basis and test function plus two.
      return Ord(10);
    }

  };

  class newton_bilinear_form_1_1 : public WeakForm::MatrixFormVol
  {
  public:
    newton_bilinear_form_1_1(int i, int j, double time_step, double Le) : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_SYM), 
    time_step(time_step), Le(Le) { }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      Func<Real>* domegady = ext->fn[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * (  1.5 * u->val[i] * v->val[i] / time_step
                          +  (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]) / Le
                          + domegady->val[i] * u->val[i] * v->val[i] );
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      // Returning the sum of the degrees of the basis and test function plus two.
      return Ord(10);
    }

  double time_step, Le;

  };

  class newton_bilinear_form_0_0_surf : public WeakForm::MatrixFormSurf
  {
  public:
    newton_bilinear_form_0_0_surf(int i, int j, Hermes::vector<std::string> newton_boundaries, double kappa) : WeakForm::MatrixFormSurf(i, j, newton_boundaries), kappa(kappa) { }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * (kappa * u->val[i] * v->val[i]);
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      // Returning the sum of the degrees of the basis and test function plus two.
      return Ord(10);
    }

    double kappa;
  };


  class newton_linear_form_0 : public WeakForm::VectorFormVol
  {
  public:
    newton_linear_form_0(int i, double time_step) : WeakForm::VectorFormVol(i),
    time_step(time_step)  {}

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      for (int i = 0; i < n; i++)
      Func<Real>* titer = u_ext[0];
      Func<Real>* t_prev_time_1 = ext->fn[0];
      Func<Real>* t_prev_time_2 = ext->fn[1];
      Func<Real>* omega = ext->fn[2];
      for (int i = 0; i < n; i++)
        result += wt[i] * ( (3.0 * titer->val[i] - 4.0 * t_prev_time_1->val[i] 
                             + t_prev_time_2->val[i]) * v->val[i] / (2.0 * time_step) +
                            (titer->dx[i] * v->dx[i] + titer->dy[i] * v->dy[i]) -
                            omega->val[i] * v->val[i]);
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const {
      // Returning the sum of the degrees of the test function and solution plus two.
      return Ord(10);
    }

  private:
    double time_step;
  };

  class newton_linear_form_1 : public WeakForm::VectorFormVol
  {
  public:
    newton_linear_form_1(int i, double time_step, double Le) : WeakForm::VectorFormVol(i),
    time_step(time_step), Le(Le)  {}

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      Func<Real>* c_prev_newton = u_ext[1];
      Func<Real>* c_prev_time_1 = ext->fn[0];
      Func<Real>* c_prev_time_2 = ext->fn[1];
      Func<Real>* omega = ext->fn[2];
      for (int i = 0; i < n; i++)
        result += wt[i] * ( (3.0 * c_prev_newton->val[i] - 4.0 * c_prev_time_1->val[i] + c_prev_time_2->val[i])
                        * v->val[i] / (2.0 * time_step) +
                        (c_prev_newton->dx[i] * v->dx[i] + c_prev_newton->dy[i] * v->dy[i]) / Le +
                        omega->val[i] * v->val[i]);
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const {
      // Returning the sum of the degrees of the test function and solution plus two.
      return Ord(10);
    }

  private:
    double time_step, Le;
  };


  class newton_linear_form_0_surf : public WeakForm::VectorFormSurf
  {
  public:
    newton_linear_form_0_surf(int i, Hermes::vector<std::string> newton_boundaries, double kappa) : WeakForm::VectorFormSurf(i, newton_boundaries), kappa(kappa)  {}

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      Func<Real>* t_prev_newton = u_ext[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * (kappa * t_prev_newton->val[i] * v->val[i]);
      return result;  
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const {
      // Returning the sum of the degrees of the test function and solution plus two.
      return Ord(10);
    }

  private:
    double kappa;
  };
};
