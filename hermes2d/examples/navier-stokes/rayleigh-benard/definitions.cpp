#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

class WeakFormRayleighBenard : public WeakForm
{
public:
  WeakFormRayleighBenard(double Pr) : WeakForm(4), Pr(Pr) {
    BilinearFormSymVel* sym_form_0 = new BilinearFormSymVel(0, 0, Pr);
    add_matrix_form(sym_form_0);
    BilinearFormSymVel* sym_form_1 = new BilinearFormSymVel(1, 1, Pr);
    add_matrix_form(sym_form_1);

    BilinearFormNonsymVel_0_0* nonsym_vel_form_0_0 = new BilinearFormNonsymVel_0_0(0, 0);
    add_matrix_form(nonsym_vel_form_0_0);
    BilinearFormNonsymVel_0_1* nonsym_vel_form_0_1 = new BilinearFormNonsymVel_0_1(0, 1);
    add_matrix_form(nonsym_vel_form_0_1);
    BilinearFormNonsymVel_1_0* nonsym_vel_form_1_0 = new BilinearFormNonsymVel_1_0(1, 0);
    add_matrix_form(nonsym_vel_form_1_0);
    BilinearFormNonsymVel_1_1* nonsym_vel_form_1_1 = new BilinearFormNonsymVel_1_1(1, 1);
    add_matrix_form(nonsym_vel_form_1_1);

    BilinearFormNonsymXVelPressure* nonsym_velx_pressure_form = new BilinearFormNonsymXVelPressure(0, 2);
    add_matrix_form(nonsym_velx_pressure_form);

    BilinearFormNonsymYVelPressure* nonsym_vely_pressure_form = new BilinearFormNonsymYVelPressure(1, 2);
    add_matrix_form(nonsym_vely_pressure_form);
    
    VectorFormNS_0* F_0 = new VectorFormNS_0(0, Pr);
    VectorFormNS_1* F_1 = new VectorFormNS_1(1, Pr);
    VectorFormNS_2* F_2 = new VectorFormNS_2(2);
  };

  class BilinearFormSymVel : public WeakForm::MatrixFormVol
  {
  public:
    BilinearFormSymVel(int i, int j, double Pr) 
            : WeakForm::MatrixFormVol(i, j, HERMES_SYM), Pr(Pr) {
      adapt_eval = false;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = int_grad_u_grad_v<double, scalar>(n, wt, u, v) / Pr;
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const {
      Ord result = int_grad_u_grad_v<Ord, Ord>(n, wt, u, v) / Pr;
      return result;
    }
  protected:
    double Pr;
  };


  class BilinearFormNonsymVel_0_0 : public WeakForm::MatrixFormVol
  {
  public:
    BilinearFormNonsymVel_0_0(int i, int j) : WeakForm::MatrixFormVol(i, j, HERMES_NONSYM) {
      adapt_eval = false;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      Func<scalar>* xvel_prev_newton = u_ext[0];
      Func<scalar>* yvel_prev_newton = u_ext[1];
      for (int i = 0; i < n; i++)
        result += wt[i] * ((xvel_prev_newton->val[i] * u->dx[i] + yvel_prev_newton->val[i]
                            * u->dy[i]) * v->val[i] + u->val[i] * v->val[i] * xvel_prev_newton->dx[i]);
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const {
    Ord result = 0;
      Func<Ord>* xvel_prev_newton = u_ext[0];
      Func<Ord>* yvel_prev_newton = u_ext[1];
      for (int i = 0; i < n; i++)
        result += wt[i] * ((xvel_prev_newton->val[i] * u->dx[i] + yvel_prev_newton->val[i]
                            * u->dy[i]) * v->val[i] + u->val[i] * v->val[i] * xvel_prev_newton->dx[i]);
      return result;
    }
  };


  class BilinearFormNonsymVel_0_1 : public WeakForm::MatrixFormVol
  {
  public:
    BilinearFormNonsymVel_0_1(int i, int j) : WeakForm::MatrixFormVol(i, j, HERMES_NONSYM) {
      adapt_eval = false;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      Func<scalar>* xvel_prev_newton = u_ext[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->val[i] * v->val[i] * xvel_prev_newton->dy[i]);

      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const {
      Ord result = 0;

      Func<Ord>* xvel_prev_newton = u_ext[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->val[i] * v->val[i] * xvel_prev_newton->dy[i]);
 
      return result;
    }
  };


  class BilinearFormNonsymVel_1_0 : public WeakForm::MatrixFormVol
  {
  public:
    BilinearFormNonsymVel_1_0(int i, int j) : WeakForm::MatrixFormVol(i, j, HERMES_NONSYM) {
      adapt_eval = false;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      Func<scalar>* yvel_prev_newton = u_ext[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->val[i] * v->val[i] * yvel_prev_newton->dx[i]);
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const {
      Ord result = 0;
      Func<Ord>* yvel_prev_newton = u_ext[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->val[i] * v->val[i] * yvel_prev_newton->dx[i]);
      return result;
    }
  };


  class BilinearFormNonsymVel_1_1 : public WeakForm::MatrixFormVol
  {
  public:
    BilinearFormNonsymVel_1_1(int i, int j) : WeakForm::MatrixFormVol(i, j, HERMES_NONSYM) {
      adapt_eval = false;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      Func<scalar>* xvel_prev_newton = u_ext[0];
      Func<scalar>* yvel_prev_newton = u_ext[1];
      for (int i = 0; i < n; i++)
        result += wt[i] * ((xvel_prev_newton->val[i] * u->dx[i] + yvel_prev_newton->val[i] * u->dy[i]) * v->val[i] + u->val[i]
                           * v->val[i] * yvel_prev_newton->dy[i]);
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const {
      Ord result = 0;
      Func<Ord>* xvel_prev_newton = u_ext[0];
      Func<Ord>* yvel_prev_newton = u_ext[1];
      for (int i = 0; i < n; i++)
        result += wt[i] * ((xvel_prev_newton->val[i] * u->dx[i] + yvel_prev_newton->val[i] * u->dy[i]) * v->val[i] + u->val[i]
                           * v->val[i] * yvel_prev_newton->dy[i]);
      return result;
    }
  };


  class BilinearFormNonsymXVelPressure : public WeakForm::MatrixFormVol
  {
  public:
    BilinearFormNonsymXVelPressure(int i, int j) : WeakForm::MatrixFormVol(i, j, HERMES_ANTISYM) {
      adapt_eval = false;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      return - int_u_dvdx<double, scalar>(n, wt, u, v);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const {
      return - int_u_dvdx<Ord, Ord>(n, wt, u, v);
    }
  };


  class BilinearFormNonsymYVelPressure : public WeakForm::MatrixFormVol
  {
  public:
    BilinearFormNonsymYVelPressure(int i, int j) : WeakForm::MatrixFormVol(i, j, HERMES_ANTISYM) {
      adapt_eval = false;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      return - int_u_dvdy<double, scalar>(n, wt, u, v);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const {
      return - int_u_dvdy<Ord, Ord>(n, wt, u, v);
    }
  };


  class VectorFormNS_0 : public WeakForm::VectorFormVol
  {
  public:
    VectorFormNS_0(int i, double Pr) : WeakForm::VectorFormVol(i), Pr(Pr) {
      adapt_eval = false;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, 
                         ExtData<scalar> *ext) const {
      scalar result = 0;
      Func<scalar>* xvel_prev_time = ext->fn[0];  
      Func<scalar>* yvel_prev_time = ext->fn[1];
      Func<scalar>* xvel_prev_newton = u_ext[0];  
      Func<scalar>* yvel_prev_newton = u_ext[1];  
      Func<scalar>* p_prev_newton = u_ext[2];
      for (int i = 0; i < n; i++)
        result += wt[i] * ((xvel_prev_newton->dx[i] * v->dx[i] + xvel_prev_newton->dy[i] * v->dy[i]) / Pr 
                          - (p_prev_newton->val[i] * v->dx[i]));
      for (int i = 0; i < n; i++)
        result += wt[i] * (
                          + ((xvel_prev_newton->val[i] * xvel_prev_newton->dx[i] 
                          + yvel_prev_newton->val[i] * xvel_prev_newton->dy[i]) * v->val[i]));
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      Ord result = 0;
      Func<Ord>* xvel_prev_time = ext->fn[0];  
      Func<Ord>* yvel_prev_time = ext->fn[1];
      Func<Ord>* xvel_prev_newton = u_ext[0];  
      Func<Ord>* yvel_prev_newton = u_ext[1];  
      Func<Ord>* p_prev_newton = u_ext[2];
      for (int i = 0; i < n; i++)
        result += wt[i] * ((xvel_prev_newton->dx[i] * v->dx[i] + xvel_prev_newton->dy[i] * v->dy[i]) / Pr 
                            - (p_prev_newton->val[i] * v->dx[i]));
      for (int i = 0; i < n; i++)
        result += wt[i] * (
                          ((xvel_prev_newton->val[i] * xvel_prev_newton->dx[i] 
                          + yvel_prev_newton->val[i] * xvel_prev_newton->dy[i]) * v->val[i]));
      return result;
    }
  protected:
    double Pr;
  };

  class VectorFormNS_1 : public WeakForm::VectorFormVol
  {
  public:
    VectorFormNS_1(int i, double Pr) 
      : WeakForm::VectorFormVol(i), Pr(Pr) {
      adapt_eval = false;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, 
                         ExtData<scalar> *ext) const {
      scalar result = 0;
      Func<scalar>* xvel_prev_time = ext->fn[0];  
      Func<scalar>* yvel_prev_time = ext->fn[1];
      Func<scalar>* xvel_prev_newton = u_ext[0];  
      Func<scalar>* yvel_prev_newton = u_ext[1];  
      Func<scalar>* p_prev_newton = u_ext[2];
      for (int i = 0; i < n; i++)
        result += wt[i] * ((yvel_prev_newton->dx[i] * v->dx[i] + yvel_prev_newton->dy[i] * v->dy[i]) / Pr 
                          - (p_prev_newton->val[i] * v->dy[i]));
      for (int i = 0; i < n; i++)
        result += wt[i] * (
                          ((xvel_prev_newton->val[i] * yvel_prev_newton->dx[i] 
                          + yvel_prev_newton->val[i] * yvel_prev_newton->dy[i]) * v->val[i]));
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      Ord result = 0;
      Func<Ord>* xvel_prev_time = ext->fn[0];  
      Func<Ord>* yvel_prev_time = ext->fn[1];
      Func<Ord>* xvel_prev_newton = u_ext[0];  
      Func<Ord>* yvel_prev_newton = u_ext[1];  
      Func<Ord>* p_prev_newton = u_ext[2];
      for (int i = 0; i < n; i++)
        result += wt[i] * ((xvel_prev_newton->dx[i] * v->dx[i] + xvel_prev_newton->dy[i] * v->dy[i]) / Pr 
                  - (p_prev_newton->val[i] * v->dx[i]));
      for (int i = 0; i < n; i++)
        result += wt[i] * (
                          ((xvel_prev_newton->val[i] * xvel_prev_newton->dx[i] 
                          + yvel_prev_newton->val[i] * xvel_prev_newton->dy[i]) * v->val[i]));
      return result;
    }
  protected:
    double Pr;
  };

  class VectorFormNS_2 : public WeakForm::VectorFormVol
  {
  public:
    VectorFormNS_2(int i) : WeakForm::VectorFormVol(i) {
      adapt_eval = false;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, 
                         ExtData<scalar> *ext) const {
      scalar result = 0;
      Func<scalar>* xvel_prev_newton = u_ext[0];  
      Func<scalar>* yvel_prev_newton = u_ext[1];  

      for (int i = 0; i < n; i++)
        result += wt[i] * (xvel_prev_newton->dx[i] * v->val[i] + yvel_prev_newton->dy[i] * v->val[i]);
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      Ord result = 0;
      Func<Ord>* xvel_prev_newton = u_ext[0];  
      Func<Ord>* yvel_prev_newton = u_ext[1];  

      for (int i = 0; i < n; i++)
        result += wt[i] * (xvel_prev_newton->dx[i] * v->val[i] + yvel_prev_newton->dy[i] * v->val[i]);
      return result;
    }
  };

protected:
  double Pr;
};

