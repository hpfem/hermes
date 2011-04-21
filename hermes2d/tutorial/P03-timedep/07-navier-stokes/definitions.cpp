#include "weakform/weakform.h"
#include "integrals/h1.h"
#include "boundaryconditions/essential_bcs.h"

class WeakFormNSSimpleLinearization : public WeakForm
{
public:
  WeakFormNSSimpleLinearization(bool Stokes, double Reynolds, double time_step, Solution* x_vel_previous_time, 
                                Solution* y_vel_previous_time) : WeakForm(3), Stokes(Stokes), 
                                Reynolds(Reynolds), time_step(time_step), x_vel_previous_time(x_vel_previous_time), 
                                y_vel_previous_time(y_vel_previous_time) {
    BilinearFormSymVel* sym_form_0 = new BilinearFormSymVel(0, 0, Stokes, Reynolds, time_step);
    add_matrix_form(sym_form_0);
    BilinearFormSymVel* sym_form_1 = new BilinearFormSymVel(1, 1, Stokes, Reynolds, time_step);
    add_matrix_form(sym_form_1);

    BilinearFormNonsymVel* nonsym_vel_form_0 = new BilinearFormNonsymVel(0, 0, Stokes);
    nonsym_vel_form_0->ext = Hermes::vector<MeshFunction*>(x_vel_previous_time, y_vel_previous_time);
    add_matrix_form(nonsym_vel_form_0);
    BilinearFormNonsymVel* nonsym_vel_form_1 = new BilinearFormNonsymVel(1, 1, Stokes);
    nonsym_vel_form_1->ext = Hermes::vector<MeshFunction*>(x_vel_previous_time, y_vel_previous_time);
    add_matrix_form(nonsym_vel_form_1);

    BilinearFormNonsymXVelPressure* nonsym_velx_pressure_form = new BilinearFormNonsymXVelPressure(0, 2);
    add_matrix_form(nonsym_velx_pressure_form);

    BilinearFormNonsymYVelPressure* nonsym_vely_pressure_form = new BilinearFormNonsymYVelPressure(1, 2);
    add_matrix_form(nonsym_vely_pressure_form);
    
    VectorFormVolVel* vector_vel_form_x = new VectorFormVolVel(0, Stokes, time_step);
    
    Hermes::vector<MeshFunction *> ext_vel_x;
    ext_vel_x.push_back(x_vel_previous_time);

    vector_vel_form_x->ext = ext_vel_x;

    VectorFormVolVel* vector_vel_form_y = new VectorFormVolVel(1, Stokes, time_step);

    Hermes::vector<MeshFunction *> ext_vel_y;
    ext_vel_y.push_back(y_vel_previous_time);
    
    vector_vel_form_y->ext = ext_vel_y;
  };

  class BilinearFormSymVel : public WeakForm::MatrixFormVol
  {
  public:
    BilinearFormSymVel(int i, int j, bool Stokes, double Reynolds, double time_step) 
      : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_SYM), Stokes(Stokes), Reynolds(Reynolds), time_step(time_step) {
      adapt_eval = false;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = int_grad_u_grad_v<double, scalar>(n, wt, u, v) / Reynolds;
      if(!Stokes)
        result += int_u_v<double, scalar>(n, wt, u, v) / time_step;
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext)
    {
      Ord result = int_grad_u_grad_v<Ord, Ord>(n, wt, u, v) / Reynolds;
      if(!Stokes)
        result += int_u_v<Ord, Ord>(n, wt, u, v) / time_step;
      return result;
    }
  protected:
    bool Stokes;
    double Reynolds;
    double time_step;
  };


  class BilinearFormNonsymVel : public WeakForm::MatrixFormVol
  {
  public:
    BilinearFormNonsymVel(int i, int j, bool Stokes) 
      : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_NONSYM), Stokes(Stokes) {
        adapt_eval = false;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      if(!Stokes) {
        Func<scalar>* xvel_prev_time = ext->fn[0];
        Func<scalar>* yvel_prev_time = ext->fn[1];
        result = int_w_nabla_u_v<double, scalar>(n, wt, xvel_prev_time, yvel_prev_time, u, v);
      }
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const {
      Ord result = 0;
      if(!Stokes) {
        Func<Ord>* xvel_prev_time = ext->fn[0];
        Func<Ord>* yvel_prev_time = ext->fn[1];
        result = int_w_nabla_u_v<Ord, Ord>(n, wt, xvel_prev_time, yvel_prev_time, u, v);
      }
      return result;
    }
  protected:
    bool Stokes;
  };


  class BilinearFormNonsymXVelPressure : public WeakForm::MatrixFormVol
  {
  public:
    BilinearFormNonsymXVelPressure(int i, int j) : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_ANTISYM) {
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
    BilinearFormNonsymYVelPressure(int i, int j) : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_ANTISYM) {
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


  class VectorFormVolVel : public WeakForm::VectorFormVol
  {
  public:
    VectorFormVolVel(int i, bool Stokes, double time_step) 
          : WeakForm::VectorFormVol(i), Stokes(Stokes), time_step(time_step) {
      adapt_eval = false;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, 
                         ExtData<scalar> *ext) const {
      scalar result = 0;
      if(!Stokes) {
        Func<scalar>* vel_prev_time = ext->fn[0]; // this form is used with both velocity components
        result = int_u_v<double, scalar>(n, wt, vel_prev_time, v) / time_step;
      }
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      Ord result = 0;
      if(!Stokes) {
        Func<Ord>* vel_prev_time = ext->fn[0]; // this form is used with both velocity components
        result = int_u_v<Ord, Ord>(n, wt, vel_prev_time, v) / time_step;
      }
      return result;
    }
  protected:
    bool Stokes;
    double time_step;
  };

protected:
  bool Stokes;
  double Reynolds;
  double time_step;
  Solution* x_vel_previous_time;
  Solution* y_vel_previous_time;
};

class WeakFormNSNewton : public WeakForm
{
public:
  WeakFormNSNewton(bool Stokes, double Reynolds, double time_step, Solution* x_vel_previous_time, 
                   Solution* y_vel_previous_time) : WeakForm(3), Stokes(Stokes), 
                   Reynolds(Reynolds), time_step(time_step), x_vel_previous_time(x_vel_previous_time), 
                   y_vel_previous_time(y_vel_previous_time) {
    BilinearFormSymVel* sym_form_0 = new BilinearFormSymVel(0, 0, Stokes, Reynolds, time_step);
    add_matrix_form(sym_form_0);
    BilinearFormSymVel* sym_form_1 = new BilinearFormSymVel(1, 1, Stokes, Reynolds, time_step);
    add_matrix_form(sym_form_1);

    BilinearFormNonsymVel_0_0* nonsym_vel_form_0_0 = new BilinearFormNonsymVel_0_0(0, 0, Stokes);
    add_matrix_form(nonsym_vel_form_0_0);
    BilinearFormNonsymVel_0_1* nonsym_vel_form_0_1 = new BilinearFormNonsymVel_0_1(0, 1, Stokes);
    add_matrix_form(nonsym_vel_form_0_1);
    BilinearFormNonsymVel_1_0* nonsym_vel_form_1_0 = new BilinearFormNonsymVel_1_0(1, 0, Stokes);
    add_matrix_form(nonsym_vel_form_1_0);
    BilinearFormNonsymVel_1_1* nonsym_vel_form_1_1 = new BilinearFormNonsymVel_1_1(1, 1, Stokes);
    add_matrix_form(nonsym_vel_form_1_1);

    BilinearFormNonsymXVelPressure* nonsym_velx_pressure_form = new BilinearFormNonsymXVelPressure(0, 2);
    add_matrix_form(nonsym_velx_pressure_form);

    BilinearFormNonsymYVelPressure* nonsym_vely_pressure_form = new BilinearFormNonsymYVelPressure(1, 2);
    add_matrix_form(nonsym_vely_pressure_form);
    
    VectorFormNS_0* F_0 = new VectorFormNS_0(0, Stokes, Reynolds, time_step);
    F_0->ext = Hermes::vector<MeshFunction*>(x_vel_previous_time, y_vel_previous_time);
    add_vector_form(F_0);
    VectorFormNS_1* F_1 = new VectorFormNS_1(1, Stokes, Reynolds, time_step);
    F_1->ext = Hermes::vector<MeshFunction*>(x_vel_previous_time, y_vel_previous_time);
    add_vector_form(F_1);
    VectorFormNS_2* F_2 = new VectorFormNS_2(2);
    add_vector_form(F_2);
  };

  class BilinearFormSymVel : public WeakForm::MatrixFormVol
  {
  public:
    BilinearFormSymVel(int i, int j, bool Stokes, double Reynolds, double time_step) 
            : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_SYM), Stokes(Stokes), 
                       Reynolds(Reynolds), time_step(time_step) {
      adapt_eval = false;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = int_grad_u_grad_v<double, scalar>(n, wt, u, v) / Reynolds;
      if(!Stokes)
        result += int_u_v<double, scalar>(n, wt, u, v) / time_step;
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const {
      Ord result = int_grad_u_grad_v<Ord, Ord>(n, wt, u, v) / Reynolds;
      if(!Stokes)
        result += int_u_v<Ord, Ord>(n, wt, u, v) / time_step;
      return result;
    }
  protected:
    bool Stokes;
    double Reynolds;
    double time_step;
  };


  class BilinearFormNonsymVel_0_0 : public WeakForm::MatrixFormVol
  {
  public:
    BilinearFormNonsymVel_0_0(int i, int j, bool Stokes) : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_NONSYM), Stokes(Stokes) {
      adapt_eval = false;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      if(!Stokes) {
        Func<scalar>* xvel_prev_newton = u_ext[0];
        Func<scalar>* yvel_prev_newton = u_ext[1];
        for (int i = 0; i < n; i++)
          result += wt[i] * ((xvel_prev_newton->val[i] * u->dx[i] + yvel_prev_newton->val[i]
                              * u->dy[i]) * v->val[i] + u->val[i] * v->val[i] * xvel_prev_newton->dx[i]);
      }
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const {
      Ord result = 0;
      if(!Stokes) {
        Func<Ord>* xvel_prev_newton = u_ext[0];
        Func<Ord>* yvel_prev_newton = u_ext[1];
        for (int i = 0; i < n; i++)
          result += wt[i] * ((xvel_prev_newton->val[i] * u->dx[i] + yvel_prev_newton->val[i]
                              * u->dy[i]) * v->val[i] + u->val[i] * v->val[i] * xvel_prev_newton->dx[i]);
      }
      return result;
    }
  protected:
    bool Stokes;
  };


  class BilinearFormNonsymVel_0_1 : public WeakForm::MatrixFormVol
  {
  public:
    BilinearFormNonsymVel_0_1(int i, int j, bool Stokes) : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_NONSYM), Stokes(Stokes) {
      adapt_eval = false;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      if(!Stokes) {
        Func<scalar>* xvel_prev_newton = u_ext[0];
        for (int i = 0; i < n; i++)
          result += wt[i] * (u->val[i] * v->val[i] * xvel_prev_newton->dy[i]);
      }
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const {
      Ord result = 0;
      if(!Stokes) {
        Func<Ord>* xvel_prev_newton = u_ext[0];
        for (int i = 0; i < n; i++)
          result += wt[i] * (u->val[i] * v->val[i] * xvel_prev_newton->dy[i]);
      }
      return result;
    }
  protected:
    bool Stokes;
  };


  class BilinearFormNonsymVel_1_0 : public WeakForm::MatrixFormVol
  {
  public:
    BilinearFormNonsymVel_1_0(int i, int j, bool Stokes) : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_NONSYM), Stokes(Stokes) {
      adapt_eval = false;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      if(!Stokes) {
        Func<scalar>* yvel_prev_newton = u_ext[0];
        for (int i = 0; i < n; i++)
          result += wt[i] * (u->val[i] * v->val[i] * yvel_prev_newton->dx[i]);
      }
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const {
      Ord result = 0;
      if(!Stokes) {
        Func<Ord>* yvel_prev_newton = u_ext[0];
        for (int i = 0; i < n; i++)
          result += wt[i] * (u->val[i] * v->val[i] * yvel_prev_newton->dx[i]);
      }
      return result;
    }
  protected:
    bool Stokes;
  };


  class BilinearFormNonsymVel_1_1 : public WeakForm::MatrixFormVol
  {
  public:
    BilinearFormNonsymVel_1_1(int i, int j, bool Stokes) : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_NONSYM), Stokes(Stokes) {
      adapt_eval = false;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      if(!Stokes) {
        Func<scalar>* xvel_prev_newton = u_ext[0];
        Func<scalar>* yvel_prev_newton = u_ext[1];
        for (int i = 0; i < n; i++)
          result += wt[i] * ((xvel_prev_newton->val[i] * u->dx[i] + yvel_prev_newton->val[i] * u->dy[i]) * v->val[i] + u->val[i]
                             * v->val[i] * yvel_prev_newton->dy[i]);
      }
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const {
      Ord result = 0;
      if(!Stokes) {
        Func<Ord>* xvel_prev_newton = u_ext[0];
        Func<Ord>* yvel_prev_newton = u_ext[1];
        for (int i = 0; i < n; i++)
          result += wt[i] * ((xvel_prev_newton->val[i] * u->dx[i] + yvel_prev_newton->val[i] * u->dy[i]) * v->val[i] + u->val[i]
                             * v->val[i] * yvel_prev_newton->dy[i]);
      }
      return result;
    }
  protected:
    bool Stokes;
  };


  class BilinearFormNonsymXVelPressure : public WeakForm::MatrixFormVol
  {
  public:
    BilinearFormNonsymXVelPressure(int i, int j) : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_ANTISYM) {
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
    BilinearFormNonsymYVelPressure(int i, int j) : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_ANTISYM) {
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
    VectorFormNS_0(int i, bool Stokes, double Reynolds, double time_step) : WeakForm::VectorFormVol(i), 
                   Stokes(Stokes), Reynolds(Reynolds), time_step(time_step) {
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
        result += wt[i] * ((xvel_prev_newton->dx[i] * v->dx[i] + xvel_prev_newton->dy[i] * v->dy[i]) / Reynolds 
                          - (p_prev_newton->val[i] * v->dx[i]));
      if(!Stokes)
        for (int i = 0; i < n; i++)
          result += wt[i] * (((xvel_prev_newton->val[i] - xvel_prev_time->val[i]) * v->val[i] / time_step )
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
        result += wt[i] * ((xvel_prev_newton->dx[i] * v->dx[i] + xvel_prev_newton->dy[i] * v->dy[i]) / Reynolds 
                            - (p_prev_newton->val[i] * v->dx[i]));
      if(!Stokes)
        for (int i = 0; i < n; i++)
          result += wt[i] * (((xvel_prev_newton->val[i] - xvel_prev_time->val[i]) * v->val[i] / time_step)
                            + ((xvel_prev_newton->val[i] * xvel_prev_newton->dx[i] 
                            + yvel_prev_newton->val[i] * xvel_prev_newton->dy[i]) * v->val[i]));
      return result;
    }
  protected:
    bool Stokes;
    double Reynolds;
    double time_step;
  };

  class VectorFormNS_1 : public WeakForm::VectorFormVol
  {
  public:
    VectorFormNS_1(int i, bool Stokes, double Reynolds, double time_step) 
        : WeakForm::VectorFormVol(i), Stokes(Stokes), Reynolds(Reynolds), time_step(time_step) {
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
        result += wt[i] * ((yvel_prev_newton->dx[i] * v->dx[i] + yvel_prev_newton->dy[i] * v->dy[i]) / Reynolds 
                          - (p_prev_newton->val[i] * v->dy[i]));
      if(!Stokes)
        for (int i = 0; i < n; i++)
          result += wt[i] * (((yvel_prev_newton->val[i] - yvel_prev_time->val[i]) * v->val[i] / time_step )
                            + ((xvel_prev_newton->val[i] * yvel_prev_newton->dx[i] 
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
        result += wt[i] * ((xvel_prev_newton->dx[i] * v->dx[i] + xvel_prev_newton->dy[i] * v->dy[i]) / Reynolds 
                  - (p_prev_newton->val[i] * v->dx[i]));
      if(!Stokes)
        for (int i = 0; i < n; i++)
          result += wt[i] * (((xvel_prev_newton->val[i] - xvel_prev_time->val[i]) * v->val[i] / time_step )
                            + ((xvel_prev_newton->val[i] * xvel_prev_newton->dx[i] 
                            + yvel_prev_newton->val[i] * xvel_prev_newton->dy[i]) * v->val[i]));
      return result;
    }
  protected:
    bool Stokes;
    double Reynolds;
    double time_step;
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
  bool Stokes;
  double Reynolds;
  double time_step;
  Solution* x_vel_previous_time;
  Solution* y_vel_previous_time;
};

class EssentialBCNonConst : public EssentialBoundaryCondition
{
public:
  EssentialBCNonConst(Hermes::vector<std::string> markers, double vel_inlet, double H, double startup_time) : 
    EssentialBoundaryCondition(markers), startup_time(startup_time), vel_inlet(vel_inlet), H(H)  {};
  EssentialBCNonConst(std::string marker, double vel_inlet, double H, double startup_time) : 
    EssentialBoundaryCondition(Hermes::vector<std::string>()), startup_time(startup_time), vel_inlet(vel_inlet), H(H)  {
    markers.push_back(marker);
  };
  
  ~EssentialBCNonConst() {};

  virtual EssentialBCValueType get_value_type() const { 
    return BC_FUNCTION; 
  };

  virtual scalar value(double x, double y, double n_x, double n_y, double t_x, double t_y) const {
    double val_y = vel_inlet * y*(H-y) / (H/2.)/(H/2.);
    if (current_time <= startup_time) 
      return val_y * current_time/startup_time;
    else 
      return val_y;
  };

protected:
  double startup_time;
  double vel_inlet;
  double H;
};
