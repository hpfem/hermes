#include "hermes2d.h"

using namespace WeakFormsH1;

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

    // Pressure term in the first velocity equation.
    add_matrix_form(new BilinearFormNonsymXVelPressure(0, 2));
    // Pressure term in the second velocity equation.
    add_matrix_form(new BilinearFormNonsymYVelPressure(1, 2));
    
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
      : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_SYM), Stokes(Stokes), Reynolds(Reynolds), time_step(time_step) { }

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
      : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_NONSYM), Stokes(Stokes) { }

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
    BilinearFormNonsymXVelPressure(int i, int j) : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_ANTISYM) { }

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
    BilinearFormNonsymYVelPressure(int i, int j) : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_ANTISYM) { }

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
          : WeakForm::VectorFormVol(i), Stokes(Stokes), time_step(time_step) { }

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
    /* Jacobian terms - first velocity equation */
    // Time derivative in the first velocity equation
    // and Laplacian divided by Re in the first velocity equation.
    add_matrix_form(new BilinearFormSymVel(0, 0, Stokes, Reynolds, time_step));
    // First part of the convective term in the first velocity equation.
    add_matrix_form(new BilinearFormNonsymVel_0_0(0, 0, Stokes));
    // Second part of the convective term in the first velocity equation.
    add_matrix_form(new BilinearFormNonsymVel_0_1(0, 1, Stokes));
    // Pressure term in the first velocity equation.
    add_matrix_form(new BilinearFormNonsymXVelPressure(0, 2));

    /* Jacobian terms - second velocity equation, continuity equation */
    // Time derivative in the second velocity equation
    // and Laplacian divided by Re in the second velocity equation.
    add_matrix_form(new BilinearFormSymVel(1, 1, Stokes, Reynolds, time_step));
    // First part of the convective term in the second velocity equation.
    add_matrix_form(new BilinearFormNonsymVel_1_0(1, 0, Stokes));
    // Second part of the convective term in the second velocity equation.
    add_matrix_form(new BilinearFormNonsymVel_1_1(1, 1, Stokes));
    // Pressure term in the second velocity equation.
    add_matrix_form(new BilinearFormNonsymYVelPressure(1, 2));

    
    /* Residual - volumetric */
    // First velocity equation.
    VectorFormNS_0* F_0 = new VectorFormNS_0(0, Stokes, Reynolds, time_step);
    F_0->ext.push_back(x_vel_previous_time);
    add_vector_form(F_0);
    // Second velocity equation.
    VectorFormNS_1* F_1 = new VectorFormNS_1(1, Stokes, Reynolds, time_step);
    F_1->ext.push_back(y_vel_previous_time);
    add_vector_form(F_1);
    // Continuity equation.
    VectorFormNS_2* F_2 = new VectorFormNS_2(2);
    add_vector_form(F_2);
  };

  class BilinearFormSymVel : public WeakForm::MatrixFormVol
  {
  public:
    BilinearFormSymVel(int i, int j, bool Stokes, double Reynolds, double time_step) 
            : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_SYM), Stokes(Stokes), 
                        Reynolds(Reynolds), time_step(time_step) { }

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
    BilinearFormNonsymVel_0_0(int i, int j, bool Stokes) 
      : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_NONSYM), Stokes(Stokes) { }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      if(!Stokes) {
        Func<scalar>* xvel_prev_newton = u_ext[0];
        Func<scalar>* yvel_prev_newton = u_ext[1];
        for (int i = 0; i < n; i++)
          result += wt[i] * (xvel_prev_newton->val[i] * u->dx[i] + yvel_prev_newton->val[i] * u->dy[i] 
                             + u->val[i] * xvel_prev_newton->dx[i]) * v->val[i];
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
    BilinearFormNonsymVel_0_1(int i, int j, bool Stokes) 
      : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_NONSYM), Stokes(Stokes) { }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      if(!Stokes) {
        Func<scalar>* xvel_prev_newton = u_ext[0];
        for (int i = 0; i < n; i++)
          result += wt[i] * u->val[i] * xvel_prev_newton->dy[i] * v->val[i];
      }
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const {
      Ord result = 0;
      if(!Stokes) {
        Func<Ord>* xvel_prev_newton = u_ext[0];
        for (int i = 0; i < n; i++)
          result += wt[i] * u->val[i] * xvel_prev_newton->dy[i] * v->val[i] ;
      }
      return result;
    }
  protected:
    bool Stokes;
  };
  class BilinearFormNonsymVel_1_0 : public WeakForm::MatrixFormVol
  {
  public:
    BilinearFormNonsymVel_1_0(int i, int j, bool Stokes) 
      : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_NONSYM), Stokes(Stokes) { }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      if(!Stokes) {
        Func<scalar>* yvel_prev_newton = u_ext[1];
        for (int i = 0; i < n; i++)
          result += wt[i] * u->val[i] * yvel_prev_newton->dx[i] * v->val[i];
      }
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const {
      Ord result = 0;
      if(!Stokes) {
        Func<Ord>* yvel_prev_newton = u_ext[1];
        for (int i = 0; i < n; i++)
          result += wt[i] * u->val[i] * yvel_prev_newton->dx[i] * v->val[i];
      }
      return result;
    }
  protected:
    bool Stokes;
  };
  class BilinearFormNonsymVel_1_1 : public WeakForm::MatrixFormVol
  {
  public:
    BilinearFormNonsymVel_1_1(int i, int j, bool Stokes) 
      : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_NONSYM), Stokes(Stokes) { }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      if(!Stokes) {
        Func<scalar>* xvel_prev_newton = u_ext[0];
        Func<scalar>* yvel_prev_newton = u_ext[1];
        for (int i = 0; i < n; i++)
          result += wt[i] * (  xvel_prev_newton->val[i] * u->dx[i] 
                             + yvel_prev_newton->val[i] * u->dy[i] 
                             + u->val[i] * yvel_prev_newton->dy[i]) * v->val[i];
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
          result += wt[i] * (  xvel_prev_newton->val[i] * u->dx[i] 
                             + yvel_prev_newton->val[i] * u->dy[i] 
                             + u->val[i] * yvel_prev_newton->dy[i]) * v->val[i];
      }
      return result;
    }
  protected:
    bool Stokes;
  };
  class BilinearFormNonsymXVelPressure : public WeakForm::MatrixFormVol
  {
  public:
    BilinearFormNonsymXVelPressure(int i, int j) : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_ANTISYM) { }

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
    BilinearFormNonsymYVelPressure(int i, int j) : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_ANTISYM) { }

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
                   Stokes(Stokes), Reynolds(Reynolds), time_step(time_step) { }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, 
                         ExtData<scalar> *ext) const {
      scalar result = 0;
      Func<scalar>* xvel_prev_time = ext->fn[0];  
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
        : WeakForm::VectorFormVol(i), Stokes(Stokes), Reynolds(Reynolds), time_step(time_step) { }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, 
                         ExtData<scalar> *ext) const {
      scalar result = 0;
      Func<scalar>* yvel_prev_time = ext->fn[0];
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
      Func<Ord>* yvel_prev_time = ext->fn[0];
      Func<Ord>* xvel_prev_newton = u_ext[0];  
      Func<Ord>* yvel_prev_newton = u_ext[1];  
      Func<Ord>* p_prev_newton = u_ext[2];
      for (int i = 0; i < n; i++)
        result += wt[i] * ((xvel_prev_newton->dx[i] * v->dx[i] + xvel_prev_newton->dy[i] * v->dy[i]) / Reynolds 
                  - (p_prev_newton->val[i] * v->dx[i]));
      if(!Stokes)
        for (int i = 0; i < n; i++)
          result += wt[i] * (((yvel_prev_newton->val[i] - yvel_prev_time->val[i]) * v->val[i] / time_step )
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
    VectorFormNS_2(int i) : WeakForm::VectorFormVol(i) { }

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

/* Essential boundary conditions */

// Time-dependent surface x-velocity of inner circle.
class EssentialBCNonConstX : public EssentialBoundaryCondition
{
public:
  EssentialBCNonConstX(Hermes::vector<std::string> markers, double vel, double startup_time) : 
    EssentialBoundaryCondition(markers), startup_time(startup_time), vel(vel)  {};
  EssentialBCNonConstX(std::string marker, double vel, double startup_time) : 
    EssentialBoundaryCondition(Hermes::vector<std::string>()), startup_time(startup_time), vel(vel)  {
    markers.push_back(marker);
  };
  
  ~EssentialBCNonConstX() {};

  virtual EssentialBCValueType get_value_type() const { 
    return BC_FUNCTION; 
  };

  virtual scalar value(double x, double y, double n_x, double n_y, double t_x, double t_y) const {
    double velocity;
    if (current_time <= startup_time) velocity = vel * current_time/startup_time;
    else velocity = vel;
    double alpha = atan2(x, y);
    double xvel = velocity*cos(alpha);
    return xvel; 
  };

  protected:
    double startup_time;
    double vel;
};

// Time-dependent surface y-velocity of inner circle.
class EssentialBCNonConstY : public EssentialBoundaryCondition
{
public:
  EssentialBCNonConstY(Hermes::vector<std::string> markers, double vel, double startup_time) : 
    EssentialBoundaryCondition(markers), startup_time(startup_time), vel(vel)  {};
  EssentialBCNonConstY(std::string marker, double vel, double startup_time) : 
    EssentialBoundaryCondition(Hermes::vector<std::string>()), startup_time(startup_time), vel(vel)  {
    markers.push_back(marker);
  };
  
  ~EssentialBCNonConstY() {};

  virtual EssentialBCValueType get_value_type() const { 
    return BC_FUNCTION; 
  };

  virtual scalar value(double x, double y, double n_x, double n_y, double t_x, double t_y) const {
    double velocity;
    if (current_time <= startup_time) velocity = vel * current_time/startup_time;
    else velocity = vel;
    double alpha = atan2(x, y);
    double yvel = -velocity*sin(alpha);
    return yvel; 
  };

  protected:
    double startup_time;
    double vel;
};
