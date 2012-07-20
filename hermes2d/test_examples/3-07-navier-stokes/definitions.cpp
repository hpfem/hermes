#include "hermes2d.h"

class WeakFormNSSimpleLinearization : public WeakForm<double>
{
public:
  WeakFormNSSimpleLinearization(bool Stokes, double Reynolds, double time_step, Solution<double>* x_vel_previous_time, Solution<double>* y_vel_previous_time) : WeakForm<double>(3), Stokes(Stokes), 
    Reynolds(Reynolds), time_step(time_step), x_vel_previous_time(x_vel_previous_time), y_vel_previous_time(y_vel_previous_time) {
      BilinearFormSymVel* sym_form_0 = new BilinearFormSymVel(0, 0, Stokes, Reynolds, time_step);
      add_matrix_form(sym_form_0);
      BilinearFormSymVel* sym_form_1 = new BilinearFormSymVel(1, 1, Stokes, Reynolds, time_step);
      add_matrix_form(sym_form_1);

      BilinearFormUnSymVel* unsym_vel_form_0 = new BilinearFormUnSymVel(0, 0, Stokes);
      unsym_vel_form_0->ext = Hermes::vector<MeshFunction<double>*>(x_vel_previous_time, y_vel_previous_time);
      add_matrix_form(unsym_vel_form_0);
      BilinearFormUnSymVel* unsym_vel_form_1 = new BilinearFormUnSymVel(1, 1, Stokes);
      unsym_vel_form_1->ext = Hermes::vector<MeshFunction<double>*>(x_vel_previous_time, y_vel_previous_time);
      add_matrix_form(unsym_vel_form_1);

      BilinearFormUnSymXVelPressure* unsym_velx_pressure_form = new BilinearFormUnSymXVelPressure(0, 2);
      add_matrix_form(unsym_velx_pressure_form);

      BilinearFormUnSymYVelPressure* unsym_vely_pressure_form = new BilinearFormUnSymYVelPressure(1, 2);
      add_matrix_form(unsym_vely_pressure_form);

      VectorFormVolVel* vector_vel_form_x = new VectorFormVolVel(0, Stokes, time_step);

      Hermes::vector<MeshFunction<double> *> ext_vel_x;
      ext_vel_x.push_back(x_vel_previous_time);

      vector_vel_form_x->ext = ext_vel_x;

      VectorFormVolVel* vector_vel_form_y = new VectorFormVolVel(1, Stokes, time_step);

      Hermes::vector<MeshFunction<double> *> ext_vel_y;
      ext_vel_y.push_back(y_vel_previous_time);

      vector_vel_form_y->ext = ext_vel_y;
  };

  class BilinearFormSymVel : public MatrixFormVol<double>
  {
  public:
    BilinearFormSymVel(int i, int j, bool Stokes, double Reynolds, double time_step) : MatrixFormVol<double>(i, j), Stokes(Stokes), 
      Reynolds(Reynolds), time_step(time_step) {
        sym = HERMES_SYM;
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const{
      double result = int_grad_u_grad_v<double, double>(n, wt, u, v) / Reynolds;
      if(!Stokes)
        result += int_u_v<double, double>(n, wt, u, v) / time_step;
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const
    {
      Ord result = int_grad_u_grad_v<Ord, Ord>(n, wt, u, v) / Reynolds;
      if(!Stokes)
        result += int_u_v<Ord, Ord>(n, wt, u, v) / time_step;
      return result;
    }

    MatrixFormVol<double>* clone()
    {
      return new BilinearFormSymVel(this->i, this->j, this->Stokes, this->Reynolds, this->time_step);
    }
  protected:
    // Members.
    bool Stokes;
    double Reynolds;
    double time_step;
  };


  class BilinearFormUnSymVel : public MatrixFormVol<double>
  {
  public:
    BilinearFormUnSymVel(int i, int j, bool Stokes) : MatrixFormVol<double>(i, j), Stokes(Stokes) {
      sym = HERMES_NONSYM;
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const{
      double result = 0;
      if(!Stokes) {
        Func<double>* xvel_prev_time = ext->fn[0];
        Func<double>* yvel_prev_time = ext->fn[1];
        result = int_w_nabla_u_v<double, double>(n, wt, xvel_prev_time, yvel_prev_time, u, v);
      }
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      Ord result = Ord(0);
      if(!Stokes) {
        Func<Ord>* xvel_prev_time = ext->fn[0];
        Func<Ord>* yvel_prev_time = ext->fn[1];
        result = int_w_nabla_u_v<Ord, Ord>(n, wt, xvel_prev_time, yvel_prev_time, u, v);
      }
      return result;
    }

    MatrixFormVol<double>* clone()
    {
      return new BilinearFormUnSymVel(this->i, this->j, this->Stokes);
    }
  protected:
    // Members.
    bool Stokes;
  };


  class BilinearFormUnSymXVelPressure : public MatrixFormVol<double>
  {
  public:
    BilinearFormUnSymXVelPressure(int i, int j) : MatrixFormVol<double>(i, j) {
      sym = HERMES_ANTISYM;
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const{
      return - int_u_dvdx<double, double>(n, wt, u, v);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return - int_u_dvdx<Ord, Ord>(n, wt, u, v);
    }
    
    MatrixFormVol<double>* clone()
    {
      return new BilinearFormUnSymXVelPressure(this->i, this->j);
    }
  };


  class BilinearFormUnSymYVelPressure : public MatrixFormVol<double>
  {
  public:
    BilinearFormUnSymYVelPressure(int i, int j) : MatrixFormVol<double>(i, j) {
      sym = HERMES_ANTISYM;
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const{
      return - int_u_dvdy<double, double>(n, wt, u, v);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const{
      return - int_u_dvdy<Ord, Ord>(n, wt, u, v);
    }

    MatrixFormVol<double>* clone()
    {
      return new BilinearFormUnSymYVelPressure(this->i, this->j);
    }
  };


  class VectorFormVolVel : public VectorFormVol<double>
  {
  public:
    VectorFormVolVel(int i, bool Stokes, double time_step) : VectorFormVol<double>(i), Stokes(Stokes), time_step(time_step) {
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext) const{
      double result = 0;
      if(!Stokes) {
        Func<double>* vel_prev_time = ext->fn[0]; // this form is used with both velocity components
        result = int_u_v<double, double>(n, wt, vel_prev_time, v) / time_step;
      }
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      Ord result = Ord(0);
      if(!Stokes) {
        Func<Ord>* vel_prev_time = ext->fn[0]; // this form is used with both velocity components
        result = int_u_v<Ord, Ord>(n, wt, vel_prev_time, v) / time_step;
      }
      return result;
    }
    
    VectorFormVol<double>* clone()
    {
      return new VectorFormVolVel(this->i, this->Stokes, this->time_step);
    }
  protected:
    // Members.
    bool Stokes;
    double time_step;
  };

protected:
  // Members.
  bool Stokes;
  double Reynolds;
  double time_step;
  Solution<double>* x_vel_previous_time;
  Solution<double>* y_vel_previous_time;
};

class WeakFormNSNewton : public WeakForm<double>
{
public:
  WeakFormNSNewton(bool Stokes, double Reynolds, double time_step, Solution<double>* x_vel_previous_time, Solution<double>* y_vel_previous_time) : WeakForm<double>(3), Stokes(Stokes), 
    Reynolds(Reynolds), x_vel_previous_time(x_vel_previous_time), y_vel_previous_time(y_vel_previous_time) 
  {
    this->current_time_step = time_step;

    BilinearFormSymVel* sym_form_0 = new BilinearFormSymVel(0, 0, Stokes, Reynolds, time_step);
    add_matrix_form(sym_form_0);
    BilinearFormSymVel* sym_form_1 = new BilinearFormSymVel(1, 1, Stokes, Reynolds, time_step);
    add_matrix_form(sym_form_1);

    BilinearFormUnSymVel_0_0* unsym_vel_form_0_0 = new BilinearFormUnSymVel_0_0(0, 0, Stokes);
    add_matrix_form(unsym_vel_form_0_0);
    BilinearFormUnSymVel_0_1* unsym_vel_form_0_1 = new BilinearFormUnSymVel_0_1(0, 1, Stokes);
    add_matrix_form(unsym_vel_form_0_1);
    BilinearFormUnSymVel_1_0* unsym_vel_form_1_0 = new BilinearFormUnSymVel_1_0(1, 0, Stokes);
    add_matrix_form(unsym_vel_form_1_0);
    BilinearFormUnSymVel_1_1* unsym_vel_form_1_1 = new BilinearFormUnSymVel_1_1(1, 1, Stokes);
    add_matrix_form(unsym_vel_form_1_1);

    BilinearFormUnSymXVelPressure* unsym_velx_pressure_form = new BilinearFormUnSymXVelPressure(0, 2);
    add_matrix_form(unsym_velx_pressure_form);

    BilinearFormUnSymYVelPressure* unsym_vely_pressure_form = new BilinearFormUnSymYVelPressure(1, 2);
    add_matrix_form(unsym_vely_pressure_form);

    VectorFormNS_0* F_0 = new VectorFormNS_0(0, Stokes, Reynolds, time_step);
    F_0->ext = Hermes::vector<MeshFunction<double>*>(x_vel_previous_time, y_vel_previous_time);
    add_vector_form(F_0);
    VectorFormNS_1* F_1 = new VectorFormNS_1(1, Stokes, Reynolds, time_step);
    F_1->ext = Hermes::vector<MeshFunction<double>*>(x_vel_previous_time, y_vel_previous_time);
    add_vector_form(F_1);
    VectorFormNS_2* F_2 = new VectorFormNS_2(2);
    add_vector_form(F_2);
  };

  class BilinearFormSymVel : public MatrixFormVol<double>
  {
  public:
    BilinearFormSymVel(int i, int j, bool Stokes, double Reynolds, double time_step) : MatrixFormVol<double>(i, j), Stokes(Stokes), 
      Reynolds(Reynolds), time_step(time_step) {
        sym = HERMES_SYM;
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const{
      double result = int_grad_u_grad_v<double, double>(n, wt, u, v) / Reynolds;
      if(!Stokes)
        result += int_u_v<double, double>(n, wt, u, v) / time_step;
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const
    {
      Ord result = int_grad_u_grad_v<Ord, Ord>(n, wt, u, v) / Reynolds;
      if(!Stokes)
        result += int_u_v<Ord, Ord>(n, wt, u, v) / time_step;
      return result;
    }

    MatrixFormVol<double>* clone()
    {
      return new BilinearFormSymVel(this->i, this->j, this->Stokes, this->Reynolds, this->time_step);
    }
  protected:
    // Members.
    bool Stokes;
    double Reynolds;
    double time_step;
  };


  class BilinearFormUnSymVel_0_0 : public MatrixFormVol<double>
  {
  public:
    BilinearFormUnSymVel_0_0(int i, int j, bool Stokes) : MatrixFormVol<double>(i, j), Stokes(Stokes) {
      sym = HERMES_NONSYM;
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const{
      double result = 0;
      if(!Stokes) {
        Func<double>* xvel_prev_newton = u_ext[0];
        Func<double>* yvel_prev_newton = u_ext[1];
        for (int i = 0; i < n; i++)
          result += wt[i] * ((xvel_prev_newton->val[i] * u->dx[i] + yvel_prev_newton->val[i]
        * u->dy[i]) * v->val[i] + u->val[i] * v->val[i] * xvel_prev_newton->dx[i]);
      }
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const{
      Ord result = Ord(0);
      if(!Stokes) {
        Func<Ord>* xvel_prev_newton = u_ext[0];
        Func<Ord>* yvel_prev_newton = u_ext[1];
        for (int i = 0; i < n; i++)
          result += wt[i] * ((xvel_prev_newton->val[i] * u->dx[i] + yvel_prev_newton->val[i]
        * u->dy[i]) * v->val[i] + u->val[i] * v->val[i] * xvel_prev_newton->dx[i]);
      }
      return result;
    }

    MatrixFormVol<double>* clone()
    {
      return new BilinearFormUnSymVel_0_0(this->i, this->j, this->Stokes);
    }
  protected:
    // Members.
    bool Stokes;
  };


  class BilinearFormUnSymVel_0_1 : public MatrixFormVol<double>
  {
  public:
    BilinearFormUnSymVel_0_1(int i, int j, bool Stokes) : MatrixFormVol<double>(i, j), Stokes(Stokes) {
      sym = HERMES_NONSYM;
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const{
      double result = 0;
      if(!Stokes) {
        Func<double>* xvel_prev_newton = u_ext[0];
        for (int i = 0; i < n; i++)
          result += wt[i] * (u->val[i] * v->val[i] * xvel_prev_newton->dy[i]);
      }
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      Ord result = Ord(0);
      if(!Stokes) {
        Func<Ord>* xvel_prev_newton = u_ext[0];
        for (int i = 0; i < n; i++)
          result += wt[i] * (u->val[i] * v->val[i] * xvel_prev_newton->dy[i]);
      }
      return result;
    }
    
    MatrixFormVol<double>* clone()
    {
      return new BilinearFormUnSymVel_0_1(this->i, this->j, this->Stokes);
    }
  protected:
    // Members.
    bool Stokes;
  };


  class BilinearFormUnSymVel_1_0 : public MatrixFormVol<double>
  {
  public:
    BilinearFormUnSymVel_1_0(int i, int j, bool Stokes) : MatrixFormVol<double>(i, j), Stokes(Stokes) {
      sym = HERMES_NONSYM;
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const{
      double result = 0;
      if(!Stokes) {
        Func<double>* yvel_prev_newton = u_ext[1];
        for (int i = 0; i < n; i++)
          result += wt[i] * (u->val[i] * v->val[i] * yvel_prev_newton->dx[i]);
      }
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const{
      Ord result = Ord(0);
      if(!Stokes) {
        Func<Ord>* yvel_prev_newton = u_ext[1];
        for (int i = 0; i < n; i++)
          result += wt[i] * (u->val[i] * v->val[i] * yvel_prev_newton->dx[i]);
      }
      return result;
    }
    
    MatrixFormVol<double>* clone()
    {
      return new BilinearFormUnSymVel_1_0(this->i, this->j, this->Stokes);
    }
  protected:
    // Members.
    bool Stokes;
  };


  class BilinearFormUnSymVel_1_1 : public MatrixFormVol<double>
  {
  public:
    BilinearFormUnSymVel_1_1(int i, int j, bool Stokes) : MatrixFormVol<double>(i, j), Stokes(Stokes) {
      sym = HERMES_NONSYM;
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const{
      double result = 0;
      if(!Stokes) {
        Func<double>* xvel_prev_newton = u_ext[0];
        Func<double>* yvel_prev_newton = u_ext[1];
        for (int i = 0; i < n; i++)
          result += wt[i] * ((xvel_prev_newton->val[i] * u->dx[i] + yvel_prev_newton->val[i] * u->dy[i]) * v->val[i] + u->val[i]
        * v->val[i] * yvel_prev_newton->dy[i]);
      }
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      Ord result = Ord(0);
      if(!Stokes) {
        Func<Ord>* xvel_prev_newton = u_ext[0];
        Func<Ord>* yvel_prev_newton = u_ext[1];
        for (int i = 0; i < n; i++)
          result += wt[i] * ((xvel_prev_newton->val[i] * u->dx[i] + yvel_prev_newton->val[i] * u->dy[i]) * v->val[i] + u->val[i]
        * v->val[i] * yvel_prev_newton->dy[i]);
      }
      return result;
    }
    
    MatrixFormVol<double>* clone()
    {
      return new BilinearFormUnSymVel_1_1(this->i, this->j, this->Stokes);
    }
  protected:
    // Members.
    bool Stokes;
  };


  class BilinearFormUnSymXVelPressure : public MatrixFormVol<double>
  {
  public:
    BilinearFormUnSymXVelPressure(int i, int j) : MatrixFormVol<double>(i, j) {
      sym = HERMES_ANTISYM;
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const{
      return - int_u_dvdx<double, double>(n, wt, u, v);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return - int_u_dvdx<Ord, Ord>(n, wt, u, v);
    }
    
    MatrixFormVol<double>* clone()
    {
      return new BilinearFormUnSymXVelPressure(this->i, this->j);
    }
  };


  class BilinearFormUnSymYVelPressure : public MatrixFormVol<double>
  {
  public:
    BilinearFormUnSymYVelPressure(int i, int j) : MatrixFormVol<double>(i, j) {
      sym = HERMES_ANTISYM;
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const{
      return - int_u_dvdy<double, double>(n, wt, u, v);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return - int_u_dvdy<Ord, Ord>(n, wt, u, v);
    }
    
    MatrixFormVol<double>* clone()
    {
      return new BilinearFormUnSymYVelPressure(this->i, this->j);
    }
  };


  class VectorFormNS_0 : public VectorFormVol<double>
  {
  public:
    VectorFormNS_0(int i, bool Stokes, double Reynolds, double time_step) : VectorFormVol<double>(i), Stokes(Stokes), Reynolds(Reynolds), time_step(time_step) {
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext) const{
      double result = 0;
      Func<double>* xvel_prev_time = ext->fn[0];  
      Func<double>* yvel_prev_time = ext->fn[1];
      Func<double>* xvel_prev_newton = u_ext[0];  
      Func<double>* yvel_prev_newton = u_ext[1];  
      Func<double>* p_prev_newton = u_ext[2];
      for (int i = 0; i < n; i++)
        result += wt[i] * ((xvel_prev_newton->dx[i] * v->dx[i] + xvel_prev_newton->dy[i] * v->dy[i]) / Reynolds - (p_prev_newton->val[i] * v->dx[i]));
      if(!Stokes)
        for (int i = 0; i < n; i++)
          result += wt[i] * (((xvel_prev_newton->val[i] - xvel_prev_time->val[i]) * v->val[i] / time_step )
          + ((xvel_prev_newton->val[i] * xvel_prev_newton->dx[i] + yvel_prev_newton->val[i] * xvel_prev_newton->dy[i]) * v->val[i]));
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      Ord result = Ord(0);
      Func<Ord>* xvel_prev_time = ext->fn[0];  
      Func<Ord>* yvel_prev_time = ext->fn[1];
      Func<Ord>* xvel_prev_newton = u_ext[0];  
      Func<Ord>* yvel_prev_newton = u_ext[1];  
      Func<Ord>* p_prev_newton = u_ext[2];
      for (int i = 0; i < n; i++)
        result += wt[i] * ((xvel_prev_newton->dx[i] * v->dx[i] + xvel_prev_newton->dy[i] * v->dy[i]) / Reynolds - (p_prev_newton->val[i] * v->dx[i]));
      if(!Stokes)
        for (int i = 0; i < n; i++)
          result += wt[i] * (((xvel_prev_newton->val[i] - xvel_prev_time->val[i]) * v->val[i] / time_step)
          + ((xvel_prev_newton->val[i] * xvel_prev_newton->dx[i] + yvel_prev_newton->val[i] * xvel_prev_newton->dy[i]) * v->val[i]));
      return result;
    }
    
    VectorFormVol<double>* clone()
    {
      return new VectorFormNS_0(this->i, this->Stokes, this->Reynolds, this->time_step);
    }
  protected:
    // Members.
    bool Stokes;
    double Reynolds;
    double time_step;
  };

  class VectorFormNS_1 : public VectorFormVol<double>
  {
  public:
    VectorFormNS_1(int i, bool Stokes, double Reynolds, double time_step) : VectorFormVol<double>(i), Stokes(Stokes), Reynolds(Reynolds), time_step(time_step) {
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext) const{
      double result = 0;
      Func<double>* xvel_prev_time = ext->fn[0];  
      Func<double>* yvel_prev_time = ext->fn[1];
      Func<double>* xvel_prev_newton = u_ext[0];  
      Func<double>* yvel_prev_newton = u_ext[1];  
      Func<double>* p_prev_newton = u_ext[2];
      for (int i = 0; i < n; i++)
        result += wt[i] * ((yvel_prev_newton->dx[i] * v->dx[i] + yvel_prev_newton->dy[i] * v->dy[i]) / Reynolds - (p_prev_newton->val[i] * v->dy[i]));
      if(!Stokes)
        for (int i = 0; i < n; i++)
          result += wt[i] * (((yvel_prev_newton->val[i] - yvel_prev_time->val[i]) * v->val[i] / time_step )
          + ((xvel_prev_newton->val[i] * yvel_prev_newton->dx[i] + yvel_prev_newton->val[i] * yvel_prev_newton->dy[i]) * v->val[i]));
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const{
      Ord result = Ord(0);
      Func<Ord>* xvel_prev_time = ext->fn[0];  
      Func<Ord>* yvel_prev_time = ext->fn[1];
      Func<Ord>* xvel_prev_newton = u_ext[0];  
      Func<Ord>* yvel_prev_newton = u_ext[1];  
      Func<Ord>* p_prev_newton = u_ext[2];
      for (int i = 0; i < n; i++)
        result += wt[i] * ((xvel_prev_newton->dx[i] * v->dx[i] + xvel_prev_newton->dy[i] * v->dy[i]) / Reynolds - (p_prev_newton->val[i] * v->dx[i]));
      if(!Stokes)
        for (int i = 0; i < n; i++)
          result += wt[i] * (((xvel_prev_newton->val[i] - xvel_prev_time->val[i]) * v->val[i] / time_step )
          + ((xvel_prev_newton->val[i] * xvel_prev_newton->dx[i] + yvel_prev_newton->val[i] * xvel_prev_newton->dy[i]) * v->val[i]));
      return result;
    }
    
    VectorFormVol<double>* clone()
    {
      return new VectorFormNS_1(this->i, this->Stokes, this->Reynolds, this->time_step);
    }
  protected:
    // Members.
    bool Stokes;
    double Reynolds;
    double time_step;
  };

  class VectorFormNS_2 : public VectorFormVol<double>
  {
  public:
    VectorFormNS_2(int i) : VectorFormVol<double>(i) {
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext) const{
      double result = 0;
      Func<double>* xvel_prev_newton = u_ext[0];  
      Func<double>* yvel_prev_newton = u_ext[1];  

      for (int i = 0; i < n; i++)
        result += wt[i] * (xvel_prev_newton->dx[i] * v->val[i] + yvel_prev_newton->dy[i] * v->val[i]);
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      Ord result = Ord(0);
      Func<Ord>* xvel_prev_newton = u_ext[0];  
      Func<Ord>* yvel_prev_newton = u_ext[1];  

      for (int i = 0; i < n; i++)
        result += wt[i] * (xvel_prev_newton->dx[i] * v->val[i] + yvel_prev_newton->dy[i] * v->val[i]);
      return result;
    }
    
    VectorFormVol<double>* clone()
    {
      return new VectorFormNS_2(this->i);
    }
  };

protected:
  // Members.
  bool Stokes;
  double Reynolds;
  Solution<double>* x_vel_previous_time;
  Solution<double>* y_vel_previous_time;
};

class EssentialBCNonConst : public EssentialBoundaryCondition<double>
{
public:
  EssentialBCNonConst(Hermes::vector<std::string> markers, double vel_inlet, double H, double startup_time) : 
      EssentialBoundaryCondition<double>(markers), vel_inlet(vel_inlet), H(H), startup_time(startup_time) {};
      EssentialBCNonConst(std::string marker, double vel_inlet, double H, double startup_time) : 
      EssentialBoundaryCondition<double>(Hermes::vector<std::string>()), vel_inlet(vel_inlet), H(H), startup_time(startup_time) {
        markers.push_back(marker);
      };

      ~EssentialBCNonConst() {};

      virtual EssentialBCValueType get_value_type() const { 
        return BC_FUNCTION; 
      };

      virtual double value(double x, double y, double n_x, double n_y, double t_x, double t_y) const {
        double val_y = vel_inlet * y*(H-y) / (H/2.)/(H/2.);
        if (current_time <= startup_time) 
          return val_y * current_time/startup_time;
        else 
          return val_y;
      };

protected:
  // Members.
  double vel_inlet;
  double H;
  double startup_time;
};
