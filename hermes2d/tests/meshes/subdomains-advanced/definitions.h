#include "hermes2d.h"
#include "runge_kutta.h"

/* Weak forms */

using namespace Hermes;
using namespace Hermes::Hermes2D;

class CustomWeakFormHeatAndFlow : public WeakForm<double>
{
public:
  CustomWeakFormHeatAndFlow(bool Stokes, double Reynolds, double time_step, Solution<double>* x_vel_previous_time, Solution<double>* y_vel_previous_time, Solution<double>* T_prev_time);

  class BilinearFormTime: public MatrixFormVol<double>
  {
  public:
    BilinearFormTime(int i, int j, double time_step) : MatrixFormVol<double>(i, j), time_step(time_step) {
        sym = HERMES_SYM;
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const
    {
      double result = int_u_v<double, double>(n, wt, u, v) / time_step;
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const
    {
      Ord result = int_u_v<Ord, Ord>(n, wt, u, v) / time_step;
      return result;
    }
  protected:
    // Members.
    double time_step;
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
  protected:
    // Members.
    bool Stokes;
    double Reynolds;
    double time_step;
  };

  class CustomJacobianAdvection : public MatrixFormVol<double>
  {
  public:
    CustomJacobianAdvection(int i, int j, std::string area) : MatrixFormVol<double>(i, j, area) {}

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const{
      double result = 0;
      Func<double>* xvel_prev_newton = u_ext[0];
      Func<double>* yvel_prev_newton = u_ext[1];
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->dx[i] * v->val[i] * xvel_prev_newton->val[i] + u->dy[i] * v->val[i] * yvel_prev_newton->val[i]);
    return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const
    {
      Ord result = Ord(0);
      Func<Ord>* xvel_prev_newton = u_ext[0];
      Func<Ord>* yvel_prev_newton = u_ext[1];
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->dx[i] * v->val[i] * xvel_prev_newton->val[i] + u->dy[i] * v->val[i] * yvel_prev_newton->val[i]);
      return result;
    }
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
        Func<double>* yvel_prev_newton = u_ext[0];
        for (int i = 0; i < n; i++)
          result += wt[i] * (u->val[i] * v->val[i] * yvel_prev_newton->dx[i]);
      }
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const{
      Ord result = Ord(0);
      if(!Stokes) {
        Func<Ord>* yvel_prev_newton = u_ext[0];
        for (int i = 0; i < n; i++)
          result += wt[i] * (u->val[i] * v->val[i] * yvel_prev_newton->dx[i]);
      }
      return result;
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
  };

  class VectorFormTime: public VectorFormVol<double>
  {
  public:
    VectorFormTime(int i, double time_step) : VectorFormVol<double>(i), time_step(time_step) {}

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext) const
    {
      Func<double>* func_prev_time = ext->fn[0];
      double result = (int_u_v<double, double>(n, wt, u_ext[3], v) - int_u_v<double, double>(n, wt, func_prev_time, v)) / time_step;
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const
    {
      Func<Ord>* func_prev_time = ext->fn[0];
      Ord result = (int_u_v<Ord, Ord>(n, wt, u_ext[3], v) - int_u_v<Ord, Ord>(n, wt, func_prev_time, v)) / time_step;
      return result;
    }
  protected:
    // Members.
    double time_step;
  };

  class CustomResidualAdvection : public VectorFormVol<double>
  {
  public:
    CustomResidualAdvection(int i, std::string area) : VectorFormVol<double>(i, area) {}

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext) const{
      double result = 0;
      Func<double>* xvel_prev_newton = u_ext[0];
      Func<double>* yvel_prev_newton = u_ext[1];
      for (int i = 0; i < n; i++)
        result += wt[i] * (u_ext[3]->dx[i] * v->val[i] * xvel_prev_newton->val[i] + u_ext[3]->dy[i] * v->val[i] * yvel_prev_newton->val[i]);
    return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const
    {
      Ord result = Ord(0);
      Func<Ord>* xvel_prev_newton = u_ext[0];
      Func<Ord>* yvel_prev_newton = u_ext[1];
      for (int i = 0; i < n; i++)
        result += wt[i] * (u_ext[3]->dx[i] * v->val[i] * xvel_prev_newton->val[i] + u_ext[3]->dy[i] * v->val[i] * yvel_prev_newton->val[i]);
      return result;
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
          + ((xvel_prev_newton->val[i] * yvel_prev_newton->dx[i] + yvel_prev_newton->val[i] * yvel_prev_newton->dy[i]) * v->val[i]));
      return result;
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
  };

protected:
  // Members.
  bool Stokes;
  double Reynolds;
  double time_step;
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

class ZeroInitialCondition : public ExactSolutionScalar<double>
{
public:
  ZeroInitialCondition(Mesh* mesh) : ExactSolutionScalar<double>(mesh) {};

  virtual double value (double x, double y) const {
    return 0.0; 
  };

  virtual void derivatives (double x, double y, double& dx, double& dy) const {
    dx = 0;
    dy = 0;
  };

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(0);
  }
};