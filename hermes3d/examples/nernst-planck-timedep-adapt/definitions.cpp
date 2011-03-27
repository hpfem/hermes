#include <hermes3d.h>
#include <norm.h>
#include <norm.cpp>
#include <vector>

#define PID_DEFAULT_TOLERANCE 0.25
#define DEFAULT_STEP 0.1

/* timestep controller */

class PidTimestepController {

public:
  PidTimestepController(double final_time, bool pid_on = true,
      double default_step = DEFAULT_STEP, double tolerance = PID_DEFAULT_TOLERANCE) {
    this->delta = tolerance;
    this->final_time = final_time;
    this->time = 0;
    this->step_number = 0;
    timestep = new double;
    (*timestep) = default_step;
    this->pid = pid_on;
    finished = false;
  };

  int get_timestep_number() {return step_number;};
  double get_time() {return time;};

  // true if next time step can be run, false if the time step must be re-run with smaller time step.
  bool end_step(Hermes::vector<Solution*> solutions, Hermes::vector<Solution *> prev_solutions);
  void begin_step();
  bool has_next();

  // reference to the current calculated time step
  double *timestep;

private:
  bool pid;
  double delta;
  double final_time;
  double time;
  int step_number;
  std::vector<double> err_vector;
  bool finished;

  // PID parameters
  const static double kp;
  const static double kl;
  const static double kD;

};

  const double PidTimestepController::kp = 0.075;
  const double PidTimestepController::kl = 0.175;
  const double PidTimestepController::kD = 0.01;

// Usage: do{ begin_step() calculations .... end_step(..);} while(has_next());

void PidTimestepController::begin_step() {

  if ((time + (*timestep)) >= final_time) {
    info("Time step would exceed the final time... reducing");
    (*timestep) = final_time - time;
    info("The last time step: %g", *timestep);
    finished = true;
  }
  time += (*timestep);
  step_number++;
  info("begin_step processed, new step number: %i and cumulative time: %g", step_number, time);
}

bool PidTimestepController::end_step(Hermes::vector<Solution *> solutions,
    Hermes::vector<Solution *> prev_solutions) {

  if (pid) {
    info("Running PID calculations...");

    unsigned int neq = solutions.size();
    if (neq == 0) {
      return true;
    }
    if (prev_solutions.empty()) {
      return true;
    }
    if (neq != prev_solutions.size()) {
      error_function("Inconsistent parameters in PidTimestepController::next(...)");
    }
    double max_rel_error = 0.0;

    for (unsigned int i = 0; i < neq; i++) {
      double abs_error = calc_error(error_fn_h1, solutions[i], prev_solutions[i]);

      double norm = calc_norm(norm_fn_h1, solutions[i]);
      max_rel_error = (abs_error/norm > max_rel_error) ? (abs_error/norm) : max_rel_error;

      info("Solution[%i] abs error %g and norm %g and the largest relative error %g",
          i, abs_error, norm, max_rel_error);
    }

    err_vector.push_back(max_rel_error);

    if (err_vector.size() > 2 && max_rel_error <= delta) {
      int size = err_vector.size();
      info("Error vector sufficient for adapting...");
      double t_coeff = pow(err_vector.at(size - 2)/err_vector.at(size-1),kp)
          * pow(0.25/err_vector.at(size - 1), kl)
          * pow(err_vector.at(size - 2)*err_vector.at(size - 2)/(err_vector.at(size -1)*err_vector.at(size-3)), kD);
       info("Coefficient %g", t_coeff);
       (*timestep) = (*timestep)*t_coeff;
       info("New time step: %g", *timestep);
    }
  } // end pid

  return true;

}

bool PidTimestepController::has_next() {
  return !finished;
}

/* Definition of residiual vectors */

template<class Real, class Scalar>
Scalar Fc_euler(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
	Scalar result = 0;
	Func<Scalar>* C_prev_time = ext->fn[0];
	Func<Scalar>* C_prev_newton = u_ext[0];
	Func<Scalar>* phi_prev_newton = u_ext[1];

	for (int i = 0; i < n; i++) {
		result += wt[i] * ((C_prev_newton->val[i] - C_prev_time->val[i]) * v->val[i] / (*TAU) +
				D * (C_prev_newton->dx[i] * v->dx[i] + C_prev_newton->dy[i] * v->dy[i] +
				    C_prev_newton->dz[i] * v->dz[i]) +
				K * C_prev_newton->val[i] * (phi_prev_newton->dx[i] * v->dx[i] + phi_prev_newton->dy[i] * v->dy[i] +
				    phi_prev_newton->dz[i] * v->dz[i]));
	}
	return result;
}

template<class Real, class Scalar>
Scalar Fphi_euler(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
	Scalar result = 0;
	Func<Scalar>* C_prev_newton = u_ext[0];
	Func<Scalar>* phi_prev_newton = u_ext[1];

	for (int i = 0; i < n; i++) {

	  result += wt[i] * ((phi_prev_newton->dx[i] * v->dx[i] + phi_prev_newton->dy[i] * v->dy[i] +
	      phi_prev_newton->dz[i] * v->dz[i]) +
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
				D * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i] + u->dz[i] * v->dz[i]) +
				K * u->val[i] * (phi_prev_newton->dx[i] * v->dx[i] + phi_prev_newton->dy[i] * v->dy[i] +
				    phi_prev_newton->dz[i] * v->dz[i]));
	}
	return result;
}

//matrix 0_1
template<class Real, class Scalar>
Scalar J_euler_DFcDYphi(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
	Scalar result = 0;
	Func<Scalar>* C_prev_newton = u_ext[0];
	for (int i = 0; i < n; i++) {
		result += wt[i] * K * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i] + u->dz[i] * v->dz[i]) * C_prev_newton->val[i];
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
		result += wt[i] * ( u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i] + u->dz[i] * v->dz[i]);
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
        0.5 * D * (C_prev_newton->dx[i] * v->dx[i] + C_prev_newton->dy[i] * v->dy[i] +
            C_prev_newton->dz[i] * v->dz[i]) +
        0.5 * D * (C_prev_time->dx[i] * v->dx[i] + C_prev_time->dy[i] * v->dy[i] +
            C_prev_time->dz[i] * v->dz[i]) +
        0.5 * K * C_prev_newton->val[i] * (phi_prev_newton->dx[i] * v->dx[i] + phi_prev_newton->dy[i] * v->dy[i] +
            phi_prev_newton->dz[i] * v->dz[i]) +
        0.5 * K * C_prev_time->val[i] * (phi_prev_time->dx[i] * v->dx[i] + phi_prev_time->dy[i] * v->dy[i] +
            phi_prev_time->dz[i] * v->dz[i]));
  }
  return result;
}

template<class Real, class Scalar>
Scalar Fphi_cranic(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
    Scalar result = 0;
    Func<Scalar>* C_prev_newton = u_ext[0];
    Func<Scalar>* phi_prev_newton = u_ext[1];
    for (int i = 0; i < n; i++) {
      result += wt[i] * ((phi_prev_newton->dx[i] * v->dx[i] + phi_prev_newton->dy[i] * v->dy[i]+
          phi_prev_newton->dz[i] * v->dz[i]) +
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
        0.5 * D * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i] + u->dz[i] * v->dz[i]) +
        0.5 * K * u->val[i] * (phi_prev_newton->dx[i] * v->dx[i] + phi_prev_newton->dy[i] * v->dy[i] +
            phi_prev_newton->dz[i] * v->dz[i]));
  }
  return result;
}

template<class Real, class Scalar>
Scalar J_cranic_DFcDYphi(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
  Scalar result = 0;
  Func<Scalar>* C_prev_newton = u_ext[0];
	for (int i = 0; i < n; i++) {
	  result += wt[i] * (0.5 * K * C_prev_newton->val[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i] +
	      u->dz[i] * v->dz[i]));
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
    result += wt[i] * ( u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i] + u->dz[i] * v->dz[i]);
  }
  return result;
}

