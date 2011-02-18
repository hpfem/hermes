#include <hermes3d.h>
#include <norm.h>
#include <norm.cpp>
#include <vector>

#define PID_DEFAULT_TOLERANCE 0.25
#define DEFAULT_STEP 0.1

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


