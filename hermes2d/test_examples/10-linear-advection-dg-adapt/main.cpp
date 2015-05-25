#include "definitions.h"

/*

Flux

\begin{eqnarray}
\int_{\partial K} \mathbf{f}\left(u\right) v \cdot \mathbf{n} & = & 	\int_{\partial K} \left(\frac{1}{2}{u^2},\ \frac{1}{2}{u^2}\right) v \cdot \mathbf{n} \\ & \approx &
\int_{\partial K} \frac{1}{2}\left(\frac{1}{2}{u_R}^2 + \frac{1}{2}{u_L}^2,\ \frac{1}{2}{u_R}^2 + \frac{1}{2}{u_L}^2\right) v \cdot \mathbf{n} \\ & - &
\int_{\partial K} \max \left\{|\nabla \left({u_R}^2\right)|, |\nabla \left({u_L}^2\right)|\right\} \left(u_R - u_L\right)\ v
\end{eqnarray}

*/

// Number of initial uniform mesh refinements.
const int INIT_REF = 3;
// Initial polynomial degrees of mesh elements in vertical and horizontal directions.
int P_INIT = 0;

// Final time
const double T_initial = 0.;
const double T_final = 5.;

// Time step
const double T_step = 2.5e-6;

bool START_FROM_ZERO = false;
int EVERY_NTH_STEP = 20000;

int main(int argc, char* args[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("square-triangular.mesh", mesh);

  // Perform initial mesh refinement.
  for (int i = 0; i < INIT_REF; i++)
    mesh->refine_all_elements();

  // Create an L2 space.
  SpaceSharedPtr<double> space(new L2Space<double>(mesh, P_INIT));

  MeshFunctionSharedPtr<double> prev_sln(START_FROM_ZERO ? (ExactSolution<double>*)new ZeroSolution<double>(mesh) : (ExactSolution<double>*)new CustomExactSolution(mesh));
  MeshFunctionSharedPtr<double> sln(new Solution<double>);

  // Initialize the weak formulation.
  WeakFormSharedPtr<double> wf(new CustomWeakForm(prev_sln, mesh));
  wf->set_current_time_step(T_step);

  ScalarView view1("Solution", new WinGeom(900, 0, 450, 350));
  view1.fix_scale_width(60);

  // Initialize linear solver.
  Hermes::Hermes2D::LinearSolver<double> linear_solver(wf, space);
  linear_solver.set_verbose_output(false);
  //linear_solver.output_matrix();
  //linear_solver.output_rhs();

  Linearizer lin(LinearizerOutputType::FileExport);

  double t = T_initial;
  int time_step_counter = 0;
  for (; t < T_final;)
  {
    try
    {
      t += T_step;
      if (time_step_counter++ % EVERY_NTH_STEP == 0)
        std::cout << "Iteration: " << time_step_counter << ", Current time : " << t << std::endl;

      linear_solver.set_time(t);

      linear_solver.solve();

      Solution<double>::vector_to_solution(linear_solver.get_sln_vector(), space, prev_sln);

      if (time_step_counter % EVERY_NTH_STEP == 0)
      {
        view1.show(prev_sln);
        std::stringstream ss;
        ss << "Solution-";
        ss << time_step_counter << "(t=" << t << "s).tcp";
        lin.save_solution_tecplot(prev_sln, ss.str().c_str(), "Sln");
        ss.clear();
        ss << "Solution-";
        ss << time_step_counter << "(t=" << t << "s).vtk";
        lin.save_solution_vtk(prev_sln, ss.str().c_str(), "Sln");
      }
    }
    catch (Exceptions::Exception& e)
    {
      std::cout << e.info();
    }
    catch (std::exception& e)
    {
      std::cout << e.what();
    }
  }
  
  // Wait for keyboard or mouse input.
  View::wait();
  return 0;
}