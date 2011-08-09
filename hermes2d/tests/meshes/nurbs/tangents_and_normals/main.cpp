#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
// This example makes sure that normal and tangential vectors to circular
// arcs nurbs are calculated correctly.
//
// NOTE: This test strongly depends on the concrete geometry (domain.mesh).
// When changing the geometry, also change the interior of the value() method
// of the CustomEssentialBCNonConst classes below.
//
// Domain: A square inscribed into the unit circle, with
// two curved edges - see file domain.mesh.

const int P_INIT = 6;                             // Initial polynomial degree of all elements.
const double TOLERANCE = 1e-1;                    // Tolerance for the match of normal and tangential
                                                  // vectors.
const int INIT_REF_NUM = 1;                       // Number of initial uniform mesh refinements.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// This variable is global since it could not be changed being local data in the const method value() below.
bool SUCCESS = true;

// Custom class for essential boundary conditions that
// tests the normal and tangential vectors on the right edge.
class CustomEssentialBCNonConstRight : public EssentialBoundaryCondition {
public:
  CustomEssentialBCNonConstRight(std::string marker)
           : EssentialBoundaryCondition(marker) { }

  ~CustomEssentialBCNonConstRight() {};

  inline EssentialBCValueType get_value_type() const {
    return EssentialBoundaryCondition::BC_FUNCTION;
  }

  // This function also checks correctness of normal and tangential vectors for
  // each boundary point that comes here.
  virtual double value(double x, double y, double n_x, double n_y, double t_x, double t_y) const
  {
    //printf("right x, y = %g, %g\n", x, y);

    // Points lie on the unit circle -> the normal vector is the point itself.
    double n_x_test = x;
    double n_y_test = y;
    double t_x_test = -y;
    double t_y_test = x;
    if (std::abs(n_x_test - n_x) > TOLERANCE || std::abs(n_y_test - n_y) > TOLERANCE ||
        std::abs(t_x_test - t_x) > TOLERANCE || std::abs(t_y_test - t_y) > TOLERANCE) {
      printf("Bdy point: [%g, %g], Hermes normal: (%g, %g), correct_normal: (%g, %g)\n",
             x, y, n_x, n_y, n_x_test, n_y_test);
      printf("                     Hermes tangent: (%g, %g), correct_tangent: (%g, %g)\n",
             t_x, t_y, t_x_test, t_y_test);
      SUCCESS = false;
    }

    return 0.0;
  }
};

// Custom class for essential boundary conditions that
// tests the normal and tangential vectors on the left edge.
class CustomEssentialBCNonConstLeft : public EssentialBoundaryCondition<double> {
public:
  CustomEssentialBCNonConstLeft(std::string marker)
           : EssentialBoundaryCondition<double>(marker) { }

  ~CustomEssentialBCNonConstLeft() {};

  inline EssentialBCValueType get_value_type() const {
    return EssentialBoundaryCondition::BC_FUNCTION;
  }

  // This function also checks correctness of normal and tangential vectors for
  // each boundary point that comes here.
  virtual double value(double x, double y, double n_x, double n_y, double t_x, double t_y) const
  {
    //printf("left x, y = %g, %g\n", x, y);

    // Normal vector is the point itself, moved by 1.0 to the right.
    double n_x_test = -(x + sqrt(2.0));
    double n_y_test = -y;
    double t_x_test = y;
    double t_y_test = -(x + sqrt(2.0));
    if (std::abs(n_x_test - n_x) > TOLERANCE || std::abs(n_y_test - n_y) > TOLERANCE ||
        std::abs(t_x_test - t_x) > TOLERANCE || std::abs(t_y_test - t_y) > TOLERANCE) {
      printf("Bdy point: [%g, %g], Hermes normal: (%g, %g), correct_normal: (%g, %g)\n",
             x, y, n_x, n_y, n_x_test, n_y_test);
      printf("                                  Hermes tangent: (%g, %g), correct_tangent: (%g, %g)\n",
             t_x, t_y, t_x_test, t_y_test);
      SUCCESS = false;
    }

    return 0.0;
  }
};

// Weak forms
class CustomWeakFormPoisson : public WeakForm
{
public:
  CustomWeakFormPoisson(double rhs_const) : WeakForm(1)
  {
    // Jacobian.
    add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion());

    // Residual.
    add_vector_form(new WeakFormsH1::DefaultResidualDiffusion());
    add_vector_form(new WeakFormsH1::DefaultVectorFormVol(0, HERMES_ANY, new HermesFunction(-rhs_const)));
  };
};

int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Hermes2D hermes2d;

  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", &mesh);

  // Initial uniform mesh refinements.
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Show the mesh.
  //MeshView mv("Mesh", new WinGeom(0, 0, 400, 400));
  //mv.show(&mesh);
  //View::wait(HERMES_WAIT_KEYPRESS);

  // Initialize boundary conditions.
  DefaultEssentialBCConst bc_bottom("Bottom", 0.0);
  DefaultEssentialBCConst bc_top("Top", 0.0);
  CustomEssentialBCNonConstRight bc_right("Right");
  CustomEssentialBCNonConstLeft bc_left("Left");
  EssentialBCs bcs(Hermes::vector<EssentialBoundaryCondition*>(&bc_bottom, &bc_top, &bc_left, &bc_right));

  // Initialize the weak formulation.
  CustomWeakFormPoisson wf(1.0);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bcs, P_INIT);
  int ndof = space.get_num_dofs();
  info("ndof = %d", ndof);

  // Initialize the FE problem.
  DiscreteProblem dp(&wf, &space);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Initialize the solution.
  Solution sln;

  // Assemble the stiffness matrix and right-hand side vector.
  //info("Assembling the stiffness matrix and right-hand side vector.");
  //dp.assemble(matrix, rhs);

  // Solve the linear system and if successful, obtain the solution.
  //info("Solving the matrix problem.");
  //if(solver->solve()) Solution::vector_to_solution(solver->get_solution(), &space, &sln);
  //else error ("Matrix solver failed.\n");

  //doubleView view("Solution", new WinGeom(0, 0, 440, 350));
  //view.show(&sln, HERMES_EPS_HIGH);
  //View::wait();

  // Clean up.
  delete solver;
  delete matrix;
  delete rhs;

  // Test
  if (SUCCESS == true) {
    printf("Success!\n");
    return TEST_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return TEST_FAILURE;
  }

  return 0;
}

