#define HERMES_REPORT_INFO
#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

// This test makes sure that NURBS work correctly.

const char* mesh_file_1 = "domain-1.mesh";          // One control point.
const char* mesh_file_2 = "domain-2.mesh";          // Two control points.
const char* mesh_file_3 = "domain-3.mesh";          // Three control points.

// The following parameters can be also changed:

const int P_INIT = 3;                             // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 2;                       // Number of initial uniform mesh refinements.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters.
const double const_f = 1.0;

int main(int argc, char* argv[])
{
  // Check number of command-line parameters.
  if(argc < 2)
    throw Hermes::Exceptions::Exception("Not enough parameters.");

  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  if(strcasecmp(argv[1], "1") == 0)
    mloader.load(mesh_file_1, &mesh);
  if(strcasecmp(argv[1], "2") == 0)
    mloader.load(mesh_file_2, &mesh);
  if(strcasecmp(argv[1], "3") == 0)
    mloader.load(mesh_file_3, &mesh);

  // Perform initial mesh refinements (optional).
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize boundary conditions.
  DefaultEssentialBCConst<double> bc_essential("Bdy", 0.0);
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space<double> space(&mesh, &bcs, P_INIT);
  int ndof = Space<double>::get_num_dofs(&space);

  // Initialize the weak formulation.
  WeakFormsH1::DefaultWeakFormPoisson<double> wf(HERMES_ANY, new Hermes1DFunction<double>(1.0), new Hermes2DFunction<double>(-const_f));

  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, &space);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix<double>* matrix = create_matrix<double>();
  Vector<double>* rhs = create_vector<double>();
  LinearMatrixSolver<double>* solver = create_linear_solver<double>(matrix, rhs);

  // Initial coefficient vector for the Newton's method.
  double* coeff_vec = new double[ndof];
  memset(coeff_vec, 0, ndof*sizeof(double));

  // Perform Newton's iteration.
  Hermes::Hermes2D::Solution<double> sln;
  Hermes::Hermes2D::NewtonSolver<double> newton(&dp);
  try{
    newton.solve(coeff_vec);
  }
  catch(Hermes::Exceptions::Exception& e)
  {
    e.printMsg();
  }
  Hermes::Hermes2D::Solution<double>::vector_to_solution(newton.get_sln_vector(), &space, &sln);

  // Clean up.
  delete solver;
  delete matrix;
  delete rhs;

  double coor_x[4] = {0.3, 0.7, 1.3, 1.7};
  double coor_y = 0.5;

  double value[4] = {0.102569, 0.167907, 0.174203, 0.109630};
  if(strcasecmp(argv[1], "2") == 0)
  {
    value[0] = 0.062896;
    value[1] = 0.096658;
    value[2] = 0.114445;
    value[3] = 0.081221;
  }

  if(strcasecmp(argv[1], "3") == 0)
  {
    value[0] = 0.048752;
    value[1] = 0.028585;
    value[2] = 0.028585;
    value[3] = 0.048752;
  }

  bool success = true;
  for (int i = 0; i < 4; i++)
  {
    if(Hermes::abs(value[i] - sln.get_pt_value(coor_x[i], coor_y)) > 1E-6)
      success = false;
  }
  if(success)
  {
    printf("Success!\n");
    return 0;
  }
  else
  {
    printf("Failure!\n");
    return -1;
  }
}