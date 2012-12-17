#define HERMES_REPORT_ALL
#include "definitions.h"

// This example shows how to solve a simple PDE that describes stationary
// heat transfer in an object consisting of two materials (aluminum and
// copper). The object is heated by constant volumetric heat sources
// (generated, for example, by a DC electric current). The temperature
// on the boundary is fixed. We will learn how to:
//
//   - load the mesh,
//   - perform initial refinements,
//   - create a H1 space over the mesh,
//   - define weak formulation,
//   - initialize matrix solver,
//   - assemble and solve the matrix system,
//   - output the solution and element orders in VTK format
//     (to be visualized, e.g., using Paraview),
//   - visualize the solution using Hermes' native OpenGL-based functionality.
//
// PDE: Poisson equation -div(LAMBDA grad u) - VOLUME_HEAT_SRC = 0.
//
// Boundary conditions: Dirichlet u(x, y) = FIXED_BDY_TEMP on the boundary.
//
// Geometry: L-Shape domain (see file domain.mesh).
//
// The following parameters can be changed:

const bool HERMES_VISUALIZATION = true;           // Set to "false" to suppress Hermes OpenGL visualization.
const bool VTK_VISUALIZATION = false;              // Set to "true" to enable VTK output.
const int P_INIT = 2;                             // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 1;                       // Number of initial uniform mesh refinements.

// Problem parameters.
const double LAMBDA_AL = 236.0;            // Thermal cond. of Al for temperatures around 20 deg Celsius.
const double LAMBDA_CU = 386.0;            // Thermal cond. of Cu for temperatures around 20 deg Celsius.
const double VOLUME_HEAT_SRC = 1.0;          // Volume heat sources generated (for example) by electric current.
const double FIXED_BDY_TEMP = 31.0;        // Fixed temperature on the boundary.

// Globals.
Hermes::Hermes2D::H1Space<double>* space;
Hermes::Hermes2D::Mesh* mesh;

void loadCache(int* i_indices_A, int* j_indices_A, int* i_indices_rhs, double** A_values, double* rhs_values);
void loadProblemData();
void calculateCache(CustomWeakFormPoisson& wf);
void calculateResult(CustomWeakFormPoisson& wf, int* i_indices_A, int* j_indices_A, int* i_indices_rhs, double** A_values, double* rhs_values);
double handleDirichlet(CustomWeakFormPoisson& wf, int shape_indexBasis, int shape_indexDirichlet, Element* e);

int main(int argc, char* argv[])
{
  // Initialize the weak formulation.
  CustomWeakFormPoisson wf("Aluminum", new Hermes::Hermes1DFunction<double>(LAMBDA_AL), "Copper",
    new Hermes::Hermes1DFunction<double>(LAMBDA_CU), new Hermes::Hermes2DFunction<double>(VOLUME_HEAT_SRC));
  
  calculateCache(wf);

  loadProblemData();

  int max_index = space->get_shapeset()->get_max_index(HERMES_MODE_QUAD) + 1;

  int* i_indices_A = new int[max_index*max_index];
  int* j_indices_A = new int[max_index*max_index];
  int* i_indices_rhs = new int[max_index];
  double** A_values = new double*[max_index];
  for(unsigned int i = 0; i < max_index; i++)
    A_values[i] = new double[max_index];
  double* rhs_values = new double[max_index];

  loadCache(i_indices_A, j_indices_A, i_indices_rhs, A_values, rhs_values);

  calculateResult(wf, i_indices_A, j_indices_A, i_indices_rhs, A_values, rhs_values);
}

void loadCache(int* i_indices_A, int* j_indices_A, int* i_indices_rhs, double** A_values, double* rhs_values)
{
  int counter = 0;
  std::ifstream matrixFormIn("matrix");
  std::ifstream rhsFormIn("rhs");

  double valueTemp;
  while(matrixFormIn.good())
  {
    matrixFormIn >> (int)i_indices_A[counter] >> (int)j_indices_A[counter] >> valueTemp;
    A_values[i_indices_A[counter]][j_indices_A[counter]] = valueTemp;
    counter++;
  }

  counter = 0;
  while(rhsFormIn.good())
  {
    rhsFormIn >> i_indices_rhs[counter] >> valueTemp;
    rhs_values[i_indices_rhs[counter]] = valueTemp;
    counter++;
  }
}

void  loadProblemData()
{
  // Load the mesh.
  Hermes::Hermes2D::MeshReaderH2DXML mloader;
  mesh = new Mesh();
  mloader.load("domain.xml", mesh);

  for(int i = 0; i < INIT_REF_NUM; i++)
    mesh->refine_all_elements();

  // Create a H1 space with default shapeset.
  Hermes::Hermes2D::DefaultEssentialBCConst<double>* bc_essential = new Hermes::Hermes2D::DefaultEssentialBCConst<double>("Bnd", FIXED_BDY_TEMP);
  Hermes::Hermes2D::EssentialBCs<double>* bcs = new Hermes::Hermes2D::EssentialBCs<double>(bc_essential);

  // Create an H1 space with default shapeset.
  space = new Hermes::Hermes2D::H1Space<double>(mesh, bcs, P_INIT);
}

void calculateCache(CustomWeakFormPoisson& wf)
{
  // Load the mesh.
  Hermes::Hermes2D::MeshReaderH2DXML mloader;
  mesh = new Mesh();
  mloader.load("domainCacheComputation.xml", mesh);

  // Create an H1 space with default shapeset.
  space = new Hermes::Hermes2D::H1Space<double>(mesh, 10);

  // Initialize linear solver.
  Hermes::Hermes2D::LinearSolver<double> linear_solver(&wf, space);

  wf.set_global_integration_order(24);

  // Solve the linear problem.
  linear_solver.solve();

  int ndof = space->get_num_dofs();

  std::ofstream matrixFormOut("matrix");
  std::ofstream rhsFormOut("rhs");

  AsmList<double> al;
  space->get_element_assembly_list(mesh->get_element(0), &al);

  for(int i = 0; i < al.get_cnt(); i++)
  {
    if(i > 0)
      rhsFormOut << std::endl;
    for(int j = 0; j < al.get_cnt(); j++)
    {
      if(i > 0 || j > 0)
        matrixFormOut << std::endl;
      matrixFormOut << al.get_idx()[i] << ' ' << al.get_idx()[j] << ' ' << linear_solver.get_jacobian()->get(i, j) / al.get_coef()[i] / al.get_coef()[j];
    }
    rhsFormOut << al.get_idx()[i] << ' ' << linear_solver.get_residual()->get(i) / al.get_coef()[i];
  }

  matrixFormOut.close();
  rhsFormOut.close();
   
  return;
}

void calculateResult(CustomWeakFormPoisson& wf, int* i_indices_A, int* j_indices_A, int* i_indices_rhs, double** A_values, double* rhs_values)
{
  // Utilities.
  int ndof = space->get_num_dofs();
  Element* e;

  // Initialize the solution.
  Hermes::Hermes2D::Solution<double> sln;

  /// The solution vector.
  double* sln_vector;

  /// Jacobian.
  SparseMatrix<double>* jacobian = create_matrix<double>();
  jacobian->prealloc(ndof);
  for_all_active_elements(e, mesh)
  {
    AsmList<double> al;
    space->get_element_assembly_list(e, &al);
    for(int i = 0; i < al.get_cnt(); i++)
    {
      if(al.get_dof()[i] < 0)
        continue;
      for(int j = 0; j < al.get_cnt(); j++)
      {
        if(al.get_dof()[j] < 0)
          continue;
        jacobian->pre_add_ij(al.get_dof()[i], al.get_dof()[j]);
        jacobian->pre_add_ij(al.get_dof()[j], al.get_dof()[i]);
      }
    }
  }
  jacobian->alloc();

  /// Residual.
  Vector<double>* residual = create_vector<double>();
  residual->alloc(ndof);
  
  /// Linear solver.
  LinearMatrixSolver<double>* matrix_solver = create_linear_solver<double>(jacobian, residual);

  for_all_active_elements(e, mesh)
  {
    RefMap refmap;
    refmap.set_active_element(e);
    double inv_ref_map_determinant = refmap.get_const_jacobian();
    AsmList<double> al;
    space->get_element_assembly_list(e, &al);
    for(int i = 0; i < al.get_cnt(); i++)
    {
      if(al.get_dof()[i] < 0)
        continue;
      for(int j = 0; j < al.get_cnt(); j++)
      {
        if(al.get_dof()[j] >= 0)
        {
          jacobian->add(al.get_dof()[i], al.get_dof()[j], A_values[al.get_idx()[i]][al.get_idx()[j]] * al.get_coef()[i] * al.get_coef()[j]);
          jacobian->add(al.get_dof()[j], al.get_dof()[i], A_values[al.get_idx()[i]][al.get_idx()[j]] * al.get_coef()[i] * al.get_coef()[j]);
        }
        else
          residual->add(al.get_dof()[i], -handleDirichlet(wf, al.get_idx()[i], al.get_idx()[j], e) * al.get_coef()[i] * al.get_coef()[j] * 2);

      }
      residual->add(al.get_dof()[i], rhs_values[al.get_idx()[i]] * al.get_coef()[i] * inv_ref_map_determinant);
    }
  }

  matrix_solver->solve();

  sln_vector = matrix_solver->get_sln_vector();

  // Translate the solution vector into the previously initialized Solution.
  Hermes::Hermes2D::Solution<double>::vector_to_solution(sln_vector, space, &sln, true);

  // VTK output.
  if(VTK_VISUALIZATION)
  {
    // Output solution in VTK format.
    Hermes::Hermes2D::Views::Linearizer lin;
    bool mode_3D = false;
    lin.save_solution_vtk(&sln, "sln.vtk", "Temperature", mode_3D, 1, Hermes::Hermes2D::Views::HERMES_EPS_LOW);

    // Output mesh and element orders in VTK format.
    Hermes::Hermes2D::Views::Orderizer ord;
    ord.save_mesh_vtk(space, "mesh.vtk");
    ord.save_orders_vtk(space, "ord.vtk");
  }

  // Visualize the solution.
  Hermes::Hermes2D::Views::ScalarView viewS("Solution", new Hermes::Hermes2D::Views::WinGeom(50, 50, 1000, 800));

  if(HERMES_VISUALIZATION)
  {
    viewS.show(&sln, Hermes::Hermes2D::Views::HERMES_EPS_LOW);
    viewS.wait_for_close();
  }
  
  return;
}

double handleDirichlet(CustomWeakFormPoisson& wf, int shape_indexBasis, int shape_indexDirichlet, Element* e)
{
  RefMap refmap;
  refmap.set_active_element(e);
  Geom<double>* geometry;
  double* jacobian_x_weights;
  int n_quadrature_points = DiscreteProblem<double>::init_geometry_points(&refmap, 24, geometry, jacobian_x_weights);
  PrecalcShapeset pssBasis(space->get_shapeset());
  pssBasis.set_active_element(e);
  PrecalcShapeset pssDirichlet(space->get_shapeset());
  pssDirichlet.set_active_element(e);
  pssBasis.set_active_shape(shape_indexBasis);
  pssDirichlet.set_active_shape(shape_indexDirichlet);
  Func<double>* fnBasis = init_fn(&pssBasis, &refmap, 24);
  Func<double>* fnDirichlet = init_fn(&pssDirichlet, &refmap, 24);
  
  return wf.get_mfvol()[0]->value(n_quadrature_points, jacobian_x_weights, NULL, fnBasis, fnDirichlet, geometry, NULL);
}