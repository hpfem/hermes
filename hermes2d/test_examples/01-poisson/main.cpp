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
const bool BASE_VISUALIZATION = true;              // Set to "true" to enable base functions output.
const int P_INIT = 6;                             // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 6;                       // Number of initial uniform mesh refinements.

// Problem parameters.
const double VOLUME_HEAT_SRC = 1.0;          // Volume heat sources generated (for example) by electric current.
const double FIXED_BDY_TEMP = 0.0;					// Fixed temperature on the boundary.

Hermes::Hermes2D::ElementMode2D elementMode = HERMES_MODE_QUAD;
int form_i = 0;
int form_matrix_i = 0;
bool complexMesh = false;

/* form_matrix_i indices:
0 - no multiplying
1 - 0, 0
2 - 0, 1
3 - 1, 0
4 - 1, 1
5 - 1 * 1
6 - 1 * 2
7 - 1 * 3
8 - 1 * 4
9 - 2 * 1
10 - 2 * 2
11 - 2 * 3
12 - 2 * 4
13 - 3 * 1
14 - 3 * 2
15 - 3 * 3
16 - 3 * 4
17 - 4 * 1
18 - 4 * 2
19 - 4 * 3
20 - 4 * 4
*/

double matrixFunction(int np, double* wt, Func<double>* u, Func<double>* v)
{
	double result = 0;
  for (int i = 0; i < np; i++)
    switch (form_i)
    {
      // Vol.
      case 0:
        switch(form_matrix_i)
        {
        case 0:
          result += wt[i] * (u->val[i] * v->val[i]);
        break;
        }
      break;
      // DX.
      case 1:
        switch(form_matrix_i)
        {
        case 5:
          result += wt[i] * (u->dx[i] * v->dx[i]);
        break;
        case 6:
          result += wt[i] * (u->dx[i] * v->dy[i]);
        break;
        case 9:
          result += wt[i] * (u->dy[i] * v->dx[i]);
        break;
        case 10:
          result += wt[i] * (u->dy[i] * v->dy[i]);
        break;
        }
      break;
      // DY.
      case 2:
        switch(form_matrix_i)
        {
        case 15:
          result += wt[i] * (u->dx[i] * v->dx[i]);
        break;
        case 16:
          result += wt[i] * (u->dx[i] * v->dy[i]);
        break;
        case 19:
          result += wt[i] * (u->dy[i] * v->dx[i]);
        break;
        case 20:
          result += wt[i] * (u->dy[i] * v->dy[i]);
        break;
        }
      break;
      // DuDxValV.
      case 3:
        switch(form_matrix_i)
        {
        case 1:
          result += wt[i] * (u->dx[i] * v->val[i]);
        break;
        case 2:
          result += wt[i] * (u->dy[i] * v->val[i]);
        break;
        }
      break;
      // DuDyValV.
      case 4:
        switch(form_matrix_i)
        {
        case 3:
          result += wt[i] * (u->dx[i] * v->val[i]);
        break;
        case 4:
          result += wt[i] * (u->dy[i] * v->val[i]);
        break;
        }
      break;
      // ValUDvDx.
      case 5:
        switch(form_matrix_i)
        {
        case 1:
          result += wt[i] * (u->val[i] * v->dx[i]);
        break;
        case 2:
          result += wt[i] * (u->val[i] * v->dy[i]);
        break;
        }
      break;
      // ValUDvDy.
      case 6:
        switch(form_matrix_i)
        {
        case 3:
          result += wt[i] * (u->val[i] * v->dx[i]);
        break;
        case 4:
          result += wt[i] * (u->val[i] * v->dy[i]);
        break;
        }
      break;
    }
	return result;
}

double rhsValue = VOLUME_HEAT_SRC;

double rhsFunction(int np, double* wt, Func<double>* v)
{
	double result = 0;
  for (int i = 0; i < np; i++)
  switch (form_i)
    {
      // Vol.
      case 0:
      switch(form_matrix_i)
        {
        case 0:
          result += wt[i] * (v->val[i]);
        break;
        }
      break;
      // DX.
      case 1:
      switch(form_matrix_i)
        {
        case 1:
          result += wt[i] * (v->dx[i]);
        break;
        case 2:
          result += wt[i] * (v->dy[i]);
        break;
        }
      break;
      // DY.
      case 2:
      switch(form_matrix_i)
        {
        case 3:
          result += wt[i] * (v->dx[i]);
        break;
        case 4:
          result += wt[i] * (v->dy[i]);
        break;
        }
      break;
      default:
      result += 0.0;
    }
	return result;
}

// Globals.
Hermes::Hermes2D::H1Space<double>* space;
int max_index;
Hermes::Hermes2D::Mesh* mesh;

// Storage.
double** A_values;
double* rhs_values;

// Cache loading.
void loadCache();

// Problem loading.
void loadProblemData();

// Cache calculation.
void calculateCache(CustomWeakFormPoissonCacheCalculation& wf, Shapeset* shapeset);

void calculateResultWithConstantForms(CustomWeakFormPoisson& wf);

void calculateResultAssembling(CustomWeakFormPoissonCacheCalculation& wf);

int main(int argc, char* argv[])
{
  try
  {
    // Initialize the weak formulation.
    CustomWeakFormPoissonCacheCalculation wf;
    CustomWeakFormPoisson wf1;

    //Hermes2DApi.set_integral_param_value(numThreads, 1);
    
    /*
    for(int i = 0; i < 2; i++)
      for(int j = 0; j < 6; j++)
        {
          elementMode = i == 0 ? HERMES_MODE_TRIANGLE : HERMES_MODE_QUAD;
          form_i = j;
          calculateCache(wf, new H1Shapeset);
        }

    for(int i = 0; i < 2; i++)
      for(int j = 0; j < 6; j++)
        {
          elementMode = i == 0 ? HERMES_MODE_TRIANGLE : HERMES_MODE_QUAD;
          form_i = j;
          calculateCache(wf, new L2Shapeset);
        }
        */

    loadProblemData();
    calculateResultWithConstantForms(wf1);
    calculateResultAssembling(wf);
  }
  catch(std::exception& e)
  {
    std::cout << e.what();
  }
}

void loadCache()
{
	std::stringstream ssMatrix;
	std::stringstream ssRhs;
	ssMatrix << "matrixCache";
	ssRhs << "rhsCache";
	if(elementMode == HERMES_MODE_TRIANGLE)
	{
		ssMatrix << "Triangle";
		ssRhs << "Triangle";
	}
	else
	{
		ssMatrix << "Quad";
		ssRhs << "Quad";
	}

	std::ifstream matrixFormIn;
  matrixFormIn.open(ssMatrix.str().c_str());
	std::ifstream rhsFormIn;
  rhsFormIn.open(ssRhs.str().c_str());

	int index_i, index_j;
	int counter = 0;
  double valueTemp;
  while(matrixFormIn.good())
  {
    matrixFormIn >> index_i >> index_j >> valueTemp;
    A_values[index_i][index_j] = valueTemp;
    counter++;
  }

  counter = 0;
  while(rhsFormIn.good())
  {
    rhsFormIn >> index_i >> valueTemp;
    rhs_values[index_i] = valueTemp;
    counter++;
  }
}

void loadProblemData()
{
  // Load the mesh.
  Hermes::Hermes2D::MeshReaderH2DXML mXMLloader;
  Hermes::Hermes2D::MeshReaderH2D mloader;
  mesh = new Mesh();
	std::stringstream ss;
	ss << "domain";
  if(complexMesh)
  {
	  if(elementMode == HERMES_MODE_TRIANGLE)
		  ss << "-tri";
	  else
		  ss << "-quad";
	  ss << ".mesh";
    mloader.load(ss.str().c_str(), mesh);
  }
  else
  {
    if(elementMode == HERMES_MODE_TRIANGLE)
		  ss << "Triangle";
	  else
		  ss << "Quad";
	  ss << ".xml";
    try
    {

      mXMLloader.load(ss.str().c_str(), mesh);
    }
    catch(std::exception& e)
    {
      std::cout << e.what();
    } 
  }
	
	for(int i = 0; i < INIT_REF_NUM; i++)
    mesh->refine_all_elements();

  // Create a H1 space with default shapeset.
  Hermes::Hermes2D::DefaultEssentialBCConst<double>* bc_essential = new Hermes::Hermes2D::DefaultEssentialBCConst<double>(HERMES_ANY, FIXED_BDY_TEMP);
  Hermes::Hermes2D::EssentialBCs<double>* bcs = new Hermes::Hermes2D::EssentialBCs<double>(bc_essential);

  // Create an H1 space with default shapeset.
  space = new Hermes::Hermes2D::H1Space<double>(mesh, bcs, P_INIT);
}

void calculateCache(CustomWeakFormPoissonCacheCalculation& wf, Shapeset* shapeset)
{
  // Load the mesh.
  Hermes::Hermes2D::MeshReaderH2DXML mloader;
  mesh = new Mesh();
	std::stringstream ss;
	ss << "domainCacheComputation";
	if(elementMode == HERMES_MODE_TRIANGLE)
		ss << "Triangle";
	else
		ss << "Quad";
	ss << ".xml";

	try
	{
		mloader.load(ss.str().c_str(), mesh);
	}
	catch(std::exception& e)
	{
		std::cout << e.what();
	}

	std::stringstream ssMatrix;
	std::stringstream ssRhs;
	
  if(form_i == 0)
  {
    ssMatrix << "MatrixFormVol";
    ssRhs << "VectorFormVol";
  }
  if(form_i == 1)
  {
    ssMatrix << "MatrixFormDx";
	  ssRhs << "VectorFormDx";
  }
  if(form_i == 2)
  {
    ssMatrix << "MatrixFormDy";
	  ssRhs << "VectorFormDy";
  }
  if(form_i == 3)
  {
    ssMatrix << "MatrixFormDuDxValV";
	  ssRhs << "xxx";
  }
  if(form_i == 4)
  {
    ssMatrix << "MatrixFormDuDyValV";
	  ssRhs << "xxx";
  }
  if(form_i == 5)
  {
    ssMatrix << "MatrixFormValUDvDx";
	  ssRhs << "xxx";
  }
  if(form_i == 6)
  {
    ssMatrix << "MatrixFormValUDvDy";
	  ssRhs << "xxx";
  }

	if(elementMode == HERMES_MODE_TRIANGLE)
	{
		ssMatrix << "Triangle";
		ssRhs << "Triangle";
	}
	else
	{
		ssMatrix << "Quad";
		ssRhs << "Quad";
	}
  if(shapeset->get_space_type() == HERMES_H1_SPACE)
  {
		ssMatrix << ".h1h1";
		ssRhs << ".h1";
	}
	else
	{
		ssMatrix << ".l2l2";
		ssRhs << ".l2";
	}

	std::ofstream matrixFormOut;
  matrixFormOut.open(ssMatrix.str().c_str());
  matrixFormOut.precision(20);
	std::ofstream rhsFormOut;
  rhsFormOut.open(ssRhs.str().c_str());
  rhsFormOut.precision(20);

  for(int mf_i = 0; mf_i < 21; mf_i++)
    switch (form_i)
    {
      // Vol.
      case 0:
        switch(mf_i)
        {
        case 0:
          matrixFormOut << 'a';
        default:
          matrixFormOut << 'n';
        break;
        }
      break;
      // DX.
      case 1:
        switch(mf_i)
        {
        case 5:
          matrixFormOut << 'a';
        break;
        case 6:
          matrixFormOut << 'a';
        break;
        case 9:
          matrixFormOut << 'a';
        break;
        case 10:
          matrixFormOut << 'a';
        break;
        default:
          matrixFormOut << 'n';
        }
      break;
      // DY.
      case 2:
        switch(mf_i)
        {
        case 15:
          matrixFormOut << 'a';
        break;
        case 16:
          matrixFormOut << 'a';
        break;
        case 19:
          matrixFormOut << 'a';
        break;
        case 20:
          matrixFormOut << 'a';
        break;
        default:
          matrixFormOut << 'n';
        }
      break;
      // DuDxValV.
      case 3:
        switch(mf_i)
        {
        case 1:
          matrixFormOut << 'a';
        break;
        case 2:
          matrixFormOut << 'a';
        break;
        default:
          matrixFormOut << 'n';
        }
      break;
      // DuDyValV.
      case 4:
        switch(form_matrix_i)
        {
        case 3:
          matrixFormOut << 'a';
        break;
        case 4:
          matrixFormOut << 'a';
        break;
        default:
          matrixFormOut << 'n';
        }
      break;
      // ValUDvDx.
      case 5:
        switch(form_matrix_i)
        {
        case 1:
          matrixFormOut << 'a';
        break;
        case 2:
          matrixFormOut << 'a';
        break;
        default:
          matrixFormOut << 'n';
        }
      break;
      // ValUDvDy.
      case 6:
        switch(form_matrix_i)
        {
        case 3:
          matrixFormOut << 'a';
        break;
        case 4:
          matrixFormOut << 'a';
        break;
        default:
          matrixFormOut << 'n';
        }
      break;
    }

  for(int mf_i = 0; mf_i < 5; mf_i++)
    switch (form_i)
    {
      // Vol.
      case 0:
        switch(mf_i)
        {
        case 0:
          rhsFormOut << 'a';
        break;
        default:
          rhsFormOut << 'n';
        }
      break;
      // DX.
      case 1:
      switch(mf_i)
        {
        case 1:
          rhsFormOut << 'a';
        break;
        case 2:
          rhsFormOut << 'a';
        break;
        default:
          rhsFormOut << 'n';
        }
      break;
      // DY.
      case 2:
      switch(mf_i)
        {
        case 3:
          rhsFormOut << 'a';
        break;
        case 4:
          rhsFormOut << 'a';
        break;
        default:
          rhsFormOut << 'n';
        }
      break;
    }

  matrixFormOut << std::endl;
  rhsFormOut << std::endl;

  max_index = shapeset->get_max_index(elementMode) + 1;
	PrecalcShapeset uP(shapeset);
	PrecalcShapeset vP(shapeset);
	uP.set_active_element(mesh->get_element(0));
	vP.set_active_element(mesh->get_element(0));
  
	int np = g_quad_2d_std.get_num_points(elementMode == HERMES_MODE_TRIANGLE ? 20 : 24, elementMode);
	double3* pt = g_quad_2d_std.get_points(elementMode == HERMES_MODE_TRIANGLE ? 20 : 24, elementMode);
	double* jwt = new double[np];
	for(int i = 0; i < np; i++)
		jwt[i] = pt[i][2];

	RefMap rm;
	rm.set_active_element(mesh->get_element(0));
	
	for(int i = 0; i < max_index; i++)
  {
    if(i > 0)
      rhsFormOut << std::endl;
		vP.set_active_shape(i);
		Func<double>* v = init_fn(&vP, &rm, (elementMode == HERMES_MODE_TRIANGLE ? 20 : 24));

    for(int j = 0; j < max_index; j++)
    {
      if(i > 0 || j > 0)
        matrixFormOut << std::endl;
			uP.set_active_shape(j);
			Func<double>* u = init_fn(&uP, &rm, (elementMode == HERMES_MODE_TRIANGLE ? 20 : 24));
      matrixFormOut << i << ' ' << j;
      for(form_matrix_i = 0; form_matrix_i < 21; form_matrix_i++)
        matrixFormOut << ' ' << matrixFunction(np, jwt, u, v);
    }
    rhsFormOut << i;
    for(form_matrix_i = 0; form_matrix_i < 5; form_matrix_i++)
        rhsFormOut << ' ' << rhsFunction(np, jwt, v);
  }

  matrixFormOut.close();
  rhsFormOut.close();

  return;
}

void calculateResultWithConstantForms(CustomWeakFormPoisson& wf)
{
  Views::BaseView<double> m;
  m.show(space);
  m.wait_for_close();


  // Utilities.
  int ndof = space->get_num_dofs();

  // Initialize the solution.
  Hermes::Hermes2D::Solution<double> sln;

  /// The solution vector.
  double* sln_vector;

	Hermes::Mixins::TimeMeasurable time;
	time.tick();

  SparseMatrix<double>* jacobian = create_matrix<double>();
  Vector<double>* residual = create_vector<double>();
	DiscreteProblemLinear<double> dp(&wf, space);
	dp.assemble(jacobian, residual);

	time.tick();
	std::cout << (std::string)"Ndofs: " << ndof << '.' << std::endl;
	std::cout << (std::string)"Assembling with constant forms: " << time.last() << '.' << std::endl;

  LinearMatrixSolver<double>* matrix_solver = create_linear_solver<double>(jacobian, residual);
	FILE* matrixFile = fopen("matrix", "w");
	FILE* rhsFile = fopen("rhs", "w");
	jacobian->dump(matrixFile, "A");
	residual->dump(rhsFile, "b");
  fclose(matrixFile);
  fclose(rhsFile);
	time.tick();
  try
  {
	  matrix_solver->solve();
  }
  catch(std::exception& e)
  {
    std::cout << e.what();
  }
	time.tick();
	std::cout << (std::string)"Solving: " << time.last() << '.' << std::endl;
  sln_vector = matrix_solver->get_sln_vector();

  // Translate the solution vector into the previously initialized Solution.
  Hermes::Hermes2D::Solution<double>::vector_to_solution(sln_vector, space, &sln, false);

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
    viewS.show(&sln);
    viewS.wait_for_close();
  }

  return;
}

void calculateResultAssembling(CustomWeakFormPoissonCacheCalculation& wf)
{
  // Utilities.
  int ndof = space->get_num_dofs();

  // Initialize the solution.
  Hermes::Hermes2D::Solution<double> sln;

  /// The solution vector.
  double* sln_vector;

	Hermes::Mixins::TimeMeasurable time;
	time.tick();

  SparseMatrix<double>* jacobian = create_matrix<double>();
  Vector<double>* residual = create_vector<double>();
	DiscreteProblemLinear<double> dp(&wf, space);
	dp.assemble(jacobian, residual);

	time.tick();
	std::cout << (std::string)"Ndofs: " << ndof << '.' << std::endl;
	std::cout << (std::string)"Assembling without constant forms: " << time.last() << '.' << std::endl;

  LinearMatrixSolver<double>* matrix_solver = create_linear_solver<double>(jacobian, residual);
	FILE* matrixFile = fopen("matrix", "w");
	FILE* rhsFile = fopen("rhs", "w");
	jacobian->dump(matrixFile, "A");
	residual->dump(rhsFile, "b");
  fclose(matrixFile);
  fclose(rhsFile);
  try
  {
	  matrix_solver->solve();
  }
  catch(std::exception& e)
  {
    std::cout << e.what();
  }
  sln_vector = matrix_solver->get_sln_vector();

  // Translate the solution vector into the previously initialized Solution.
  Hermes::Hermes2D::Solution<double>::vector_to_solution(sln_vector, space, &sln, false);

  return;
}