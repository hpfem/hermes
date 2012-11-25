#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

// This example explains how to use the multimesh adaptive hp-FEM,
// where different physical fields (or solution components) can be
// approximated using different meshes and equipped with mutually
// independent adaptivity mechanisms. For the tutorial purposes,
// we manufactured an exact solution for a simplified version of
// the FitzHugh-Nagumo equation. This equation, in its full form,
// is a prominent example of activator-inhibitor systems in two-component
// reaction-diffusion equations, It describes a prototype of an
// excitable system (e.g., a neuron).
//
// PDE: Linearized FitzHugh-Nagumo equation
//      -d_u^2 \Delta u - f(u) + \sigma v - g1 = 0,
//      -d_v^2 \Delta v - u + v - g2 = 0.
// In the original equation, f(u) = \lambda u - u^3 - \kappa. For
// simplicity, here we just take f(u) = u.
//
// Domain: Square (-1,1)^2.
//
// BC: Both solution components are zero on the boundary.
//
// Exact solution: The functions g1 and g2 were calculated so that
//                 the exact solution is:
//        u(x,y) = U(x)*U(y) where U(t) = Hermes::cos(M_PI*t/2)
//        v(x,y) = V(x)V(y) where V(t) = 1 - (exp(K*t)+exp(-K*t))/(exp(K) + exp(-K))
// Note: V(t) is the exact solution of the 1D singularly perturbed equation
//       -u'' + K*K*u = K*K in (-1, 1) with zero Dirichlet BC.
//
// The following parameters can be changed: In particular, compare hp- and
// h-adaptivity via the CAND_LIST option, and compare the multi-mesh vs.
// single-mesh using the MULTI parameter.

// Initial polynomial degree for u.
const int P_INIT_U = 2;                           
// Initial polynomial degree for v.
const int P_INIT_V = 1;                           
// Number of initial boundary refinements
const int INIT_REF_BDY = 5;                       
// MULTI = true  ... use multi-mesh,
// MULTI = false ... use single-mesh.
// Note: In the single mesh option, the meshes are
// forced to be geometrically the same but the
// polynomial degrees can still vary.
const bool MULTI = true;                          
// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies.
const double THRESHOLD = 1.3;                     
// Adaptive strategy:
// STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
//   error is processed. If more elements have similar errors, refine
//   all to keep the mesh symmetric.
// STRATEGY = 1 ... refine all elements whose error is larger
//   than THRESHOLD times maximum element error.
// STRATEGY = 2 ... refine all elements whose error is larger
//   than THRESHOLD.
// More adaptive strategies can be created in adapt_ortho_h1.cpp.
const int STRATEGY = 0;                           
// Predefined list of element refinement candidates. Possible values are
// H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
// H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
const CandList CAND_LIST = H2D_HP_ANISO;              
// Maximum allowed level of hanging nodes:
// MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
// MESH_REGULARITY = 1 ... at most one-level hanging nodes,
// MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
// Note that regular meshes are not supported, this is due to
// their notoriously bad performance.
const int MESH_REGULARITY = -1;                   
// This parameter influences the selection of
// candidates in hp-adaptivity. Default value is 1.0. 
const double CONV_EXP = 1;                        
// Stopping criterion for adaptivity.
const double ERR_STOP = 5.91348;
// Adaptivity process stops when the number of degrees of freedom grows over
// this limit. This is mainly to prevent h-adaptivity to go on forever.
const int NDOF_STOP = 60000;
// Newton's method.
double NEWTON_TOL_FINE = 1e-0;
int NEWTON_MAX_ITER = 10;
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  

// Problem parameters.
const double D_u = 1;
const double D_v = 1;
const double SIGMA = 1;
const double LAMBDA = 1;
const double KAPPA = 1;
const double K = 100.;

int main(int argc, char* argv[])
{
  // Time measurement.
  Hermes::Mixins::TimeMeasurable cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh u_mesh, v_mesh;
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", &u_mesh);
  if (MULTI == false) u_mesh.refine_towards_boundary("Bdy", INIT_REF_BDY);

  // Create initial mesh (master mesh).
  v_mesh.copy(&u_mesh);

  // Initial mesh refinements in the v_mesh towards the boundary.
  if (MULTI == true) v_mesh.refine_towards_boundary("Bdy", INIT_REF_BDY);

  // Set exact solutions.
  ExactSolutionFitzHughNagumo1 exact_u(&u_mesh);
  ExactSolutionFitzHughNagumo2 exact_v(&v_mesh, K);

  // Define right-hand sides.
  CustomRightHandSide1 g1(K, D_u, SIGMA);
  CustomRightHandSide2 g2(K, D_v);

  // Initialize the weak formulation.
  CustomWeakForm wf(&g1, &g2);
  
  // Initialize boundary conditions
  DefaultEssentialBCConst<double> bc_u("Bdy", 0.0);
  EssentialBCs<double> bcs_u(&bc_u);
  DefaultEssentialBCConst<double> bc_v("Bdy", 0.0);
  EssentialBCs<double> bcs_v(&bc_v);

  // Create H1 spaces with default shapeset for both displacement components.
  H1Space<double> u_space(&u_mesh, &bcs_u, P_INIT_U);
  H1Space<double> v_space(MULTI ? &v_mesh : &u_mesh, &bcs_v, P_INIT_V);


  // Initialize coarse and reference mesh solutions.
  Solution<double> u_sln, v_sln, u_ref_sln, v_ref_sln;

  // Initialize refinement selector.
  H1ProjBasedSelector<double> selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  NewtonSolver<double> newton;

  // Adaptivity loop:
  int as = 1;
  bool done = false;
  do
  {
    // Construct globally refined reference mesh and setup reference space.
    Mesh::ReferenceMeshCreator u_ref_mesh_creator(&u_mesh);
    Mesh* u_ref_mesh = u_ref_mesh_creator.create_ref_mesh();
    Mesh::ReferenceMeshCreator v_ref_mesh_creator(&v_mesh);
    Mesh* v_ref_mesh = v_ref_mesh_creator.create_ref_mesh();
    Space<double>::ReferenceSpaceCreator u_ref_space_creator(&u_space, u_ref_mesh);
    Space<double>* u_ref_space = u_ref_space_creator.create_ref_space();
    Space<double>::ReferenceSpaceCreator v_ref_space_creator(&v_space, v_ref_mesh);
    Space<double>* v_ref_space = v_ref_space_creator.create_ref_space();

    Hermes::vector<const Space<double> *> ref_spaces_const(u_ref_space, v_ref_space);

    int ndof_ref = Space<double>::get_num_dofs(ref_spaces_const);

    // Perform Newton's iteration.
    try
    {
      newton.set_spaces(ref_spaces_const);

      newton.set_weak_formulation(&wf);

      newton.set_newton_tol(1e-1);

      newton.solve();
    }
    catch(Hermes::Exceptions::Exception& e)
    {
      std::cout << e.what();
    }
    catch(std::exception& e)
    {
      std::cout << e.what();
    }

    // Translate the resulting coefficient vector into the instance of Solution.
    Solution<double>::vector_to_solutions(newton.get_sln_vector(), ref_spaces_const,
                                          Hermes::vector<Solution<double> *>(&u_ref_sln, &v_ref_sln));

    // Project the fine mesh solution onto the coarse mesh.
    OGProjection<double> ogProjection; ogProjection.project_global(Hermes::vector<const Space<double> *>(&u_space, &v_space),
                                                                   Hermes::vector<Solution<double> *>(&u_ref_sln, &v_ref_sln),
                                                                   Hermes::vector<Solution<double> *>(&u_sln, &v_sln));

    // Calculate element errors.
    Adapt<double>* adaptivity = new Adapt<double>(Hermes::vector<Space<double> *>(&u_space, &v_space));

    // Calculate error estimate for each solution component and the total error estimate.
    Hermes::vector<double> err_est_rel;
    double err_est_rel_total = adaptivity->calc_err_est(Hermes::vector<Solution<double> *>(&u_sln, &v_sln),
                                                        Hermes::vector<Solution<double> *>(&u_ref_sln, &v_ref_sln),
                                                        &err_est_rel) * 100;

    // If err_est too large, adapt the mesh.
    if (err_est_rel_total < ERR_STOP)
      done = true;
    else
    {
      done = adaptivity->adapt(Hermes::vector<RefinementSelectors::Selector<double> *>(&selector, &selector),
                               THRESHOLD, STRATEGY, MESH_REGULARITY);
    }
    if (Space<double>::get_num_dofs(Hermes::vector<const Space<double> *>(&u_space, &v_space)) >= NDOF_STOP) done = true;

    // Clean up.
    delete adaptivity;

    // Increase counter.
    as++;
  }
  while (done == false);

	if(std::abs(u_ref_sln.get_pt_value(-0.98, -0.98)->val[0]- 0.000986633) > 1e-4) 
	{
		printf("Failure!\n");
		return -1;
	}
	if(std::abs(v_ref_sln.get_pt_value(-0.98, -0.98)->val[0]- 0.747675) > 1e-4) 
	{
		printf("Failure!\n");
		return -1;
	}
	if(std::abs(u_ref_sln.get_pt_value(-0.98, 0.98)->val[0]- 0.000986633) > 1e-4) 
	{
		printf("Failure!\n");
		return -1;
	}
	if(std::abs(v_ref_sln.get_pt_value(-0.98, 0.98)->val[0]- 0.747675) > 1e-4) 
	{
		printf("Failure!\n");
		return -1;
	}

	if(std::abs(u_ref_sln.get_pt_value(-0.98, -0.98)->dx[0]- 0.0493155) > 1e-4) 
	{
		printf("Failure!\n");
		return -1;
	}
	if(std::abs(v_ref_sln.get_pt_value(-0.98, -0.98)->dx[0]- 11.7667) > 1e-4) 
	{
		printf("Failure!\n");
		return -1;
	}
	if(std::abs(u_ref_sln.get_pt_value(-0.98, 0.98)->dx[0]- 0.0493155) > 1e-4) 
	{
		printf("Failure!\n");
		return -1;
	}
	if(std::abs(v_ref_sln.get_pt_value(-0.98, 0.98)->dx[0]- 11.7667) > 1e-4) 
	{
		printf("Failure!\n");
		return -1;
	}

	if(std::abs(u_ref_sln.get_pt_value(-0.98, -0.98)->dy[0]- 0.0493155) > 1e-4) 
	{
		printf("Failure!\n");
		return -1;
	}
	if(std::abs(v_ref_sln.get_pt_value(-0.98, -0.98)->dy[0]- 11.7667) > 1e-4) 
	{
		printf("Failure!\n");
		return -1;
	}
	if(std::abs(u_ref_sln.get_pt_value(-0.98, 0.98)->dy[0] + 0.0493155) > 1e-4) 
	{
		printf("Failure!\n");
		return -1;
	}
	if(std::abs(v_ref_sln.get_pt_value(-0.98, 0.98)->dy[0] + 11.7667) > 1e-4) 
	{
		printf("Failure!\n");
		return -1;
	}

	printf("Success!\n");
	return 0;
}
