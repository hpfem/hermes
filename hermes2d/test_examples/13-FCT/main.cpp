#define HERMES_REPORT_ALL
#include "definitions.h"
#include "lumped_projection.h"
#include "hp_adapt.h"
#include "highOrder.h"
#include "lowOrder.h"
#include "fct.h"
#include "reg_estimator.h"
#include "prev_solution.h"
#include <list>


using namespace RefinementSelectors;
using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

//This example solves adaptively a linear convection equation.
//It describes a counterclockwise rotation about the center of the domain.
//The initial data is given as a slotted-cylinder, a hump and a cone.
//After each full revolution the exact solution coincides with the initial data.

// 1. Step: Low-order solution u^L :   (M_L/tau -theta(K+D)) u^L =   (M_L/tau + (1-theta)(K+D)) u^n
// 2. Step: antidiffusive fluxes:     f_ij = (M_c)_ij (dt_u_L(i)- dt_u_L(j)) + D_ij (u_L(i)- u_L(j)); 
//          limiting:    f_i = sum_(j!=i) alpha_ij f_ij
// 3. Step: explicit correction:  M_L u^(n+1) = M_L u^L + tau * f 



const int INIT_REF_NUM =3;                   // Number of initial refinements.
const int P_INIT = 1;       						// Initial polynomial degree.
const int P_MAX = 3; 										//Maximal polynomial degree.

const double time_step = 1e-3;                           // Time step.
const double T_FINAL = 2*PI;                       // Time interval length.

const double EPS_smooth = 1e-5;   		//constant for the smoothness indicator (a<b => a+eps<=b)
const double theta = 0.5;   			 // theta-scheme for time (theta =0 -> explizit, theta=1 -> implizit)

MatrixSolverType matrix_solver = SOLVER_UMFPACK; 


// Adaptivity
const double THRESHOLD_UNREF = 0.001; 							// Unrefinement: error of all sons is smaller than THRESHOLD_UNREF times maximum element error
const int UNREF_FREQ = 1;                         // Every UNREF_FREQth time step the mesh is derefined.
const int UNREF_METHOD = 1;                       // 1... mesh reset to basemesh and poly degrees to P_INIT.   
// 2... one ref. layer shaved off, poly degrees reset to P_INIT.
// 3... one ref. layer shaved off, poly degrees decreased by one. 
const double THRESHOLD = 0.75;                      // This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;                           // Adaptive strategy:
// STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
//   error is processed. If more elements have similar errors, refine
//   all to keep the mesh symmetric.
// STRATEGY = 1 ... refine all elements whose error is larger
//   than THRESHOLD times maximum element error.
// STRATEGY = 2 ... refine all elements whose error is larger
//   than THRESHOLD.
// More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO;          // Predefined list of element refinement candidates. Possible values are
// H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
// H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
// See the User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
// MESH_REGULARITY = -1 ... arbitrary level hanging nodes (default),
// MESH_REGULARITY = 1 ... at most one-level hanging nodes,
// MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
// Note that regular meshes are not supported, this is due to
// their notoriously bad performance.
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
// cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 0.1;                      // Stopping criterion for adaptivity (rel. error tolerance between the
// fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
// over this limit. This is to prevent h-adaptivity to go on forever.

const int ADAPSTEP_MAX = 5;												// max. numbers of adaptivity steps


//Visualization
const bool HERMES_VISUALIZATION = true;           // Set to "false" to suppress Hermes OpenGL visualization.
const bool VTK_VISUALIZATION = false;              // Set to "true" to enable VTK output.
const int VTK_FREQ = 1000;													//Every VTK_FREQth time step the solution is saved as VTK output.

// Boundary markers.
const std::string BDY_IN = "inlet";
const std::string BDY_OUT = "no_bdry";


int main(int argc, char* argv[])
{  
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh), basemesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", basemesh);

  /*  MeshView meshview("mesh", new WinGeom(0, 0, 500, 400));
  meshview.show(basemesh);
  View::wait();*/

  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_REF_NUM; i++)
    basemesh->refine_all_elements();
  mesh->copy(basemesh);

  // Initialize boundary conditions.
  DefaultEssentialBCConst<double>  bc_essential(BDY_IN, 0.0);
  EssentialBCs<double>  bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, &bcs, P_INIT));	

  // Initialize solution of lower & higher order
  MeshFunctionSharedPtr<double>  low_sln(new Solution<double>()), ref_sln(new Solution<double>()), high_sln(new Solution<double>()), sln(new Solution<double>());
  MeshFunctionSharedPtr<double>u_prev_time(new PrevSolution);

  // Previous time level solution (initialized by the initial condition).
  MeshFunctionSharedPtr<double>initial_condition(new CustomInitialCondition(mesh));

  // Initialize the weak formulation.
  CustomWeakFormMassmatrix  massmatrix(time_step);
  CustomWeakFormConvection  convection;

  // Initialize views.
  ScalarView sview("Solution", new WinGeom(0, 500, 500, 400));
  OrderView mview("mesh", new WinGeom(0, 0, 500, 400));
  OrderView ref_mview("ref_mesh", new WinGeom(500, 0, 500, 400));
  char title[100];

  // Output solution in VTK format.
  Linearizer lin;
  Orderizer ord;
  bool mode_3D = true;

  // Create a refinement selector.
  H1ProjBasedSelector<double> selector(CAND_LIST);

  //Initialize
  CSCMatrix<double> * mass_matrix = new CSCMatrix<double> ;   //M_c/tau
  CSCMatrix<double> * conv_matrix = new CSCMatrix<double> ;   //K
  double* u_L = NULL; 
  double* u_H =NULL;
  double* ref_sln_double =NULL;

  int ref_ndof, ndof; double err_est_rel_total;
  DefaultErrorCalculator<double, HERMES_L2_NORM> error_calculator(RelativeErrorToGlobalNorm, 1);

  Adapt<double> adaptivity(space, &error_calculator);
  AdaptStoppingCriterionCumulative<double> stoppingCriterion(THRESHOLD);
  adaptivity.set_strategy(&stoppingCriterion);

  OGProjection<double> ogProjection;	
  Lumped_Projection lumpedProjection;	
  Low_Order lowOrder(theta);
  High_Order highOrd(theta);
  Flux_Correction fluxCorrection(theta);
  Regularity_Estimator regEst(EPS_smooth);

  DiscreteProblem<double> dp_mass(&massmatrix, space);
  DiscreteProblem<double> dp_convection(&convection, space); 

  // Time stepping loop:
  double current_time = 0.0; 
  int ts = 1;
  do
  { 
    Hermes::Mixins::Loggable::Static::info("Time step %d, time %3.5f", ts, current_time);

    // Periodic global derefinement. 
    if ((ts > 1 && ts % UNREF_FREQ == 0)||(space->get_num_dofs() >= NDOF_STOP)) 
    { 
      Hermes::Mixins::Loggable::Static::info("Global mesh derefinement.");
      switch (UNREF_METHOD) {
      case 1: mesh->copy(basemesh);
        space->set_uniform_order(P_INIT);
        break;
      case 2: mesh->unrefine_all_elements();
        space->set_uniform_order(P_INIT);
        break;
      case 3: mesh->unrefine_all_elements();
        space->adjust_element_order(-1, -1, P_INIT, P_INIT);
        break;
      default: Exceptions::Exception("Wrong global derefinement method.");
      }      
    }

    bool done = false; int as = 1;	
    space->assign_dofs();

    do 
    {
      Hermes::Mixins::Loggable::Static::info("Time step %i, adap_step %i", ts, as);			
      ndof = space->get_num_dofs(); 

      ndof = space->get_num_dofs();  	    		
      double* coeff_vec_smooth = new double[ndof];
      int* smooth_elem_ref;	

      //smoothness-check for projected data		
      Hermes::Mixins::Loggable::Static::info("Projecting...");
      if(ts==1)
        ogProjection.project_global(space, initial_condition, coeff_vec_smooth, HERMES_L2_NORM);		
      else
        ogProjection.project_global(space, u_prev_time, coeff_vec_smooth, HERMES_L2_NORM);			

      Hermes::Mixins::Loggable::Static::info("Calling get_smooth_elems()...");
      smooth_elem_ref = regEst.get_smooth_elems(space,coeff_vec_smooth);
      // Construct reference mesh and setup reference space->
      MeshSharedPtr ref_mesh(new Mesh);
      ref_mesh->copy(space->get_mesh());
      Space<double>::ReferenceSpaceCreator ref_space_creator(space, ref_mesh, 0);
      SpaceSharedPtr<double> ref_space = ref_space_creator.create_ref_space();

      HPAdapt* adapting = new HPAdapt(ref_space);	

      // increase p in smooth regions, h refine in non-smooth regions 
      Hermes::Mixins::Loggable::Static::info("Calling adapt_smooth()...");
      if(adapting->adapt_smooth(smooth_elem_ref, P_MAX) == false) 
        throw Exceptions::Exception("reference space couldn't be constructed");							

      delete adapting;
      delete [] coeff_vec_smooth; 	 


      ref_ndof = ref_space->get_num_dofs();

      Hermes::Mixins::Loggable::Static::info("Visualization...");
      if(HERMES_VISUALIZATION) 
      {
        sprintf(title, "Ref_Mesh: Time %3.2f,timestep %i,as=%i,", current_time,ts,as);
        ref_mview.set_title(title);
        ref_mview.show(ref_space);
      }


      dp_mass.set_space(ref_space);
      dp_convection.set_space(ref_space);
      fluxCorrection.init(ref_space);	

      double* coeff_vec = new double[ref_ndof];
      double* coeff_vec_2 = new double[ref_ndof];
      double* limited_flux = new double[ref_ndof];	


      dp_mass.assemble(mass_matrix); 										//M_c/tau
      dp_convection.assemble(conv_matrix, NULL);		//K

      //----------------------MassLumping  & Artificial Diffusion --------------------------------------------------------------------	
      CSCMatrix<double>* lumped_matrix = fluxCorrection.massLumping(mass_matrix); // M_L/tau
      CSCMatrix<double>* diffusion = fluxCorrection.artificialDiffusion(conv_matrix);	

      //-----------------Assembling of matrices ---------------------------------------------------------------------	
      lowOrder.assemble_Low_Order(conv_matrix,diffusion,lumped_matrix);	
      highOrd.assemble_High_Order(conv_matrix,mass_matrix);

      mass_matrix->multiply_with_Scalar(time_step);  // massmatrix = M_C
      lumped_matrix->multiply_with_Scalar(time_step);  // M_L


      //--------- Project the previous timestep solution on the FE space (FCT is applied )----------------			
      // coeff_vec : FCT -Projection, coeff_vec_2: L2 Projection (ogProjection)	
      if(ts==1)
        fluxCorrection.project_FCT(initial_condition, coeff_vec, coeff_vec_2,mass_matrix,lumped_matrix,time_step,&ogProjection,&lumpedProjection, &regEst);	
      else		
        fluxCorrection.project_FCT(u_prev_time, coeff_vec, coeff_vec_2,mass_matrix,lumped_matrix,time_step,&ogProjection,&lumpedProjection, &regEst);			
      //------------------------- lower order solution------------					
      u_L = lowOrder.solve_Low_Order(lumped_matrix, coeff_vec,time_step);								
      //-------------high order solution (standard galerkin) ------				
      u_H = highOrd.solve_High_Order(coeff_vec);		
      //------------------------------Assemble antidiffusive fluxes and limit these-----------------------------------	
      fluxCorrection.antidiffusiveFlux(mass_matrix,lumped_matrix,conv_matrix,diffusion,u_H, u_L,coeff_vec, limited_flux,time_step,&regEst);
      //-------------Compute final solution ---------------			
      ref_sln_double = lowOrder.explicit_Correction(limited_flux);
      Solution<double>::vector_to_solution(ref_sln_double, ref_space, ref_sln);	

      // Project the fine mesh solution onto the coarse mesh.
      ogProjection.project_global(space, ref_sln, sln, HERMES_L2_NORM); 
      // Calculate element errors and total error estimate.
      error_calculator.calculate_errors(sln, ref_sln);
      double err_est_rel_total = error_calculator.get_total_error_squared() * 100;
      // Report results.
      Hermes::Mixins::Loggable::Static::info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%", ndof,ref_ndof, err_est_rel_total);				
      // If err_est_rel too large, adapt the mesh->
      if((err_est_rel_total < ERR_STOP)||(as>=ADAPSTEP_MAX)) done = true;
      else
      {
        done = adaptivity.adapt(&selector);
        // Increase the counter of performed adaptivity steps.
        if(done == false)  as++;
      }
      if(space->get_num_dofs() >= NDOF_STOP)
        done = true;

      if(done) 
      {  
        u_prev_time->copy(ref_sln);
      }

      // Visualize the solution and mesh->
      if(HERMES_VISUALIZATION)
      {
        sprintf(title, "Ref-Loesung: Time %3.2f,timestep %i,as=%i,", current_time,ts,as);
        sview.set_title(title);
        sview.show(ref_sln);
        sprintf(title, "Mesh: Time %3.2f,timestep %i,as=%i,", current_time,ts,as);
        mview.set_title(title);
        mview.show(space);
      }

      if((VTK_VISUALIZATION) &&((done==true)&&(ts  % VTK_FREQ == 0)))
      {
        // Output solution in VTK format.
        char filename[40];
        sprintf(filename, "solution-%i.vtk", ts );
        lin.save_solution_vtk(u_prev_time, filename, "solution", mode_3D);  
        sprintf(filename, "ref_space_order-%i.vtk", ts);
        ord.save_orders_vtk(ref_space, filename);
        sprintf(filename, "ref_mesh-%i.vtk", ts );
        ord.save_mesh_vtk(ref_space, filename);       
      } 

      // Clean up.
      delete lumped_matrix; 
      delete diffusion;
      delete [] coeff_vec_2;
      delete [] coeff_vec; 
      delete [] limited_flux; 
    }
    while (done == false);

    // Update global time.
    current_time += time_step;
    // Increase time step counter
    ts++;
  }
  while (current_time < T_FINAL); 

  // Visualize the solution.
  if(VTK_VISUALIZATION) {
    lin.save_solution_vtk(u_prev_time, "end_solution.vtk", "solution", mode_3D);
    ord.save_mesh_vtk(space, "end_mesh");
    ord.save_orders_vtk(space, "end_order.vtk");
  }

  delete mass_matrix;  
  delete conv_matrix;

  // Wait for the view to be closed.
  View::wait();
  return 0;
}

