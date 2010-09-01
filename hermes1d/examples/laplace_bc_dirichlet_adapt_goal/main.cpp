#include "hermes1d.h"

// ********************************************************************

// This example uses hp-adaptivity to solve the Poisson equation 
// -u'' - f = 0 in an interval (A, B), equipped with Dirichlet 
// boundary conditions on both end points. A series of small reference 
// solutions (we call them fast trial refinements, FTR) is used both 
// to decide what elements will be refined, and how they will be refined.
// One has to define a goal of computation (quantity of interest) which 
// is any linear functional of the solution. Adaptivity can be guided either
// by the error in global norm or by the error in the quantity of interest. 
// One can use either standard Newton's method or JFNK to solve discrete 
// problems for trial refinements.

// General input:
const int N_eq = 1;
const int N_elem = 30;                   // Number of elements
//const double A = 0, B = 1;              // Domain end points
const double A = -M_PI, B = M_PI;      // Domain end points
const double EPSILON = 0.01;            // Equation parameter
const int P_init = 1;                   // Initial polynomal degree

// JFNK or classical Newton?
const double MATRIX_SOLVER_TOL = 1e-5;
const int MATRIX_SOLVER_MAXITER = 150;
const int JFNK = 0;                     // 0... classical Newton
                                        // 1... JFNK  
                                        // Only relevant for JFNK:
const double JFNK_EPSILON = 1e-4;       // Parameter in the JFNK finite difference
 
// Newton's method
const double NEWTON_TOL_COARSE = 1e-8;  // Coarse mesh
const double NEWTON_TOL_REF = 1e-8;     // Reference mesh
const int NEWTON_MAXITER = 150;

// Adaptivity
const int GOAL_ORIENTED = 0;            // 0... standard adaptivity in norm
                                        // 1... goal-oriented adaptivity 
const double X_QOI = 0.9*B;             // Value of u[0] at X_QOI is the 
                                        // Quantity of interest
const double TOL_ERR_QOI = 1e-8;        // Tolerance for the maximum FTR error
const int ADAPT_TYPE = 0;               // 0... hp-adaptivity
                                        // 1... h-adaptivity
                                        // 2... p-adaptivity
const double THRESHOLD = 0.7;           // Refined will be all elements whose FTR error
                                        // is greater than THRESHOLD*max_ftr_error
const int NORM = 1;                     // To measure errors:
                                        // 1... H1 norm
                                        // 0... L2 norm

// Boundary conditions
//double Val_dir_left = pow(A + EPSILON, 1./3.);   // Dirichlet condition left
//double Val_dir_right = pow(B + EPSILON, 1./3.);  // Dirichlet condition right
double Val_dir_left = 0;                           // Dirichlet condition left
double Val_dir_right = 0;                          // Dirichlet condition right

// Function f(x)
double f(double x) {
  return sin(x);
  //return 2;
  //return 2./9. * pow(x + EPSILON, -5./3.);
}

// Exact solution:
// When changing exact solution, do not 
// forget to update interval accordingly
const int EXACT_SOL_PROVIDED = 1;
double exact_sol(double x, double u[MAX_EQN_NUM], double dudx[MAX_EQN_NUM]) {
  u[0] = sin(x);
  dudx[0] = cos(x);
  //u[0] = 1. - x*x;
  //dudx[0] = -2.*x;
  //u[0] = pow(x + EPSILON, 1./3.);
  //dudx[0] = 1./3. * pow(x + EPSILON, -2./3.);
}

// Weak forms for Jacobi matrix and residual
#include "forms.cpp"

/******************************************************************************/

// Sample quantity of interest - solution value u[0] at 'x'.
// NOTE: quantity of interest is any linear functional of solution.
double quantity_of_interest(Mesh *mesh, double x)
{
  // Traversal of the mesh  
  Iterator *I = new Iterator(mesh);
  Element *e;
  double val[MAX_EQN_NUM];
  double der[MAX_EQN_NUM];
  while ((e = I->next_active_element()) != NULL) {
    if (e->x1 <= x && x <= e->x2) {
      e->get_solution_point(x, val, der);
      return val[0];
    }
  }

  error("computation of quantity of interest failed.");
}

/******************************************************************************/
int main() {
  // Create coarse mesh, set Dirichlet BC, enumerate 
  // basis functions
  Mesh *mesh = new Mesh(A, B, N_elem, P_init, N_eq);
  mesh->set_bc_left_dirichlet(0, Val_dir_left);
  mesh->set_bc_right_dirichlet(0, Val_dir_right);
  mesh->assign_dofs();

  // Create discrete problem on coarse mesh
  DiscreteProblem *dp = new DiscreteProblem();
  dp->add_matrix_form(0, 0, jacobian);
  dp->add_vector_form(0, residual);

  // Convergence graph wrt. the number of degrees of freedom
  // (goal-oriented adaptivity)
  GnuplotGraph graph_ftr;
  graph_ftr.set_log_y();
  graph_ftr.set_captions("Convergence History", "Degrees of Freedom", "QOI error");
  graph_ftr.add_row("QOI error - FTR (exact)", "k", "-", "o");
  graph_ftr.add_row("QOI error - FTR (est)", "k", "--");

  // Main adaptivity loop
  int adapt_iterations = 1;
  double ftr_errors[MAX_ELEM_NUM];      // This array decides what 
                                         // elements will be refined.
  ElemPtr2 ref_ftr_pairs[MAX_ELEM_NUM]; // To store element pairs from the 
                                         // FTR solution. Decides how 
                                         // elements will be hp-refined. 
  for (int i=0; i < MAX_ELEM_NUM; i++) {
    ref_ftr_pairs[i][0] = new Element();
    ref_ftr_pairs[i][1] = new Element();
  }
  while(1) {
    printf("============ Adaptivity step %d ============\n", adapt_iterations); 

    printf("N_dof = %d\n", mesh->get_n_dof());
 
    // Newton's loop on coarse mesh
    int success;
    if(JFNK == 0) {
      newton(dp, mesh, NULL, NEWTON_TOL_COARSE, NEWTON_MAXITER);
    }
    else {
      jfnk_cg(dp, mesh, MATRIX_SOLVER_TOL, MATRIX_SOLVER_MAXITER,
              JFNK_EPSILON, NEWTON_TOL_COARSE, NEWTON_MAXITER);
    }
    // For every element perform its fast trial refinement (FTR),
    // calculate the norm of the difference between the FTR
    // solution and the coarse mesh solution, and store the
    // error in the ftr_errors[] array.
    int n_elem = mesh->get_n_active_elem();
    double max_qoi_err_est = 0;
    for (int i=0; i < n_elem; i++) {

      printf("=== Starting FTR of Elem [%d]\n", i);

      // Replicate coarse mesh including solution.
      Mesh *mesh_ref_local = mesh->replicate();

      // Perform FTR of element 'i'
      mesh_ref_local->reference_refinement(i, 1);
      printf("Elem [%d]: fine mesh created (%d DOF).\n", 
             i, mesh_ref_local->assign_dofs());

      // Newton's loop on the FTR mesh
      if(JFNK == 0) {
        newton(dp, mesh_ref_local, NULL, NEWTON_TOL_COARSE, NEWTON_MAXITER);
      }
      else {
        jfnk_cg(dp, mesh_ref_local, MATRIX_SOLVER_TOL, MATRIX_SOLVER_MAXITER, 
                JFNK_EPSILON, NEWTON_TOL_REF, NEWTON_MAXITER);
      }

      // Print FTR solution (enumerated) 
      Linearizer *lxx = new Linearizer(mesh_ref_local);
      char out_filename[255];
      sprintf(out_filename, "solution_ref_%d.gp", i);
      lxx->plot_solution(out_filename);
      delete lxx;

      // Calculate FTR errors for refinement purposes
      if (GOAL_ORIENTED == 1) {
        // Use quantity of interest.
        double qoi_est = quantity_of_interest(mesh, X_QOI);
        double qoi_ref_est = quantity_of_interest(mesh_ref_local, X_QOI);
        ftr_errors[i] = fabs(qoi_ref_est - qoi_est);
      }
      else {
        // Use global norm
        double err_est_array[MAX_ELEM_NUM];
        ftr_errors[i] = calc_error_estimate(NORM, mesh, mesh_ref_local, 
                                            err_est_array);
      }

      // Calculating maximum of QOI FTR error for plotting purposes
      if (GOAL_ORIENTED == 1) {
        if (ftr_errors[i] > max_qoi_err_est) 
	  max_qoi_err_est = ftr_errors[i];
      }
      else {
        double qoi_est = quantity_of_interest(mesh, X_QOI);
        double qoi_ref_est = quantity_of_interest(mesh_ref_local, X_QOI);
        double err_est = fabs(qoi_ref_est - qoi_est);
        if (err_est > max_qoi_err_est) 
	  max_qoi_err_est = err_est;
      }

      // Copy the reference element pair for element 'i'
      // into the ref_ftr_pairs[i][] array
      Iterator *I = new Iterator(mesh);
      Iterator *I_ref = new Iterator(mesh_ref_local);
      Element *e, *e_ref;
      while (1) {
        e = I->next_active_element();
        e_ref = I_ref->next_active_element();
        if (e->id == i) {
  	  e_ref->copy_into(ref_ftr_pairs[e->id][0]);
          // coarse element 'e' was split in space
          if (e->level != e_ref->level) {
            e_ref = I_ref->next_active_element();
            e_ref->copy_into(ref_ftr_pairs[e->id][1]);
          }
          break;
        }
      }

      delete I;
      delete I_ref;
      delete mesh_ref_local;
    }  

    // Add entries to convergence graphs
    if (EXACT_SOL_PROVIDED) {
      double qoi_est = quantity_of_interest(mesh, X_QOI);
      double u[MAX_EQN_NUM], dudx[MAX_EQN_NUM];
      exact_sol(X_QOI, u, dudx);
      double err_qoi_exact = fabs(u[0] - qoi_est);
      // plotting error in quantity of interest wrt. exact value
      graph_ftr.add_values(0, mesh->get_n_dof(), err_qoi_exact);
    }
    graph_ftr.add_values(1, mesh->get_n_dof(), max_qoi_err_est);

    // Decide whether the max. FTR error in the quantity of interest 
    // is sufficiently small
    if(max_qoi_err_est < TOL_ERR_QOI) break;

    // debug
    if (adapt_iterations == 3) break;

    // Returns updated coarse mesh with the last solution on it. 
    adapt(NORM, ADAPT_TYPE, THRESHOLD, ftr_errors,
          mesh, ref_ftr_pairs);

    adapt_iterations++;
  }

  // Plot meshes, results, and errors
  adapt_plotting(mesh, ref_ftr_pairs,
                 NORM, EXACT_SOL_PROVIDED, exact_sol);

  // Save convergence graph
  graph_ftr.save("conv_dof.gp");

  printf("Done.\n");
  return 1;
}









