#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

/** \addtogroup t_cand_proj Candidate Projection (Zero Error)
*  \{
*  \brief This test tests projection of a candidate in H1 space on a quad.
*
*  For each combination of orders it creates a polynom \f$ p(x,y) = \sum_h \sum_v a_{h,v} x^h y^v \f$.
*  The test ensures that all coefficients \f$a_{h,v}\f$ are not zero. The polynom is defined
*  in the reference domain.
*  Then, it creates a mesh which consist of a single element. The order of the element is
*  \f$(h-1, v-1)\f$ where \f$h\f$ is a current horizontal order and \f$v\f$ is current vertical order.
*  Using the mesh and the polynom it calculates a solution and reference solution. The reference solution
*  is used the control the selection of candidates. Using projection based selector,
*  the test generates HP_ANISO candidates and calculates their errors.
*
*  The test succeeds if:
*   - Values of reference solution equals to values of the original function in random points.
*   - Errors of all candidates whose orders of all sons are greater or equal to order of the polynom are zero.
*
*  Output table legend:
*   - '-': Test was not done
*   - ' ': Test succesfull
*   - 'S': Values of the reference solution differs from values of the test function despite that the order of the reference solution is either the same of higher than the test function.
*   - 'C': Some candidates has non-zero error despite that order of their sons is greater or equal to the order of the test function.
*/

#include "functions.h"

/* global definitions */
#define H2D_TEST_ELEM_ID 0 ///< ID of an alement which is going to be handled.
#define H2D_TEST_ZERO_QUADS 4e-12 ///< Numerical zero: quads. Since polynoms are defined on a reference domain, some high-order polynoms might yield error a little bit higher than \f$10^{-12}\f$ in H1. In L2, high-order polynomials yields an error higher than \f$3\times 10^{-12}\f$.
#define H2D_TEST_ZERO_TRIS 6e-8 ///< Numerical zero: triangles. Since polynoms are defined on a reference domain, some high-order polynoms might yield error a little bit higher than \f$10^{-12}\f$ in H1. In L2, high-order polynomials yields an error higher than \f$6\times 2^{-8}\f$.
#define delete_not_null(__ptr) if (__ptr != NULL) delete __ptr; ///< Deletes an instance if the pointer is not NULL.

#define H2D_TEST_NOT_DONE -1 ///< Flag: Combination of orders was not tested.
#define H2D_TEST_SUCCESS 0x00 ///< Flag: Combination of orders was tested successfully.
#define H2D_CAND_FAILED 0x01  ///< Flag: Appropriate candidates have non-zero error.
#define H2D_RSLN_FAILED 0x02  ///< Flag: failure of the condition of equality between ref. solution failed and the polynom.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

/* global variables */
Mesh* mesh = NULL; ///< Mesh used by the test.
Shapeset* shapeset = NULL; ///< Shapeset used by the test.
Space<double>* space = NULL; ///< Space used by the test.
WeakForm<double>* weakform = NULL; ///< Weakform used by the test.
OptimumSelector<double>* selector = NULL; ///< Selector used by the test.
double numerical_zero = 0.0; ///< A numberical zero.

TestCase* cur_test_case = NULL; ///< Current test case: required by callbacks.

class H1TestWeakForm : public WeakForm<double>
{
public:
  H1TestWeakForm()
  {
    add_matrix_form(new BiForm(0, 0));
    add_vector_form(new LiForm(0));
  }

  class BiForm : public MatrixFormVol<double>
  {
  public:
    BiForm(int i, int j) : MatrixFormVol<double>(i, j) {};

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
      Geom<double> *e, ExtData<double> *ext) const
    {
      return int_u_v<double, double>(n, wt, u, v) + int_grad_u_grad_v<double, double>(n, wt, u, v);
    }

    Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u, Func<Hermes::Ord> *v,
      Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const
    {
      return int_u_v<Ord, Ord>(n, wt, u, v) + int_grad_u_grad_v<Ord, Ord>(n, wt, u, v);
    }
  };

  class LiForm : public VectorFormVol<double>
  {
  public:
    LiForm(int i) : VectorFormVol<double>(i) {};

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
      Geom<double> *e, ExtData<double> *ext) const
    {
      double result = 0;
      for(int i = 0; i < n; i++) {
        double x = e->x[i], y = e->y[i];
        result += wt[i] * (v->val[i] * func_val(cur_test_case->poly_matrix(), cur_test_case->quad_order(), x, y)
          + v->dx[i] * func_dx(cur_test_case->poly_matrix(), cur_test_case->quad_order(), x, y)
          + v->dy[i] * func_dy(cur_test_case->poly_matrix(), cur_test_case->quad_order(), x, y));
      }
      return result;
    }

    Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
      Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const
    {
      return Ord(H2D_GET_H_ORDER(cur_test_case->quad_order()) + H2D_GET_V_ORDER(cur_test_case->quad_order()) + 2*v->val->get_order());
    }
  };
};

class L2TestWeakForm : public WeakForm<double>
{
public:
  L2TestWeakForm()
  {
    add_matrix_form(new BiForm(0, 0));
    add_vector_form(new LiForm(0));
  }

  class BiForm : public MatrixFormVol<double>
  {
  public:
    BiForm(int i, int j) : MatrixFormVol<double>(i, j) {};

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
      Geom<double> *e, ExtData<double> *ext) const
    {
      return int_u_v<double, double>(n, wt, u, v);
    }

    Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u, Func<Hermes::Ord> *v,
      Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const
    {
      return int_u_v<Ord, Ord>(n, wt, u, v);
    }
  };

  class LiForm : public VectorFormVol<double>
  {
  public:
    LiForm(int i) : VectorFormVol<double>(i) {};

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
      Geom<double> *e, ExtData<double> *ext) const
    {
      double result = 0;
      for(int i = 0; i < n; i++) 
      {
        double x = e->x[i], y = e->y[i];
        result += wt[i] * v->val[i] * func_val(cur_test_case->poly_matrix(), cur_test_case->quad_order(), x, y);
      }
      return result;
    }

    Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
      Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const
    {
      return Ord(H2D_GET_H_ORDER(cur_test_case->quad_order()) + H2D_GET_V_ORDER(cur_test_case->quad_order()) + 2*v->val->get_order());
    }
  };
};

/// Function that is used to exact solution.
double h1_exact_func(double x, double y, double& dx, double& dy) {
  double value = func_val(cur_test_case->poly_matrix(), cur_test_case->quad_order(), x, y);
  dx = func_dx(cur_test_case->poly_matrix(), cur_test_case->quad_order(), x, y);
  dy = func_dy(cur_test_case->poly_matrix(), cur_test_case->quad_order(), x, y);
  return value;
}

/// Cleanup.
void cleanup() {
  delete_not_null(selector);
  delete_not_null(weakform);
  delete_not_null(space);
  delete_not_null(mesh);
  delete_not_null(shapeset);
}

/// Initialize mesh.
Mesh* init_mesh(bool tri) {
  Mesh* mesh = new Mesh();
  if (tri) {
    const int vertex_num = 3, tria_num = 1, quad_num = 0, marker_num = 3;
    double2 vertex_array[3] = {{-1,-1}, { 1,-1}, {-1, 1}};
    int4 tria_array[1] = {{0, 1, 2, 0}};
    int5 *quad_array = NULL;
    int3 marker_array[3] = {{0, 1, 2}, {1, 2, 1}, {2, 0, 3}};
    mesh->create(vertex_num, vertex_array, tria_num, tria_array, quad_num, quad_array, marker_num, marker_array);
  }
  else {
    const int vertex_num = 4, tria_num = 0, quad_num = 1, marker_num = 4;
    double2 vertex_array[4] = {{-1,-1}, { 1,-1}, { 1, 1}, {-1, 1}};
    int4 *tria_array = NULL;
    int5 quad_array[1] = {{0, 1, 2, 3, 0}};
    int3 marker_array[4] = {{0, 1, 3}, {1, 2, 2}, {2, 3, 4}, {3, 0, 1}};
    mesh->create(vertex_num, vertex_array, tria_num, tria_array, quad_num, quad_array, marker_num, marker_array);
  }
  return mesh;
}

/// Initialize internal data structures: H1 space.
bool init_h1(bool tri) 
{
  try 
  {
    // Shapeset, cache and mesh.
    H1Shapeset* h1_shapeset = new H1Shapeset();
    mesh = init_mesh(tri);
    shapeset = h1_shapeset;

    // Finite element space.
    space = new H1Space<double>(mesh);

    // Weakform.
    weakform = new H1TestWeakForm();

    // Prepare selector.
    selector = new H1ProjBasedSelector<double>(H2D_HP_ANISO, 1.0, H2DRS_DEFAULT_ORDER, h1_shapeset);

    return true;
  }
  catch (std::exception& e) 
  {
    info("failed: %s", e.what());
    return false;
  }
}

/// Initialize internal data structures: L2 space.
bool init_l2(bool tri) 
{
  try 
  {
    // Shapeset, cache and mesh.
    L2Shapeset* l2_shapeset = new L2Shapeset();
    mesh = init_mesh(tri);
    shapeset = l2_shapeset;

    // Finite element space.
    space = new L2Space<double>(mesh, 1);

    // Weakform.
    weakform = new L2TestWeakForm();

    // Prepare selector.
    selector = new L2ProjBasedSelector<double>(H2D_HP_ANISO, 1.0, H2DRS_DEFAULT_ORDER, l2_shapeset);

    return true;
  }
  catch (std::exception& e) 
  {
    info("failed: %s", e.what());
    return false;
  }
}

/// Prints failure matrix.
void show_fail_matrix(int** fail_matrix, const std::string& space_name, const int max_quad_order) 
{
  const int max_order_h = H2D_GET_H_ORDER(max_quad_order), max_order_v = H2D_GET_V_ORDER(max_quad_order);
#define NUMBER_W 2 ///< A size of a cell in the table.
#define NUMBER_FMT "%02d" ///< A format for writing orders.
  char buff_number[1024];

  info("!Test summary (V/H): %s", space_name.c_str());

  // Print header.
  {
    std::stringstream str;
    for(int i = 0; i < NUMBER_W; i++)
      str << ' ';
    for(int i = 0; i <= max_order_h; i++) {
      str << '|';
      sprintf(buff_number, NUMBER_FMT, i);
      str << buff_number;
    }
    info(" %s", str.str().c_str());
  }

  // Print body.
  for(int i = 0; i <= max_order_v; i++) 
  {
    // Build row head.
    std::stringstream str;
    sprintf(buff_number, NUMBER_FMT, i);
    str << buff_number;

    // Build row body.
    for(int k = 0; k <= max_order_h; k++) 
    {
      str << '|';
      if (fail_matrix[i][k] == H2D_TEST_NOT_DONE) 
      {
        for(int j = 0; j < NUMBER_W; j++)
          str << '-';
      }
      else {
        if ((fail_matrix[i][k] & H2D_RSLN_FAILED) != 0)
          str << 'S';
        else
          str << ' ';
        if ((fail_matrix[i][k] & H2D_CAND_FAILED) != 0)
          str << 'C';
        else
          str << ' ';
        for(int j = 2; j < NUMBER_W; j++)
          str << ' ';
      }
    }

    // Print row.
    info(" %s", str.str().c_str());
  }
}

/// Test reference solution whether its values matche values of the polynom in the test case.
bool test_ref_solution(Solution<double>& rsln) 
{
  verbose("Testing reference solution");

  // Check FN order.
  int test_quad_order = cur_test_case->quad_order();

  // Check function values.
  bool failed = false;
  Mesh* mesh = rsln.get_mesh();
  Element* element;
  for_all_active_elements(element, mesh) {
    // Set element.
    rsln.set_active_element(element);
    int rsln_order = rsln.get_fn_order();
    if (rsln_order < H2D_GET_H_ORDER(test_quad_order) || rsln_order < H2D_GET_V_ORDER(test_quad_order)) 
    {
      verbose(" Failed: reference solution order (%d/%d) is lower than test case order (%d/%d) at element #%d", rsln_order, rsln_order, H2D_GET_H_ORDER(test_quad_order), H2D_GET_V_ORDER(test_quad_order), element->id);
      return true;
    }

    // Set GIP order at which points will be inspected.
    int check_gip_order = 2*rsln_order;
    rsln.set_quad_order(check_gip_order);

    // Get rsln values.
    double* rsln_vals = rsln.get_fn_values(0);

    // Get physical locations.
    RefMap* ref_map = rsln.get_refmap();
    double* phys_x = ref_map->get_phys_x(check_gip_order);
    double* phys_y = ref_map->get_phys_y(check_gip_order);

    // Compare to test-case values.
    const int num_pt = rsln.get_quad_2d()->get_num_points(check_gip_order);
    for(int i = 0; i < num_pt; i++) {
      double func_value = func_val(cur_test_case->poly_matrix(), test_quad_order, phys_x[i], phys_y[i]);
      if (Hermes::abs(func_value - rsln_vals[i]) > numerical_zero) {
        verbose(" Failed: rsln and function value differs (%g)", Hermes::abs(func_value - rsln_vals[i]));
        return true;
      }
    }
  }

  return false;
}

/// Test.
bool test(bool tri, const std::string& space_name, int min_order, int max_order = H2DRS_MAX_ORDER) 
{
  bool failed = false;

  // Prepare cases.
  int min_quad_order = -1, max_quad_order = -1;
  if (tri) 
  {
    min_quad_order = H2D_MAKE_QUAD_ORDER(min_order, 0);
    max_quad_order = H2D_MAKE_QUAD_ORDER(max_order, 0);
  }
  else 
  {
    min_quad_order = H2D_MAKE_QUAD_ORDER(min_order, min_order);
    max_quad_order = H2D_MAKE_QUAD_ORDER(max_order, max_order);
  }
  OrderPermutator order_perm(min_quad_order, max_quad_order, false);

  // Prepare place for result summary.
  int end_order_v = H2D_GET_V_ORDER(order_perm.get_end_quad_order());
  int end_order_h = H2D_GET_H_ORDER(order_perm.get_end_quad_order());
  int** fail_matrix = new_matrix<int>(end_order_v+1, end_order_h+1);
  for(int i = 0; i <= end_order_v; i++)
    for(int k = 0; k <= end_order_h; k++)
      fail_matrix[i][k] = H2D_TEST_NOT_DONE;

  // Initialize the FE problem.
  /*  bool is_linear = true;
  DiscreteProblem dp(&weakform, &ref_space, is_linear);

  // Initialize matrix solver.
  SparseMatrix* mat = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, mat, rhs);
  */
  // Process cases.
  do 
  {
    int test_result = H2D_TEST_SUCCESS;

    // Prepare test case.
    TestCase test_case(order_perm.get_quad_order());
    cur_test_case = &test_case;
    verbose("!test case: %s", test_case.title().c_str());

    // Process.
    space->set_element_order(H2D_TEST_ELEM_ID, test_case.start_quad_order());

    // Assemble the stiffness matrix and right-hand side vector.
    //  info("Assembling the stiffness matrix and right-hand side vector.");
    //  dp.assemble(mat, rhs);

    // Create and solve the reference system.
    Solution<double> rsln;

    // Construct globally refined reference mesh
    // and setup reference space.
    Mesh *ref_mesh = new Mesh();
    ref_mesh->copy(space->get_mesh());
    ref_mesh->refine_all_elements();
    int order_increase = 1;
    Space<double>* ref_space = space->dup(ref_mesh, order_increase);

    // Solve the reference problem.
    // solve_linear(ref_space, weakform, matrix_solver, &rsln);

    // Initialize the FE problem.
    bool is_linear = true;
    DiscreteProblem<double> dp(weakform, ref_space);

    // Initialize matrix solver.
    SparseMatrix<double>* mat = create_matrix<double>(matrix_solver);
    Vector<double>* rhs = create_vector<double>(matrix_solver);
    LinearSolver<double>* solver = create_linear_solver<double>(matrix_solver, mat, rhs);

    // Assemble the stiffness matrix and right-hand side vector.
    info("Assembling the stiffness matrix and right-hand side vector.");
    dp.assemble(mat, rhs);

    info("Solving the matrix problem.");
    if(solver->solve())
      Solution<double>::vector_to_solution(solver->get_solution(), ref_space, &rsln);
    else
      error ("Matrix solver failed.\n");



    // Check projected functions.
    if (test_ref_solution(rsln))
      test_result |= H2D_RSLN_FAILED;

    // Select candidate.
    ElementToRefine refinement;
    Element* e = mesh->get_element(H2D_TEST_ELEM_ID);
    int order = space->get_element_order(H2D_TEST_ELEM_ID);
    selector->select_refinement(e, order, &rsln, refinement);

    // Check candidates.
    verbose("Testing candidates");
    const std::vector<OptimumSelector<double>::Cand>& candidates = selector->get_candidates();
    std::vector<OptimumSelector<double>::Cand>::const_iterator cand = candidates.begin();
    while (cand != candidates.end()) 
    {
      if (cur_test_case->should_match(*cand) && Hermes::abs(cand->error) > numerical_zero) 
        test_result |= H2D_CAND_FAILED;
      cand++;
    }
    if (test_result == H2D_TEST_SUCCESS)
      verbose("Test success");
    else 
    {
      verbose("Test failed!");
      failed = true;
    }
    fail_matrix[order_perm.get_order_v()][order_perm.get_order_h()] = test_result;
  } while(order_perm.next());

  // Print result matrix.
  show_fail_matrix(fail_matrix, space_name, order_perm.get_end_quad_order());

  // Clenup.
  delete[] fail_matrix;

  return !failed;
}

/// Test entry-point.
int main(int argc, char* argv[]) 
{
  // Check input parameters.
  bool tri = false;
  if (argc > 1 && strcmp(argv[1], "-tri") == 0)
    tri = true;
  bool test_success = true;

  // Set zero error.
  if (tri)
    numerical_zero = H2D_TEST_ZERO_TRIS;
  else
    numerical_zero = H2D_TEST_ZERO_QUADS;

  // H1.
  info("! Space: H1");
  if (init_h1(tri))
    test_success &= test(tri, "H1", 2, H2DRS_MAX_ORDER);
  cleanup();
  if (!test_success)
    goto quit;

  // L2.
  info("! Space: L2");
  if (init_l2(tri))
    // test_success &= test(tri, "L2", 2, 3);
    test_success &= test(tri, "L2", 2, H2DRS_MAX_ORDER);
  cleanup();
  if (!test_success)
    goto quit;

quit:
  if (test_success)
  {
    info("!Test: Success");
    return TEST_SUCCESS;
  }
  else {
    info("!Test: Failed!");
    return TEST_FAILURE;
  }
}