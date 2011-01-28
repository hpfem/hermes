#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#include "config.h"
//#include <getopt.h>
#include <hermes3d.h>

// The following parameters can be changed:
const char* iterative_method = "bicgstab";        // Name of the iterative method employed by AztecOO (ignored
                                                  // by the other solvers). 
                                                  // Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* preconditioner = "jacobi";            // Name of the preconditioner employed by AztecOO (ignored by
                                                  // the other solvers). 
                                                  // Possibilities: none, jacobi, neumann, least-squares, or a
                                                  // preconditioner from IFPACK (see solver/aztecoo.h).
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters.
//#define XM_YN_ZO
#define XM_YN_ZO_2
//#define X2_Y2_Z2
//#define X3_Y3_Z3

int m = 2, n = 2, o = 2;

template<typename S, typename T>
T fnc(S x, S y, S z) {
#ifdef XM_YN_ZO
  return pow(x, m) * pow(y, n) * pow(z, o) + pow(x, 2) * pow(y, 3) - pow(x, 3) * z + pow(z, 4);
#elif defined XM_YN_ZO_2
  return pow(x, m) * pow(y, n) * pow(z, o) + pow(x, 2) * pow(y, 3) - pow(x, 2) * z + pow(z, 4);
#elif defined X2_Y2_Z2
  return x*x + y*y + z*z;
#elif defined X3_Y3_Z3
  return x*x*x + y*y*y + z*z*z;
#endif
}

template<typename S, typename T>
T dfnc(S x, S y, S z) {
#ifdef XM_YN_ZO
  T ddxx = m * (m - 1) * pow(x, m - 2) * pow(y, n) * pow(z, o) + 2 * pow(y, 3) - 6 * x * z;
  T ddyy = n * (n - 1) * pow(x, m) * pow(y, n - 2) * pow(z, o) + 6 * pow(x, 2) * y;
  T ddzz = o * (o - 1) * pow(x, m) * pow(y, n) * pow(z, o - 2) + 12 * pow(z, 2);
  return -(ddxx + ddyy + ddzz);

#elif defined XM_YN_ZO_2
  T ddxx = m*(m-1) * pow(x, m-2) * pow(y, n) * pow(z, o) + 2 * pow(y, 3) - 2 * z;
  T ddyy = n*(n-1) * pow(x, m) * pow(y, n-2) * pow(z, o) + 6 * pow(x, 2) * y;
  T ddzz = o*(o-1) * pow(x, m) * pow(y, n) * pow(z, o-2) + 12 * pow(z, 2);
  return -(ddxx + ddyy + ddzz);
#elif defined X2_Y2_Z2
  return -6.0;
#elif defined X3_Y3_Z3
  return -6.0 * x - 6.0 * y - 6.0 * z;
#endif
}

// Exact solution.
double exact_solution(double x, double y, double z, double &dx, double &dy, double &dz) {
#ifdef XM_YN_ZO
  dx = m * pow(x, m - 1) * pow(y, n) * pow(z, o) + 2 * x * pow(y, 3) - 3 * pow(x, 2) * z;
  dy = n * pow(x, m) * pow(y, n - 1) * pow(z, o) + 3 * pow(x, 2) * pow(y, 2);
  dz = o * pow(x, m) * pow(y, n) * pow(z, o - 1) - pow(x, 3) + 4 * pow(z, 3);
#elif defined XM_YN_ZO_2
  dx = m * pow(x, m-1) * pow(y, n) * pow(z, o) + 2 * x * pow(y, 3) - 2 * x * z;
  dy = n * pow(x, m) * pow(y, n-1) * pow(z, o) + 3 * pow(x, 2) * pow(y, 2);
  dz = o * pow(x, m) * pow(y, n) * pow(z, o-1) - pow(x, 2) + 4 * pow(z, 3);
#elif defined X2_Y2_Z2
  dx = 2 * x;
  dy = 2 * y;
  dz = 2 * z;
#elif defined X3_Y3_Z3
  dx = 3 * x*x;
  dy = 3 * y*y;
  dz = 3 * z*z;
#endif

  return fnc<double, double>(x, y, z);
}

// Boundary condition types.
BCType bc_types(int marker) 
{
#ifdef DIRICHLET
  return BC_ESSENTIAL;
#elif defined NEWTON
  return BC_NATURAL;
#endif
}

// Dirichlet boundary conditions.
scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z) {
#ifdef DIRICHLET
  return fnc<double, scalar>(x, y, z);
#else
  return 0;
#endif
}

template<typename Real, typename Scalar>
Scalar bilinear_form(int np, double *jwt, Func<Scalar> *u_ext[], Func<Real> *fu, Func<Real> *fv, Geom<Real> *e, ExtData<Scalar> *ud) {
  return int_grad_u_grad_v<Real, Scalar>(np, jwt, fu, fv, e);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_surf(int np, double *jwt, Func<Scalar> *u_ext[], Func<Real> *fu, Func<Real> *fv, Geom<Real> *e, ExtData<Scalar> *ud) {
  return int_u_v<Real, Scalar>(np, jwt, fu, fv, e);
}

template<typename Real, typename Scalar>
Scalar linear_form(int np, double *jwt, Func<Scalar> *u_ext[], Func<Real> *fv, Geom<Real> *e, ExtData<Scalar> *ud) {
  return int_F_v<Real, Scalar>(np, jwt, dfnc, fv, e);
}

template<typename Real, typename Scalar>
Scalar linear_form_surf(int np, double *jwt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ud) {
  Scalar result = 0;
#ifdef XM_YN_ZO
  for (int i = 0; i < np; i++) {
    Scalar dx = m * pow(e->x[i], m - 1) * pow(e->y[i], n) * pow(e->z[i], o) + 2 * e->x[i] * pow(e->y[i], 3) - 3 * pow(e->x[i], 2) * e->z[i];
    Scalar dy = n * pow(e->x[i], m) * pow(e->y[i], n - 1) * pow(e->z[i], o) + 3 * pow(e->x[i], 2) * pow(e->y[i], 2);
    Scalar dz = o * pow(e->x[i], m) * pow(e->y[i], n) * pow(e->z[i], o - 1) - pow(e->x[i], 3) + 4 * pow(e->z[i], 3);
    result += jwt[i] * (v->fn[i] * (dx * e->nx[i] + dy * e->ny[i] + dz * e->nz[i] + fnc<Real, Scalar>(e->x[i], e->y[i], e->z[i])));
  }
#elif defined XM_YN_ZO_2
  for (int i = 0; i < np; i++) {
    Scalar dx = m * pow(e->x[i], m-1) * pow(e->y[i], n) * pow(e->z[i], o) + 2 * e->x[i] * pow(e->y[i], 3) - 2 * e->x[i] * e->z[i];
    Scalar dy = n * pow(e->x[i], m) * pow(e->y[i], n-1) * pow(e->z[i], o) + 3 * pow(e->x[i], 2) * pow(e->y[i], 2);
    Scalar dz = o * pow(e->x[i], m) * pow(e->y[i], n) * pow(e->z[i], o-1) - pow(e->x[i], 2) + 4 * pow(e->z[i], 3);
    result += jwt[i] * (v->val[i] * (dx * e->nx[i] + dy * e->ny[i] + dz * e->nz[i] + fnc<Real, Scalar>(e->x[i], e->y[i], e->z[i])));
  }
#elif defined X2_Y2_Z2
  for (int i = 0; i < np; i++) {
    Scalar dx = 2 * e->x[i];
    Scalar dy = 2 * e->y[i];
    Scalar dz = 2 * e->z[i];
    result += jwt[i] * (v->fn[i] * (dx * e->nx[i] + dy * e->ny[i] + dz * e->nz[i] + fnc<Real, Scalar>(e->x[i], e->y[i], e->z[i])));
  }
#elif defined X3_Y3_Z3
  for (int i = 0; i < np; i++) {
    Scalar dx = 3 * e->x[i] * e->x[i];
    Scalar dy = 3 * e->y[i] * e->y[i];
    Scalar dz = 3 * e->z[i] * e->z[i];
    result += jwt[i] * (v->fn[i] * (dx * e->nx[i] + dy * e->ny[i] + dz * e->nz[i] + fnc<Real, Scalar>(e->x[i], e->y[i], e->z[i])));
  }
#endif
  return result;
}

// Helpers.

int parse_reft(char *str) {
  if (strcasecmp(str, "x") == 0) return H3D_REFT_HEX_X;
  else if (strcasecmp(str, "y") == 0) return H3D_REFT_HEX_Y;
  else if (strcasecmp(str, "z") == 0) return H3D_REFT_HEX_Z;
  else if (strcasecmp(str, "xy") == 0 || strcasecmp(str, "yx") == 0) return H3D_H3D_REFT_HEX_XY;
  else if (strcasecmp(str, "xz") == 0 || strcasecmp(str, "zx") == 0) return H3D_H3D_REFT_HEX_XZ;
  else if (strcasecmp(str, "yz") == 0 || strcasecmp(str, "zy") == 0) return H3D_H3D_REFT_HEX_YZ;
  else if (strcasecmp(str, "xyz") == 0) return H3D_H3D_H3D_REFT_HEX_XYZ;
  else return H3D_REFT_HEX_NONE;
}


// Maximum level of refinement considered (= 0 .. # of ref)
#ifdef DEV_TESTS
  #define MAX_LEVEL          4
#else
  #define MAX_LEVEL          2
#endif

#define NUM_RULES            ((MAX_LEVEL + 1) * (MAX_LEVEL + 1))

// Special quadrature used in this test to define set of points, where
// continuity is tested. Quadrature is used in order to use RefMap abilities
// calculate physical coordinates of points from the reference domain.
//
// points are chosen in this way:
// - only face points are used
// - there are more levels of points
// - 0-level consist of 16 points chosen symmetricaly
// - other levels (or orders) simulate division of adjacent elements and map
//   points in such way, that some points from divided and some from nondivided
//   face will still match in the physical domain
// - this is of course only possible up to some level of division, controled by
//   constant LEVEL. for example, for hex1 :
//   LEVEL = 1  divisions 0 x 1 y or 0 x 1 y 3 z are ok, but 0 x 1 y 3 y not
//   LEVEL = 2  division 0 x 1 y 3 y is ok, but 0 x 1 y 3 y 5 y not
// - if LEVEL is not sufficient, there will be some faces, that will not be tested,
//   because no points from the face will match to points from the constraining face.
class ContQuad : public Quad3D {
public:
  ContQuad() {
    _F_  
    int my_np_1d[MAX_LEVEL + 1];
    double my_tables_1d[MAX_LEVEL + 1][1000];

    my_np_1d[0] = 4;
    my_tables_1d[0][0] = -0.71;  my_tables_1d[0][1] = -0.59;
    my_tables_1d[0][2] =  0.59;  my_tables_1d[0][3] =  0.71;

    for(int order = 1; order <= MAX_LEVEL; order++) {
      my_np_1d[order] = 2 * my_np_1d[order - 1];
      for(int i = 0; i < my_np_1d[order - 1]; i++) {
        my_tables_1d[order][i] = (my_tables_1d[order - 1][i] - 1.) / 2.;
        my_tables_1d[order][i + my_np_1d[order - 1]] = (my_tables_1d[order - 1][i] + 1.) / 2.;
      }
    }
    for (unsigned int order = 0; order < NUM_RULES; order++) {
      int ord1 = order / (MAX_LEVEL + 1);
      int ord2 = order % (MAX_LEVEL + 1);
      (*np_face)[order] = my_np_1d[ord1] * my_np_1d[ord2];
      for (int face = 0; face < Hex::NUM_FACES; face++) {
        if((*face_tables)[face] == NULL)
          (*face_tables)[face] = new std::map<unsigned int, QuadPt3D *>;
        (*(*face_tables)[face])[order] = new QuadPt3D[(*np_face)[order]];
      }

      for (int k = 0, n = 0; k < my_np_1d[ord1]; k++) {
        for (int l = 0; l < my_np_1d[ord2]; l++, n++) {
          assert(n < (*np_face)[order]);
          (*(*face_tables)[0])[order][n].x = -1;
          (*(*face_tables)[0])[order][n].y = my_tables_1d[ord1][k];
          (*(*face_tables)[0])[order][n].z = my_tables_1d[ord2][l];
        }
      }

      for (int k = 0, n = 0; k < my_np_1d[ord1]; k++) {
        for (int l = 0; l < my_np_1d[ord2]; l++, n++) {
          assert(n < (*np_face)[order]);
          (*(*face_tables)[1])[order][n].x = 1;
          (*(*face_tables)[1])[order][n].y = my_tables_1d[ord1][k];
          (*(*face_tables)[1])[order][n].z = my_tables_1d[ord2][l];
        }
      }

      for (int k = 0, n = 0; k < my_np_1d[ord1]; k++) {
        for (int l = 0; l < my_np_1d[ord2]; l++, n++) {
          assert(n < (*np_face)[order]);
          (*(*face_tables)[2])[order][n].x = my_tables_1d[ord1][k];
          (*(*face_tables)[2])[order][n].y = -1;
          (*(*face_tables)[2])[order][n].z = my_tables_1d[ord2][l];
        }
      }

      for (int k = 0, n = 0; k < my_np_1d[ord1]; k++) {
        for (int l = 0; l < my_np_1d[ord2]; l++, n++) {
          assert(n < (*np_face)[order]);
          (*(*face_tables)[3])[order][n].x = my_tables_1d[ord1][k];
          (*(*face_tables)[3])[order][n].y = 1;
          (*(*face_tables)[3])[order][n].z = my_tables_1d[ord2][l];
        }
      }

      for (int k = 0, n = 0; k < my_np_1d[ord1]; k++) {
        for (int l = 0; l < my_np_1d[ord2]; l++, n++) {
          assert(n < (*np_face)[order]);
          (*(*face_tables)[4])[order][n].x = my_tables_1d[ord1][k];
          (*(*face_tables)[4])[order][n].y = my_tables_1d[ord2][l];
          (*(*face_tables)[4])[order][n].z = -1;
        }
      }

      for (int k = 0, n = 0; k < my_np_1d[ord1]; k++) {
        for (int l = 0; l < my_np_1d[ord2]; l++, n++) {
          assert(n < (*np_face)[order]);
          (*(*face_tables)[5])[order][n].x = my_tables_1d[ord1][k];
          (*(*face_tables)[5])[order][n].y = my_tables_1d[ord2][l];
          (*(*face_tables)[5])[order][n].z = 1;
        }
      }
    }
  }

  ~ContQuad() {
    _F_
    for(std::map<unsigned int, std::map<unsigned int, QuadPt3D *>*>::iterator it = face_tables->begin(); it != face_tables->end(); it++) {
      for(std::map<unsigned int, QuadPt3D *>::iterator it_inner = it->second->begin(); it_inner != it->second->end(); it_inner++)
        delete [] it_inner->second;
      delete it->second;
    }
  }
};

struct Point {
  double ref_x, ref_y, ref_z;
  double phys_x, phys_y, phys_z;
  unsigned int elm_idx;
  Facet::Key fac_idx;

  Point(int idx, Facet::Key f_idx, double rx, double ry, double rz, double px, double py, double pz) {
    elm_idx = idx; fac_idx = f_idx;
    ref_x = rx;  ref_y = ry;  ref_z = rz;
    phys_x = px; phys_y = py; phys_z = pz;
  }
};

typedef
  int (*compfn)(const void*, const void*);

int compare(Point **pt1, Point **pt2) {
  _F_
  double val1 = 1000000. * (*pt1)->phys_x + 1000. * (*pt1)->phys_y + (*pt1)->phys_z;
  double val2 = 1000000. * (*pt2)->phys_x + 1000. * (*pt2)->phys_y + (*pt2)->phys_z;

  if (val1 < val2) return -1;
  else if (val1 > val2) return 1;
  else return 0;
}

// The error should be smaller than this epsilon.  
const double EPS = 1e-10;
// For the testing of continuity, the jump can not be higher than this.
const double TOLERANCE = 1e-10;

bool equal(Point *pt1, Point *pt2) {
  _F_
  if (fabs(pt1->phys_x - pt2->phys_x) > TOLERANCE) return false;
  if (fabs(pt1->phys_y - pt2->phys_y) > TOLERANCE) return false;
  if (fabs(pt1->phys_z - pt2->phys_z) > TOLERANCE) return false;
  return true;
}

int main(int argc, char **args) 
{
  // Test variable.
  int success_test = 1;

  if (argc < 2) error("Not enough parameters.");

  // Load the mesh.
  Mesh mesh;
  H3DReader mloader;
  if (!mloader.load(args[1], &mesh)) error("Loading mesh file '%s'.", args[1]);

  // Apply refinements according to the command line parameters passed.
  for (int i = 2; i < argc; i += 2) {
    int elem_id, reft_id;
    sscanf(args[i], "%d", &elem_id);
    reft_id = parse_reft(args[i + 1]);
    mesh.refine_element(elem_id + 1, reft_id);
  }

  Ord3 order(2, 3, 4);

  // Initialize the space.
  H1Space space(&mesh, bc_types, essential_bc_values, order);

  
  // Initialize the weak formulation.
  WeakForm wf(1);
#ifdef DIRICHLET
  wf.add_matrix_form(0, 0, bilinear_form<double, scalar>, bilinear_form<Ord, Ord>, HERMES_SYM);
  wf.add_vector_form(0, linear_form<double, scalar>, linear_form<Ord, Ord>);
#elif defined NEWTON
  wf.add_matrix_form(0, 0, bilinear_form<double, scalar>, bilinear_form<Ord, Ord>, HERMES_SYM);
  wf.add_matrix_form_surf(0, 0, bilinear_form_surf<double, scalar>, bilinear_form_surf<Ord, Ord>);
  wf.add_vector_form(0, linear_form<double, scalar>, linear_form<Ord, Ord>);
  wf.add_vector_form_surf(0, linear_form_surf<double, scalar>, linear_form_surf<Ord, Ord>);
#endif

  // Initialize the FE problem.
  bool is_linear = true;
  DiscreteProblem dp(&wf, &space, is_linear);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);
  
  // Initialize the preconditioner in the case of SOLVER_AZTECOO.
  if (matrix_solver == SOLVER_AZTECOO) 
  {
    ((AztecOOSolver*) solver)->set_solver(iterative_method);
    ((AztecOOSolver*) solver)->set_precond(preconditioner);
    // Using default iteration parameters (see solver/aztecoo.h).
  }

  // Assemble the linear problem.
  info("Assembling (ndof: %d).", Space::get_num_dofs(&space));
  dp.assemble(matrix, rhs);
    
  // Solve the linear system. If successful, obtain the solution.
  info("Solving.");
  Solution sln(&mesh);
  if(solver->solve()) Solution::vector_to_solution(solver->get_solution(), &space, &sln);
  else error ("Matrix solver failed.\n");

  ExactSolution ex_sln(&mesh, exact_solution);

  // Calculate exact error.
  info("Calculating exact error.");
  Adapt *adaptivity = new Adapt(&space, HERMES_H1_NORM);
  bool solutions_for_adapt = false;
  double err_exact = adaptivity->calc_err_exact(&sln, &ex_sln, solutions_for_adapt, HERMES_TOTAL_ERROR_ABS);

  if (err_exact > EPS)
    // Calculated solution is not precise enough.
    success_test = 0;

  // Special code for this test starts here.
  ContQuad my_quad;
  RefMap ref_map(&mesh);

  int num_points = 0;
  for (int order = 0; order < NUM_RULES; order++)
    for(std::map<unsigned int, Element*>::iterator it = mesh.elements.begin(); it != mesh.elements.end(); it++)
      if (it->second->used && it->second->active)
        for (int iface = 0; iface < Hex::NUM_FACES; iface++)
          num_points += my_quad.get_face_num_points(iface, order);

  Point **points = new Point *[num_points];
  int ipt = 0;

  // Find points.
  for (int order = 0; order < NUM_RULES; order++) {
    for(std::map<unsigned int, Element*>::iterator it = mesh.elements.begin(); it != mesh.elements.end(); it++)
      if (it->second->used && it->second->active) {
        Element *e = mesh.elements[it->first];
        ref_map.set_active_element(e);
        for (int iface = 0; iface < Hex::NUM_FACES; iface++) {
          Facet::Key fac_idx = mesh.get_facet_id(e, iface);

          QuadPt3D *quad_pts = my_quad.get_face_points(iface, order);
          int np = my_quad.get_face_num_points(iface, order);
          double *phys_x = ref_map.get_phys_x(np, quad_pts);
          double *phys_y = ref_map.get_phys_y(np, quad_pts);
          double *phys_z = ref_map.get_phys_z(np, quad_pts);

          // For each face and each integration point store reference and physical coordinates.
          for (int pt = 0; pt < np; pt++) {
            points[ipt++] = new Point(it->first, fac_idx,
              quad_pts[pt].x, quad_pts[pt].y, quad_pts[pt].z,
              phys_x[pt], phys_y[pt], phys_z[pt]);
          }
        }
      }
  }

  // Sort points according to first phys_x, then phys_y and phys_z
  // it means, that two points, with almost identical physical coordinates
  // (even though from different elements) will be next to each other in the array.
  qsort((void *) points, num_points, sizeof(Point*), (compfn)compare);

  int *pairs = new int [num_points];
  int num_pairs = 0;

  // Choose those indicies, that correspond to pairs with identical physical coordinates
  // and store them in field pairs.
  for (int i = 0; i < num_points - 1; i++) {
    if (equal(points[i], points[i+1])) {
      pairs[num_pairs++] = i;
    }
  }

  // Check, whether we tested points from all inner active facets
  // this is done only for testing of correctness of the test itself.
  int nonchecked_faces = 0;
  for(std::map<Facet::Key, Facet*>::iterator it = mesh.facets.begin(); it != mesh.facets.end(); it++) {
    bool ok = false;
    Facet *fac = it->second;
    if (fac->type == Facet::OUTER) continue;
    if (!(fac->ractive || fac->lactive)) continue;
    for (int i = 0; i < num_pairs; i++) {
      if ((points[pairs[i]]->fac_idx == it->first) || (points[pairs[i] + 1]->fac_idx == it->first)) {
        ok = true;
        break;
      }
    }

    if (!ok) nonchecked_faces++;
  }


  // Loop over all basis functions.
  int ndofs = Space::get_num_dofs(&space);
  for (int dof = 0; dof < ndofs; dof++) {
    info("processing dof %d...", dof);

    // Prepare solution corresponding to basis function with dof dof.
    double *sln_vector = new double[ndofs];
    memset(sln_vector, 0, ndofs * sizeof(double));
    sln_vector[dof] = 1.0;
    Solution::vector_to_solution(sln_vector, &space, &sln);

    double max_difference = 0.;
    double max_pt_x, max_pt_y, max_pt_z, max_val_1, max_val_2;
    unsigned int max_elm_1, max_elm_2;

    // Loop over all pairs of points, that correspond to one point in the physical domain.
    for (int pair = 0; pair < num_pairs; pair++) {
      int i = pairs[pair];
      Element *e1 = mesh.elements[points[i]->elm_idx];
      sln.set_active_element(e1);
      double val1 = sln.get_pt_value(points[i]->ref_x, points[i]->ref_y, points[i]->ref_z);

      Element *e2 = mesh.elements[points[i + 1]->elm_idx];
      sln.set_active_element(e2);
      double val2 = sln.get_pt_value(points[i + 1]->ref_x, points[i + 1]->ref_y, points[i + 1]->ref_z);

      // Value in this point should be the same, no matter from which element we go.
      double difference = fabs(val1 - val2);
      if (difference > max_difference) {
        max_difference = difference;
        max_pt_x = points[i]->phys_x;
        max_pt_y = points[i]->phys_y;
        max_pt_z = points[i]->phys_z;
        max_val_1 = val1;
        max_val_2 = val2;
        max_elm_1 = points[i]->elm_idx;
        max_elm_2 = points[i + 1]->elm_idx;
      }
    }

    if (max_difference > TOLERANCE) {
      info("failed\nbase fn %d NOT continuous between elements %ld and %ld @ (% lf, % lf, % lf), max difference %g (%.15g <-> %.15g)",
           dof, max_elm_1, max_elm_2 , max_pt_x, max_pt_y, max_pt_z, max_difference, max_val_1, max_val_2);
      success_test = 0;
    }
    else
      info("ok");
  }

  for (int i = 0; i < num_points; i++) delete points[i];
  delete [] points;
  delete [] pairs;

  info("continuity tested in %d points and %d inner faces with at least one active adjacent element were not tested.", num_pairs, nonchecked_faces);

  if (success_test) {
    info("Success!");
    return ERR_SUCCESS;
  }
  else {
    info("Failure!");
    return ERR_FAILURE;
  }
}

