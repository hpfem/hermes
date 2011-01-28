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
//#define X2_Y2_Z2
#define FN4

const double alpha = 1.0;

// Exact solution.
scalar3 &exact_solution(double x, double y, double z, scalar3 &dx, scalar3 &dy, scalar3 &dz) {
  _F_
  static scalar3 val;

#ifdef FN4
  dx[0] = (1 - y*y)*(1 - z*z)*(z - 6*x*x);
  dx[1] = 2*(1 - x*x)*(1 - z*z) - 2*x*(1 - z*z)*(y*y*y + 2*x);
  dx[2] = -2*x*(1 - y*y)*(z*z - 3*x*y*z) - 3*y*z*(1 - x*x)*(1 - y*y);

  dy[0] = -2*y*(1 - z*z)*(1 - 2*x*x*x + x*z);
  dy[1] = 3*y*y*(1 - x*x)*(1 - z*z);
  dy[2] = -2*y*(1 - x*x)*(z*z - 3*x*y*z) - 3*x*z*(1 - x*x)*(1 - y*y);

  dz[0] = x*(1 - y*y)*(1 - z*z) - 2*z*(1 - y*y)*(1 - 2*x*x*x + x*z);
  dz[1] = -2*z*(1 - x*x)*(y*y*y + 2*x);
  dz[2] = (1 - x*x)*(1 - y*y)*(2*z - 3*x*y);

  val[0] = (1-y*y) * (1-z*z) * (x*z - 2*x*x*x + 1);
  val[1] = (1-x*x) * (1-z*z) * (y*y*y + 2*x);
  val[2] = (1-x*x) * (1-y*y) * (z*z - 3*x*y*z);
#elif defined X2_Y2_Z2
  dx[0] = 0;
  dx[1] = 0;
  dx[2] = -2 * x * (1 - y*y) * (1 - z*z);

  dy[0] = 0;
  dy[1] = 0;
  dy[2] = -2 * (1 - x*x) * y * (1 - z*z);

  dz[0] = 0;
  dz[1] = 0;
  dz[2] = -2 * (1 - x*x) * (1 - y*y) * z;

  val[0] = 0;
  val[1] = 0;
  val[2] = (1 - x*x) * (1 - y*y) * (1 - z*z);
#endif

  return val;
}

template<typename S, typename T>
void exact_sln(S x, S y, S z, T (&val)[3], T (&dx)[3], T (&dy)[3], T (&dz)[3]) {
  _F_
#ifdef FN4
  val[0] = (1-y*y) * (1-z*z) * (x*z - 2*x*x*x + 1);
  val[1] = (1-x*x) * (1-z*z) * (y*y*y + 2*x);
  val[2] = (1-x*x) * (1-y*y) * (z*z - 3*x*y*z);

  dx[0] = (1 - y*y)*(1 - z*z)*(z - 6*x*x);
  dx[1] = 2*(1 - x*x)*(1 - z*z) - 2*x*(1 - z*z)*(y*y*y + 2*x);
  dx[2] = -2*x*(1 - y*y)*(z*z - 3*x*y*z) - 3*y*z*(1 - x*x)*(1 - y*y);

  dy[0] = -2*y*(1 - z*z)*(1 - 2*x*x*x + x*z);
  dy[1] = 3*y*y*(1 - x*x)*(1 - z*z);
  dy[2] = -2*y*(1 - x*x)*(z*z - 3*x*y*z) - 3*x*z*(1 - x*x)*(1 - y*y);

  dz[0] = x*(1 - y*y)*(1 - z*z) - 2*z*(1 - y*y)*(1 - 2*x*x*x + x*z);
  dz[1] = -2*z*(1 - x*x)*(y*y*y + 2*x);
  dz[2] = (1 - x*x)*(1 - y*y)*(2*z - 3*x*y);
#elif defined X2_Y2_Z2
  dx[0] = 0;
  dx[1] = 0;
  dx[2] = -2 * x * (1 - y*y) * (1 - z*z);

  dy[0] = 0;
  dy[1] = 0;
  dy[2] = -2 * (1 - x*x) * y * (1 - z*z);

  dz[0] = 0;
  dz[1] = 0;
  dz[2] = -2 * (1 - x*x) * (1 - y*y) * z;

  val[0] = 0;
  val[1] = 0;
  val[2] = (1 - x*x) * (1 - y*y) * (1 - z*z);
#endif
}

template<typename S, typename T>
void f(S x, S y, S z, T (&val)[3]) {
  _F_
  T ev[3], dx[3], dy[3], dz[3];
  exact_sln(x, y, z, ev, dx, dy, dz);

#ifdef FN4
  T curlpart[3] = {
    2*(1 - y*y)*(1 - 2*x*x*x + x*z) + 2*(1 - z*z)*(1 - 2*x*x*x + x*z) - 6*x*y*y*(1 - z*z) - 3*y*(1 - x*x)*(1 - y*y) - 2*x*(1 - y*y)*(2*z - 3*x*y) + 4*x*z*(1 - y*y),
    2*(1 - x*x)*(y*y*y + 2*x) + 2*(1 - z*z)*(y*y*y + 2*x) + 8*x*(1 - z*z) - 3*x*(1 - x*x)*(1 - y*y) - 2*y*(1 - x*x)*(2*z - 3*x*y) - 2*y*(1 - z*z)*(z - 6*x*x),
    (1 - y*y)*(1 - z*z) + 2*(1 - x*x)*(z*z - 3*x*y*z) + 2*(1 - y*y)*(z*z - 3*x*y*z) - 6*z*y*y*(1 - x*x) - 2*z*(1 - y*y)*(z - 6*x*x) - 12*x*y*z*(1 - x*x) - 12*x*y*z*(1 - y*y)
  };
#elif defined X2_Y2_Z2
  // \nabla x \nabla x val
  T curlpart[3] = {
    4 * x * (1 - y*y) * z,
    4 * (1 - x*x) * y * z,
    2 * (1 - y*y) * (1 - z*z) + 2 * (1 - x*x) * (1 - z*z)
  };
#endif

  val[0] = curlpart[0] - alpha * ev[0];
  val[1] = curlpart[1] - alpha * ev[1];
  val[2] = curlpart[2] - alpha * ev[2];
}

// Boundary condition types.
BCType bc_types(int marker) 
{
  return BC_ESSENTIAL;
}

template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *data) {
  return
    hcurl_int_curl_u_curl_v<Real, Scalar>(n, wt, u, v, e) -
    alpha * hcurl_int_u_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Geom<Real> *e, ExtData<Scalar> *data) {
  return hcurl_int_F_v<Real, Scalar>(n, wt, f<Real, Scalar>, u, e);
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
    delete face_tables;
  }
};

struct Point {
  double ref_x, ref_y, ref_z;
  double phys_x, phys_y, phys_z;
  unsigned int elm_idx, fac_idx;

  Point(int idx, int f_idx, double rx, double ry, double rz, double px, double py, double pz) {
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

  // Initialize the space.
#if defined FN4
  Ord3 order(4, 4, 4);
#elif defined X2_Y2_Z2
  Ord3 order(2, 2, 2);
#else
  Ord3 order(2, 2, 2);
#endif
  HcurlSpace space(&mesh, bc_types, NULL, order);

  

  // Initialize the weak formulation.
  WeakForm wf(1);
  wf.add_matrix_form(0, 0, bilinear_form<double, scalar>, bilinear_form<Ord, Ord>, HERMES_SYM);
  wf.add_vector_form(0, linear_form<double, scalar>, linear_form<Ord, Ord>);

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
  Adapt *adaptivity = new Adapt(&space, HERMES_HCURL_NORM);
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
          for (int fac_idx = 0; fac_idx < Hex::NUM_FACES; fac_idx++) {
            QuadPt3D *quad_pts = my_quad.get_face_points(fac_idx, order);
            int np = my_quad.get_face_num_points(fac_idx, order);
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
    qsort((void *) points, num_points, sizeof(Point *), (compfn) compare);

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
      for (int i = 0; i < num_pairs - 1; i++) {
        Facet::Key fac_idx1 = mesh.get_facet_id(mesh.elements[points[pairs[i]]->elm_idx], points[pairs[i]]->fac_idx);
        Facet::Key fac_idx2 = mesh.get_facet_id(mesh.elements[points[pairs[i + 1]]->elm_idx], points[pairs[i + 1]]->fac_idx);
        if ((fac_idx1 == it->first) || (fac_idx2 == it->first)) {
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

      double max_difference[3] = { 0., 0.0, 0.0 };
      double max_pt_x, max_pt_y, max_pt_z, max_val_1[3], max_val_2[3];
      unsigned int max_elm_1, max_elm_2;

      RefMap *rm;
      double *nx, *ny, *nz;

      // Loop over all pairs of points, that correspond to one point in the physical domain.
      for(int pair = 0; pair < num_pairs; pair++) {
        int i = pairs[pair];

        Element *e1 = mesh.elements[points[i]->elm_idx];
        sln.set_active_element(e1);

        rm = sln.get_refmap();
        // FIXME: !!!
        QuadPt3D pt1(points[i]->ref_x, points[i]->ref_y, points[i]->ref_z, 1.0);
        rm->calc_face_normal(points[i]->fac_idx, 1, &pt1, nx, ny, nz);

        // TODO: improve me!
        double3 val1 = {
          sln.get_pt_value(points[i]->ref_x, points[i]->ref_y, points[i]->ref_z, 0),
          sln.get_pt_value(points[i]->ref_x, points[i]->ref_y, points[i]->ref_z, 1),
          sln.get_pt_value(points[i]->ref_x, points[i]->ref_y, points[i]->ref_z, 2),
        };

        double3 tp1;
        calc_tan_proj<double, double>(nx[0], ny[0], nz[0], val1, tp1);

        delete [] nx;
        delete [] ny;
        delete [] nz;


        Element *e2 = mesh.elements[points[i + 1]->elm_idx];
        sln.set_active_element(e2);

        rm = sln.get_refmap();
        QuadPt3D pt2(points[i]->ref_x, points[i]->ref_y, points[i]->ref_z, 1.0);
        rm->calc_face_normal(points[i]->fac_idx, 1, &pt2, nx, ny, nz);

        // TODO: improve me!
        double val2[3] = {
          sln.get_pt_value(points[i + 1]->ref_x, points[i + 1]->ref_y, points[i + 1]->ref_z, 0),
          sln.get_pt_value(points[i + 1]->ref_x, points[i + 1]->ref_y, points[i + 1]->ref_z, 1),
          sln.get_pt_value(points[i + 1]->ref_x, points[i + 1]->ref_y, points[i + 1]->ref_z, 2),
        };

        double3 tp2;
        calc_tan_proj(nx[0], ny[0], nz[0], val2, tp2);

        delete [] nx;
        delete [] ny;
        delete [] nz;

        //value in this point should be the same, no matter from which element we go
        double difference[3] = {
          fabs(tp1[0] - tp2[0]),
          fabs(tp1[1] - tp2[1]),
          fabs(tp1[2] - tp2[2])
        };


        double norm = sqrt(sqr(difference[0]) + sqr(difference[1] + sqr(difference[2])));
        double md_norm = sqrt(sqr(max_difference[0]) + sqr(max_difference[1] + sqr(max_difference[2])));
        if (norm > md_norm) {
          max_difference[0] = difference[0];
          max_difference[1] = difference[1];
          max_difference[2] = difference[2];

          max_pt_x = points[i]->phys_x;
          max_pt_y = points[i]->phys_y;
          max_pt_z = points[i]->phys_z;
          memcpy(max_val_1, val1, 3 * sizeof(double));
          memcpy(max_val_2, val2, 3 * sizeof(double));
          max_elm_1 = points[i]->elm_idx;
          max_elm_2 = points[i + 1]->elm_idx;
        }
      }

      double md_norm = sqrt(sqr(max_difference[0]) + sqr(max_difference[1] + sqr(max_difference[2])));
      if (md_norm > TOLERANCE) {
        info("base fn %d NOT continuous between elements %ld and %ld @ (% lf, % lf, % lf), "
          "max difference [%g, %g, %g] ([%.15g, %.15g, %.15g] <-> [%.15g, %.15g, %.15g]).",
           dof, max_elm_1, max_elm_2 , max_pt_x, max_pt_y, max_pt_z, max_difference[0], max_difference[1], max_difference[2],
           max_val_1[0], max_val_1[1], max_val_1[2], max_val_2[0], max_val_2[1], max_val_2[2]);
        success_test = 0;
      }
      else
        info("ok");
    }

    for (int i = 0; i < num_points; i++) delete points[i];
    delete [] pairs;
    delete [] points;

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

