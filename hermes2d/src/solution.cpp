// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#include "h2d_common.h"
#include "solution.h"
#include "../../hermes_common/matrix.h"
#include "precalc.h"
#include "refmap.h"
#include "auto_local_array.h"

//// MeshFunction //////////////////////////////////////////////////////////////////////////////////

MeshFunction::MeshFunction()
            : ScalarFunction()
{
  refmap = new RefMap;
  mesh = NULL;
  element = NULL;
}

MeshFunction::MeshFunction(Mesh *mesh) :
	ScalarFunction()
{
	this->mesh = mesh;
	this->refmap = new RefMap;
	// FIXME - this was in H3D: MEM_CHECK(this->refmap);
	this->element = NULL;		// this comes with Transformable
}

MeshFunction::~MeshFunction()
{
  delete refmap;
  if(overflow_nodes != NULL) {
    for(std::map<unsigned int, Node*>::iterator it = overflow_nodes->begin(); it != overflow_nodes->end(); it++)
      ::free(it->second);
    delete overflow_nodes;
  }
}


void MeshFunction::set_quad_2d(Quad2D* quad_2d)
{
  ScalarFunction::set_quad_2d(quad_2d);
  refmap->set_quad_2d(quad_2d);
}


void MeshFunction::set_active_element(Element* e)
{
  element = e;
  mode = e->get_mode();
  refmap->set_active_element(e);
  reset_transform();
}

void MeshFunction::handle_overflow_idx()
{
  if(overflow_nodes != NULL) {
    for(std::map<unsigned int, Node*>::iterator it = overflow_nodes->begin(); it != overflow_nodes->end(); it++)
      ::free(it->second);
    delete overflow_nodes;
  }
  nodes = new std::map<unsigned int, Node*>;
  overflow_nodes = nodes;
}

//// Quad2DCheb ////////////////////////////////////////////////////////////////////////////////////

static double3* cheb_tab_tri[11];
static double3* cheb_tab_quad[11];
static int      cheb_np_tri[11];
static int      cheb_np_quad[11];

static double3** cheb_tab[2] = { cheb_tab_tri, cheb_tab_quad };
static int*      cheb_np[2]  = { cheb_np_tri,  cheb_np_quad  };

/// Quad2DCheb is a special "quadrature" consisting of product Chebyshev
/// points on the reference triangle and quad. It is used for expressing
/// the solution on an element as a linear combination of monomials.
///
static class Quad2DCheb : public Quad2D
{
public:

  Quad2DCheb()
  {
    mode = HERMES_MODE_TRIANGLE;
    max_order[0]  = max_order[1]  = 10;
    num_tables[0] = num_tables[1] = 11;
    tables = cheb_tab;
    np = cheb_np;

    tables[0][0] = tables[1][0] = NULL;
    np[0][0] = np[1][0] = 0;

    int i, j, k, n, m;
    double3* pt;
    for (mode = 0; mode <= 1; mode++)
    {
      for (k = 0; k <= 10; k++)
      {
        np[mode][k] = n = mode ? sqr(k+1) : (k+1)*(k+2)/2;
        tables[mode][k] = pt = new double3[n];

        for (i = k, m = 0; i >= 0; i--)
          for (j = k; j >= (mode ? 0 : k-i); j--, m++) {
            pt[m][0] = k ? cos(j * M_PI / k) : 1.0;
            pt[m][1] = k ? cos(i * M_PI / k) : 1.0;
            pt[m][2] = 1.0;
          }
      }
    }
  };

  ~Quad2DCheb()
  {
    for (int mode = 0; mode <= 1; mode++)
      for (int k = 1; k <= 10; k++)
        delete[] tables[mode][k];
  }

  virtual void dummy_fn() {}

} g_quad_2d_cheb;


//// Solution //////////////////////////////////////////////////////////////////////////////////////

//  The higher-order solution on elements is best calculated not as a linear  combination
//  of shape functions (the usual approach), but as a linear combination of monomials.
//  This has the advantage that no shape function table calculation and look-ups are
//  necessary (except for the conversion of the solution coefficients), which means that
//  visualization and multi-mesh calculations are much faster (all the push_transforms
//  and table searches take the most time when evaluating the solution).
//
//  The linear combination of monomials can be calculated using the Horner's scheme, which
//  requires the same number of multiplications as the calculation of the linear combination
//  of shape functions. However, sub-element transforms are trivial and cheap. Moreover,
//  after the solution on all elements is expressed as a combination of monomials, the
//  Space can be forgotten. This is comfortable for the user, since the Solution class acts
//  as a self-contained unit, internally containing just a copy of the mesh and a table of
//  monomial coefficients. It is also very straight-forward to obtain all derivatives of
//  a solution defined in this way. Finally, it is possible to store the solution on the
//  disk easily (no need to store the Space, which is difficult).
//
//  The following is an example of the set of monomials for a cubic quad and a cubic triangle.
//  (Note that these are actually the definitions of the polynomial spaces on these elements.)
//
//    [ x^3*y^3  x^2*y^3  x*y^3  y^3 ]       [                    y^3 ]
//    [ x^3*y^2  x^2*y^2  x*y^2  y^2 ]       [             x*y^2  y^2 ]
//    [ x^3*y    x^2*y    x*y    y   ]       [      x^2*y  x*y    y   ]
//    [ x^3      x^2      x      1   ]       [ x^3  x^2    x      1   ]
//
//  (The number of monomials is (n+1)^2 for quads and (n+1)*(n+2)/2 for triangles, where
//   'n' is the polynomial degree.)
//

void Solution::init()
{
  memset(tables, 0, sizeof(tables));
  memset(elems,  0, sizeof(elems));
  memset(oldest, 0, sizeof(oldest));
  transform = true;
  sln_type = HERMES_UNDEF;
  own_mesh = false;
  num_components = 0;
  e_last = NULL;
  exact_mult = 1.0;
  
  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++)
      tables[i][j] = new std::map<uint64_t, std::map<unsigned int, Node*>*>;
  
  mono_coefs = NULL;
  elem_coefs[0] = elem_coefs[1] = NULL;
  elem_orders = NULL;
  dxdy_buffer = NULL;
  num_coefs = num_elems = 0;
  num_dofs = -1;

  set_quad_2d(&g_quad_2d_std);
}

Solution::Solution()
        : MeshFunction()
{
  space_type = HERMES_INVALID_SPACE;
  this->init();
}

Solution::Solution(Mesh *mesh) : MeshFunction(mesh) 
{
  space_type = HERMES_INVALID_SPACE;
  this->init();
  this->mesh = mesh;
  this->own_mesh = false;
}

Solution::Solution(Mesh *mesh, ExactFunction exactfn) : MeshFunction(mesh) 
{
  space_type = HERMES_INVALID_SPACE;
  this->init();
  this->mesh = mesh;
  this->own_mesh = false;
  this->set_exact(mesh, exactfn);
}

Solution::Solution(Space* s, Vector* coeff_vec) : MeshFunction(s->get_mesh()) 
{
  space_type = s->get_type();
  this->init();
  this->mesh = s->get_mesh();
  this->own_mesh = false;
  Solution::vector_to_solution(coeff_vec, s, this);
}

Solution::Solution(Space* s, scalar* coeff_vec) : MeshFunction(s->get_mesh()) 
{
  space_type = s->get_type();
  this->init();
  this->mesh = s->get_mesh();
  this->own_mesh = false;
  Solution::vector_to_solution(coeff_vec, s, this);
}

void Solution::assign(Solution* sln)
{
  if (sln->sln_type == HERMES_UNDEF) error("Solution being assigned is uninitialized.");
  if (sln->sln_type != HERMES_SLN) { copy(sln); return; }

  free();

  mesh = sln->mesh;
  own_mesh = sln->own_mesh;
  sln->own_mesh = false;

  mono_coefs = sln->mono_coefs;        sln->mono_coefs = NULL;
  elem_coefs[0] = sln->elem_coefs[0];  sln->elem_coefs[0] = NULL;
  elem_coefs[1] = sln->elem_coefs[1];  sln->elem_coefs[1] = NULL;
  elem_orders = sln->elem_orders;      sln->elem_orders = NULL;
  dxdy_buffer = sln->dxdy_buffer;      sln->dxdy_buffer = NULL;
  num_coefs = sln->num_coefs;          sln->num_coefs = 0;
  num_elems = sln->num_elems;          sln->num_elems = 0;

  sln_type = sln->sln_type;
  space_type = sln->get_space_type();
  num_components = sln->num_components;

  sln->sln_type = HERMES_UNDEF;
  memset(sln->tables, 0, sizeof(sln->tables));
}


void Solution::copy(const Solution* sln)
{
  if (sln->sln_type == HERMES_UNDEF) error("Solution being copied is uninitialized.");

  free();

  mesh = new Mesh;
  //printf("Copying mesh from Solution and setting own_mesh = true.\n");
  mesh->copy(sln->mesh);
  own_mesh = true;

  sln_type = sln->sln_type;
  space_type = sln->get_space_type();
  num_components = sln->num_components;
  num_dofs = sln->num_dofs;

  if (sln->sln_type == HERMES_SLN) // standard solution: copy coefficient arrays
  {
    num_coefs = sln->num_coefs;
    num_elems = sln->num_elems;

    mono_coefs = new scalar[num_coefs];
    memcpy(mono_coefs, sln->mono_coefs, sizeof(scalar) * num_coefs);

    for (int l = 0; l < num_components; l++) {
      elem_coefs[l] = new int[num_elems];
      memcpy(elem_coefs[l], sln->elem_coefs[l], sizeof(int) * num_elems);
    }

    elem_orders = new int[num_elems];
    memcpy(elem_orders, sln->elem_orders, sizeof(int) * num_elems);

    init_dxdy_buffer();
  }
  else // exact, const
  {
    exactfn1 = sln->exactfn1;
    exactfn2 = sln->exactfn2;
    cnst[0] = sln->cnst[0];
    cnst[1] = sln->cnst[1];
  }
}


void Solution::free_tables()
{
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      if(tables[i][j] != NULL) {
	std::map<uint64_t, std::map<unsigned int, Node*>*>::iterator it;
	for (it = tables[i][j]->begin(); it != tables[i][j]->end(); it++) {
	  std::map<unsigned int, Node*>::iterator it_inner;
	    for (it_inner = it->second->begin(); it_inner != it->second->end(); it_inner++)
	      ::free(it_inner->second);
	  it->second->clear();
	}
      }
}


void Solution::free()
{
  if (mono_coefs  != NULL) { delete [] mono_coefs;   mono_coefs = NULL;  }
  if (elem_orders != NULL) { delete [] elem_orders;  elem_orders = NULL; }
  if (dxdy_buffer != NULL) { delete [] dxdy_buffer;  dxdy_buffer = NULL; }

  for (int i = 0; i < num_components; i++)
    if (elem_coefs[i] != NULL)
      { delete [] elem_coefs[i];  elem_coefs[i] = NULL; }

  if (own_mesh == true && mesh != NULL)
  {
    //printf("Deleting mesh in Solution (own_mesh == true).\n");
    delete mesh;
    own_mesh = false;
  }

  e_last = NULL;

  free_tables();
}


Solution::~Solution()
{
  free();
  space_type = HERMES_INVALID_SPACE;
}


//// set_coeff_vector ///////////////////////////////////////////////////////////////////////////////

static struct mono_lu_init
{
public:

  // this is a set of LU-decomposed matrices shared by all Solutions
  double** mat[2][11];
  int* perm[2][11];

  mono_lu_init()
  {
    memset(mat, 0, sizeof(mat));
  }

  ~mono_lu_init()
  {
    for (int m = 0; m <= 1; m++)
      for (int i = 0; i <= 10; i++)
        if (mat[m][i] != NULL) {
          delete [] mat[m][i];
          delete [] perm[m][i];
        }
  }
}
mono_lu;


double** Solution::calc_mono_matrix(int o, int*& perm)
{
  int i, j, k, l, m, row;
  double x, y, xn, yn;
  int n = mode ? sqr(o+1) : (o+1)*(o+2)/2;

  // loop through all chebyshev points
  double** mat = new_matrix<double>(n, n);
  for (k = o, row = 0; k >= 0; k--) {
    y = o ? cos(k * M_PI / o) : 1.0;
    for (l = o; l >= (mode ? 0 : o-k); l--, row++) {
      x = o ? cos(l * M_PI / o) : 1.0;

      // each row of the matrix contains all the monomials x^i*y^j
      for (i = 0, yn = 1.0, m = n-1;  i <= o;  i++, yn *= y)
        for (j = (mode ? 0 : i), xn = 1.0;  j <= o;  j++, xn *= x, m--)
          mat[row][m] = xn * yn;
    }
  }

  double d;
  perm = new int[n];
  ludcmp(mat, n, perm, &d);
  return mat;
}

// using coefficient vector
void Solution::set_coeff_vector(Space* space, Vector* vec, bool add_dir_lift)
{
    // sanity check
    if (space == NULL) error("Space == NULL in Solutin::set_coeff_vector().");
    
    space_type = space->get_type();
    scalar* coeffs = new scalar [vec->length()];
    vec->extract(coeffs);
    // debug
    //printf("coeffs:\n");
    //for (int i=0; i<9; i++) printf("%g ", coeffs[i]);
    //printf("\n");
    this->set_coeff_vector(space, coeffs, add_dir_lift);
    delete [] coeffs;
}

// using coefficient array (no pss)
void Solution::set_coeff_vector(Space* space, scalar* coeffs, bool add_dir_lift)
{
    // sanity check
    if (space == NULL) error("Space == NULL in Solutin::set_coeff_vector().");
    int ndof = Space::get_num_dofs(space);

    // initialize precalc shapeset using the space's shapeset
    Shapeset *shapeset = space->get_shapeset();
    if (space->get_shapeset() == NULL) error("Space->shapeset == NULL in Solution::set_coeff_vector().");
    PrecalcShapeset *pss = new PrecalcShapeset(shapeset);
    if (pss == NULL) error("PrecalcShapeset could not be allocated in Solution::set_coeff_vector().");

    set_coeff_vector(space, pss, coeffs, add_dir_lift);

    delete pss;
}

// using pss and coefficient array
void Solution::set_coeff_vector(Space* space, PrecalcShapeset* pss, scalar* coeffs, bool add_dir_lift)
{
  int o;

  // some sanity checks
  if (space == NULL) error("Space == NULL in Solution::set_coeff_vector().");
  if (space->get_mesh() == NULL) error("Mesh == NULL in Solution::set_coeff_vector().");
  if (pss == NULL) error("PrecalcShapeset == NULL in Solution::set_coeff_vector().");
  if (coeffs == NULL) error("Coefficient vector == NULL in Solution::set_coeff_vector().");
  if (!space->is_up_to_date())
    error("Provided 'space' is not up to date.");
  if (space->get_shapeset() != pss->get_shapeset())
    error("Provided 'space' and 'pss' must have the same shapesets.");
  int ndof = Space::get_num_dofs(space);

  free();
  
  space_type = space->get_type();

  num_components = pss->get_num_components();
  sln_type = HERMES_SLN;
  num_dofs = Space::get_num_dofs(space);

  // copy the mesh   TODO: share meshes between solutions // WHAT???
  mesh = space->get_mesh();

  // allocate the coefficient arrays
  num_elems = mesh->get_max_element_id();
  if(elem_orders != NULL)
    delete [] elem_orders;
  elem_orders = new int[num_elems];
  memset(elem_orders, 0, sizeof(int) * num_elems);
  for (int l = 0; l < num_components; l++) {
    if(elem_coefs[l] != NULL)
      delete [] elem_coefs[l];
    elem_coefs[l] = new int[num_elems];
    memset(elem_coefs[l], 0, sizeof(int) * num_elems);
  }

  // obtain element orders, allocate mono_coefs
  Element* e;
  num_coefs = 0;
  for_all_active_elements(e, mesh)
  {
    mode = e->get_mode();
    o = space->get_element_order(e->id);
    o = std::max(H2D_GET_H_ORDER(o), H2D_GET_V_ORDER(o));
    for (unsigned int k = 0; k < e->nvert; k++) {
      int eo = space->get_edge_order(e, k);
      if (eo > o) o = eo;
    }

    // Hcurl: actual order of functions is one higher than element order
    if ((space->get_shapeset())->get_num_components() == 2) o++;

    num_coefs += mode ? sqr(o+1) : (o+1)*(o+2)/2;
    elem_orders[e->id] = o;
  }
  num_coefs *= num_components;
  if(mono_coefs != NULL)
    delete [] mono_coefs;
  mono_coefs = new scalar[num_coefs];

  // express the solution on elements as a linear combination of monomials
  Quad2D* quad = &g_quad_2d_cheb;
  pss->set_quad_2d(quad);
  scalar* mono = mono_coefs;
  for_all_active_elements(e, mesh)
  {
    mode = e->get_mode();
    quad->set_mode(mode);
    o = elem_orders[e->id];
    int np = quad->get_num_points(o);

    AsmList al;
    space->get_element_assembly_list(e, &al);
    pss->set_active_element(e);

    for (int l = 0; l < num_components; l++)
    {
      // obtain solution values for the current element
      scalar* val = mono;
      elem_coefs[l][e->id] = (int) (mono - mono_coefs);
      memset(val, 0, sizeof(scalar)*np);
      for (int k = 0; k < al.cnt; k++)
      {
        pss->set_active_shape(al.idx[k]);
        pss->set_quad_order(o, H2D_FN_VAL);
        int dof = al.dof[k];
        double dir_lift_coeff = add_dir_lift ? 1.0 : 0.0;
        scalar coef = al.coef[k] * (dof >= 0 ? coeffs[dof] : dir_lift_coeff);
        double* shape = pss->get_fn_values(l);
        for (int i = 0; i < np; i++)
          val[i] += shape[i] * coef;
      }
      mono += np;

      // solve for the monomial coefficients
      if (mono_lu.mat[mode][o] == NULL)
        mono_lu.mat[mode][o] = calc_mono_matrix(o, mono_lu.perm[mode][o]);
      lubksb(mono_lu.mat[mode][o], np, mono_lu.perm[mode][o], val);
    }
  }

  if(mesh == NULL) error("mesh == NULL.\n");
  init_dxdy_buffer();
}


//// set_exact etc. ////////////////////////////////////////////////////////////////////////////////

void Solution::set_exact(Mesh* mesh, ExactFunction exactfn)
{
  free();
  this->mesh = mesh;
  exactfn1 = exactfn;
  num_components = 1;
  sln_type = HERMES_EXACT;
  exact_mult = 1.0;
  num_dofs = -1;
}


void Solution::set_exact(Mesh* mesh, ExactFunction2 exactfn)
{
  free();
  this->mesh = mesh;
  exactfn2 = exactfn;
  num_components = 2;
  sln_type = HERMES_EXACT;
  exact_mult = 1.0;
  num_dofs = -1;
}


void Solution::set_const(Mesh* mesh, scalar c)
{
  free();
  this->mesh = mesh;
  cnst[0] = c;
  cnst[1] = 0.0;
  num_components = 1;
  sln_type = HERMES_CONST;
  num_dofs = -1;
}


void Solution::set_const(Mesh* mesh, scalar c0, scalar c1)
{
  free();
  this->mesh = mesh;
  cnst[0] = c0;
  cnst[1] = c1;
  num_components = 2;
  sln_type = HERMES_CONST;
  num_dofs = -1;
}


void Solution::set_zero(Mesh* mesh)
{
  set_const(mesh, 0.0);
}


void Solution::set_zero_2(Mesh* mesh)
{
  set_const(mesh, 0.0, 0.0);
}

void Solution::vector_to_solutions(scalar* solution_vector, 
                                   Hermes::Tuple<Space*> spaces, 
                                   Hermes::Tuple<Solution*> solutions, 
                                   Hermes::Tuple<bool> add_dir_lift)
{
  assert(spaces.size() == solutions.size());
  for(unsigned int i = 0; i < solutions.size(); i++)
    if(add_dir_lift == Hermes::Tuple<bool>())
      solutions[i]->set_coeff_vector(spaces[i], solution_vector, true);
    else
      solutions[i]->set_coeff_vector(spaces[i], solution_vector,
              add_dir_lift.at(i));
  return;
}

void Solution::vector_to_solution(scalar* solution_vector, Space* space, 
                                  Solution* solution, bool add_dir_lift)
{
  Solution::vector_to_solutions(solution_vector, Hermes::Tuple<Space*>(space), 
                                Hermes::Tuple<Solution*>(solution), 
                                Hermes::Tuple<bool>(add_dir_lift));
}

void Solution::vector_to_solutions(Vector* solution_vector, Hermes::Tuple<Space*> spaces, 
                                   Hermes::Tuple<Solution*> solutions, 
                                   Hermes::Tuple<bool> add_dir_lift)
{
  assert(spaces.size() == solutions.size());
  for(unsigned int i = 0; i < solutions.size(); i++)
    if(add_dir_lift == Hermes::Tuple<bool>())
      solutions[i]->set_coeff_vector(spaces[i], solution_vector, true);
    else
      solutions[i]->set_coeff_vector(spaces[i], solution_vector,
              add_dir_lift.at(i));
  return;
}

void Solution::vector_to_solution(Vector* solution_vector, Space* space, 
                                  Solution* solution, bool add_dir_lift)
{
  Solution::vector_to_solutions(solution_vector, Hermes::Tuple<Space*>(space), 
                                Hermes::Tuple<Solution*>(solution), 
                                Hermes::Tuple<bool>(add_dir_lift));
}

void Solution::vector_to_solutions(scalar* solution_vector, Hermes::Tuple<Space*> spaces, 
                                   Hermes::Tuple<Solution*> solutions, 
                                   Hermes::Tuple<PrecalcShapeset *> pss, 
                                   Hermes::Tuple<bool> add_dir_lift)
{
  assert(spaces.size() == solutions.size());
  for(unsigned int i = 0; i < solutions.size(); i++)
    if(add_dir_lift == Hermes::Tuple<bool>())
      solutions[i]->set_coeff_vector(spaces[i], pss[i], solution_vector, true);
    else
      solutions[i]->set_coeff_vector(spaces[i], pss[i], solution_vector,
              add_dir_lift.at(i));
  return;
}

void Solution::vector_to_solution(scalar* solution_vector, Space* space, Solution* solution, 
                                  PrecalcShapeset* pss, bool add_dir_lift)
{
  Solution::vector_to_solutions(solution_vector, Hermes::Tuple<Space*>(space), 
                                Hermes::Tuple<Solution*>(solution), 
                                Hermes::Tuple<PrecalcShapeset *>(pss), 
                                Hermes::Tuple<bool>(add_dir_lift));
}


void Solution::set_dirichlet_lift(Space* space, PrecalcShapeset* pss)
{
  space_type = space->get_type();
  int ndof = Space::get_num_dofs(space);
  scalar *temp = new scalar[ndof];
  memset(temp, 0, sizeof(scalar)*ndof);
  this->set_coeff_vector(space, pss, temp, true);
  delete [] temp;
}

void Solution::enable_transform(bool enable)
{
  if (transform != enable) free_tables();
  transform = enable;
}


void Solution::multiply(scalar coef)
{
  if (sln_type == HERMES_SLN)
  {
    for (int i = 0; i < num_coefs; i++)
      mono_coefs[i] *= coef;
  }
  else if (sln_type == HERMES_CONST)
  {
    cnst[0] *= coef;
    cnst[1] *= coef;
  }
  else if (sln_type == HERMES_EXACT)
  {
    exact_mult *= coef;
  }
  else
    error("Uninitialized solution.");
}


//// set_active_element ////////////////////////////////////////////////////////////////////////////

// differentiates the mono coefs by x
static void make_dx_coefs(int mode, int o, scalar* mono, scalar* result)
{
  int i, j, k;
  for (i = 0; i <= o; i++) {
    *result++ = 0.0;
    k = mode ? o : i;
    for (j = 0; j < k; j++)
      *result++ = (scalar) (k-j) * mono[j];
    mono += k+1;
  }
}

// differentiates the mono coefs by y
static void make_dy_coefs(int mode, int o, scalar* mono, scalar* result)
{
  int i, j;
  if (mode) {
    for (j = 0; j <= o; j++)
      *result++ = 0.0;
    for (i = 0; i < o; i++)
      for (j = 0; j <= o; j++)
        *result++ = (scalar) (o-i) * (*mono++);
  }
  else {
    for (i = 0; i <= o; i++) {
      *result++ = 0.0;
      for (j = 0; j < i; j++)
        *result++ = (scalar) (o+1-i) * (*mono++);
    }
  }
}

void Solution::init_dxdy_buffer()
{
  if (dxdy_buffer != NULL) {
    delete [] dxdy_buffer;
    dxdy_buffer = NULL;
  }
  dxdy_buffer = new scalar[num_components * 5 * sqr(11)];
}


void Solution::set_active_element(Element* e)
{
  // if (e == element) return; // FIXME
  if (!e->active) error("Cannot select inactive element. Wrong mesh?");
  MeshFunction::set_active_element(e);

  // try finding an existing table for e
  for (cur_elem = 0; cur_elem < 4; cur_elem++)
    if (elems[cur_quad][cur_elem] == e)
      break;

  // if not found, free the oldest one and use its slot
  if (cur_elem >= 4)
  {
    if (tables[cur_quad][oldest[cur_quad]] != NULL) {
      std::map<uint64_t, std::map<unsigned int, Node*>*>::iterator it;
      for (it = tables[cur_quad][oldest[cur_quad]]->begin(); it != tables[cur_quad][oldest[cur_quad]]->end(); it++) {
	std::map<unsigned int, Node*>::iterator it_inner;
	  for (it_inner = it->second->begin(); it_inner != it->second->end(); it_inner++)
	    ::free(it_inner->second);
	it->second->clear();
      }
    }

    cur_elem = oldest[cur_quad];
    if (++oldest[cur_quad] >= 4)
      oldest[cur_quad] = 0;

    elems[cur_quad][cur_elem] = e;
  }

  if (sln_type == HERMES_SLN)
  {
    int o = order = elem_orders[element->id];
    int n = mode ? sqr(o+1) : (o+1)*(o+2)/2;

    for (int i = 0, m = 0; i < num_components; i++)
    {
      scalar* mono = mono_coefs + elem_coefs[i][e->id];
      dxdy_coefs[i][0] = mono;

      make_dx_coefs(mode, o, mono, dxdy_coefs[i][1] = dxdy_buffer+m);  m += n;
      make_dy_coefs(mode, o, mono, dxdy_coefs[i][2] = dxdy_buffer+m);  m += n;
      make_dx_coefs(mode, o, dxdy_coefs[i][1], dxdy_coefs[i][3] = dxdy_buffer+m);  m += n;
      make_dy_coefs(mode, o, dxdy_coefs[i][2], dxdy_coefs[i][4] = dxdy_buffer+m);  m += n;
      make_dx_coefs(mode, o, dxdy_coefs[i][2], dxdy_coefs[i][5] = dxdy_buffer+m);  m += n;
    }
  }
  else if (sln_type == HERMES_EXACT)
  {
    order = 20; // fixme
  }
  else if (sln_type == HERMES_CONST)
  {
    order = 0;
  }
  else
    error("Uninitialized solution.");

  sub_tables = tables[cur_quad][cur_elem];

  update_nodes_ptr();
}


//// precalculate //////////////////////////////////////////////////////////////////////////////////

// sets all elements of y[] to num
static inline void set_vec_num(int n, scalar* y, scalar num)
{
  for (int i = 0; i < n; i++)
    y[i] = num;
}

// y = y .* x + num
static inline void vec_x_vec_p_num(int n, scalar* y, scalar* x, scalar num)
{
  for (int i = 0; i < n; i++)
    y[i] = y[i]*x[i] + num;
}

// y = y .* x + z
static inline void vec_x_vec_p_vec(int n, scalar* y, scalar* x, scalar* z)
{
  for (int i = 0; i < n; i++)
    y[i] = y[i]*x[i] + z[i];
}


static const int H2D_GRAD = H2D_FN_DX_0 | H2D_FN_DY_0;
static const int H2D_SECOND = H2D_FN_DXX_0 | H2D_FN_DXY_0 | H2D_FN_DYY_0;
static const int H2D_CURL = H2D_FN_DX | H2D_FN_DY;


void Solution::transform_values(int order, Node* node, int newmask, int oldmask, int np)
{
  double2x2 *mat, *m;
  double3x2 *mat2, *mm;
  int i, mstep = 0;

  // H1 space
  if (space_type == HERMES_H1_SPACE)
  {
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
    if (((newmask & H2D_SECOND) == H2D_SECOND && (oldmask & H2D_SECOND) != H2D_SECOND))
    {
      update_refmap();
      mat = refmap->get_inv_ref_map(order);
      mat2 = refmap->get_second_ref_map(order);
      for (i = 0, m = mat, mm = mat2; i < np; i++, m++, mm++)
      {
        scalar vx = node->values[0][1][i];
        scalar vy = node->values[0][2][i];
        scalar vxx = node->values[0][3][i];
        scalar vyy = node->values[0][4][i];
        scalar vxy = node->values[0][5][i];

        node->values[0][3][i] = sqr((*m)[0][0])*vxx + 2*(*m)[0][1]*(*m)[0][0]*vxy + sqr((*m)[0][1])*vyy + (*mm)[0][0]*vx + (*mm)[0][1]*vy;   // dxx
        node->values[0][4][i] = sqr((*m)[1][0])*vxx + 2*(*m)[1][1]*(*m)[1][0]*vxy + sqr((*m)[1][1])*vyy + (*mm)[2][0]*vx + (*mm)[2][1]*vy;   // dyy
        node->values[0][5][i] = (*m)[0][0]*(*m)[1][0]*vxx + ((*m)[0][0]*(*m)[1][1]+(*m)[1][0]*(*m)[0][1])*vxy + (*m)[0][1]*(*m)[1][1]*vyy + (*mm)[1][0]*vx + (*mm)[1][1]*vy;   //dxy
      }
    }
#endif
    if ((newmask & H2D_GRAD) == H2D_GRAD && (oldmask & H2D_GRAD) != H2D_GRAD)
    {
      update_refmap();
      mat = refmap->get_const_inv_ref_map();
      if (!refmap->is_jacobian_const()) { mat = refmap->get_inv_ref_map(order); mstep = 1; }

      for (i = 0, m = mat; i < np; i++, m += mstep)
      {
        scalar vx = node->values[0][1][i];
        scalar vy = node->values[0][2][i];
        node->values[0][1][i] = (*m)[0][0]*vx + (*m)[0][1]*vy;
        node->values[0][2][i] = (*m)[1][0]*vx + (*m)[1][1]*vy;
      }
    }
  }

  // Hcurl space
  else if (space_type == HERMES_HCURL_SPACE)
  {
    bool trans_val = false, trans_curl = false;
    if ((newmask & H2D_FN_VAL) == H2D_FN_VAL && (oldmask & H2D_FN_VAL) != H2D_FN_VAL) trans_val  = true;
    if ((newmask &   H2D_CURL) ==   H2D_CURL && (oldmask &   H2D_CURL) !=   H2D_CURL) trans_curl = true;

    if (trans_val || trans_curl)
    {
      update_refmap();
      mat = refmap->get_const_inv_ref_map();
      if (!refmap->is_jacobian_const()) { mat = refmap->get_inv_ref_map(order); mstep = 1; }

      for (i = 0, m = mat; i < np; i++, m += mstep)
      {
        if (trans_val)
        {
          scalar vx = node->values[0][0][i];
          scalar vy = node->values[1][0][i];
          node->values[0][0][i] = (*m)[0][0]*vx + (*m)[0][1]*vy;
          node->values[1][0][i] = (*m)[1][0]*vx + (*m)[1][1]*vy;
        }
        if (trans_curl)
        {
          scalar e0x = node->values[0][1][i], e0y = node->values[0][2][i];
          scalar e1x = node->values[1][1][i], e1y = node->values[1][2][i];
          node->values[1][1][i] = (*m)[0][0]*((*m)[1][0]*e0x + (*m)[1][1]*e1x) + (*m)[0][1]*((*m)[1][0]*e0y + (*m)[1][1]*e1y);
          node->values[0][2][i] = (*m)[1][0]*((*m)[0][0]*e0x + (*m)[0][1]*e1x) + (*m)[1][1]*((*m)[0][0]*e0y + (*m)[0][1]*e1y);
        }
      }
    }
  }

  // Hdiv space
  else if (space_type == HERMES_HDIV_SPACE)
  {
    if ((newmask & H2D_FN_VAL) == H2D_FN_VAL && (oldmask & H2D_FN_VAL) != H2D_FN_VAL)
    {
      update_refmap();
      mat = refmap->get_const_inv_ref_map();
      if (!refmap->is_jacobian_const()) { mat = refmap->get_inv_ref_map(order); mstep = 1; }

      for (i = 0, m = mat; i < np; i++, m += mstep)
      {
        scalar vx = node->values[0][0][i];
        scalar vy = node->values[1][0][i];
        node->values[0][0][i] =   (*m)[1][1]*vx - (*m)[1][0]*vy;
        node->values[1][0][i] = - (*m)[0][1]*vx + (*m)[0][0]*vy;
      }
    }
  }
}

int Solution::get_edge_fn_order(int edge, Space* space, Element* e)
{
  if (e == NULL) e = element;
  
  if (sln_type == HERMES_SLN && space != NULL) {
    return space->get_edge_order(e, edge); 
  } else {
    return ScalarFunction::get_edge_fn_order(edge);
  }
}


void Solution::precalculate(int order, int mask)
{
  int i, j, k, l;
  Node* node;
  Quad2D* quad = quads[cur_quad];
  quad->set_mode(mode);
  H2D_CHECK_ORDER(quad, order);
  int np = quad->get_num_points(order);

  if (sln_type == HERMES_SLN)
  {
    // if we are required to transform vectors, we must precalculate both their components
    const int H2D_GRAD = H2D_FN_DX_0 | H2D_FN_DY_0;
    const int H2D_SECOND = H2D_FN_DXX_0 | H2D_FN_DXY_0 | H2D_FN_DYY_0;
    const int H2D_CURL = H2D_FN_DX | H2D_FN_DY; // sic
    if (transform)
    {
      if (num_components == 1)                                            // H1 space or L2 space
      {
        if ((mask & H2D_FN_DX_0)  || (mask & H2D_FN_DY_0))  mask |= H2D_GRAD;
        if ((mask & H2D_FN_DXX_0)  || (mask & H2D_FN_DXY_0) || (mask & H2D_FN_DYY_0))  mask |= H2D_SECOND;
      }
      else if (space_type == HERMES_HCURL_SPACE)                                           // Hcurl space
        { if ((mask & H2D_FN_VAL_0) || (mask & H2D_FN_VAL_1)) mask |= H2D_FN_VAL;
          if ((mask & H2D_FN_DX_1)  || (mask & H2D_FN_DY_0))  mask |= H2D_CURL; }
      else                                                                // Hdiv space
        { if ((mask & H2D_FN_VAL_0) || (mask & H2D_FN_VAL_1)) mask |= H2D_FN_VAL; }
    }

    int oldmask = (cur_node != NULL) ? cur_node->mask : 0;
    int newmask = mask | oldmask;
    node = new_node(newmask, np);

    // transform integration points by the current matrix
    AUTOLA_OR(scalar, x, np); AUTOLA_OR(scalar, y, np); AUTOLA_OR(scalar, tx, np);
    double3* pt = quad->get_points(order);
    for (i = 0; i < np; i++)
    {
      x[i] = pt[i][0] * ctm->m[0] + ctm->t[0];
      y[i] = pt[i][1] * ctm->m[1] + ctm->t[1];
    }

    // obtain the solution values, this is the core of the whole module
    int o = elem_orders[element->id];
    for (l = 0; l < num_components; l++)
    {
      for (k = 0; k < 6; k++)
      {
        if (newmask & idx2mask[k][l])
        {
          scalar* result = node->values[l][k];
          if (oldmask & idx2mask[k][l])
          {
            // copy the old table if we have it already
            memcpy(result, cur_node->values[l][k], np * sizeof(scalar));
          }
          else
          {
            // calculate the solution values using Horner's scheme
            scalar* mono = dxdy_coefs[l][k];
            for (i = 0; i <= o; i++)
            {
              set_vec_num(np, tx, *mono++);
              for (j = 1; j <= (mode ? o : i); j++)
                vec_x_vec_p_num(np, tx, x, *mono++);

              if (!i) memcpy(result, tx, sizeof(scalar)*np);
                 else vec_x_vec_p_vec(np, result, y, tx);
            }
          }
        }
      }
    }

    // transform gradient or vector solution, if required
    if (transform)
      transform_values(order, node, newmask, oldmask, np);
  }
  else if (sln_type == HERMES_EXACT)
  {
    if (mask & ~H2D_FN_DEFAULT)
      error("Cannot obtain second derivatives of an exact solution.");
    node = new_node(mask = H2D_FN_DEFAULT, np);

    update_refmap();
    double* x = refmap->get_phys_x(order);
    double* y = refmap->get_phys_y(order);

    // evaluate the exact solution
    if (num_components == 1)
    {
      // untransform values
      if (!transform)
      {
        double2x2 *mat, *m;
        int mstep = 0;
        mat = refmap->get_const_inv_ref_map();
        if (!refmap->is_jacobian_const()) { mat = refmap->get_inv_ref_map(order); mstep = 1; }

        for (i = 0, m = mat; i < np; i++, m += mstep)
        {
          double jac = (*m)[0][0] *  (*m)[1][1] - (*m)[1][0] *  (*m)[0][1];
          scalar val, dx = 0.0, dy = 0.0;
          val = exactfn1(x[i], y[i], dx, dy);
          node->values[0][0][i] = val * exact_mult;
          node->values[0][1][i] = (  (*m)[1][1]*dx - (*m)[0][1]*dy) / jac * exact_mult;
          node->values[0][2][i] = (- (*m)[1][0]*dx + (*m)[0][0]*dy) / jac * exact_mult;
        }
      }
      else
      {
        for (i = 0; i < np; i++)
        {
          scalar val, dx = 0.0, dy = 0.0;
          val = exactfn1(x[i], y[i], dx, dy);
          node->values[0][0][i] = val * exact_mult;
          node->values[0][1][i] = dx * exact_mult;
          node->values[0][2][i] = dy * exact_mult;
        }
      }
    }
    else
    {
      for (i = 0; i < np; i++)
      {
        scalar2 dx = { 0.0, 0.0 }, dy = { 0.0, 0.0 };
        scalar2& val = exactfn2(x[i], y[i], dx, dy);
        for (j = 0; j < 2; j++) {
          node->values[j][0][i] = val[j] * exact_mult;
          node->values[j][1][i] = dx[j] * exact_mult;
          node->values[j][2][i] = dy[j] * exact_mult;
        }
      }
    }
  }
  else if (sln_type == HERMES_CONST)
  {
    if (mask & ~H2D_FN_DEFAULT)
      error("Second derivatives of a constant solution not implemented.");
    node = new_node(mask = H2D_FN_DEFAULT, np);

    for (j = 0; j < num_components; j++)
      for (i = 0; i < np; i++)
      {
        node->values[j][0][i] = cnst[j];
        node->values[j][1][i] = 0.0;
        node->values[j][2][i] = 0.0;
      }
  }
  else
  {
    error("Cannot obtain values -- uninitialized solution. The solution was either "
          "not calculated yet or you used the assignment operator which destroys "
          "the solution on its right-hand side.");
  }

  if((*nodes)[order] != NULL) {
    assert((*nodes)[order] == cur_node);
    ::free((*nodes)[order]);
  }
  (*nodes)[order] = node;
  cur_node = node;
}


//// save & load ///////////////////////////////////////////////////////////////////////////////////

void Solution::save(const char* filename, bool compress)
{
  int i;

  if (sln_type == HERMES_EXACT) error("Exact solution cannot be saved to a file.");
  if (sln_type == HERMES_CONST)  error("Constant solution cannot be saved to a file.");
  if (sln_type == HERMES_UNDEF) error("Cannot save -- uninitialized solution.");

  // open the stream
  std::string fname = filename;
  if (compress) fname += ".gz";
  FILE* f = fopen(fname.c_str(), "wb");
  if (f == NULL) error("Could not open %s for writing.", filename);

  if (compress)
  {
    fclose(f);
    std::stringstream cmdline;
    cmdline << "gzip > " << filename << ".gz";
    f = popen(cmdline.str().c_str(), "w");
    if (f == NULL) error("Could not create compressed stream (command line: %s).", cmdline.str().c_str());
  }

  // write header
  hermes_fwrite("H2DS\001\000\000\000", 1, 8, f);
  int ssize = sizeof(scalar);
  hermes_fwrite(&ssize, sizeof(int), 1, f);
  hermes_fwrite(&num_components, sizeof(int), 1, f);
  hermes_fwrite(&num_elems, sizeof(int), 1, f);
  hermes_fwrite(&num_coefs, sizeof(int), 1, f);

  // write monomial coefficients
  hermes_fwrite(mono_coefs, sizeof(scalar), num_coefs, f);

  // write element orders
  char* temp_orders = new char[num_elems];
  for (i = 0; i < num_elems; i++) {
    temp_orders[i] = elem_orders[i];
  }
  hermes_fwrite(temp_orders, sizeof(char), num_elems, f);
  delete [] temp_orders;

  // write element coef table
  for (i = 0; i < num_components; i++)
    hermes_fwrite(elem_coefs[i], sizeof(int), num_elems, f);

  // write the mesh
  mesh->save_raw(f);

  if (compress) pclose(f); else fclose(f);
}


void Solution::load(const char* filename)
{
  int i;

  free();
  sln_type = HERMES_SLN;

  int len = strlen(filename);
  bool compressed = (len > 3 && !strcmp(filename + len - 3, ".gz"));

  // open the stream
  FILE* f = fopen(filename, "rb");
  if (f == NULL) error("Could not open %s", filename);

  if (compressed)
  {
    fclose(f);
    std::stringstream cmdline;
    cmdline << "gunzip < " << filename << ".gz";
    f = popen(cmdline.str().c_str(), "r");
    if (f == NULL) error("Could not read from compressed stream (command line: %s).", cmdline.str().c_str());
  }

  // load header
  struct {
    char magic[4];
    int  ver, ss, nc, ne, nf;
  } hdr;
  hermes_fread(&hdr, sizeof(hdr), 1, f);

  // some checks
  if (hdr.magic[0] != 'H' || hdr.magic[1] != '2' || hdr.magic[2] != 'D' || hdr.magic[3] != 'S')
    error("Not a Hermes2D solution file.");
  if (hdr.ver > 1)
    error("Unsupported file version.");

  // load monomial coefficients
  num_coefs = hdr.nf;
  if (hdr.ss == sizeof(double))
  {
    double* temp = new double[num_coefs];
    hermes_fread(temp, sizeof(double), num_coefs, f);

    #ifndef H2D_COMPLEX
      mono_coefs = temp;
    #else
      mono_coefs = new scalar[num_coefs];
      for (i = 0; i < num_coefs; i++)
        mono_coefs[i] = temp[i];
      delete [] temp;
    #endif
  }
  else if (hdr.ss == 2*sizeof(double))
  {
    #ifndef H2D_COMPLEX
      warn("Ignoring imaginary part of the complex solution since this is not H2D_COMPLEX code.");
      scalar* temp = new double[num_coefs*2];
      hermes_fread(temp, sizeof(scalar), num_coefs*2, f);
      mono_coefs = new double[num_coefs];
      for (i = 0; i < num_coefs; i++)
        mono_coefs[i] = temp[2*i];
      delete [] temp;

    #else
      mono_coefs = new scalar[num_coefs];;
      hermes_fread(mono_coefs, sizeof(scalar), num_coefs, f);
    #endif
  }
  else
    error("Corrupt solution file.");

  // load element orders
  num_elems = hdr.ne;
  char* temp_orders = new char[num_elems];
  hermes_fread(temp_orders, sizeof(char), num_elems, f);
  elem_orders = new int[num_elems];
  for (i = 0; i < num_elems; i++)
    elem_orders[i] = temp_orders[i];
  delete [] temp_orders;

  // load element coef table
  num_components = hdr.nc;
  for (i = 0; i < num_components; i++)
  {
    elem_coefs[i] = new int[num_elems];
    hermes_fread(elem_coefs[i], sizeof(int), num_elems, f);
  }

  // load the mesh
  mesh = new Mesh;
  mesh->load_raw(f);
  //printf("Loading mesh from file and setting own_mesh = true.\n");
  own_mesh = true;

  if (compressed) pclose(f); else fclose(f);

  init_dxdy_buffer();
}


//// getting solution values in arbitrary points ///////////////////////////////////////////////////////////////

scalar Solution::get_ref_value(Element* e, double xi1, double xi2, int component, int item)
{
  set_active_element(e);

  int o = elem_orders[e->id];
  scalar* mono = dxdy_coefs[component][item];
  scalar result = 0.0;
  int k = 0;
  for (int i = 0; i <= o; i++)
  {
    scalar row = mono[k++];
    for (int j = 0; j < (mode ? o : i); j++)
      row = row * xi1 + mono[k++];
    result = result * xi2 + row;
  }
  return result;
}


static inline bool is_in_ref_domain(Element* e, double xi1, double xi2)
{
  const double TOL = 1e-11;
  if (e->is_triangle())
    return (xi1 + xi2 <= TOL) && (xi1 + 1.0 >= -TOL) && (xi2 + 1.0 >= -TOL);
  else
    return (xi1 - 1.0 <= TOL) && (xi1 + 1.0 >= -TOL) && (xi2 - 1.0 <= TOL) && (xi2 + 1.0 >= -TOL);
}


scalar Solution::get_ref_value_transformed(Element* e, double xi1, double xi2, int a, int b)
{

  if (num_components == 1)
  {
    if (b == 0)
      return get_ref_value(e, xi1, xi2, a, b);
    if (b == 1 || b == 2)
    {
      double2x2 m;
      double xx, yy;
      refmap->inv_ref_map_at_point(xi1, xi2, xx, yy, m);
      scalar dx = get_ref_value(e_last = e, xi1, xi2, a, 1);
      scalar dy = get_ref_value(e, xi1, xi2, a, 2);
      if (b == 1) return m[0][0]*dx + m[0][1]*dy; // H2D_FN_DX
      if (b == 2) return m[1][0]*dx + m[1][1]*dy; // H2D_FN_DY
    }
    else
      error("Getting second derivatives of the solution: Not implemented yet.");
  }
  else // vector solution
  {
    if (b == 0)
    {
      double2x2 m;
      double xx, yy;
      refmap->inv_ref_map_at_point(xi1, xi2, xx, yy, m);
      scalar vx = get_ref_value(e, xi1, xi2, 0, 0);
      scalar vy = get_ref_value(e, xi1, xi2, 1, 0);
      if (a == 0) return m[0][0]*vx + m[0][1]*vy; // H2D_FN_VAL_0
      if (a == 1) return m[1][0]*vx + m[1][1]*vy; // H2D_FN_VAL_1
    }
    else
      error("Getting derivatives of the vector solution: Not implemented yet.");
  }
  error("internal error: reached end of non-void function");
  return 0;
}

scalar Solution::get_pt_value(double x, double y, int item)
{
  double xi1, xi2;

  int a = 0, b = 0, mask = item; // a = component, b = val, dx, dy, dxx, dyy, dxy
  if (num_components == 1) mask = mask & H2D_FN_COMPONENT_0;
  if ((mask & (mask - 1)) != 0) error("'item' is invalid. ");
  if (mask >= 0x40) { a = 1; mask >>= 6; }
  while (!(mask & 1)) { mask >>= 1; b++; }

  if (sln_type == HERMES_EXACT)
  {
    if (num_components == 1)
    {
      scalar val, dx = 0.0, dy = 0.0;
      val = exactfn1(x, y, dx, dy);
      if (b == 0) return val;
      if (b == 1) return dx;
      if (b == 2) return dy;
    }
    else
    {
      scalar2 dx = {0.0, 0.0}, dy = {0.0, 0.0};
      scalar2& val = exactfn2(x, y, dx, dy);
      if (b == 0) return val[a];
      if (b == 1) return dx[a];
      if (b == 2) return dy[a];
    }
    error("Cannot obtain second derivatives of an exact solution.");
  }
  else if (sln_type == HERMES_CONST)
  {
    if (b == 0) return cnst[a];
    return 0.0;
  }
  else if (sln_type == HERMES_UNDEF)
  {
    error("Cannot obtain values -- uninitialized solution. The solution was either "
          "not calculated yet or you used the assignment operator which destroys "
          "the solution on its right-hand side.");
  }

  // try the last visited element and its neighbours
  if (e_last != NULL)
  {
    Element* elem[5];
    elem[0] = e_last;
    for (unsigned int i = 1; i <= e_last->nvert; i++)
      elem[i] = e_last->get_neighbor(i-1);

    for (unsigned int i = 0; i <= e_last->nvert; i++)
      if (elem[i] != NULL)
      {
        refmap->set_active_element(elem[i]);
        refmap->untransform(elem[i], x, y, xi1, xi2);
        if (is_in_ref_domain(elem[i], xi1, xi2))
        {
          e_last = elem[i];
          return get_ref_value_transformed(elem[i], xi1, xi2, a, b);
        }
      }
  }

  // go through all elements
  Element *e;
  for_all_active_elements(e, mesh)
  {
    refmap->set_active_element(e);
    refmap->untransform(e, x, y, xi1, xi2);
    if (is_in_ref_domain(e, xi1, xi2))
    {
      e_last = e;
      return get_ref_value_transformed(e, xi1, xi2, a, b);
    }
  }

  warn("Point (%g, %g) does not lie in any element.", x, y);
  return NAN;
}

