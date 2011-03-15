// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes3D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes3D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes3D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#include "h3d_common.h"
#include "solution.h"
#include "function.cpp" // non-inline template members
#include "quadcheb.h"
#include "../../hermes_common/error.h"
#include "../../hermes_common/matrix.h"
#include "../../hermes_common/callstack.h"

//// MeshFunction //////////////////////////////////////////////////////////////////////////////////

static int g_mfn_seq = 0;

MeshFunction::MeshFunction(Mesh *mesh) :
	ScalarFunction()
{
	_F_
	this->mesh = mesh;
	this->refmap = new RefMap(mesh);
	MEM_CHECK(this->refmap);
	this->element = NULL;		// this comes with Transformable
	this->seq = g_mfn_seq++;
	this->noinc = false;
}

MeshFunction::~MeshFunction() {
	_F_
	delete refmap;
}

void MeshFunction::set_active_element(Element *e) {
	_F_
	element = e;
	mode = e->get_mode();
	refmap->set_active_element(e);
	reset_transform();
}

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
//                                                         [ x^2*y^2*z^2  x*y^2*z^2  y^2*z^2 ]
//                           [ x^2*y^2*z  x*y^2*z  y^2*z ] [ x^2*y*z^2    x*y*z^2    y*z^2   ]
//   [ x^2*y^2  x*y^2  y^2 ] [ x^2*y*z    x*y*z    y*z   ] [ x^2*z^2      x*z^2      z^2     ]
//   [ x^2*y    x*y    y   ] [ x^2*z      x*z      z     ]
//   [ x^2      x      1   ]
//
//  The number of monomials is (n+1)^3 for hexahedra and (n+1)*(n+2)*(n+3)/6 for tetrahedra,
//  'n' is the polynomial degree.
//

Solution::Solution(Mesh *mesh) : MeshFunction(mesh) {
	_F_

	transform = true;
	type = HERMES_UNDEF;
	num_components = 0;
	mono_coefs = NULL;
	elem_coefs[0] = elem_coefs[1] = elem_coefs[2] = NULL;
	elem_orders = NULL;
	dxdydz_buffer = NULL;
	num_coefs = num_elems = 0;
	num_dofs = -1;
}

Solution::~Solution() {
	_F_
	free();
}

void Solution::free() {
	_F_
	free_cur_node();

	if (mono_coefs  != NULL)   { delete [] mono_coefs;    mono_coefs = NULL;  }
	if (elem_orders != NULL)   { delete [] elem_orders;   elem_orders = NULL; }
	if (dxdydz_buffer != NULL) { delete [] dxdydz_buffer; dxdydz_buffer = NULL; }

	for (int i = 0; i < num_components; i++)
		if (elem_coefs[i] != NULL) {
			delete [] elem_coefs[i];
			elem_coefs[i] = NULL;
		}
}

void Solution::assign(Solution *sln) {
	_F_

	if (sln->type == HERMES_UNDEF) EXIT("Solution being assigned is uninitialized.");
	if (sln->type != HERMES_SLN) { copy(sln); return; }

	free();

	mesh = sln->mesh;

	mono_coefs = sln->mono_coefs;        sln->mono_coefs = NULL;
	elem_coefs[0] = sln->elem_coefs[0];  sln->elem_coefs[0] = NULL;
	elem_coefs[1] = sln->elem_coefs[1];  sln->elem_coefs[1] = NULL;
	elem_coefs[2] = sln->elem_coefs[2];  sln->elem_coefs[2] = NULL;
	elem_orders = sln->elem_orders;      sln->elem_orders = NULL;
	dxdydz_buffer = sln->dxdydz_buffer;  sln->dxdydz_buffer = NULL;
	num_coefs = sln->num_coefs;          sln->num_coefs = 0;
	num_elems = sln->num_elems;          sln->num_elems = 0;

	sptype = sln->sptype;
	type = sln->type;
	num_components = sln->num_components;
	seq = sln->seq;

	sln->type = HERMES_UNDEF;
}

void Solution::copy(const Solution *sln) {
	_F_

	if (sln->type == HERMES_UNDEF) EXIT("Solution being copied is uninitialized.");

	free();

//	mesh = new Mesh;
//	mesh->copy(sln->mesh);
//	own_mesh = true;
	mesh = sln->mesh;

	type = sln->type;
	num_components = sln->num_components;

	if (sln->type == HERMES_SLN) { // standard solution: copy coefficient arrays
		num_coefs = sln->num_coefs;
		num_elems = sln->num_elems;

		mono_coefs = new scalar[num_coefs];
		memcpy(mono_coefs, sln->mono_coefs, sizeof(scalar) * num_coefs);

		for (int l = 0; l < num_components; l++) {
			elem_coefs[l] = new int[num_elems + 1];
			memcpy(elem_coefs[l], sln->elem_coefs[l], sizeof(int) * (num_elems + 1));
		}

	    elem_orders = new Ord3[num_elems + 1];
	    memcpy(elem_orders, sln->elem_orders, sizeof(Ord3) * (num_elems + 1));

		init_dxdydz_buffer();
	}
	else {
		// exact, const
		exact_fn = sln->exact_fn;
		exact_vec_fn = sln->exact_vec_fn;
		cnst[0] = sln->cnst[0];
		cnst[1] = sln->cnst[1];
		cnst[2] = sln->cnst[2];
	}
	seq = sln->seq;
}

void Solution::set_exact(exact_fn_t exactfn) {
	_F_
	free();
	this->mesh = mesh;
	exact_fn = exactfn;
	num_components = 1;
	type = HERMES_EXACT;
	num_dofs = -1;
	seq = g_mfn_seq++;
}

void Solution::set_exact(exact_vec_fn_t exactfn) {
	_F_
	free();
	this->mesh = mesh;
	exact_vec_fn = exactfn;
	num_components = 3;
	type = HERMES_EXACT;
	num_dofs = -1;
	seq = g_mfn_seq++;
}

void Solution::set_const(scalar c) {
	_F_
	free();
	this->mesh = mesh;
	cnst[0] = c;
	cnst[1] = cnst[2] = 0.0;
	num_components = 1;
	type = HERMES_CONST;
	num_dofs = -1;
	seq = g_mfn_seq++;
}

void Solution::set_const(scalar c0, scalar c1, scalar c2) {
	_F_
	free();
	this->mesh = mesh;
	cnst[0] = c0;
	cnst[1] = c1;
	cnst[2] = c2;
	num_components = 3;
	type = HERMES_CONST;
	num_dofs = -1;
	seq = g_mfn_seq++;
}

void Solution::set_zero() {
	_F_
	set_const(0.0);
}

void Solution::set_zero_3() {
	_F_
	set_const(0.0, 0.0, 0.0);
}

// differentiates the mono coefs by x
static void make_dx_coefs(int mode, Ord3 ord, scalar *mono, scalar *result) {
	int i, j, k;

	switch (mode) {
		case HERMES_MODE_TET:
			for (i = 0; i <= ord.order; i++)
				for (j = 0; j <= i; j++) {
					*result++ = 0.0;
					for (k = 0; k < j; k++)
						*result++ = (scalar) (j - k) * mono[k];
					mono += j + 1;
				}
			break;

		case HERMES_MODE_HEX:
			for (i = 0; i <= ord.z; i++) {
				for (j = 0; j <= ord.y; j++) {
					*result++ = 0.0;
					for (k = 0; k < ord.x; k++)
						*result++ = (scalar) (ord.x - k) * mono[k];
					mono += ord.x + 1;
				}
			}
			break;

		case HERMES_MODE_PRISM:
			EXIT(HERMES_ERR_NOT_IMPLEMENTED);

		default:
			EXIT(HERMES_ERR_UNKNOWN_MODE);
	}

}

// differentiates the mono coefs by y
static void make_dy_coefs(int mode, Ord3 ord, scalar *mono, scalar *result) {
	int i, j, k;

	switch (mode) {
		case HERMES_MODE_TET:
			for (i = 0; i <= ord.order; i++) {
				for (j = 0; j <= i; j++) {
					*result++ = 0.0;
					for (k = 0; k < j; k++)
						*result++ = (scalar) (i + 1 - j) * (*mono++);
				}
				mono += i + 1;
			}
			break;

		case HERMES_MODE_HEX:
			for (i = 0; i <= ord.z; i++) {
				for (j = 0; j <= ord.x; j++)
					*result++ = 0.0;
				for (j = 0; j < ord.y; j++)
					for (k = 0; k <= ord.x; k++)
						*result++ = (scalar) (ord.y - j) * (*mono++);
				mono += ord.x + 1;
			}
			break;

		case HERMES_MODE_PRISM:
			EXIT(HERMES_ERR_NOT_IMPLEMENTED);

		default:
			EXIT(HERMES_ERR_UNKNOWN_MODE);
	}
}

// differentiates the mono coefs by z
static void make_dz_coefs(int mode, Ord3 ord, scalar *mono, scalar *result) {
	int i, j, k;

	switch (mode) {
		case HERMES_MODE_TET:
			for (i = 0; i <= ord.order; i++)
				for (j = 0; j <= i; j++) {
					*result++ = 0.0;
					for (k = 0; k < j; k++)
						*result++ = (scalar) (ord.order + 1 - i) * (*mono++);
				}
			break;

		case HERMES_MODE_HEX:
			for (j = 0; j <= ord.y; j++)
				for (k = 0; k <= ord.x; k++)
					*result++ = 0.0;

			for (i = 0; i < ord.z; i++)
				for (j = 0; j <= ord.y; j++)
					for (k = 0; k <= ord.x; k++)
						*result++ = (scalar) (ord.z - i) * (*mono++);
			break;

		case HERMES_MODE_PRISM:
			EXIT(HERMES_ERR_NOT_IMPLEMENTED);

		default:
			EXIT(HERMES_ERR_UNKNOWN_MODE);
	}
}

void Solution::init_dxdydz_buffer() {
	if (dxdydz_buffer != NULL) delete [] dxdydz_buffer;
	dxdydz_buffer = new scalar[num_components * 5 * (int) pow(11., 3)];		// FIXME: 5?
}

void Solution::set_active_element(Element *e) {
	_F_
	MeshFunction::set_active_element(e);

	free_cur_node();

	mode = e->get_mode();
	if (type == HERMES_SLN) {
		order = elem_orders[element->id];
		int np;
		switch (mode) {
			case HERMES_MODE_TET: np = (order.order + 1) * (order.order + 2) * (order.order + 3) / 6; break;
			case HERMES_MODE_HEX: np = (order.x + 1) * (order.y + 1) * (order.z + 1); break;
			default: EXIT(HERMES_ERR_NOT_IMPLEMENTED); break;
		}

		for (int i = 0, m = 0; i < num_components; i++) {
			scalar *mono = mono_coefs + elem_coefs[i][e->id];
			dxdydz_coefs[i][FN] = mono;

			make_dx_coefs(mode, order, mono, dxdydz_coefs[i][DX] = dxdydz_buffer + m);  m += np;
			make_dy_coefs(mode, order, mono, dxdydz_coefs[i][DY] = dxdydz_buffer + m);  m += np;
			make_dz_coefs(mode, order, mono, dxdydz_coefs[i][DZ] = dxdydz_buffer + m);  m += np;
		}
	}
	else if (type == HERMES_EXACT) {
		switch (mode) {
			case HERMES_MODE_TET: order = Ord3(H3D_MAX_QUAD_ORDER_TETRA); break;
			case HERMES_MODE_HEX: order = Ord3(H3D_MAX_QUAD_ORDER, H3D_MAX_QUAD_ORDER, H3D_MAX_QUAD_ORDER); break;
			default: EXIT(HERMES_ERR_NOT_IMPLEMENTED); break;
		}
	}
	else if (type == HERMES_CONST) {
    switch (mode) {
			case HERMES_MODE_TET: order = Ord3(0); break;
			case HERMES_MODE_HEX: order = Ord3(0, 0, 0); break;
			default: EXIT(HERMES_ERR_NOT_IMPLEMENTED); break;
		}
	}
	else
		EXIT("Uninitialized solution.");
}

static struct mono_lu_init {
public:
	// this is a set of LU-decomposed matrices shared by all Solutions
	std::map<unsigned int, double **> mat[3];
	std::map<unsigned int, int *> perm[3];

	mono_lu_init() {
	}

	~mono_lu_init() {
		for (int m = 0; m <= 2; m++) {
			for(std::map<unsigned int, double**>::iterator it = mat[m].begin(); it != mat[m].end(); it++)
        delete [] it->second;
			for(std::map<unsigned int, int*>::iterator it = perm[m].begin(); it != perm[m].end(); it++)
        delete [] it->second;
    }
	}
} mono_lu;

void calc_mono_matrix(const Ord3 &ord, mono_lu_init &mono) {
	int i, j, k, p, q, r, m, row;
	double x, y, z, xn, yn, zn;

	int np;
	double **mat;
	switch (ord.type) {
		case HERMES_MODE_TET:
			np = (ord.order + 1) * (ord.order + 2) * (ord.order + 3) / 6;
			mat = new_matrix<double>(np, np);
			for (p = ord.order, row = 0; p >= 0; p--) {
				z = ord.order ? cos(p * M_PI / ord.order) : 1.0;

				for (q = ord.order; q >= (int)ord.order - p; q--) {
					y = ord.order ? cos(q * M_PI / ord.order) : 1.0;

					for (r = ord.order; r >= (int)ord.order - q + (int)ord.order - p; r--, row++) {
						x = ord.order ? cos(r * M_PI / ord.order) : 1.0;

						// each row of the matrix contains all the monomials x^i * y^j * z^k
						for (i = 0, zn = 1.0, m = np - 1; i <= ord.order; i++, zn *= z)
							for (j = i, yn = 1.0; j <= ord.order; j++, yn *= y)
								for (k = j, xn = 1.0; k <= ord.order; k++, xn *= x, m--)
									mat[row][m] = xn * yn * zn;
					}
				}
			}
			break;

		case HERMES_MODE_HEX:
			np = (ord.x + 1) * (ord.y + 1) * (ord.z + 1);
			mat = new_matrix<double>(np, np);
			for (p = ord.z, row = 0; p >= 0; p--) {
				z = ord.z ? cos(p * M_PI / ord.z) : 1.0;

				for (q = ord.y; q >= 0; q--) {
					y = ord.y ? cos(q * M_PI / ord.y) : 1.0;

					for (r = ord.x; r >= 0; r--, row++) {
						x = ord.x ? cos(r * M_PI / ord.x) : 1.0;

						// each row of the matrix contains all the monomials x^i * y^j * z^k
						for (i = 0, zn = 1.0, m = np - 1; i <= ord.z; i++, zn *= z)
							for (j = 0, yn = 1.0; j <= ord.y; j++, yn *= y)
								for (k = 0, xn = 1.0; k <= ord.x; k++, xn *= x, m--)
									mat[row][m] = xn * yn * zn;
					}
				}
			}
			break;

		default:
			EXIT(HERMES_ERR_NOT_IMPLEMENTED);
			break;
	}


	double d;
	int *perm = new int [np];
	ludcmp(mat, np, perm, &d);

	mono.mat[ord.type][ord.get_idx()] = mat;
	mono.perm[ord.type][ord.get_idx()] = perm;
}

#ifdef WITH_TETRA
	static QuadChebTetra		quad_cheb_tetra;
	#define H3D_QUAD_CHEB_TETRA		&quad_cheb_tetra
#else
	#define H3D_QUAD_CHEB_TETRA		NULL
#endif

#ifdef WITH_HEX
	static QuadChebHex			quad_cheb_hex;
	#define H3D_QUAD_CHEB_HEX		&quad_cheb_hex
#else
	#define H3D_QUAD_CHEB_HEX		NULL
#endif

#define H3D_QUAD_CHEB_PRISM		NULL

static Quad3D *cheb_quad[] = { H3D_QUAD_CHEB_TETRA, H3D_QUAD_CHEB_HEX, H3D_QUAD_CHEB_PRISM };

void Solution::set_coeff_vector(Space *space, scalar *vec, double dir) {
	_F_

	free();

	sptype = space->get_type();
	Shapeset *ss = space->get_shapeset();
	num_components = ss->get_num_components();
	type = HERMES_SLN;
	num_dofs = space->get_dof_count();

	// allocate the coefficient arrays
	num_elems = mesh->get_max_element_id();
	elem_orders = new Ord3 [num_elems + 1];
	for (int l = 0; l < num_components; l++) {
		elem_coefs[l] = new int [num_elems + 1];
		memset(elem_coefs[l], 0, sizeof(int) * num_elems);
	}

	// obtain element orders, allocate mono_coefs
	num_coefs = 0;
	for(std::map<unsigned int, Element*>::iterator it = mesh->elements.begin(); it != mesh->elements.end(); it++)
		if (it->second->used && it->second->active) {
		  Element *e = mesh->elements[it->first];
		  int mode = e->get_mode();
		  Quad3D *quad = cheb_quad[mode];
		  Ord3 ord = space->get_element_order(e->id);
		  // FIXME: this is not very nice, could we handle this in a better (=more general) way
		  if (space->get_type() == HERMES_HCURL_SPACE) ord += Ord3(1, 1, 1);		// FIXME: tetras need Ord3(1)

		  num_coefs += quad->get_num_points(ord);
		  elem_orders[e->id] = ord;
	  }
	num_coefs *= num_components;
	mono_coefs = new scalar[num_coefs];

	ShapeFunction shfn(ss);
	// express the solution on elements as a linear combination of monomials
	scalar *mono = mono_coefs;
	for(std::map<unsigned int, Element*>::iterator it = mesh->elements.begin(); it != mesh->elements.end(); it++)
		if (it->second->used && it->second->active) {
      Element *e = mesh->elements[it->first];
		  int mode = e->get_mode();
		  Quad3D *quad = cheb_quad[mode];

		  Ord3 ord = elem_orders[e->id];
		  int np = quad->get_num_points(ord);
		  QuadPt3D *pt = quad->get_points(ord);

		  AsmList al;
		  space->get_element_assembly_list(e, &al);
		  shfn.set_active_element(e);

		  for (int l = 0; l < num_components; l++) {
			  // obtain solution values for the current element
			  scalar *val = mono;
			  elem_coefs[l][e->id] = (int) (mono - mono_coefs);
			  memset(val, 0, sizeof(scalar) * np);
			  for (int k = 0; k < al.cnt; k++) {
				  shfn.set_active_shape(al.idx[k]);
				  shfn.precalculate(np, pt, FN_VAL);
				  int dof = al.dof[k];
				  scalar coef = al.coef[k] * (dof >= 0 ? vec[dof] : dir);
				  double *shape = shfn.get_fn_values(l);
				  for (int i = 0; i < np; i++)
					  val[i] += shape[i] * coef;
			  }
			  mono += np;

			  // solve for the monomial coefficients
			  if (mono_lu.mat[mode].find(ord.get_idx()) == mono_lu.mat[mode].end())
				  calc_mono_matrix(ord, mono_lu);
			  lubksb(mono_lu.mat[mode][ord.get_idx()], np, mono_lu.perm[mode][ord.get_idx()], val);
		  }
	  }

	init_dxdydz_buffer();
	seq = g_mfn_seq++;
}

// sets all elements of y[] to num

void Solution::vector_to_solutions(scalar* solution_vector, Hermes::vector<Space*> spaces, Hermes::vector<Solution*> solutions, Hermes::vector<double> dir)
{
  assert(spaces.size() == solutions.size());
  for(unsigned int i = 0; i < solutions.size(); i++)
    if(dir == Hermes::vector<double>())
      solutions[i]->set_coeff_vector(spaces[i], solution_vector);
    else
      solutions[i]->set_coeff_vector(spaces[i], solution_vector, dir[i]);
  return;
}

void Solution::vector_to_solution(scalar* solution_vector, Space* space, Solution* solution, double dir)
{
  Hermes::vector<Space*> spaces;
  spaces.push_back(space);
  Hermes::vector<Solution*> solutions;
  solutions.push_back(solution);
  Hermes::vector<double> dirs;
  dirs.push_back(dir);
  Solution::vector_to_solutions(solution_vector, spaces, solutions, dirs);
}


static inline void set_vec_num(int n, scalar *y, scalar num) {
	for (int i = 0; i < n; i++)
		y[i] = num;
}

// y = y .* x + num
static inline void vec_x_vec_p_num(int n, scalar *y, scalar *x, scalar num) {
	for (int i = 0; i < n; i++)
		y[i] = y[i] * x[i] + num;
}

// y = y .* x + z
static inline void vec_x_vec_p_vec(int n, scalar *y, scalar *x, scalar *z) 
{
  for (int i = 0; i < n; i++) y[i] = y[i] * x[i] + z[i];
}

void Solution::precalculate(const int np, const QuadPt3D *pt, int mask) 
{
  _F_
  switch (this->type) {
    case HERMES_SLN: precalculate_fe(np, pt, mask); break;
    case HERMES_EXACT: precalculate_exact(np, pt, mask); break;
    case HERMES_CONST: precalculate_const(np, pt, mask); break;

    default: EXIT("Unknown solution type in Solution::precalculate().");
  }
}

void Solution::precalculate_fe(const int np, const QuadPt3D *pt, int mask) 
{
  _F_

  // if we are required to transform vectors, we must precalculate both their components
  const int GRAD = FN_DX_0 | FN_DY_0 | FN_DZ_0;
  const int CURL = FN_DX | FN_DY | FN_DZ;
  if (transform) {
    if (num_components == 1) {
      if ((mask & FN_DX_0) || (mask & FN_DY_0) || (mask & FN_DZ_0))  mask |= GRAD;
    }
    else {
      if ((mask & FN_VAL_0) || (mask & FN_VAL_1) || (mask & FN_VAL_2)) mask |= FN_VAL;
      if ((mask & FN_DX_0) || (mask & FN_DX_1) || (mask & FN_DX_2) ||
	  (mask & FN_DY_0) || (mask & FN_DY_1) || (mask & FN_DY_2) ||
	  (mask & FN_DZ_0) || (mask & FN_DZ_1) || (mask & FN_DZ_2))
	   mask |= CURL;
    }
  }

  int newmask = mask;
  Node *node = new_node(newmask, np);

  // transform integration points by the current matrix
  scalar * x = new scalar[np];
  scalar * y = new scalar[np];
  scalar * z = new scalar[np];
  scalar * tx = new scalar[np];
  scalar * ty = new scalar[np];
  for (int i = 0; i < np; i++) {
    x[i] = pt[i].x * ctm->m[0] + ctm->t[0];
    y[i] = pt[i].y * ctm->m[1] + ctm->t[1];
    z[i] = pt[i].z * ctm->m[2] + ctm->t[2];
  }

  // obtain the solution values, this is the core of the whole module
  Ord3 ord = elem_orders[element->id];
  for (int l = 0; l < num_components; l++) {
    for (int v = 0; v < 6; v++) {
      if (newmask & idx2mask[v][l]) {
	scalar *result = node->values[l][v];
	// calculate the solution values using Horner's scheme
        scalar *mono = dxdydz_coefs[l][v];
	switch (mode) {
	  case HERMES_MODE_TET:
	  for (int k = 0; k <= ord.order; k++) {					// z
	    for (int i = 0; i <= k; i++) {				// y
	      set_vec_num(np, tx, *mono++);
	      for (int j = 1; j <= i; j++)			// x
	        vec_x_vec_p_num(np, tx, x, *mono++);

	      if (i == 0) memcpy(ty, tx, sizeof(scalar) * np);
	      else vec_x_vec_p_vec(np, ty, y, tx);
	    }

	    if (k == 0) memcpy(result, ty, sizeof(scalar) * np);
	    else vec_x_vec_p_vec(np, result, z, ty);
          }
          break;

	  case HERMES_MODE_HEX:
	    for (int k = 0; k <= ord.z; k++) {					// z
	      for (int i = 0; i <= ord.y; i++) {				// y
	        set_vec_num(np, tx, *mono++);
		for (int j = 1; j <= ord.x; j++)			// x
		  vec_x_vec_p_num(np, tx, x, *mono++);

		if (i == 0) memcpy(ty, tx, sizeof(scalar) * np);
		else vec_x_vec_p_vec(np, ty, y, tx);
	      }

	      if (k == 0) memcpy(result, ty, sizeof(scalar) * np);
	      else vec_x_vec_p_vec(np, result, z, ty);
	    }
	    break;
	}
      }
    }
  }

  // transform gradient or vector solution, if required
  if (transform) {
    bool trans1 = false, trans2 = false;
    scalar *tab1, *tab2, *tab3;
    scalar *tabx[3], *taby[3], *tabz[3];
    if (num_components == 1 && (newmask & GRAD) == GRAD) { // && (oldmask & GRAD) != GRAD) {
      trans1 = true;
      tab1 = node->values[0][DX];
      tab2 = node->values[0][DY];
      tab3 = node->values[0][DZ];
    }
    else if (num_components == 3 && (newmask & FN_VAL) == FN_VAL) { // && (oldmask & FN_VAL) != FN_VAL) {
      trans1 = true;
      tab1 = node->values[0][FN];
      tab2 = node->values[1][FN];
      tab3 = node->values[2][FN];
    }

    if (num_components == 3 && (newmask & CURL) == CURL) { // && (oldmask & CURL) != CURL) {
      trans2 = true;
      for (int i = 0; i < 3; i++) {
	tabx[i] = node->values[i][DX];
	taby[i] = node->values[i][DY];
	tabz[i] = node->values[i][DZ];
      }
    }

    double3x3 *mat = NULL, *m;
    if (trans1 || trans2) mat = refmap->get_inv_ref_map(np, pt);

    // transformation of derivatives in H1 or transformation of values in Hcurl
    int i;
    if (trans1) {
      for (i = 0, m = mat; i < np; i++, m++) {
        scalar vx = tab1[i], vy = tab2[i], vz = tab3[i];
        tab1[i] = (*m)[0][0]*vx + (*m)[0][1]*vy + (*m)[0][2]*vz;
        tab2[i] = (*m)[1][0]*vx + (*m)[1][1]*vy + (*m)[1][2]*vz;
        tab3[i] = (*m)[2][0]*vx + (*m)[2][1]*vy + (*m)[2][2]*vz;
      }
    }

    // transformation of derivatives in Hcurl
    if (trans2) {
      for (i = 0, m = mat; i < np; i++, m++) {
	scalar vhx[3], vhy[3], vhz[3];
	for(int c = 0; c < 3; c++) {
	  vhx[c] = (*m)[0][0] * tabx[c][i] + (*m)[0][1] * taby[c][i] + (*m)[0][2] * tabz[c][i];
	  vhy[c] = (*m)[1][0] * tabx[c][i] + (*m)[1][1] * taby[c][i] + (*m)[1][2] * tabz[c][i];
	  vhz[c] = (*m)[2][0] * tabx[c][i] + (*m)[2][1] * taby[c][i] + (*m)[2][2] * tabz[c][i];
	}

	for(int c = 0; c < 3; c++) {
	  node->values[c][DX][i] = (*m)[c][0] * vhx[0] + (*m)[c][1] * vhx[1] + (*m)[c][2] * vhx[2];
	  node->values[c][DY][i] = (*m)[c][0] * vhy[0] + (*m)[c][1] * vhy[1] + (*m)[c][2] * vhy[2];
	  node->values[c][DZ][i] = (*m)[c][0] * vhz[0] + (*m)[c][1] * vhz[1] + (*m)[c][2] * vhz[2];
	}
      }
    }

    delete [] mat;
  }

  delete [] x;
  delete [] y;
  delete [] z;
  delete [] tx;
  delete [] ty;
  replace_cur_node(node);
}

void Solution::precalculate_exact(const int np, const QuadPt3D *pt, int mask) 
{
  _F_

  mask = FN_DEFAULT;
  Node *node = new_node(mask, np);

  // transform points from ref. domain to physical one
  double *x = refmap->get_phys_x(np, pt);
  double *y = refmap->get_phys_y(np, pt);
  double *z = refmap->get_phys_z(np, pt);

  // evaluate the exact solution
  if (num_components == 1) {
    if (transform) {
       for (int i = 0; i < np; i++) {
	 scalar val, dx = 0.0, dy = 0.0, dz = 0.0;
	 val = exact_fn(x[i], y[i], z[i], dx, dy, dz);
	 node->values[0][FN][i] = val;
	 node->values[0][DX][i] = dx;
	 node->values[0][DY][i] = dy;
	 node->values[0][DZ][i] = dz;
       }
     }
     else {
       // untransform values
       double3x3 *mat = NULL, *m;
       mat = refmap->get_ref_map(np, pt);

       int i;
       for (i = 0, m = mat; i < np; i++, m++) {
         scalar val, dx = 0.0, dy = 0.0, dz = 0.0;
         val = exact_fn(x[i], y[i], z[i], dx, dy, dz);

         node->values[0][FN][i] = val;
         node->values[0][DX][i] = ((*m)[0][0]*dx + (*m)[0][1]*dy + (*m)[0][2]*dz);
         node->values[0][DY][i] = ((*m)[1][0]*dx + (*m)[1][1]*dy + (*m)[1][2]*dz);
         node->values[0][DZ][i] = ((*m)[2][0]*dx + (*m)[2][1]*dy + (*m)[2][2]*dz);
       }

       delete [] mat;
     }
   }
   else if (num_components == 3) {
     // TODO: untransform Hcurl and vector-valued functions
     assert(transform == true);
     for (int i = 0; i < np; i++) {
       scalar3 dx(0.0, 0.0, 0.0);
       scalar3 dy(0.0, 0.0, 0.0);
       scalar3 dz(0.0, 0.0, 0.0);
       scalar3 fn = exact_vec_fn(x[i], y[i], z[i], dx, dy, dz);
       for (int j = 0; j < num_components; j++) {
	 node->values[j][FN][i] = fn[j];
	 node->values[j][DX][i] = dx[j];
	 node->values[j][DY][i] = dy[j];
	 node->values[j][DZ][i] = dz[j];
       }
     }
   }
   else
   EXIT("Invalid number of components.");

   replace_cur_node(node);

   delete [] x;
   delete [] y;
   delete [] z;
}

void Solution::precalculate_const(const int np, const QuadPt3D *pt, int mask) 
{
  _F_
  mask = FN_DEFAULT;
  Node *node = new_node(mask, np);

  assert(num_components == 1 || num_components == 3);
  for (int i = 0; i < np; i++) {
    for (int j = 0; j < num_components; j++) {
      node->values[j][FN][i] = cnst[j];
      node->values[j][DX][i] = 0.0;
      node->values[j][DY][i] = 0.0;
      node->values[j][DZ][i] = 0.0;
    }
  }

  replace_cur_node(node);
}

void Solution::enable_transform(bool enable) 
{
  _F_
  transform = enable;
}

Ord3 Solution::get_order()
{
  _F_
  switch (element->get_mode()) {
    case HERMES_MODE_HEX:
      switch (type) {
	case HERMES_SLN: return elem_orders[element->id];
	case HERMES_EXACT: return Ord3(10, 10, 10);
        case HERMES_CONST: return Ord3(0, 0, 0);
	default: EXIT("Internal error in Solution::get_order() - A.");
      }
      break;

    case HERMES_MODE_TET:
      switch (type) {
	case HERMES_SLN: return elem_orders[element->id];
	case HERMES_EXACT: return Ord3(10);
        case HERMES_CONST: return Ord3(0);
	default: EXIT("Internal error in Solution::get_order() - A.");
      }
      break;

    default: EXIT(HERMES_ERR_NOT_IMPLEMENTED);
    break;
  }

  return Ord3(0);
}
