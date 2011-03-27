/* Problem specification (core geometry, material properties, initial FE space) */

// Geometry definition
const int N_CONF = 1;						// Number of core configurations
const int N_ASSY = 7;						// Number of assemblies in each core configuration

// Configuration of assemblies in testing cores
int core_conf[N_CONF][N_ASSY] = {
			{ 0, 1, 2, 1, 2, 2, 1 }
		};
int TEST_CONF = 0;		
		
enum Cell_type { pin_set=0 };
const int N_CELL_TYPES = 1;
const int N_CELLS_PER_ASSY = 1;	// Number of cells within each assembly (1 if heterogeneous assemblies are not considered)
const int N_ASSY_TYPES = 3;			// Number of assemblies with different cell composition

// Cell configuration of the assemblies
Cell_type assy_geo[N_ASSY_TYPES][N_CELLS_PER_ASSY] = { 
						{ pin_set },
						{ pin_set },
						{ pin_set }	
					};

// Width of each cell type [cm]					
double cell_widths[N_CELL_TYPES] = { 100 };	 

// Materials definition
enum Material_type { mat1=0, mat2, mat3 };													
const int N_MAT = 3;			  		// Number of different materials
const int N_GRP = 2;			  		// Number of energy groups in multigroup approximation

// Material-composition of the core
Material_type assy_mat[N_ASSY_TYPES][N_CELLS_PER_ASSY] =	{ 
								{ mat1 },
								{ mat2 },
								{ mat3 }				
							};																
															
// Physical properties of each material type (order as in the enumeration above)
static double D[N_MAT][N_GRP] = 	
	{ 
		{ 1.2, 0.4 },
		{ 1.2, 0.4 },
		{ 1.2, 0.4 } 
	};															// diffusion coefficient
static double Sr[N_MAT][N_GRP] = 	
	{ 
		{ 0.03, 0.10 },
		{ 0.03, 0.20 },
		{ 0.03, 0.25 }
	};															// removal cross-section ( \Sigma_{r,g} = \Sigma_{a,g} + 
																	//  sum_{g'\in{all other groups}} \Sigma_{s,g->g'} )
static double nSf[N_MAT][N_GRP] = 	
	{ 
		{ 0.0050, 0.1 },
		{ 0.0075, 0.1 },
		{ 0.0075, 0.1 }
	};															// fission cross-section	
static double chi[N_GRP] = 
	{ 1.0, 0.0 };										// fission spectrum
static double Sgg[N_MAT][N_GRP][N_GRP] = 
	{
		{ 
			{ 0.00, 0.00 }, 
			{ 0.02, 0.00 }  
		},
		{ 
			{ 0.000, 0.00 }, 
			{ 0.015, 0.00 }  
		},
		{ 
			{ 0.000, 0.00 }, 
			{ 0.015, 0.00 }  
		}
	};															// scattering matrix: 1st index : material
																 	// 										2nd index : target group
																 	//										3rd index :	source group

static double Vsrc[N_MAT][N_GRP] =
	{
		{ 0.0, 0.0 },
		{ 1.5, 0.0 },
		{	1.8, 0.0 }
	};															// external sources

// Initial FE mesh																		
int cell_N_subdiv[N_CELL_TYPES] = { 10 };	// Equidistant subdivisions for each cell type
int cell_P_init[N_CELL_TYPES] = { 3 }; 		// Initial polynomial degree for each cell type

// Boundary conditions
double flux_left_surf[N_GRP] = { 0.0, 0.0 };	// zero-flux on the left
double flux_right_surf[N_GRP] = { 0.0, 0.0 };	// zero-flux on the right

/* Common functions for neutronics problems (requires variable declarations from "Problem specification") */

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) < (b)) ? (b) : (a))

// Structure that translates Space input data specified in "neutronics_problem_def.cpp"
// to format understood by the "Space" object from Hermes
struct SpaceData
{
	double *interfaces;
	int *poly_orders, *material_markers, *subdivisions;
	int N_macroel;
	
	SpaceData(bool report=false)
	{
		// One macroelement is defined for each cell 
		N_macroel = N_ASSY * N_CELLS_PER_ASSY;
		
		interfaces = new double [N_macroel+1]; 	// coordinates of macroelement interfaces [cm] 
		interfaces[0] = 0;
		poly_orders = new int [N_macroel];			// initial poly degrees of macroel
		material_markers = new int [N_macroel]; // material markers of macroel
		subdivisions = new int [N_macroel];			// equidistant subdivision of macroel
	
		if (report)
			printf("\nMacroelements: (a,b), [material, subdivisions, poly_orders]\n");

		int macroel_id = 0;		
		for (int a = 0; a < N_ASSY; a++)
		{
			int assy_type = core_conf[TEST_CONF][a];
		
			for (int c = 0; c < N_CELLS_PER_ASSY; c++)
			{
				Cell_type 		cell_type = assy_geo[assy_type][c];
				Material_type	mtrl_type = assy_mat[assy_type][c];
			
				material_markers[macroel_id] = mtrl_type;																		
				subdivisions		[macroel_id] = cell_N_subdiv[cell_type];  									
				poly_orders			[macroel_id] = cell_P_init[cell_type];											

				interfaces[macroel_id+1] = interfaces[macroel_id] + cell_widths[cell_type];	
				
				if (report)
					printf("(%.2f,%.2f)\t: [%d,%d,%d]\n", interfaces[macroel_id], interfaces[macroel_id+1],
																						material_markers[macroel_id], subdivisions[macroel_id], 
																						poly_orders[macroel_id]);
				macroel_id++;																																
			}
		}
	}
	
	~SpaceData()
	{
		delete [] interfaces;
		delete [] poly_orders;
		delete [] material_markers;
		delete [] subdivisions;
	}
};
	
// Calculate the reaction rate specified by cross-section \Sigma ('xsec'),
// summed (or 'condensed') over all energy groups (or 'spectrum'): 
//	\int_S \sum_{g=1}^{N_GRP} \Sigma^g(x) u^g(x) dV,
// where 'S' is the set defined as the intersection of element "e" with the 
// global integration interval ("x1", "x2"). Set S is assumed to be non-empty.
double calc_elem_reaction_rate_full_spectrum(Element *e, const double xsec[N_MAT][N_GRP], double x1, double x2)
{
	// S = [a,b]
	double a = MAX(x1,e->x1), b = MIN(x2,e->x2);
	
  // numerical quadrature in element 'e'
  int order = e->p; // this is enough since cross-sections are constant in elements
  double phys_x[MAX_QUAD_PTS_NUM];
  double phys_weights[MAX_QUAD_PTS_NUM];

  int pts_num;
  create_phys_element_quadrature(a, b, order, phys_x, phys_weights, &pts_num); 

	// variables to store solution components in each quadrature point
  double val_phys_i[N_GRP], der_phys_i[N_GRP];
  
  int m = e->marker;	// material type
  
  double rr = 0;
	for (int i=0; i < pts_num; i++) {    
    // get all solution components in the i-th quadrature point
		e->get_solution_point(phys_x[i], val_phys_i, der_phys_i);
		
		// calculate the reaction rate
    double val = 0;
		for (int g=0; g < N_GRP; g++)
			val += xsec[m][g] * val_phys_i[g];

		// add contribution to the integral		
		rr += val * phys_weights[i];
		
		// check degeneracy of the integration interval
		if (fabs(a-b) < 1e-14) break;
	}
  return rr;
}

// Calculate sum over all groups of neutron flux integrated over the set S, which is defined
// as the intersection of element "e" with the global integration interval ("x1", "x2").
// Set S is assumed to be non-empty
double calc_elem_flux_full_spectrum(Element *e, double x1, double x2)
{
	// S = [a,b]
	double a = MAX(x1,e->x1), b = MIN(x2,e->x2);

  // numerical quadrature in element 'e'
  int order = e->p; // this is enough since cross-sections are constant in elements
  double phys_x[MAX_QUAD_PTS_NUM];
  double phys_weights[MAX_QUAD_PTS_NUM];
  int pts_num;
  create_phys_element_quadrature(a, b, order, phys_x, phys_weights, &pts_num); 

	// variables to store solution components in each quadrature point
  double val_phys_i[N_GRP];
  double der_phys_i[N_GRP];
  
  int m = e->marker;	// material type
  
  double f = 0;
	for (int i=0; i < pts_num; i++) {    
    // get all solution components in the i-th quadrature point
		e->get_solution_point(phys_x[i], val_phys_i, der_phys_i);
		
		// calculate the flux spectrum
    double val = 0;
		for (int g=0; g < N_GRP; g++)
			val += val_phys_i[g];
		
		// add contribution to the integral	
		f += val * phys_weights[i];
		
		// check degeneracy of the integration interval
		if (fabs(a-b) < 1e-14) break;
	}
	 
  return f;
}

// Calculate integral of neutron flux in group "g" over the set S, which is defined
// as the intersection of element "e" with the global integration interval ("x1", "x2").
// Set S is assumed to be non-empty
double calc_elem_flux(Element *e, int g, double x1, double x2)
{
	// S = [a,b]
	double a = MAX(x1,e->x1), b = MIN(x2,e->x2);

  // numerical quadrature in element 'e'
  int order = e->p; // order with which the solution in element "e" has been approximated
  double phys_x[MAX_QUAD_PTS_NUM];
  double phys_weights[MAX_QUAD_PTS_NUM];
  int pts_num;
  create_phys_element_quadrature(a, b, order, phys_x, phys_weights, &pts_num); 

	// variables to store solution components in each quadrature point
  double val_phys_i[N_GRP];
  double der_phys_i[N_GRP];
  
  double f = 0;
	for (int i=0; i < pts_num; i++) {    
    // get all solution components in the i-th quadrature point
		e->get_solution_point(phys_x[i], val_phys_i, der_phys_i);
		// add contribution to the integral
		f += val_phys_i[g] * phys_weights[i];
		// check degeneracy of the integration interval
		if (fabs(a-b) < 1e-14) break;
	}
	 
  return f;
}

// Calculate sum over all groups of reaction rates specified by cross-section "xsec",
// integrated between points "x1" and "x2"
double calc_total_reaction_rate(Space* Space, const double xsec[N_MAT][N_GRP], double x1, double x2)
{
	double rr = 0;
  Iterator *I = new Iterator(Space);
  Element *e;
  while ((e = I->next_active_element()) != NULL) {
  	// iterate until an element with at least a part lying within the specified integration range is encountered 
  	if (e->x2 < x1 || e->x1 > x2)	continue;
  	// add contribution of the element to the total reaction rate
		rr += calc_elem_reaction_rate_full_spectrum(e, xsec, x1, x2);
  }
  delete I;
  return rr;
}

// Calculate sum over all groups of neutron flux integrated between points "x1" and "x2"
double calc_total_flux(Space* Space, double x1, double x2)
{
	double rr = 0;
  Iterator *I = new Iterator(Space);
  Element *e;
  while ((e = I->next_active_element()) != NULL) {
  	if (e->x2 < x1 || e->x1 > x2) continue;
		rr += calc_elem_flux_full_spectrum(e, x1, x2);
  }
  delete I;
  return rr;
}

// Calculate integral of group-"g" flux between points "x1" and "x2"
double calc_integrated_flux(Space* Space, int g, double x1, double x2)
{
	double rr = 0;
  Iterator *I = new Iterator(Space);
  Element *e;
  while ((e = I->next_active_element()) != NULL) {
  	if (e->x2 < x1 || e->x1 > x2) continue;
		rr += calc_elem_flux(e, g, x1, x2);
  }
  delete I;
  return rr;
}	

// Calculates neutron flux and neutron current in all groups at point "x"
void get_solution_at_point(Space *Space, double x, double flux[N_GRP], double J[N_GRP] )
{
	Iterator *I = new Iterator(Space);
  Element *e;
  int m;
  
  while ((e = I->next_active_element()) != NULL) {
  	if (x < e->x1 || x > e->x2) continue;
  	e->get_solution_point(x, flux, J);
  	m = e->marker;
  }
 
 	for (int g = 0; g < N_GRP; g++)
 		J[g] *= -D[m][g];
}

/* Right-hand side */

// Construct the right-hand side of the "group"-th equation of the multigroup
// system at node "node" in current element. Contribution from both fission and 
// in-scattering is treated implicitly.
double group_sources(int group, double flux[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], int node, int material)
{
	double fission_src = 0.0, scatter_src = 0.0, ext_src = 0.0;
	int last_newton = 0;  // solution indices
	
	for (int g = 0; g<N_GRP; g++) {
		// fission sources
		fission_src += chi[group]*nSf[material][g] * flux[last_newton][g][node];
		scatter_src += Sgg[material][group][g] * flux[last_newton][g][node];
	}
	ext_src += Vsrc[material][group];
	
	return fission_src + scatter_src + ext_src;
}

// bilinear forms for the Jacobi matrix 
// num...number of Gauss points in element
// x[]...Gauss points
// weights[]...Gauss weights for points in x[]
// u...basis function
// v...test function
// u_prev...previous solution


/* JACOBIAN */

/* 0 0 */

double jacobian_mat1_0_0(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                void *user_data)
{
  Material_type m = mat1;	// material type (enumerated in neutronics_problem_def.cpp)		
  int comp = 0;    				// solution component (energy group)  
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += (D[m][comp] * dudx[i] * dvdx[i] + (Sr[m][comp] - nSf[m][comp]) * u[i] * v[i]) * weights[i];
  }
  return val;
}
double jacobian_mat2_0_0(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                void *user_data)
{
  Material_type m = mat2;	// material type (enumerated in neutronics_problem_def.cpp)		
  int comp = 0;    				// solution component (energy group)  
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += (D[m][comp] * dudx[i] * dvdx[i] + (Sr[m][comp] - nSf[m][comp]) * u[i] * v[i]) * weights[i];
  }
  return val;
}
double jacobian_mat3_0_0(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                void *user_data)
{
  Material_type m = mat3;	// material type (enumerated in neutronics_problem_def.cpp)		
  int comp = 0;    				// solution component (energy group)  
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += (D[m][comp] * dudx[i] * dvdx[i] + (Sr[m][comp] - nSf[m][comp]) * u[i] * v[i]) * weights[i];
  }
  return val;
}

/* 0 1 */

double jacobian_mat1_0_1(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                void *user_data)
{
  Material_type m = mat1;	// material type (enumerated in neutronics_problem_def.cpp)		
  int comp = 1;    				// solution component (energy group)  
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += -(nSf[m][comp] + Sgg[m][0][1]) * u[i] * v[i] * weights[i];
  }
  return val;
}
double jacobian_mat2_0_1(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                void *user_data)
{
  Material_type m = mat2;	// material type (enumerated in neutronics_problem_def.cpp)		
  int comp = 1;    				// solution component (energy group)  
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += -(nSf[m][comp] + Sgg[m][0][1]) * u[i] * v[i] * weights[i];
  }
  return val;
}
double jacobian_mat3_0_1(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                void *user_data)
{
  Material_type m = mat3;	// material type (enumerated in neutronics_problem_def.cpp)		
  int comp = 1;    				// solution component (energy group)  
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += -(nSf[m][comp] + Sgg[m][0][1]) * u[i] * v[i] * weights[i];
  }
  return val;
}

/* 1 0 */

double jacobian_mat1_1_0(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                void *user_data)
{
  Material_type m = mat1;	// material type (enumerated in neutronics_problem_def.cpp)		
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += -Sgg[m][1][0] * u[i] * v[i] * weights[i];
  }

  return val;
}
double jacobian_mat2_1_0(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                void *user_data)
{
  Material_type m = mat2;	// material type (enumerated in neutronics_problem_def.cpp)		
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += -Sgg[m][1][0] * u[i] * v[i] * weights[i];
  }

  return val;
}
double jacobian_mat3_1_0(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                void *user_data)
{
  Material_type m = mat3;	// material type (enumerated in neutronics_problem_def.cpp)		
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += -Sgg[m][1][0] * u[i] * v[i] * weights[i];
  }

  return val;
}

/* 1 1 */

double jacobian_mat1_1_1(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                void *user_data)
{
  Material_type m = mat1;	// material type (enumerated in neutronics_problem_def.cpp)		
  int comp = 1;    				// solution component (energy group)  
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += (D[m][comp] * dudx[i] * dvdx[i] + Sr[m][comp] * u[i] * v[i]) * weights[i];
  }
  return val;
}
double jacobian_mat2_1_1(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                void *user_data)
{
  Material_type m = mat2;	// material type (enumerated in neutronics_problem_def.cpp)		
  int comp = 1;    				// solution component (energy group)  
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += (D[m][comp] * dudx[i] * dvdx[i] + Sr[m][comp] * u[i] * v[i]) * weights[i];
  }
  return val;
}
double jacobian_mat3_1_1(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                void *user_data)
{
  Material_type m = mat3;	// material type (enumerated in neutronics_problem_def.cpp)		
  int comp = 1;    				// solution component (energy group)  
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += (D[m][comp] * dudx[i] * dvdx[i] + Sr[m][comp] * u[i] * v[i]) * weights[i];
  }
  return val;
}



/* RESIDUALS */

/* 1st group equations */

double residual_mat1_0(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double *v, double *dvdx, void *user_data)
{
	Material_type m = mat1;		// material type
	int comp = 0;    					// solution component (energy group)
	int last_newton = 0;			// solution indices
 
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += ( D[m][comp] * du_prevdx[last_newton][comp][i] * dvdx[i] + Sr[m][comp] * u_prev[last_newton][comp][i] * v[i]
    				- group_sources(comp, u_prev, i, m) * v[i] ) * weights[i];
  }
  return val;
}
double residual_mat2_0(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double *v, double *dvdx, void *user_data)
{
	Material_type m = mat2;		// material type
	int comp = 0;    					// solution component (energy group)
	int last_newton = 0;			// solution indices
 
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += ( D[m][comp] * du_prevdx[last_newton][comp][i] * dvdx[i] + Sr[m][comp] * u_prev[last_newton][comp][i] * v[i]
    				- group_sources(comp, u_prev, i, m) * v[i] ) * weights[i];
  }
  return val;
}
double residual_mat3_0(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double *v, double *dvdx, void *user_data)
{
	Material_type m = mat3;		// material type
	int comp = 0;    					// solution component (energy group)
	int last_newton = 0;			// solution indices
 
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += ( D[m][comp] * du_prevdx[last_newton][comp][i] * dvdx[i] + Sr[m][comp] * u_prev[last_newton][comp][i] * v[i]
    				- group_sources(comp, u_prev, i, m) * v[i] ) * weights[i];
  }
  return val;
}

/* 2nd group equations */

double residual_mat1_1(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double *v, double *dvdx, void *user_data)
{
  Material_type m = mat1;									// material type
  int comp = 1;    												// solution component (energy group)
  int last_newton = 0;   // solution indices
 
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += ( D[m][comp] * du_prevdx[last_newton][comp][i] * dvdx[i] + Sr[m][comp]*u_prev[last_newton][comp][i] * v[i]
    				- group_sources(comp, u_prev, i, m) * v[i] ) * weights[i];
  }
  return val;
}
double residual_mat2_1(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double *v, double *dvdx, void *user_data)
{
  Material_type m = mat2;									// material type
  int comp = 1;    												// solution component (energy group)
  int last_newton = 0;   // solution indices
 
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += ( D[m][comp] * du_prevdx[last_newton][comp][i] * dvdx[i] + Sr[m][comp]*u_prev[last_newton][comp][i] * v[i]
    				- group_sources(comp, u_prev, i, m) * v[i] ) * weights[i];
  }
  return val;
}
double residual_mat3_1(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double *v, double *dvdx, void *user_data)
{
  Material_type m = mat3;									// material type
  int comp = 1;    												// solution component (energy group)
  int last_newton = 0;   // solution indices
 
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += ( D[m][comp] * du_prevdx[last_newton][comp][i] * dvdx[i] + Sr[m][comp]*u_prev[last_newton][comp][i] * v[i]
    				- group_sources(comp, u_prev, i, m) * v[i] ) * weights[i];
  }
  return val;
}
