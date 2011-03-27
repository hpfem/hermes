/* Problem specification (core geometry, material properties, initial FE space) */

// Geometry definition
const int N_CONF = 1;						// Number of core configurations
const int N_ASSY = 2;						// Number of assemblies in each core configuration

// Configuration of assemblies in testing cores
int core_conf[N_CONF][N_ASSY] = {
			{ 0, 1 }
		};
int TEST_CONF = 0;		
		
enum Cell_type { pin_group=0, reflector };
const int N_CELL_TYPES = 2;
const int N_CELLS_PER_ASSY = 1;	// Number of cells within each assembly (1 if heterogeneous assemblies are not considered)
const int N_ASSY_TYPES = 2;			// Number of assemblies with different cell composition

// Cell configuration of the assemblies
Cell_type assy_geo[N_ASSY_TYPES][N_CELLS_PER_ASSY] = { 
						{ pin_group },			// 1st assembly
						{ reflector } 			// reflector
					};

// Width of each cell type [cm]					
double cell_widths[N_CELL_TYPES] = { 40, 30 };	 

// Materials definition
enum Material_type { fuel=0, water };													
const int N_MAT = 2;			  		// Number of different materials
const int N_GRP = 2;			  		// Number of energy groups in multigroup approximation

// Material-composition of the core
Material_type assy_mat[N_ASSY_TYPES][N_CELLS_PER_ASSY] =	{ 
								{ fuel },						// 1st assembly
							 	{ water }						// reflector
							};																

// Physical properties of each material type (order as in the enumeration above)
static double D[N_MAT][N_GRP] = 	
	{ 
		{ 1.1, 0.89 },		// assembly
		{ 1.1, 0.89 } 		// reflector
	};															// diffusion coefficient
static double Sr[N_MAT][N_GRP] = 	
	{ 
		{ 0.0034, 0.00818 },	
		{ 0.0034, 0.00040 }  
	};															// removal cross-section ( \Sigma_{r,g} = \Sigma_{a,g} + 
																	//  sum_{g'\in{all other groups}} \Sigma_{s,g->g'} )
static double nSf[N_MAT][N_GRP] = 
	{ 
		{ 0.0, 0.0102382},
		{ 0.0, 0.0			}	
	};															// fission-yield cross section (\nu \Sigma_f)
static double chi[N_GRP] = 
	{ 1.0, 0.0 	};									// fission spectrum 
static double Sgg[N_MAT][N_GRP][N_GRP] = 
	{
		{ 
			{ 0.0, 		0.0 }, 
			{ 0.0034, 0.0 }  
		},								// assembly
		{ 
			{ 0.0, 		0.0 }, 
			{ 0.0034, 0.0 }  
		}									// reflector
	};															// scattering matrix: 1st index : material
																 	// 										2nd index : target group
																 	//										3rd index :	source group

static double Vsrc[N_GRP][N_MAT] = 
	{{}};														// no external sources

// Initial FE mesh																		
int cell_N_subdiv[N_CELL_TYPES] = { 4, 3 };	// Equidistant subdivisions for each cell type
int cell_P_init[N_CELL_TYPES] = { 3, 3 }; 		// Initial polynomial degree for each cell type

// Boundary conditions
double current_left_surf[N_GRP] = { 0.0, 0.0 };	// total reflection on the left (zero Neumann)
double flux_right_surf[N_GRP] = { 0.0, 0.0 };		// zero-flux on the right (zero Dirichlet)


/* Common functions for neutronics problems (requires variable declarations from "Problem specification") */

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) < (b)) ? (b) : (a))

// Structure that translates space input data specified in "Problem specification"
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
double calc_total_reaction_rate(Space* space, const double xsec[N_MAT][N_GRP], double x1, double x2)
{
	double rr = 0;
  Iterator *I = new Iterator(space);
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
double calc_total_flux(Space* space, double x1, double x2)
{
	double rr = 0;
  Iterator *I = new Iterator(space);
  Element *e;
  while ((e = I->next_active_element()) != NULL) {
  	if (e->x2 < x1 || e->x1 > x2) continue;
		rr += calc_elem_flux_full_spectrum(e, x1, x2);
  }
  delete I;
  return rr;
}

// Calculate integral of group-"g" flux between points "x1" and "x2"
double calc_integrated_flux(Space* space, int g, double x1, double x2)
{
	double rr = 0;
  Iterator *I = new Iterator(space);
  Element *e;
  while ((e = I->next_active_element()) != NULL) {
  	if (e->x2 < x1 || e->x1 > x2) continue;
		rr += calc_elem_flux(e, g, x1, x2);
  }
  delete I;
  return rr;
}	

// Calculates neutron flux and neutron current in all groups at point "x"
void get_solution_at_point(Space *space, double x, double flux[N_GRP], double J[N_GRP] )
{
	Iterator *I = new Iterator(space);
  Element *e;
  int m;
  
  while ((e = I->next_active_element()) != NULL) {
  	if (x < e->x1 || x > e->x2) continue;
  	e->get_solution_point(x, flux, J);
  	m = e->marker;
  }
 
 	for (int g = 0; g < N_GRP; g++)
 		J[g] *= -D[m][g];
        delete I;
}

/* Right-hand side */

// Construct the right-hand side of the "group"-th equation of the multigroup
// system at node "node" in current element. Contribution from fission sources 
// is treated explictly (using the solution from previous source iteration), 
// while scattering is treated implicitly. Previous flux iterate is assumed to 
// be normalized to one fission neutron
double group_sources(int group, double flux[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], int node, int material)
{
	double fission_src = 0.0, scatter_src = 0.0;
	int last_newton = 0, last_global = 1;   // solution indices
	
	for (int g = 0; g<N_GRP; g++) {
		fission_src += chi[group]*nSf[material][g] * flux[last_global][g][node];
		scatter_src += Sgg[material][group][g] * flux[last_newton][g][node];
	}
		
	return fission_src + scatter_src;
}

// bilinear form for the Jacobi matrix 
// num...number of Gauss points in element
// x[]...Gauss points
// weights[]...Gauss weights for points in x[]
// u...basis function
// v...test function
// u_prev...previous solution

/* JACOBIAN */

/* 0 0 */

double jacobian_fuel_0_0(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                void *user_data)
{
  Material_type m = fuel;	// material type (enumerated in neutronics_problem_def.cpp)		
  int comp = 0;    				// solution component (energy group)  
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += (D[m][comp] * dudx[i] * dvdx[i] + Sr[m][comp] * u[i] * v[i]) * weights[i];
    
  }
  return val;
}
double jacobian_water_0_0(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                void *user_data)
{
  Material_type m = water;	// material type (enumerated in neutronics_problem_def.cpp)		
  int comp = 0;    					// solution component (energy group)  
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += (D[m][comp] * dudx[i] * dvdx[i] + Sr[m][comp] * u[i] * v[i]) * weights[i];
  }
  return val;
}

/* 0 1 */

double jacobian_fuel_0_1(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                void *user_data)
{
  Material_type m = fuel;	// material type (enumerated in neutronics_problem_def.cpp)		
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += -Sgg[m][0][1] * u[i] * v[i] * weights[i];
  }
  return val;
}
double jacobian_water_0_1(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                void *user_data)
{
  Material_type m = water;	// material type (enumerated in neutronics_problem_def.cpp)		
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += -Sgg[m][0][1] * u[i] * v[i] * weights[i];
  }
  return val;
}

/* 1 0 */

double jacobian_fuel_1_0(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                void *user_data)
{
  Material_type m = fuel;	// material type (enumerated in neutronics_problem_def.cpp)		
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += -Sgg[m][1][0] * u[i] * v[i] * weights[i];
  }
  return val;
}
double jacobian_water_1_0(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                void *user_data)
{
  Material_type m = water;	// material type (enumerated in neutronics_problem_def.cpp)		
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += -Sgg[m][1][0] * u[i] * v[i] * weights[i];
  }
  return val;
}

/* 1 1 */

double jacobian_fuel_1_1(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                void *user_data)
{
  Material_type m = fuel;	// material type (enumerated in neutronics_problem_def.cpp)		
  int comp = 1;    				// solution component (energy group)  
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += (D[m][comp] * dudx[i] * dvdx[i] + Sr[m][comp] * u[i] * v[i]) * weights[i];
  }
  return val;
}
double jacobian_water_1_1(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                void *user_data)
{
  Material_type m = water;	// material type (enumerated in neutronics_problem_def.cpp)		
  int comp = 1;    					// solution component (energy group)  
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += (D[m][comp] * dudx[i] * dvdx[i] + Sr[m][comp] * u[i] * v[i]) * weights[i];
  }
  return val;
}



/* RESIDUALS */


/* 1st group equations */

double residual_fuel_0(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double *v, double *dvdx, void *user_data)
{
  Material_type m = fuel;									// material type
  int comp = 0;    												// solution component (energy group)
  int last_newton = 0, last_global = 1;   // solution indices
 
  double val = 0;

  for(int i = 0; i<num; i++) {
    val += ( D[m][comp] * du_prevdx[last_newton][comp][i] * dvdx[i] + Sr[m][comp] * u_prev[last_newton][comp][i] * v[i]
    				- group_sources(comp, u_prev, i, m) * v[i] ) * weights[i];
  }
  return val;
}
double residual_water_0(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double *v, double *dvdx, void *user_data)
{
  Material_type m = water;								// material type
  int comp = 0;    												// solution component (energy group)
  int last_newton = 0, last_global = 1;   // solution indices
 
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += ( D[m][comp] * du_prevdx[last_newton][comp][i] * dvdx[i] + Sr[m][comp] * u_prev[last_newton][comp][i] * v[i]
    				- group_sources(comp, u_prev, i, m) * v[i] ) * weights[i];
 	}
  return val;
}


/* 2nd group equations */

double residual_fuel_1(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double *v, double *dvdx, void *user_data)
{
  Material_type m = fuel;									// material type
  int comp = 1;    												// solution component (energy group)
  int last_newton = 0, last_global = 1;   // solution indices
 
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += ( D[m][comp] * du_prevdx[last_newton][comp][i] * dvdx[i] + Sr[m][comp]*u_prev[last_newton][comp][i] * v[i]
    				- group_sources(comp, u_prev, i, m) * v[i] ) * weights[i];
  }
  return val;
}
double residual_water_1(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double *v, double *dvdx, void *user_data)
{
  Material_type m = water;								// material type
  int comp = 1;    												// solution component (energy group)
  int last_newton = 0, last_global = 1;   // solution indices
 
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += ( D[m][comp] * du_prevdx[last_newton][comp][i] * dvdx[i] + Sr[m][comp]*u_prev[last_newton][comp][i] * v[i]
    				- group_sources(comp, u_prev, i, m) * v[i] ) * weights[i];
  }
  return val;
}

/* NEUMANN B.C. */

double residual_surf_left_0(double x, double u_prev[MAX_SLN_NUM][MAX_EQN_NUM], 
        double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM], double v,
        double dvdx, void *user_data)
{
	Material_type m = fuel;	// material of the leftmost cell
  int comp = 0;    				// solution component (energy group)
  return current_left_surf[comp] * v; 
}
double residual_surf_left_1(double x, double u_prev[MAX_SLN_NUM][MAX_EQN_NUM], 
        double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM], double v,
        double dvdx, void *user_data)
{
	Material_type m = fuel; // material of the leftmost cell
  int comp = 1;    				// solution component (energy group)
  return current_left_surf[comp] * v; 
}

/* EXTRAPOLATED ZERO FLUX B.C. */

double jacobian_surf_right_0(double x, double u, double dudx,
        double v, double dvdx, double u_prev[MAX_SLN_NUM][MAX_EQN_NUM], 
        double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM], void *user_data)
{
	Material_type m = water;
	int comp = 0;
  return 0.5 * u * v;
}
double jacobian_surf_right_1(double x, double u, double dudx,
        double v, double dvdx, double u_prev[MAX_SLN_NUM][MAX_EQN_NUM], 
        double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM], void *user_data)
{
	Material_type m = water;
	int comp = 1;
  return 0.5 * u * v;
}

double residual_surf_right_0(double x, double u_prev[MAX_SLN_NUM][MAX_EQN_NUM], 
        double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM], double v,
        double dvdx, void *user_data)
{
	Material_type m = water;
  int comp = 0;    // solution component
  int last_newton = 0, last_global = 1;   // solution indices
  return 0.5 * u_prev[last_newton][comp] * v; 
}

double residual_surf_right_1(double x, double u_prev[MAX_SLN_NUM][MAX_EQN_NUM], 
        double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM], double v,
        double dvdx, void *user_data)
{
	Material_type m = water;
  int comp = 1;    // solution component
  int last_newton = 0, last_global = 1;   // solution indices
  return 0.5 * u_prev[last_newton][comp] * v; 
}
