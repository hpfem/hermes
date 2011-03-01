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


