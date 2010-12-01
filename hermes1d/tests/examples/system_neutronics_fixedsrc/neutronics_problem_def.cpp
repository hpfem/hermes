// Geometry definition
const int N_CONF = 1;						// Number of core configurations
const int N_ASSY = 1;						// Number of assemblies in each core configuration

// Configuration of assemblies in testing cores
int core_conf[N_CONF][N_ASSY] = {
			{ 0 }
		};
int TEST_CONF = 0;		
		
enum Cell_type { pin_group=0 };
const int N_CELL_TYPES = 1;
const int N_CELLS_PER_ASSY = 1;	// Number of cells within each assembly (1 if heterogeneous assemblies are not considered)
const int N_ASSY_TYPES = 1;			// Number of assemblies with different cell composition

// Cell configuration of the assemblies
Cell_type assy_geo[N_ASSY_TYPES][N_CELLS_PER_ASSY] = { 
						{ pin_group }			// 1st assembly
					};

// Width of each cell type [cm]					
double cell_widths[N_CELL_TYPES] = { 80 };	 

// Materials definition
enum Material_type { fuel=0 };													
const int N_MAT = 1;			  		// Number of different materials
const int N_GRP = 2;			  		// Number of energy groups in multigroup approximation

// Material-composition of the core
Material_type assy_mat[N_ASSY_TYPES][N_CELLS_PER_ASSY] =	{ 
								{ fuel }				// 1st assembly
							};																

// Physical properties of each material type (order as in the enumeration above)
static double D[N_MAT][N_GRP] = 	
	{ 
		{ 1.2, 0.4 } 
	};															// diffusion coefficient
static double Sr[N_MAT][N_GRP] = 	
	{ 
		{ 0.03, 0.10 }
	};															// removal cross-section ( \Sigma_{r,g} = \Sigma_{a,g} + 
																	//  sum_{g'\in{all other groups}} \Sigma_{s,g->g'} )
static double nSf[N_MAT][N_GRP] = 
	{{}};														
static double chi[N_GRP] = 
	{};															// no fission
static double Sgg[N_MAT][N_GRP][N_GRP] = 
	{
		{ 
			{ 0.00, 0.00 }, 
			{ 0.02, 0.00 }  
		}
	};															// scattering matrix: 1st index : material
																 	// 										2nd index : target group
																 	//										3rd index :	source group
static double Vsrc[N_MAT][N_GRP] = 
	{
		{ 1.5, 0.0 }
	};															// volumetric source


// Initial FE mesh																		
int cell_N_subdiv[N_CELL_TYPES] = { 10 };	// Equidistant subdivisions for each cell type
int cell_P_init[N_CELL_TYPES] = { 3 }; 		// Initial polynomial degree for each cell type

// Boundary conditions
double current_left_surf[N_GRP] = { 0.0, 0.0 };	// total reflection on the left (zero Neumann)
double flux_right_surf[N_GRP] = { 0.0, 0.0 };		// zero-flux on the right (zero Dirichlet)

