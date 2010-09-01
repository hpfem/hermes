#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) < (b)) ? (b) : (a))

// Structure that translates mesh input data specified in "neutronics_problem_def.cpp"
// to format understood by the "Mesh" object from Hermes
struct MeshData
{
	double *interfaces;
	int *poly_orders, *material_markers, *subdivisions;
	int N_macroel;
	
	MeshData(bool report=false)
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
	
	~MeshData()
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
double calc_total_reaction_rate(Mesh* mesh, const double xsec[N_MAT][N_GRP], double x1, double x2)
{
	double rr = 0;
  Iterator *I = new Iterator(mesh);
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
double calc_total_flux(Mesh* mesh, double x1, double x2)
{
	double rr = 0;
  Iterator *I = new Iterator(mesh);
  Element *e;
  while ((e = I->next_active_element()) != NULL) {
  	if (e->x2 < x1 || e->x1 > x2) continue;
		rr += calc_elem_flux_full_spectrum(e, x1, x2);
  }
  delete I;
  return rr;
}

// Calculate integral of group-"g" flux between points "x1" and "x2"
double calc_integrated_flux(Mesh* mesh, int g, double x1, double x2)
{
	double rr = 0;
  Iterator *I = new Iterator(mesh);
  Element *e;
  while ((e = I->next_active_element()) != NULL) {
  	if (e->x2 < x1 || e->x1 > x2) continue;
		rr += calc_elem_flux(e, g, x1, x2);
  }
  delete I;
  return rr;
}	

// Calculates neutron flux and neutron current in all groups at point "x"
void get_solution_at_point(Mesh *mesh, double x, double flux[N_GRP], double J[N_GRP] )
{
	Iterator *I = new Iterator(mesh);
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
