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
  return 1./D[m][comp] * current_left_surf[comp] * v; 
}
double residual_surf_left_1(double x, double u_prev[MAX_SLN_NUM][MAX_EQN_NUM], 
        double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM], double v,
        double dvdx, void *user_data)
{
	Material_type m = fuel; // material of the leftmost cell
  int comp = 1;    				// solution component (energy group)
  return 1./D[m][comp] * current_left_surf[comp] * v; 
}

/* EXTRAPOLATED ZERO FLUX B.C. */

double jacobian_surf_right_0(double x, double u, double dudx,
        double v, double dvdx, double u_prev[MAX_SLN_NUM][MAX_EQN_NUM], 
        double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM], void *user_data)
{
	Material_type m = water;
	int comp = 0;
  return 0.5 / D[m][comp] * u * v;
}
double jacobian_surf_right_1(double x, double u, double dudx,
        double v, double dvdx, double u_prev[MAX_SLN_NUM][MAX_EQN_NUM], 
        double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM], void *user_data)
{
	Material_type m = water;
	int comp = 1;
  return 0.5 / D[m][comp] * u * v;
}

double residual_surf_right_0(double x, double u_prev[MAX_SLN_NUM][MAX_EQN_NUM], 
        double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM], double v,
        double dvdx, void *user_data)
{
	Material_type m = water;
  int comp = 0;    // solution component
  int last_newton = 0, last_global = 1;   // solution indices
  return 0.5 / D[m][comp] * u_prev[last_newton][comp] * v; 
}

double residual_surf_right_1(double x, double u_prev[MAX_SLN_NUM][MAX_EQN_NUM], 
        double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM], double v,
        double dvdx, void *user_data)
{
	Material_type m = water;
  int comp = 1;    // solution component
  int last_newton = 0, last_global = 1;   // solution indices
  return 0.5 / D[m][comp] * u_prev[last_newton][comp] * v; 
}
