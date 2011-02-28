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
