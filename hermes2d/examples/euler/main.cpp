// Development version of Euler equations example. For testing the equations are solved on a simple rectangular mesh.
// Boundary conditions are specified for Inlet / Outlet / Solid wall parts of the boundary in a common way.
// Explicit Euler's method used is used for the time discretization.
// Initial conditions :  w0 = 0.0, w1 = 50.0, w2 = 0.0, w3 = 1E5.
// Prescribed boundary values correspond with the initial conditions.


const int P_INIT = 0;

double TAU = 1E-6;

double R = 300; // Gas constant

double c_v = 700; // The specific heat capacity at constant volume

double t = 0;

#include "forms.cpp"
#include "filters.cpp"
MatrixSolverType matrix_solver = SOLVER_UMFPACK;

BCType bc_types(int marker)
{
	return BC_NATURAL;
}

int main(int argc, char* argv[])
{
	Mesh mesh;
	H2DReader mloader;
	mloader.load("GAMM-channel.mesh", &mesh);

	mesh.refine_all_elements();
	mesh.refine_all_elements();
	mesh.refine_all_elements();
	mesh.refine_all_elements();

	L2Space space_rho(&mesh,P_INIT);
	L2Space space_rho_v_x(&mesh,P_INIT);
	L2Space space_rho_v_y(&mesh,P_INIT);
	L2Space space_e(&mesh,P_INIT);

	space_rho.set_bc_types(bc_types);
	space_rho_v_x.set_bc_types(bc_types);
	space_rho_v_y.set_bc_types(bc_types);
	space_e.set_bc_types(bc_types);

	/*
	BaseView bview;
	bview.show(&space_rho);
	//MeshView mview("Mesh");
	//mview.show(&mesh);
	View::wait();
	*/


	Solution sln_rho, sln_rho_v_x, sln_rho_v_y, sln_e, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e;
	sln_rho.set_exact(&mesh, ic_density);
	sln_rho_v_x.set_exact(&mesh, ic_density_vel_x);
	sln_rho_v_y.set_exact(&mesh, ic_density_vel_y);
	sln_e.set_exact(&mesh, ic_energy);
	prev_rho.set_exact(&mesh, ic_density);
	prev_rho_v_x.set_exact(&mesh, ic_density_vel_x);
	prev_rho_v_y.set_exact(&mesh, ic_density_vel_y);
	prev_e.set_exact(&mesh, ic_energy);

	WeakForm wf(4);

	// Bilinear forms coming from time discretization by explicit Euler's method.
	wf.add_matrix_form(0,0,callback(bilinear_form_0_0_time));
	wf.add_matrix_form(1,1,callback(bilinear_form_1_1_time));
	wf.add_matrix_form(2,2,callback(bilinear_form_2_2_time));
	wf.add_matrix_form(3,3,callback(bilinear_form_3_3_time));

	//Volumetric linear forms.
	// Linear forms coming from the linearization by taking the Eulerian fluxes' Jacobian matrices from the previous time step.
	// First flux.
	wf.add_matrix_form(0,1,callback(linear_form_0_1),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_rho_v_x));
	wf.add_matrix_form(1,0,callback(linear_form_1_0_first_flux),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
	wf.add_matrix_form(1,1,callback(linear_form_1_1_first_flux),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
	wf.add_matrix_form(1,2,callback(linear_form_1_2_first_flux),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
	wf.add_matrix_form(1,3,callback(linear_form_1_3_first_flux),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
	wf.add_matrix_form(2,0,callback(linear_form_2_0_first_flux),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
	wf.add_matrix_form(2,1,callback(linear_form_2_1_first_flux),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
	wf.add_matrix_form(2,2,callback(linear_form_2_2_first_flux),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
	wf.add_matrix_form(2,3,callback(linear_form_2_3_first_flux),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
	wf.add_matrix_form(3,0,callback(linear_form_3_0_first_flux),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
	wf.add_matrix_form(3,1,callback(linear_form_3_1_first_flux),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
	wf.add_matrix_form(3,2,callback(linear_form_3_2_first_flux),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
	wf.add_matrix_form(3,3,callback(linear_form_3_3_first_flux),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
	// Second flux.
	wf.add_matrix_form(0,2,callback(linear_form_0_2),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_rho_v_y));
	wf.add_matrix_form(1,0,callback(linear_form_1_0_second_flux),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
	wf.add_matrix_form(1,1,callback(linear_form_1_1_second_flux),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
	wf.add_matrix_form(1,2,callback(linear_form_1_2_second_flux),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
	wf.add_matrix_form(1,3,callback(linear_form_1_3_second_flux),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
	wf.add_matrix_form(2,0,callback(linear_form_2_0_second_flux),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
	wf.add_matrix_form(2,1,callback(linear_form_2_1_second_flux),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
	wf.add_matrix_form(2,2,callback(linear_form_2_2_second_flux),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
	wf.add_matrix_form(2,3,callback(linear_form_2_3_second_flux),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
	wf.add_matrix_form(3,0,callback(linear_form_3_0_second_flux),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
	wf.add_matrix_form(3,1,callback(linear_form_3_1_second_flux),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
	wf.add_matrix_form(3,2,callback(linear_form_3_2_second_flux),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
	wf.add_matrix_form(3,3,callback(linear_form_3_3_second_flux),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));

	// Volumetric linear forms coming from the time discretization.
	wf.add_vector_form(0,linear_form, linear_form_order, H2D_ANY, &prev_rho);
	wf.add_vector_form(1,linear_form, linear_form_order, H2D_ANY, &prev_rho_v_x);
	wf.add_vector_form(2,linear_form, linear_form_order, H2D_ANY, &prev_rho_v_y);
	wf.add_vector_form(3,linear_form, linear_form_order, H2D_ANY, &prev_e);

	// Surface linear forms - inner edges coming from the DG formulation.
	wf.add_vector_form_surf(0,linear_form_interface_0, linear_form_order, H2D_DG_INNER_EDGE, Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
	wf.add_vector_form_surf(1,linear_form_interface_1, linear_form_order, H2D_DG_INNER_EDGE, Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
	wf.add_vector_form_surf(2,linear_form_interface_2, linear_form_order, H2D_DG_INNER_EDGE, Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
	wf.add_vector_form_surf(3,linear_form_interface_3, linear_form_order, H2D_DG_INNER_EDGE, Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));

	// Surface linear forms - inlet / outlet edges.
	wf.add_vector_form_surf(0,bdy_flux_inlet_outlet_comp_0, linear_form_order, 2, Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
	wf.add_vector_form_surf(1,bdy_flux_inlet_outlet_comp_1, linear_form_order, 2, Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
	wf.add_vector_form_surf(2,bdy_flux_inlet_outlet_comp_2, linear_form_order, 2, Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
	wf.add_vector_form_surf(3,bdy_flux_inlet_outlet_comp_3, linear_form_order, 2, Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));

	// Surface linear forms - Solid wall edges.
	wf.add_vector_form_surf(0,bdy_flux_solid_wall_comp_0, linear_form_order, 1, Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
	wf.add_vector_form_surf(1,bdy_flux_solid_wall_comp_1, linear_form_order, 1, Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
	wf.add_vector_form_surf(2,bdy_flux_solid_wall_comp_2, linear_form_order, 1, Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
	wf.add_vector_form_surf(3,bdy_flux_solid_wall_comp_3, linear_form_order, 1, Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));

	LinearProblem lp(&wf, Tuple<Space*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));

	// Filters for visualization of pressure and the two components of velocity.
	SimpleFilter pressure(calc_pressure_func, Tuple<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e));
	SimpleFilter u(calc_u_func, Tuple<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e));
	SimpleFilter w(calc_w_func, Tuple<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e));

	VectorView vview("Velocity", 0, 0, 600, 300);
	ScalarView sview("Pressure", 700, 0, 600, 300);

	// Iteration number.
	int iteration = 0;
	// Select matrix solver.
	Matrix* mat; Vector* rhs; CommonSolver* solver;
	init_matrix_solver(matrix_solver, get_num_dofs(Tuple<Space*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e)), mat, rhs, solver);

	for(t = 0.0; t < TAU * 1000; t += TAU)
	{		
		iteration++;
		
		// Assemble stiffness matrix and rhs.
		lp.assemble(mat, rhs, iteration == 1 ? false : true);

		// Debugging.
		/*
		std::ofstream out("matrix");
		for(int i = 0; i < mat->get_size(); i++)
			for(int j = 0; j < mat->get_size(); j++)
				if(mat->get(i,j) != 0)
					out << '(' << i << ',' << j << ')' << ':' << mat->get(i,j) << std::endl;
		out.close();

		out.open("rhs");
			for(int j = 0; j < mat->get_size(); j++)
				if(rhs->get(j) != 0)
					out << '(' << j << ')' << ':' << rhs->get(j) << std::endl;
		out.close();
		*/
		
		// Solve the matrix problem.
		if (!solver->solve(mat, rhs)) error ("Matrix solver failed.\n");

		// Debugging continued.
		/*
		out.open("sol");
		for(int j = 0; j < mat->get_size(); j++)
				if(rhs->get(j) != 0)
					out << '(' << j << ')' << ':' << rhs->get(j) << std::endl;
		out.close();
		*/

		sln_rho.set_coeff_vector(&space_rho, rhs);
		sln_rho_v_x.set_coeff_vector(&space_rho_v_x, rhs);
		sln_rho_v_y.set_coeff_vector(&space_rho_v_y, rhs);
		sln_e.set_coeff_vector(&space_e, rhs);

		prev_rho.copy(&sln_rho);
		prev_rho_v_x.copy(&sln_rho_v_x);
		prev_rho_v_y.copy(&sln_rho_v_y);
		prev_e.copy(&sln_e);

		// Visualization.
		pressure.reinit();
		u.reinit();
		w.reinit();
		sview.show(&pressure);
		vview.show(&u,&w);
	}
	vview.close();
	return 0;
}

