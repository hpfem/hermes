#include "forms.cpp"
#include "filters.cpp"

const int p_init = 0;

double tau = 1E-5;

double R = 300; // Gas constant

double c_v = 700; // The specific heat capacity at constant volume

double t = 0;

MatrixSolverType matrix_solver = SOLVER_UMFPACK;

BCType bc_types(int marker)
{
	return BC_NATURAL;
}

int main(int argc, char* argv[])
{
	Mesh mesh;
  H2DReader mloader;
  mloader.load("quad.mesh", &mesh);

	mesh.refine_all_elements();
	mesh.refine_all_elements();
	mesh.refine_all_elements();
	mesh.refine_all_elements();

	L2Space space_d(&mesh,p_init);
  L2Space space_dv1(&mesh,p_init);
  L2Space space_dv2(&mesh,p_init);
  L2Space space_e(&mesh,p_init);

	space_d.set_bc_types(bc_types);
	space_dv1.set_bc_types(bc_types);
	space_dv2.set_bc_types(bc_types);
	space_e.set_bc_types(bc_types);

	/*
  BaseView bview;
  bview.show(&space_d);
  //MeshView mview("Mesh");
  //mview.show(&mesh);
  View::wait();
	*/
	
	Solution sln_d, sln_dv1, sln_dv2, sln_e, prev_d, prev_dv1, prev_dv2, prev_e;
	sln_d.set_exact(&mesh, ic_density);
	sln_dv1.set_exact(&mesh, ic_density_vel1);
	sln_dv2.set_exact(&mesh, ic_density_vel2);
	sln_e.set_exact(&mesh, ic_energy);
	prev_d.set_exact(&mesh, ic_density);
	prev_dv1.set_exact(&mesh, ic_density_vel1);
	prev_dv2.set_exact(&mesh, ic_density_vel2);
	prev_e.set_exact(&mesh, ic_energy);

  WeakForm wf(4);
  
	//Volumetric bilinear forms.
	wf.add_matrix_form(0,0,callback(bilinear_form_0_0));
  wf.add_matrix_form(0,1,callback(bilinear_form_0_1));
  wf.add_matrix_form(0,2,callback(bilinear_form_0_2));
  wf.add_matrix_form(1,0,callback(bilinear_form_1_0),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_d, &prev_dv1, &prev_dv2));
  wf.add_matrix_form(1,1,callback(bilinear_form_1_1),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_d, &prev_dv1, &prev_dv2));
  wf.add_matrix_form(1,2,callback(bilinear_form_1_2),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_d, &prev_dv1, &prev_dv2));
  wf.add_matrix_form(1,3,callback(bilinear_form_1_3),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_d, &prev_dv1, &prev_dv2));
  wf.add_matrix_form(2,0,callback(bilinear_form_2_0),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_d, &prev_dv1, &prev_dv2));
  wf.add_matrix_form(2,1,callback(bilinear_form_2_1),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_d, &prev_dv1, &prev_dv2));
  wf.add_matrix_form(2,2,callback(bilinear_form_2_2),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_d, &prev_dv1, &prev_dv2));
  wf.add_matrix_form(2,3,callback(bilinear_form_2_3),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_d, &prev_dv1, &prev_dv2));
  wf.add_matrix_form(3,0,callback(bilinear_form_3_0),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_d, &prev_dv1, &prev_dv2, &prev_e));
  wf.add_matrix_form(3,1,callback(bilinear_form_3_1),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_d, &prev_dv1, &prev_dv2, &prev_e));
  wf.add_matrix_form(3,2,callback(bilinear_form_3_2),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_d, &prev_dv1, &prev_dv2, &prev_e));
  wf.add_matrix_form(3,3,callback(bilinear_form_3_3),H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&prev_d, &prev_dv1, &prev_dv2));

	// Volumetric linear forms.
	wf.add_vector_form(0,linear_form, linear_form_order, H2D_ANY, &prev_d);
	wf.add_vector_form(1,linear_form, linear_form_order, H2D_ANY, &prev_dv1);
	wf.add_vector_form(2,linear_form, linear_form_order, H2D_ANY, &prev_dv2);
	wf.add_vector_form(3,linear_form, linear_form_order, H2D_ANY, &prev_e);

	// Surface linear forms - inner edges.
	wf.add_vector_form_surf(0,linear_form_interface_0, linear_form_order, H2D_DG_INNER_EDGE, Tuple<MeshFunction*>(&prev_d, &prev_dv1, &prev_dv2, &prev_e));
	wf.add_vector_form_surf(1,linear_form_interface_1, linear_form_order, H2D_DG_INNER_EDGE, Tuple<MeshFunction*>(&prev_d, &prev_dv1, &prev_dv2, &prev_e));
	wf.add_vector_form_surf(2,linear_form_interface_2, linear_form_order, H2D_DG_INNER_EDGE, Tuple<MeshFunction*>(&prev_d, &prev_dv1, &prev_dv2, &prev_e));
	wf.add_vector_form_surf(3,linear_form_interface_3, linear_form_order, H2D_DG_INNER_EDGE, Tuple<MeshFunction*>(&prev_d, &prev_dv1, &prev_dv2, &prev_e));

	// Surface linear forms - IO edges.
	wf.add_vector_form_surf(0,linear_form_IO_0, linear_form_order, 2, Tuple<MeshFunction*>(&prev_d, &prev_dv1, &prev_dv2, &prev_e));
	wf.add_vector_form_surf(1,linear_form_IO_1, linear_form_order, 2, Tuple<MeshFunction*>(&prev_d, &prev_dv1, &prev_dv2, &prev_e));
	wf.add_vector_form_surf(2,linear_form_IO_2, linear_form_order, 2, Tuple<MeshFunction*>(&prev_d, &prev_dv1, &prev_dv2, &prev_e));
	wf.add_vector_form_surf(3,linear_form_IO_3, linear_form_order, 2, Tuple<MeshFunction*>(&prev_d, &prev_dv1, &prev_dv2, &prev_e));

	// Surface linear forms - IMP edges.
	wf.add_vector_form_surf(0,linear_form_IMP_0, linear_form_order, 1, Tuple<MeshFunction*>(&prev_d, &prev_dv1, &prev_dv2, &prev_e));
	wf.add_vector_form_surf(1,linear_form_IMP_1, linear_form_order, 1, Tuple<MeshFunction*>(&prev_d, &prev_dv1, &prev_dv2, &prev_e));
	wf.add_vector_form_surf(2,linear_form_IMP_2, linear_form_order, 1, Tuple<MeshFunction*>(&prev_d, &prev_dv1, &prev_dv2, &prev_e));
	wf.add_vector_form_surf(3,linear_form_IMP_3, linear_form_order, 1, Tuple<MeshFunction*>(&prev_d, &prev_dv1, &prev_dv2, &prev_e));


	LinearProblem lp(&wf, Tuple<Space*>(&space_d, &space_dv1, &space_dv2, &space_e));

	std::ofstream out;
	
	//SimpleFilter pressure(calc_pressure_func, &sln_d, &sln_dv1, &sln_dv2, &sln_e);
	SimpleFilter u(calc_u_func, Tuple<MeshFunction*>(&sln_d, &sln_dv1, &sln_dv2, &sln_e));
  SimpleFilter w(calc_w_func, Tuple<MeshFunction*>(&sln_d, &sln_dv1, &sln_dv2, &sln_e));

	ScalarView view;
	ScalarView view2;
	VectorView vview;

	// Iteration, for screenshot saving.
	int i = 0;
  // assemble and solve the finite element problem
	for(t = 0.0; t < tau * 1000; t += tau)
	{		
		i++;

    // Select matrix solver.
    Matrix* mat; Vector* rhs; CommonSolver* solver;
    init_matrix_solver(matrix_solver, get_num_dofs(Tuple<Space*>(&space_d, &space_dv1, &space_dv2, &space_e)), mat, rhs, solver);

    // Assemble stiffness matrix and rhs.
    lp.assemble(mat, rhs, false);

    // Solve the matrix problem.
    if (!solver->solve(mat, rhs)) error ("Matrix solver failed.\n");

		prev_d.copy(&sln_d);
		prev_dv1.copy(&sln_dv1);
		prev_dv2.copy(&sln_dv2);
		prev_e.copy(&sln_e);

		// Visualization.
		//pressure.reinit();
		//view.set_title("Pressure");
		//view.show(&pressure);

		u.reinit();
		w.reinit();
		vview.set_title("Velocity");
		vview.show(&u,&w);

		View::wait();

		/* Screenshot saving.
		char * t_n = new char[10];
		itoa(i,t_n,10);
		vview.save_screenshot(t_n);
		*/

	}
	view.close();
	vview.close();

	return 0;
}

