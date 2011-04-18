#include "hermes3d.h"


void OGProjection::project_internal(Hermes::vector<Space *> spaces, WeakForm *proj_wf, scalar* target_vec, MatrixSolverType matrix_solver)
{
  _F_
  unsigned int n = spaces.size();

  // sanity checks
  if (n <= 0 || n > 10) error("Wrong number of projected functions in project_internal().");
  for (unsigned int i = 0; i < n; i++) if(spaces[i] == NULL) error("this->spaces[%d] == NULL in project_internal().", i);
  if (spaces.size() != n) error("Number of spaces must matchnumber of projected functions in project_internal().");

  // this is needed since spaces may have their DOFs enumerated only locally.
  int ndof = Space::assign_dofs(spaces);

  // Initialize DiscreteProblem.
  bool is_linear = true;
  DiscreteProblem* dp = new DiscreteProblem(proj_wf, spaces, is_linear);

  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  dp->assemble(matrix, rhs);

  // Calculate the coefficient vector.
  bool solved = solver->solve();
  scalar* coeffs;
  if (solved) 
    coeffs = solver->get_solution();

  if (target_vec != NULL) 
    for (int i=0; i<ndof; i++) target_vec[i] = coeffs[i];
    
  delete solver;
  delete matrix;
  delete rhs;
  delete dp;
  delete proj_wf;
}


void OGProjection::project_global(Hermes::vector<Space *> spaces, Hermes::vector<MeshFunction *> source_meshfns, 
                              scalar* target_vec, MatrixSolverType matrix_solver, Hermes::vector<ProjNormType> proj_norms)
{
  _F_
  int n = spaces.size();  

  // define temporary projection weak form
  WeakForm* proj_wf = new WeakForm(n);
  int found[100];
  for (int i = 0; i < 100; i++) found[i] = 0;
  for (int i = 0; i < n; i++) 
  {
    int norm;
    if (proj_norms == Hermes::vector<ProjNormType>()) norm = HERMES_DEFAULT_PROJ_NORM;
    else norm = proj_norms[i];
    if (norm == HERMES_L2_NORM) 
    {
      found[i] = 1;
      proj_wf->add_matrix_form(i, i, L2projection_biform<double, scalar>, L2projection_biform<Ord, Ord>);
      Hermes::vector<MeshFunction*> mesh_fns;
      mesh_fns.push_back(source_meshfns[i]);
      proj_wf->add_vector_form(i, L2projection_liform<double, scalar>, L2projection_liform<Ord, Ord>,
                               HERMES_ANY_INT, mesh_fns);
    }
    if (norm == HERMES_H1_NORM) 
    {
      found[i] = 1;
      proj_wf->add_matrix_form(i, i, H1projection_biform<double, scalar>, H1projection_biform<Ord, Ord>);
      Hermes::vector<MeshFunction*> mesh_fns;
      mesh_fns.push_back(source_meshfns[i]);
      proj_wf->add_vector_form(i, H1projection_liform<double, scalar>, H1projection_liform<Ord, Ord>,
                               HERMES_ANY_INT, mesh_fns);
    }
    if (norm == HERMES_H1_SEMINORM) 
    {
      found[i] = 1;
      proj_wf->add_matrix_form(i, i, H1_semi_projection_biform<double, scalar>, H1_semi_projection_biform<Ord, Ord>);
      Hermes::vector<MeshFunction*> mesh_fns;
      mesh_fns.push_back(source_meshfns[i]);
      proj_wf->add_vector_form(i, H1_semi_projection_liform<double, scalar>, H1_semi_projection_liform<Ord, Ord>,
                               HERMES_ANY_INT, mesh_fns);
    }
    if (norm == HERMES_HCURL_NORM) 
    {
      found[i] = 1;
      proj_wf->add_matrix_form(i, i, Hcurlprojection_biform<double, scalar>, Hcurlprojection_biform<Ord, Ord>);
      Hermes::vector<MeshFunction*> mesh_fns;
      mesh_fns.push_back(source_meshfns[i]);
      proj_wf->add_vector_form(i, Hcurlprojection_liform<double, scalar>, Hcurlprojection_liform<Ord, Ord>,
                               HERMES_ANY_INT, mesh_fns);
    }
  }
  for (int i=0; i < n; i++) 
  {
    if (found[i] == 0) 
    {
      printf("index of component: %d\n", i);
      error("Wrong projection norm in project_global().");
    }
  }

  project_internal(spaces, proj_wf, target_vec, matrix_solver);
}

void OGProjection::project_global(Hermes::vector<Space *> spaces, 
                              Hermes::vector<Solution*> sols_src, Hermes::vector<Solution*> sols_dest, 
                              MatrixSolverType matrix_solver, Hermes::vector<ProjNormType> proj_norms)
{
  _F_
  scalar* target_vec = new scalar[Space::get_num_dofs(spaces)];
  Hermes::vector<MeshFunction *> ref_slns_mf;
  for (unsigned int i = 0; i < sols_src.size(); i++) 
    ref_slns_mf.push_back(static_cast<MeshFunction*>(sols_src[i]));
  
  OGProjection::project_global(spaces, ref_slns_mf, target_vec, matrix_solver, proj_norms);
  
  Solution::vector_to_solutions(target_vec, spaces, sols_dest);
  
  delete [] target_vec;
}

void OGProjection::project_global(Space* space, 
                              Solution* sol_src, Solution* sol_dest, 
                              MatrixSolverType matrix_solver, ProjNormType proj_norm)
{
  Hermes::vector<Space *> spaces;
  spaces.push_back(space);
  Hermes::vector<Solution *> sols_src;
  sols_src.push_back(sol_src);
  Hermes::vector<Solution *> sols_dest;
  sols_dest.push_back(sol_dest);
  Hermes::vector<ProjNormType> proj_norms;
  proj_norms.push_back(proj_norm);
  project_global(spaces, sols_src, sols_dest, matrix_solver, proj_norms);
}
