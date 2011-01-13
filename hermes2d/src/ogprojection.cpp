#include "hermes2d.h"

void OGProjection::project_internal(Hermes::vector<Space *> spaces, WeakForm* wf, 
                                    scalar* target_vec, MatrixSolverType matrix_solver)
{
  _F_
  unsigned int n = spaces.size();

  // sanity checks
  if (n <= 0 || n > 10) error("Wrong number of projected functions in project_internal().");
  for (unsigned int i = 0; i < n; i++) if(spaces[i] == NULL) error("this->spaces[%d] == NULL in project_internal().", i);
  if (spaces.size() != n) error("Number of spaces must match number of projected functions in project_internal().");

  // this is needed since spaces may have their DOFs enumerated only locally.
  int ndof = Space::assign_dofs(spaces);

  // Initialize DiscreteProblem.
  bool is_linear = true;
  DiscreteProblem* dp = new DiscreteProblem(wf, spaces, is_linear);

  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  dp->assemble(matrix, rhs, false);

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
  delete wf;
}

void OGProjection::project_global(Hermes::vector<Space *> spaces, Hermes::vector<MeshFunction*> source_meshfns, 
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
    if (proj_norms == Hermes::vector<ProjNormType>()) {
      ESpaceType space_type = spaces[i]->get_type();
      switch (space_type) {
        case HERMES_H1_SPACE: norm = HERMES_H1_NORM; break;
        case HERMES_HCURL_SPACE: norm = HERMES_HCURL_NORM; break;
        case HERMES_HDIV_SPACE: norm = HERMES_HDIV_NORM; break;
        case HERMES_L2_SPACE: norm = HERMES_L2_NORM; break;
        default: error("Unknown space type in OGProjection::project_global().");
      }
    }
    else norm = proj_norms[i];
    if (norm == HERMES_H1_NORM) 
    {
      found[i] = 1;
      proj_wf->add_matrix_form(i, i, H1projection_biform<double, scalar>, H1projection_biform<Ord, Ord>);
      proj_wf->add_vector_form(i, H1projection_liform<double, scalar>, H1projection_liform<Ord, Ord>,
                               HERMES_ANY, source_meshfns[i]);
    }
    if (norm == HERMES_H1_SEMINORM) 
    {
      found[i] = 1;
      proj_wf->add_matrix_form(i, i, H1_semi_projection_biform<double, scalar>, H1_semi_projection_biform<Ord, Ord>);
      proj_wf->add_vector_form(i, H1_semi_projection_liform<double, scalar>, H1_semi_projection_liform<Ord, Ord>,
                               HERMES_ANY, source_meshfns[i]);
    }
    if (norm == HERMES_HCURL_NORM) 
    {
      found[i] = 1;
      proj_wf->add_matrix_form(i, i, Hcurlprojection_biform<double, scalar>, Hcurlprojection_biform<Ord, Ord>);
      proj_wf->add_vector_form(i, Hcurlprojection_liform<double, scalar>, Hcurlprojection_liform<Ord, Ord>,
                               HERMES_ANY, source_meshfns[i]);
    }
    if (norm == HERMES_HDIV_NORM) 
    {
      found[i] = 1;
      proj_wf->add_matrix_form(i, i, Hdivprojection_biform<double, scalar>, Hdivprojection_biform<Ord, Ord>);
      proj_wf->add_vector_form(i, Hdivprojection_liform<double, scalar>, Hdivprojection_liform<Ord, Ord>,
                               HERMES_ANY, source_meshfns[i]);
    }
    if (norm == HERMES_L2_NORM) 
    {
      found[i] = 1;
      proj_wf->add_matrix_form(i, i, L2projection_biform<double, scalar>, L2projection_biform<Ord, Ord>);
      proj_wf->add_vector_form(i, L2projection_liform<double, scalar>, L2projection_liform<Ord, Ord>,
                               HERMES_ANY, source_meshfns[i]);
    }
  }
  for (int i=0; i < n; i++) 
  {
    if (found[i] == 0) 
    {
      warn("index of component: %d\n", i);
      error("Wrong projection norm in project_global().");
    }
  }

  project_internal(spaces, proj_wf, target_vec, matrix_solver);
}

void OGProjection::project_global(Hermes::vector<Space *> spaces, Hermes::vector<Solution *> sols_src, 
                                  Hermes::vector<Solution *> sols_dest, MatrixSolverType matrix_solver, 
                                  Hermes::vector<ProjNormType> proj_norms)
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

void OGProjection::project_global(Hermes::vector<Space *> spaces, 
                                  Hermes::vector< std::pair<WeakForm::matrix_form_val_t, 
                                  WeakForm::matrix_form_ord_t> > proj_biforms, 
                                  Hermes::vector< std::pair<WeakForm::vector_form_val_t, 
                                  WeakForm::vector_form_ord_t> > proj_liforms, 
                                  Hermes::vector<MeshFunction*> source_meshfns, 
                                  scalar* target_vec, MatrixSolverType matrix_solver)
{
  _F_
  unsigned int n = spaces.size();
  unsigned int n_biforms = proj_biforms.size();
  if (n_biforms == 0)
    error("Please use the simpler version of project_global with the argument Hermes::vector<ProjNormType> proj_norms if you do not provide your own projection norm.");
  if (n_biforms != proj_liforms.size())
    error("Mismatched numbers of projection forms in project_global().");
  if (n != n_biforms)
    error("Mismatched numbers of projected functions and projection forms in project_global().");

  // This is needed since spaces may have their DOFs enumerated only locally
  // when they come here.
  int ndof = Space::assign_dofs(spaces);

  // Define projection weak form.
  WeakForm* proj_wf = new WeakForm(n);
  for (unsigned int i = 0; i < n; i++) {
    proj_wf->add_matrix_form(i, i, proj_biforms[i].first, proj_biforms[i].second);
    proj_wf->add_vector_form(i, proj_liforms[i].first, proj_liforms[i].second,
                    HERMES_ANY, source_meshfns[i]);
  }

  project_internal(spaces, proj_wf, target_vec, matrix_solver);
}

void OGProjection::project_global(Space *space, ExactFunction2 source_fn, scalar* target_vec, 
                                  MatrixSolverType matrix_solver)
{
  _F_
  ProjNormType norm;
  ESpaceType space_type = space->get_type();
  switch (space_type) {
    case HERMES_H1_SPACE: norm = HERMES_H1_NORM; break;
    case HERMES_HCURL_SPACE: norm = HERMES_HCURL_NORM; break;
    case HERMES_HDIV_SPACE: norm = HERMES_HDIV_NORM; break;
    case HERMES_L2_SPACE: norm = HERMES_L2_NORM; break;
    default: error("Unknown space type in OGProjection::project_global().");
  }
  // Since the projected function is vector-valued, H1 and L2 spaces are not admissible.
  if (space_type == HERMES_H1_SPACE) error("Mismatched space and projection norm in OGProjection::project_global().");
  if (space_type == HERMES_L2_SPACE) error("Mismatched space and projection norm in OGProjection::project_global().");
  Mesh *mesh = space->get_mesh();
  if (mesh == NULL) error("Mesh is NULL in project_global().");
  Solution source_sln;
  source_sln.set_exact(mesh, source_fn);
  project_global(space, (MeshFunction*)&source_sln, target_vec, matrix_solver, norm);
};

void OGProjection::project_global(Space *space, ExactFunction source_fn, scalar* target_vec, 
                                  MatrixSolverType matrix_solver)
{
  _F_
  ProjNormType norm;
  ESpaceType space_type = space->get_type();
  switch (space_type) {
    case HERMES_H1_SPACE: norm = HERMES_H1_NORM; break;
    case HERMES_HCURL_SPACE: norm = HERMES_HCURL_NORM; break;
    case HERMES_HDIV_SPACE: norm = HERMES_HDIV_NORM; break;
    case HERMES_L2_SPACE: norm = HERMES_L2_NORM; break;
    default: error("Unknown space type in OGProjection::project_global().");
  }
  // Since the projected function is scalar, Hcurl and Hdiv spaces are not admissible.
  if (space_type == HERMES_HCURL_SPACE) error("Mismatched space and projection norm in OGProjection::project_global().");
  if (space_type == HERMES_HDIV_SPACE) error("Mismatched space and projection norm in OGProjection::project_global().");

  Mesh *mesh = space->get_mesh();
  if (mesh == NULL) error("Mesh is NULL in project_global().");
  Solution source_sln;
  source_sln.set_exact(mesh, source_fn);
  project_global(space, (MeshFunction*)&source_sln, target_vec, matrix_solver, norm);
};

void OGProjection::project_local(Space *space, int proj_norm, ExactFunction source_fn, Mesh* mesh,
                   scalar* target_vec)
{
  _F_
  error("OGProjection::project_local(): not implemented.");
  /// TODO
}
