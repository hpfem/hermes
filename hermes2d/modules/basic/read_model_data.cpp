// This file is part of the module Basic.
// It reads model-dependent data.

  // Read list of material markers.
  int n_mat_markers;
  if(!Get(f, &n_mat_markers)) error("Could not read number of material markers.");
  info("n_mat_markers: %d", n_mat_markers);
  if(n_mat_markers <= 0) error("At least one material marker must be given.");
  Hermes::vector<int> mat_markers;
  for (int i = 0; i < n_mat_markers; i++) {
    int tmp;
    if(!Get(f, &tmp)) error("Could not read a material marker.");
    info("mat_marker[%d]: %d", i, tmp);
    mat_markers.push_back(tmp);
  }
  B.set_material_markers(mat_markers);

  // Read list of c1 constants.
  std::vector<double> c1_array;
  for (int i = 0; i < n_mat_markers; i++) {
    double tmp;
    if(!Get(f, &tmp)) error("Could not read a c1 constant.");
    info("c1_array[%d]: %g", i, tmp);
    c1_array.push_back(tmp);
  }
  B.set_c1_array(c1_array);

  // Read list of c2 constants.
  std::vector<double> c2_array;
  for (int i = 0; i < n_mat_markers; i++) {
    double tmp;
    if(!Get(f, &tmp)) error("Could not read a c2 constant.");
    info("c2_array[%d]: %g", i, tmp);
    c2_array.push_back(tmp);
  }
  B.set_c2_array(c1_array);

  // Read list of c3 constants.
  std::vector<double> c3_array;
  for (int i = 0; i < n_mat_markers; i++) {
    double tmp;
    if(!Get(f, &tmp)) error("Could not read a c3 constant.");
    info("c3_array[%d]: %g", i, tmp);
    c3_array.push_back(tmp);
  }
  B.set_c3_array(c1_array);

  // Read list of c4 constants.
  std::vector<double> c4_array;
  for (int i = 0; i < n_mat_markers; i++) {
    double tmp;
    if(!Get(f, &tmp)) error("Could not read a c4 constant.");
    info("c4_array[%d]: %g", i, tmp);
    c4_array.push_back(tmp);
  }
  B.set_c4_array(c1_array);

  // Read list of c5 constants.
  std::vector<double> c5_array;
  for (int i = 0; i < n_mat_markers; i++) {
    double tmp;
    if(!Get(f, &tmp)) error("Could not read a c5 constant.");
    info("c5_array[%d]: %g", i, tmp);
    c5_array.push_back(tmp);
  }
  B.set_c5_array(c1_array);

  // Read Dirichlet boundary markers.
  int n_bc_dirichlet;
  if(!Get(f, &n_bc_dirichlet)) error("Could not read number of Dirichlet boundary markers.");
  info("n_bc_dirichlet: %d", n_bc_dirichlet);
  Hermes::vector<int> bdy_markers_dirichlet;
  for (int i = 0; i < n_bc_dirichlet; i++) {
    int tmp;
    if(!Get(f, &tmp)) error("Could not read a VALUE boundary marker.");
    info("bdy_markers_dirichlet[%d]: %d", i, tmp);
    bdy_markers_dirichlet.push_back(tmp);
  }
  if (n_bc_dirichlet > 0) B.set_dirichlet_markers(bdy_markers_dirichlet);

  // Read Dirichlet boundary values.
  Hermes::vector<double> bdy_values_dirichlet;
  for (int i = 0; i < n_bc_dirichlet; i++) {
    double tmp;
    if(!Get(f, &tmp)) error("Could not read a Dirichlet boundary value.");
    info("bdy_values_dirichlet[%d]: %g", i, tmp);
    bdy_values_dirichlet.push_back(tmp);
  }
  if (n_bc_dirichlet > 0) B.set_dirichlet_values(bdy_markers_dirichlet, bdy_values_dirichlet);

  // Read Neumann boundary markers.
  int n_bc_neumann;
  if(!Get(f, &n_bc_neumann)) error("Could not read number of Neumann boundary markers.");
  info("n_bc_neumann: %d", n_bc_neumann);
  Hermes::vector<int> bdy_markers_neumann;
  for (int i = 0; i < n_bc_neumann; i++) {
    int tmp;
    if(!Get(f, &tmp)) error("Could not read a Neumann boundary marker.");
    info("bdy_markers_neumann[%d]: %d", i, tmp);
    bdy_markers_neumann.push_back(tmp);
  }
  if (n_bc_neumann > 0) B.set_neumann_markers(bdy_markers_neumann);

  // Read list of Neumann boundary values.
  Hermes::vector<double> bdy_values_neumann;
  for (int i = 0; i < n_bc_neumann; i++) {
    double tmp;
    if(!Get(f, &tmp)) error("Could not read a Neumann boundary value.");
    info("bdy_values_neumann[%d]: %g", i, tmp);
    bdy_values_neumann.push_back(tmp);
  }
  if (n_bc_neumann > 0) B.set_neumann_values(bdy_values_neumann);

  // Read Newton boundary markers.
  int n_bc_newton;
  if(!Get(f, &n_bc_newton)) error("Could not read number of Newton boundary markers.");
  info("n_bc_newton: %d", n_bc_newton);
  Hermes::vector<int> bdy_markers_newton;
  for (int i = 0; i < n_bc_newton; i++) {
    int tmp;
    if(!Get(f, &tmp)) error("Could not read a Newton boundary marker.");
    info("bdy_markers_newton[%d]: %d", i, tmp);
    bdy_markers_newton.push_back(tmp);
  }
  if (n_bc_newton > 0) B.set_newton_markers(bdy_markers_newton);

  // Read list of Newton boundary value pairs.
  Hermes::vector<double_pair> bdy_values_newton;
  for (int i = 0; i < n_bc_newton; i++) {
    double tmp1, tmp2;
    if(!Get(f, &tmp1)) error("Could not read a Newton boundary value (first in pair).");
    if(!Get(f, &tmp2)) error("Could not read a Newton boundary value (second in pair).");
    info("bdy_values_newton[%d]: %g %g", i, tmp1, tmp2);
    bdy_values_newton.push_back(std::make_pair(tmp1, tmp2));
  }
  if (n_bc_newton > 0) B.set_newton_values(bdy_values_newton);
