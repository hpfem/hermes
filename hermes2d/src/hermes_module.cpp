// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#include "hermes_module.h"

void HermesModule::add_mesh(Mesh *mesh)
{
  this->meshes.push_back(mesh);
}

Mesh *HermesModule::get_mesh(int index)
{
  return this->meshes.at(index);
}

Hermes::vector<Mesh *> HermesModule::get_meshes()
{
  return this->meshes;
}

void HermesModule::add_boundary(BoundaryData *boundary)
{
  this->boundaries.push_back(boundary);
}

void HermesModule::add_material(MaterialData *material)
{
  this->materials.push_back(material);
}

Space *HermesModule::get_space(int index)
{
  return this->spaces.at(index);
}

Hermes::vector<Space *> HermesModule::get_spaces()
{
  return this->spaces;
}

bool HermesModule::solve()
{
  Hermes2D hermes2d;
  TimePeriod cpu_time;

  bool success = false;

  /*
  RefinementSelectors::Selector* selector = NULL;
  Hermes::vector<RefinementSelectors::Selector *> selectors;

  if (this->properties()->adaptivity()->cand_list != H2D_NONE)
    selector = new RefinementSelectors::H1ProjBasedSelector(this->properties()->adaptivity()->cand_list,
                                                            properties()->adaptivity()->conv_exp,
                                                            H2DRS_DEFAULT_ORDER);
  */

  // Set essential BCs.
  this->set_essential_bcs();

  // Initialize weak forms.
  this->set_weakforms();

  // Initialize spaces.
  this->set_spaces();

  // Initialize solutions.
  for (int i = 0; i < this->properties()->solution()->num_sol; i++)
  {
      this->slns.push_back(new Solution(this->meshes.at(i)));

      //if (this->properties()->adaptivity()->cand_list != H2D_NONE)
      //    selectors.push_back(selector);
  }

  // Initialize matrix solver.
  SparseMatrix *matrix = create_matrix(this->properties()->solver()->mat_solver);
  Vector *rhs = create_vector(this->properties()->solver()->mat_solver);
  Solver *solver = create_linear_solver(this->properties()->solver()->mat_solver, matrix, rhs);

  if (this->properties()->adaptivity()->cand_list == H2D_NONE)
  {
    int ndof = Space::get_num_dofs(this->get_spaces());
    if (ndof != 0)
      info("ndof = %d", ndof);
    else
      error("ndof = %d", ndof);

    DiscreteProblem dp(wf, this->spaces);

    scalar* coeff_vec = new scalar[ndof];
    memset(coeff_vec, 0, ndof*sizeof(scalar));

    if (hermes2d.solve_newton(coeff_vec, &dp, solver, matrix, rhs, this->properties()->solver()->jacobian_changed,
            this->properties()->solver()->newton_tol, this->properties()->solver()->newton_max_iter))
      success = true;
    else
      error("Newton's iteration failed.");

    Solution::vector_to_solutions(coeff_vec, this->spaces, this->slns);
    delete [] coeff_vec;
  }
  /*
  else
  {
    for (int i = 0; i <= this->properties()->adaptivity()->max_steps; i++)
    {
      info("---- Adaptivity step %d:", i);

      Hermes::vector<Solution *> ref_slns;
      for (int j = 0; j < this->properties()->solution()->num_sol; j++)
           ref_slns.push_back(new Solution());

      Hermes::vector<Space *> ref_spaces = *Space::construct_refined_spaces(this->spaces);

      DiscreteProblem dp(this->wf, ref_spaces);

      int ndof_ref = Space::get_num_dofs(ref_spaces);
      if (ndof_ref == 0)
        error("ndof_fine = %d", ndof_ref);

      scalar* coeff_vec = new scalar[ndof_ref];
      memset(coeff_vec, 0, ndof_ref * sizeof(scalar));

      if (!hermes2d.solve_newton(coeff_vec, &dp, solver, matrix, rhs, this->properties()->solver()->jacobian_changed,
            this->properties()->solver()->newton_tol, this->properties()->solver()->newton_max_iter))
        error("Newton's iteration failed.");

      Solution::vector_to_solutions(solver->get_solution(), ref_spaces, ref_slns);

      OGProjection::project_global(this->spaces, ref_slns, this->slns, this->properties()->solver()->mat_solver);

      Adapt adaptivity(this->spaces, this->proj_norms);

      Hermes::vector<double> err_est_rel;
      double error = adaptivity.calc_err_est(this->slns, ref_slns, &err_est_rel) * 100;

      info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%",
        Space::get_num_dofs(this->spaces), ndof_ref, err_est_rel);

      if (error < this->properties()->adaptivity()->tolerance || Space::get_num_dofs(this->spaces) >= this->properties()->adaptivity()->max_dofs)
        break;

      if (i != this->properties()->adaptivity()->max_steps-1)
        adaptivity.adapt(selectors, this->properties()->adaptivity()->threshold, this->properties()->adaptivity()->strategy,
                         this->properties()->adaptivity()->regularize);

      for (unsigned int i = 0; i < ref_spaces.size(); i++)
        delete ref_spaces.at(i);
      ref_spaces.clear();

      for (unsigned int i = 0; i < ref_slns.size(); i++)
        delete ref_slns.at(i);
      ref_slns.clear();
    }
  }
  */

  delete solver;
  delete matrix;
  delete rhs;
  /*
  if (this->properties()->adaptivity()->cand_list != H2D_NONE)
  {
    delete selector;
    selectors.clear();
  }
  */

  return success;
}

Solution *HermesModule::get_solution(int index)
{
  return this->slns.at(index);
}

Hermes::vector<Solution *> HermesModule::get_solutions()
{
  return this->slns;
}
