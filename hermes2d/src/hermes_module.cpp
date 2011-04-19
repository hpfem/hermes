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

void HermesModule::set_boundary(BoundaryData *boundary) {
  if(boundary->bc_type == HERMES_ESSENTIAL)
    this->essential_boundaries.push_back(boundary);
  else if (boundary->bc_type == HERMES_NATURAL)
    this->natural_boundaries.push_back(boundary);
}

void HermesModule::set_material(MaterialData *material) {
  this->materials.push_back(material);
}

void HermesModule::solve() {
  RefinementSelectors::Selector* selector = NULL;
  Hermes::vector<RefinementSelectors::Selector *> selectors;

  if (properties()->adaptivity()->cand_list != H2D_NONE)
    selector = new RefinementSelectors::H1ProjBasedSelector(this->properties()->adaptivity()->cand_list,
                                                            properties()->adaptivity()->conv_exp,
                                                            H2DRS_DEFAULT_ORDER);


  for (int i = 0; i < properties()->solution()->num_sol; i++)
  {
      this->slns.push_back(new Solution(this->meshes.at(i)));

      if (properties()->adaptivity()->cand_list != H2D_NONE)
          selectors.push_back(selector);
  }

  set_weakform();

  SparseMatrix* matrix = create_matrix(properties()->solver()->mat_solver);
  Vector *rhs = create_vector(properties()->solver()->mat_solver);
  Solver *solver = create_linear_solver(properties()->solver()->mat_solver, matrix, rhs);

  if (properties()->adaptivity()->cand_list == H2D_NONE)
  {

    DiscreteProblem dp(this->wf, this->spaces, this->properties()->is_linear);
    dp.assemble(matrix, rhs);

    if(solver->solve())
      Solution::vector_to_solutions(solver->get_solution(), this->spaces, this->slns);
  }
  else
  {
    Hermes::vector<Solution *> ref_slns;

    for (int i = 0; i < properties()->adaptivity()->max_steps; i++)
    {
      for (int j = 0; j < properties()->solution()->num_sol; j++)
           ref_slns.push_back(new Solution());

      Hermes::vector<Space *> ref_spaces = *Space::construct_refined_spaces(this->spaces);

      DiscreteProblem dp(this->wf, ref_spaces, true); // FIXME
      dp.assemble(matrix, rhs);

      if (solver->solve())
        Solution::vector_to_solutions(solver->get_solution(), ref_spaces, ref_slns);
      else
        break;

      OGProjection::project_global(this->spaces, ref_slns, this->slns, properties()->solver()->mat_solver);

      Adapt adaptivity(this->spaces, proj_norms);

      Hermes::vector<double> err_est_rel;
      double error = adaptivity.calc_err_est(this->slns, ref_slns, &err_est_rel) * 100;

      if (error < properties()->adaptivity()->tolerance || Space::get_num_dofs(this->spaces) >= properties()->adaptivity()->max_dofs)
        break;

      if (i != properties()->adaptivity()->max_steps-1)
        adaptivity.adapt(selectors, properties()->adaptivity()->threshold, properties()->adaptivity()->strategy,
                         properties()->adaptivity()->regularize);

      for (unsigned int i = 0; i < ref_spaces.size(); i++)
      {
        delete ref_spaces.at(i);
      }

      ref_spaces.clear();

      for (unsigned int i = 0; i < ref_slns.size(); i++)
        delete ref_slns.at(i);

      ref_slns.clear();
    }
  }

  delete solver;
  delete matrix;
  delete rhs;

  if (this->properties()->adaptivity()->cand_list != H2D_NONE)
  {
    delete selector;
    selectors.clear();
  }

  this->meshes.clear();

  for (unsigned int i = 0; i < this->spaces.size(); i++)
      delete this->spaces.at(i);
  this->spaces.clear();

  for (unsigned int i = 0; i < this->slns.size(); i++)
      delete this->slns.at(i);
  this->slns.clear();
}
