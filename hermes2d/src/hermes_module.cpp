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
  if(boundary->type == ESSENTIAL)
    this->essential_boundaries.push_back(boundary);
  else if (boundary->type == NATURAL)
    this->natural_boundaries.push_back(boundary);
}

void HermesModule::set_material(MaterialData *material) {
  this->materials.push_back(material);
}

void HermesModule::refine_mesh(Mesh *mesh, const int refinement) {
  for (int i = 0; i < refinement; i++)
    this->mesh->refine_all_elements();
}

void HermesModule::solve() {
  if (properties()->mesh()->init_ref != 0)
    refine_mesh(this->mesh, properties()->mesh()->init_ref);

  RefinementSelectors::Selector *select = NULL;
  switch (properties()->adaptivity_type)
  {
    case H:
      select = new RefinementSelectors::HOnlySelector();
      break;
    case P:
      select = new RefinementSelectors::H1ProjBasedSelector(RefinementSelectors::H2D_P_ANISO,
                                                            properties()->adaptivity()->conv_exp,
                                                            H2DRS_DEFAULT_ORDER);
      break;
    case HP:
      select = new RefinementSelectors::H1ProjBasedSelector(RefinementSelectors::H2D_HP_ANISO,
                                                            properties()->adaptivity()->conv_exp,
                                                            H2DRS_DEFAULT_ORDER);
      break;
  }

  Hermes::vector<ProjNormType> proj_norm_type;
  Hermes::vector<RefinementSelectors::Selector *> selector;
  for (int i = 0; i < properties()->solution()->num_sol; i++)
  {
      space.push_back(new H1Space(this->mesh, bcs, properties()->mesh()->init_deg)); // FIXME
      sln.push_back(new Solution());

      if (properties()->adaptivity_type != NONE)
      {
          proj_norm_type.push_back(properties()->adaptivity()->proj_norm_type);
          selector.push_back(select);
      }
  }

  set_weakform();

  SparseMatrix* matrix = create_matrix(properties()->solution()->mat_solver);;
  Vector *rhs = create_vector(properties()->solution()->mat_solver);
  Solver *solver = create_linear_solver(properties()->solution()->mat_solver, matrix, rhs);

  if (properties()->adaptivity_type == NONE)
  {

    DiscreteProblem dp(this->wf, this->space, true); // FIXME
    dp.assemble(matrix, rhs);

    if(solver->solve())
      Solution::vector_to_solutions(solver->get_solution(), this->space, this->sln);
  }
  else
  {
    Hermes::vector<Solution *> ref_sln;

    for (int i = 0; i < properties()->adaptivity()->max_steps; i++)
    {
      for (int j = 0; j < properties()->solution()->num_sol; j++)
           ref_sln.push_back(new Solution());

      Hermes::vector<Space *> ref_space = *Space::construct_refined_spaces(this->space);

      DiscreteProblem dp(this->wf, ref_space, true); // FIXME
      dp.assemble(matrix, rhs);

      if (solver->solve())
        Solution::vector_to_solutions(solver->get_solution(), ref_space, ref_sln);
      else
        break;

      OGProjection::project_global(this->space, ref_sln, this->sln, properties()->solution()->mat_solver);

      Adapt adaptivity(this->space, proj_norm_type);

      Hermes::vector<double> err_est_rel;
      double error = adaptivity.calc_err_est(this->sln, ref_sln, &err_est_rel) * 100;

      if (error < properties()->adaptivity()->tolerance || Space::get_num_dofs(this->space) >= properties()->adaptivity()->max_dofs)
        break;

      if (i != properties()->adaptivity()->max_steps-1)
        adaptivity.adapt(selector, properties()->adaptivity()->threshold, properties()->adaptivity()->strategy,
                         properties()->adaptivity()->regularize);

      for (unsigned int i = 0; i < ref_space.size(); i++)
      {
        delete ref_space.at(i);
      }

      ref_space.clear();

      for (unsigned int i = 0; i < ref_sln.size(); i++)
        delete ref_sln.at(i);

      ref_sln.clear();
    }
  }

  delete solver;
  delete matrix;
  delete rhs;

  delete select;
  selector.clear();

  delete this->mesh;

  for (unsigned int i = 0; i < this->space.size(); i++)
      delete this->space.at(i);
  this->space.clear();

  for (unsigned int i = 0; i < this->sln.size(); i++)
      delete this->sln.at(i);
  this->sln.clear();
}
