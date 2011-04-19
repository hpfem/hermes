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

#include "hermes2d.h"

#ifndef __H2D_HERMES_MODULE_H
#define __H2D_HERMES_MODULE_H

using namespace RefinementSelectors;

/* BoundaryData */

class BoundaryData {
public:
  // Constant boundary conditions
  BoundaryData(std::string marker, BCType bc_type) {
    this->markers.push_back(marker);
    this->bc_type = bc_type;
  };

  BoundaryData(Hermes::vector<std::string> markers, BCType bc_type) {
    this->markers = markers;
    this->bc_type = bc_type;
  };

  virtual ~BoundaryData();

  Hermes::vector<std::string> markers;
  BCType bc_type;
};

class BoundaryDataH1 : public BoundaryData {
public:
  BoundaryDataH1(std::string marker, BCTypeH1 bc_type_h1, scalar value1, scalar value2 = 0) :
    BoundaryData(marker, (bc_type_h1 == HERMES_DIRICHLET) ? HERMES_ESSENTIAL : HERMES_NATURAL), value1(value1), value2(value2) { }

  BoundaryDataH1(Hermes::vector<std::string> markers, BCTypeH1 bc_type_h1, scalar value1, scalar value2 = 0) :
    BoundaryData(markers, (bc_type_h1 == HERMES_DIRICHLET) ? HERMES_ESSENTIAL : HERMES_NATURAL), value1(value1), value2(value2) { }

  BCTypeH1 bc_type_h1;
  scalar value1, value2;
};

/* MaterialData */

class MaterialData {
public:
  MaterialData(std::string marker) {
    this->markers.push_back(marker);
  };
  MaterialData(Hermes::vector<std::string> markers) {
    this->markers = markers;
  };

  virtual ~MaterialData();

  Hermes::vector<std::string> markers;
};

/* ModuleProperties */

struct MeshProperties {
  int init_deg;
};

struct SolverProperties {
  MatrixSolverType mat_solver;
};

struct SolutionProperties {
  int num_sol;
};

struct AdaptivityProperties {
  CandList cand_list;
  double tolerance;
  int max_dofs;
  int max_steps;
  double conv_exp;
  int threshold;
  int strategy;
  int regularize;
};

class ModuleProperties {
public:
  ModuleProperties();
  virtual ~ModuleProperties();

  GeomType geometry;

  AnalysisType analysis;

  bool is_linear;

  MeshProperties *mesh() {
    return mesh_properties;
  }

  SolverProperties *solver() {
    return solver_properties;
  }

  SolutionProperties *solution() {
    return solution_properties;
  }

  AdaptivityProperties *adaptivity() {
    return adaptivity_properties;
  }

protected:
  MeshProperties *mesh_properties;
  SolverProperties *solver_properties;
  SolutionProperties *solution_properties;
  AdaptivityProperties *adaptivity_properties;
};

/* HermesModule */

class HermesModule {
public:
  HermesModule() {
  };
  virtual ~HermesModule();

  virtual void set_boundary(BoundaryData *boundary);
  virtual void set_material(MaterialData *material);

  virtual void set_meshes(Hermes::vector<std::string> meshes) = 0; // FIXME (Input from file.)
  virtual void set_spaces(Hermes::vector<Space *> spaces) = 0;
  virtual void set_proj_norms(Hermes::vector<ProjNormType *> proj_norms) = 0;
  virtual void set_weakform() = 0;

  virtual void solve();

  ModuleProperties *properties() {
    return module_properties;
  }

protected:
  ModuleProperties *module_properties;

  Hermes::vector<Mesh *> meshes;

  Hermes::vector<BoundaryData *> natural_boundaries;
  Hermes::vector<BoundaryData *> essential_boundaries;
  Hermes::vector<MaterialData *> materials;

  EssentialBCs *bcs;
  WeakForm *wf;
  Hermes::vector<Space *> spaces;
  Hermes::vector<ProjNormType> proj_norms ;
  DiscreteProblem *dp;
  Hermes::vector<Solution *> slns;
};
#endif
