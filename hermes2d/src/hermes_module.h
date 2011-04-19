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

enum BoundaryConditionType {
  NATURAL, ESSENTIAL
};

enum BoundaryConditionTypeH1 {
  DIRICHLET, NEUMANN, NEWTON
};


enum GeometryType {
  PLANAR, AXISYM_X, AXISYM_Y
};

enum AnalysisType {
  STEADY_STATE, TRANSIENT, HARMONIC
};

enum AdaptivityType {
  NONE, H, P, HP
};

/* BoundaryData */

class BoundaryData {
public:
  // Constant boundary conditions
  BoundaryData(std::string marker, BoundaryConditionType type, scalar value) {
    this->markers.push_back(marker);
    this->type = type;
    this->value = value;
    this->constant_value = true;
  };

  BoundaryData(Hermes::vector<std::string> markers, BoundaryConditionType type, scalar value) {
    this->markers = markers;
    this->type = type;
    this->value = value;
    this->constant_value = true;
  };

  virtual ~BoundaryData();

  inline bool is_constant() {
    return constant_value;
  }

  Hermes::vector<std::string> markers;
  BoundaryConditionType type;
  scalar value;

protected:
  bool constant_value;
};

class BoundaryDataH1 : BoundaryData {
public:
  BoundaryDataH1(std::string marker, BoundaryConditionType type, scalar value, BoundaryConditionTypeH1 h1_type = DIRICHLET) :
    BoundaryData(marker, type, value), h1_type(h1_type) { }

  BoundaryDataH1(Hermes::vector<std::string> marker, BoundaryConditionType type, scalar value, BoundaryConditionTypeH1 h1_type = DIRICHLET) :
    BoundaryData(markers, type, value), h1_type(h1_type) { }

  BoundaryConditionTypeH1 h1_type;
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
  int init_ref;
  int init_deg;
};

struct SolutionProperties {
  int num_sol;
  MatrixSolverType mat_solver;
};

struct AdaptivityProperties {
  double tolerance;
  int max_dofs;
  int max_steps;
  double conv_exp;
  ProjNormType proj_norm_type;
  int threshold;
  int strategy;
  int regularize;
};

class ModuleProperties {
public:
  ModuleProperties();
  virtual ~ModuleProperties();

  // FIXME
  GeometryType geometry;
  AnalysisType analysis;
  AdaptivityType adaptivity_type;


  MeshProperties *mesh() {
    return mesh_properties;
  }

  SolutionProperties *solution() {
    return solution_properties;
  }

  AdaptivityProperties *adaptivity() {
    return adaptivity_properties;
  }

protected:
  MeshProperties *mesh_properties;
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

  virtual void set_mesh(const std::string &string) = 0;
  virtual void refine_mesh(Mesh *mesh, const int refinement);

  virtual void set_weakform() = 0;

  virtual void solve();

  ModuleProperties *properties() {
    return module_properties;
  }

protected:
  ModuleProperties *module_properties;

  Mesh *mesh;

  Hermes::vector<BoundaryData *> natural_boundaries;
  Hermes::vector<BoundaryData *> essential_boundaries;
  Hermes::vector<MaterialData *> materials;

  EssentialBCs *bcs;
  WeakForm *wf;
  Hermes::vector<Space *> space;
  DiscreteProblem *dp;
  Hermes::vector<Solution *> sln;
};
#endif
