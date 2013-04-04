/// This file is part of Hermes2D.
///
/// Hermes2D is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 2 of the License, or
/// (at your option) any later version.
///
/// Hermes2D is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY;without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with Hermes2D. If not, see <http:///www.gnu.org/licenses/>.

#ifndef __H2D_DISCRETE_PROBLEM_DG_ASSEMBLER_H
#define __H2D_DISCRETE_PROBLEM_DG_ASSEMBLER_H

#include "hermes_common.h"
#include "adapt/adapt.h"
#include "graph.h"
#include "forms.h"
#include "weakform/weakform.h"
#include "function/function.h"
#include "neighbor.h"
#include "refinement_selectors/selector.h"
#include "exceptions.h"
#include "mixins2d.h"
#include "multimesh_dg_neighbor_tree.h"
#include "discrete_problem_dg_assembly_data.h"
#include "discrete_problem/discrete_problem_helpers.h"
#include "multimesh_dg_neighbor_tree_node.h"

namespace Hermes
{
  namespace Hermes2D
  {
    class PrecalcShapeset;
    /// @ingroup inner
    /// Discrete problem class.
    ///
    /// This class does assembling into external matrix / vector structures.
    ///
    template<typename Scalar>
    class HERMES_API DiscreteProblemDGAssembler : 
      public Hermes::Mixins::TimeMeasurable, 
      public Hermes::Hermes2D::Mixins::SettableSpaces<Scalar>, 
      public Hermes::Hermes2D::Mixins::StateQueryable,
      public Hermes::Hermes2D::Mixins::DiscreteProblemWeakForm<Scalar>,
      public Hermes::Hermes2D::Mixins::DiscreteProblemRungeKutta<Scalar>,
      public Hermes::Hermes2D::Mixins::DiscreteProblemMatrixVector<Scalar>,
      public Hermes::Hermes2D::Mixins::DiscreteProblemSingleAssemblyData
    {
    public:
      /// Constructor for multiple components / equations.
      DiscreteProblemDGAssembler(WeakForm<Scalar>* wf, Hermes::vector<SpaceSharedPtr<Scalar> > spaces);

      /// Constructor for one equation.
      DiscreteProblemDGAssembler(WeakForm<Scalar>* wf, SpaceSharedPtr<Scalar> space);

      DiscreteProblemDGAssembler();

      void init();
      void init_assembling();
      virtual std::string getClassName() const { return "DiscreteProblemDGAssembler"; }

      void free();
      bool is_linear;

      /// The form will be assembled.
      bool form_to_be_assembled(MatrixForm<Scalar>* form, Traverse::State* current_state);
      bool form_to_be_assembled(MatrixFormVol<Scalar>* form, Traverse::State* current_state);
      bool form_to_be_assembled(MatrixFormSurf<Scalar>* form, Traverse::State* current_state);
      bool form_to_be_assembled(VectorForm<Scalar>* form, Traverse::State* current_state);
      bool form_to_be_assembled(VectorFormVol<Scalar>* form, Traverse::State* current_state);
      bool form_to_be_assembled(VectorFormSurf<Scalar>* form, Traverse::State* current_state);

      bool isOkay() const;

      Hermes::vector<SpaceSharedPtr<Scalar> > get_spaces() const
      {
        return this->spaces;
      }

      /// There is a matrix form set on DG_INNER_EDGE area or not.
      bool DG_matrix_forms_present;

      /// There is a vector form set on DG_INNER_EDGE area or not.
      bool DG_vector_forms_present;

      /// Sets new spaces for the instance.
      virtual void set_spaces(Hermes::vector<SpaceSharedPtr<Scalar> > spaces);
      virtual void set_space(SpaceSharedPtr<Scalar> space);

      /// Assemble DG forms.
      void assemble_one_DG_state(PrecalcShapeset** current_pss, RefMap** current_refmaps,  Solution<Scalar>** current_u_ext, AsmList<Scalar>** current_als,
        Traverse::State* current_state, Hermes::vector<MatrixFormDG<Scalar>*> current_mfDG, Hermes::vector<VectorFormDG<Scalar>*> current_vfDG, Transformable** fn, WeakForm<Scalar>* current_wf);

      /// Assemble one DG neighbor.
      void assemble_DG_one_neighbor(bool edge_processed, unsigned int neighbor_i,
        PrecalcShapeset** current_pss, RefMap** current_refmaps,  Solution<Scalar>** current_u_ext, AsmList<Scalar>** current_als,
        Traverse::State* current_state, Hermes::vector<MatrixFormDG<Scalar>*> current_mfDG, Hermes::vector<VectorFormDG<Scalar>*> current_vfDG, Transformable** fn,
        std::map<unsigned int, PrecalcShapeset *> npss, std::map<unsigned int, PrecalcShapeset *> nspss, std::map<unsigned int, RefMap *> nrefmap,
        LightArray<NeighborSearch<Scalar>*>& neighbor_searches, unsigned int min_dg_mesh_seq, WeakForm<Scalar>* current_wf);

      /// Assemble DG matrix forms.
      void assemble_DG_matrix_forms(PrecalcShapeset** current_pss, RefMap** current_refmaps, AsmList<Scalar>** current_als,
        Traverse::State* current_state, MatrixFormDG<Scalar>** current_mfDG, std::map<unsigned int, PrecalcShapeset*> npss,
        std::map<unsigned int, PrecalcShapeset*> nspss, std::map<unsigned int, RefMap*> nrefmap, LightArray<NeighborSearch<Scalar>*>& neighbor_searches);

      /// Assemble DG vector forms.
      void assemble_DG_vector_forms(RefMap** current_refmaps, AsmList<Scalar>** current_als,
        Traverse::State* current_state, VectorFormDG<Scalar>** current_vfDG, std::map<unsigned int, PrecalcShapeset*> nspss,
        std::map<unsigned int, RefMap*> nrefmap, LightArray<NeighborSearch<Scalar>*>& neighbor_searches);

      /// Initialize external functions for DG forms.
      DiscontinuousFunc<Scalar>** init_ext_fns(Hermes::vector<MeshFunctionSharedPtr<Scalar> > ext,
        LightArray<NeighborSearch<Scalar>*>& neighbor_searches,
        int order, unsigned int min_dg_mesh_seq);

      /// Initialize neighbors.
      bool init_neighbors(LightArray<NeighborSearch<Scalar>*>& neighbor_searches, Traverse::State* current_state, unsigned int min_dg_mesh_seq);

      /// Initialize the tree for traversing multimesh neighbors.
      void build_multimesh_tree(MultimeshDGNeighborTreeNode* root, LightArray<NeighborSearch<Scalar>*>& neighbor_searches);

      /// Recursive insertion function into the tree.
      void insert_into_multimesh_tree(MultimeshDGNeighborTreeNode* node, unsigned int* transformations, unsigned int transformation_count);

      /// Return a global (unified list of central element transformations representing the neighbors on the union mesh.
      Hermes::vector<Hermes::vector<unsigned int>*> get_multimesh_neighbors_transformations(MultimeshDGNeighborTreeNode* multimesh_tree);

      /// Traverse the multimesh tree. Used in the function get_multimesh_neighbors_transformations().
      void traverse_multimesh_tree(MultimeshDGNeighborTreeNode* node, Hermes::vector<Hermes::vector<unsigned int>*>& running_transformations);

      /// Update the NeighborSearch according to the multimesh tree.
      void update_neighbor_search(NeighborSearch<Scalar>* ns, MultimeshDGNeighborTreeNode* multimesh_tree);

      /// Finds a node in the multimesh tree that corresponds to the array transformations, with the length of transformation_count,
      /// starting to look for it in the MultimeshDGNeighborTreeNode node.
      MultimeshDGNeighborTreeNode* find_node(unsigned int* transformations, unsigned int transformation_count, MultimeshDGNeighborTreeNode* node);

      /// Updates the NeighborSearch ns according to the subtree of MultimeshDGNeighborTreeNode node.
      /// Returns 0 if no neighbor was deleted, -1 otherwise.
      int update_ns_subtree(NeighborSearch<Scalar>* ns, MultimeshDGNeighborTreeNode* node, unsigned int ith_neighbor);

      /// Traverse the multimesh subtree. Used in the function update_ns_subtree().
      void traverse_multimesh_subtree(MultimeshDGNeighborTreeNode* node, Hermes::vector<Hermes::vector<unsigned int>*>& running_central_transformations,
        Hermes::vector<Hermes::vector<unsigned int>*>& running_neighbor_transformations, const typename NeighborSearch<Scalar>::NeighborEdgeInfo& edge_info, const int& active_edge, const int& mode);

      /// Space instances for all equations in the system.
      Hermes::vector<SpaceSharedPtr<Scalar> > spaces;
      int spaces_size;

      template<typename T> friend class DiscreteProblem;
    };
  }
}
#endif
