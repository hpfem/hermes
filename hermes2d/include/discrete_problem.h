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

#ifndef __H2D_DISCRETE_PROBLEM_H
#define __H2D_DISCRETE_PROBLEM_H

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

namespace Hermes
{
  namespace Hermes2D
  {
    class PrecalcShapeset;

    /// @ingroup inner
    /// Multimesh neighbors traversal class.
    class NeighborNode
    {
    private:
      NeighborNode(NeighborNode* parent, unsigned int transformation);
      ~NeighborNode();
      void set_left_son(NeighborNode* left_son);
      void set_right_son(NeighborNode* right_son);
      void set_transformation(unsigned int transformation);
      NeighborNode* get_left_son();
      NeighborNode* get_right_son();
      unsigned int get_transformation();
      NeighborNode* parent;
      NeighborNode* left_son;
      NeighborNode* right_son;
      unsigned int transformation;
      template<typename Scalar> friend class DiscreteProblem;
      template<typename Scalar> friend class KellyTypeAdapt;
    };

    /// @ingroup inner
    /// Discrete problem class.
    ///
    /// This class does assembling into external matrix / vector structures.
    ///
    template<typename Scalar>
    class HERMES_API DiscreteProblem : public DiscreteProblemInterface<Scalar>, public Hermes::Mixins::TimeMeasurable, public Hermes::Hermes2D::Mixins::SettableSpaces<Scalar>
    {
    public:
      /// Constructor for multiple components / equations.
      DiscreteProblem(const WeakForm<Scalar>* wf, Hermes::vector<const Space<Scalar> *> spaces);

      /// Constructor for one equation.
      DiscreteProblem(const WeakForm<Scalar>* wf, const Space<Scalar>* space);
      
      /// Set this problem to Finite Volume (no integration order calculation).
      void set_fvm();

      /// Sets new spaces for the instance.
      virtual void set_spaces(Hermes::vector<const Space<Scalar>*> spaces);
      virtual void set_space(const Space<Scalar>* space);
      
      /// Non-parameterized constructor.
      DiscreteProblem();

      /// Destuctor.
      virtual ~DiscreteProblem();

      /// If the cache should not be used for any reason.
      inline void setDoNotUseCache() { this->doNotUseCache = true; }

      // GET functions.
      /// Get the weak forms.
      const WeakForm<Scalar>* get_weak_formulation() const;

      /// Get all spaces as a Hermes::vector.
      virtual Hermes::vector<const Space<Scalar>*> get_spaces() const;

       /// Get the number of unknowns.
      int get_num_dofs() const;

      /// Get info about presence of a matrix.
      bool is_matrix_free() const;
      
      /// set time information for time-dependent problems.
      virtual void setTime(double time);
      virtual void setTimeStep(double timeStep);

      void delete_cache();

      /// Assembling.
      /// General assembling procedure for nonlinear problems. coeff_vec is the
      /// previous Newton vector. If force_diagonal_block == true, then (zero) matrix
      /// antries are created in diagonal blocks even if corresponding matrix weak
      /// forms do not exist. This is useful if the matrix is later to be merged with
      /// a matrix that has nonzeros in these blocks. The Table serves for optional
      /// weighting of matrix blocks in systems. The parameter add_dir_lift decides
      /// whether Dirichlet lift will be added while coeff_vec is converted into
      /// Solutions.
      void assemble(Scalar* coeff_vec, SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs = NULL,
        bool force_diagonal_blocks = false, Table* block_weights = NULL);

      /// Assembling.
      /// Without the matrix.
      void assemble(Scalar* coeff_vec, Vector<Scalar>* rhs = NULL,
        bool force_diagonal_blocks = false, Table* block_weights = NULL);

      /// Light version passing NULL for the coefficient vector. External solutions
      /// are initialized with zeros.
      virtual void assemble(SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs = NULL, bool force_diagonal_blocks = false,
        Table* block_weights = NULL);

      /// Light version passing NULL for the coefficient vector. External solutions
      /// are initialized with zeros.
      /// Without the matrix.
      void assemble(Vector<Scalar>* rhs = NULL, bool force_diagonal_blocks = false,
        Table* block_weights = NULL);

    protected:
      void init_assembling(Scalar* coeff_vec, PrecalcShapeset*** pss , PrecalcShapeset*** spss, RefMap*** refmaps, Solution<Scalar>*** u_ext, AsmList<Scalar>*** als, Hermes::vector<MeshFunction<Scalar>*>& ext_functions, MeshFunction<Scalar>*** ext,
          Hermes::vector<MatrixFormVol<Scalar>*>* mfvol, Hermes::vector<MatrixFormSurf<Scalar>*>* mfsurf, Hermes::vector<VectorFormVol<Scalar>*>* vfvol, Hermes::vector<VectorFormSurf<Scalar>*>* vfsurf);

      void deinit_assembling(PrecalcShapeset*** pss , PrecalcShapeset*** spss, RefMap*** refmaps, Solution<Scalar>*** u_ext, AsmList<Scalar>*** als, Hermes::vector<MeshFunction<Scalar>*>& ext_functions, MeshFunction<Scalar>*** ext,
          Hermes::vector<MatrixFormVol<Scalar>*>* mfvol, Hermes::vector<MatrixFormSurf<Scalar>*>* mfsurf, Hermes::vector<VectorFormVol<Scalar>*>* vfvol, Hermes::vector<VectorFormSurf<Scalar>*>* vfsurf);

      /// The form will be assembled.
      bool form_to_be_assembled(MatrixForm<Scalar>* form, Traverse::State* current_state);
      bool form_to_be_assembled(MatrixFormVol<Scalar>* form, Traverse::State* current_state);
      bool form_to_be_assembled(MatrixFormSurf<Scalar>* form, Traverse::State* current_state);
      bool form_to_be_assembled(VectorForm<Scalar>* form, Traverse::State* current_state);
      bool form_to_be_assembled(VectorFormVol<Scalar>* form, Traverse::State* current_state);
      bool form_to_be_assembled(VectorFormSurf<Scalar>* form, Traverse::State* current_state);

      // Return scaling coefficient.
      double block_scaling_coeff(MatrixForm<Scalar>* form) const;

      /// Preassembling.
      /// Precalculate matrix sparse structure.
      /// If force_diagonal_block == true, then (zero) matrix
      /// antries are created in diagonal blocks even if corresponding matrix weak
      /// forms do not exist. This is useful if the matrix is later to be merged with
      /// a matrix that has nonzeros in these blocks. The Table serves for optional
      /// weighting of matrix blocks in systems.
      void create_sparse_structure();
      void create_sparse_structure(SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs = NULL);

      /// Initialize a state, returns a non-NULL Element.
      void init_state(PrecalcShapeset** current_pss, PrecalcShapeset** current_spss, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, AsmList<Scalar>** current_als, Traverse::State* current_state);
      void init_surface_state(AsmList<Scalar>** current_als, Traverse::State* current_state);

      /// Set the special handling of external functions of Runge-Kutta methods, including information how many spaces were there in the original problem.
      inline void set_RK(int original_spaces_count) { this->RungeKutta = true; RK_original_spaces_count = original_spaces_count; }

      /// Assemble one state.
      void assemble_one_state(PrecalcShapeset** current_pss, PrecalcShapeset** current_spss, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, AsmList<Scalar>** current_als, Traverse::State* current_state,
           MatrixFormVol<Scalar>** current_mfvol, MatrixFormSurf<Scalar>** current_mfsurf, VectorFormVol<Scalar>** current_vfvol, VectorFormSurf<Scalar>** current_vfsurf);

      /// Adjusts order to refmaps.
      void adjust_order_to_refmaps(Form<Scalar> *form, int& order, Hermes::Ord* o, RefMap** current_refmaps);

      /// Matrix volumetric forms - calculate the integration order.
      int calc_order_matrix_form(MatrixForm<Scalar>* mfv, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, Traverse::State* current_state);

      /// Matrix volumetric forms - assemble the form.
      virtual void assemble_matrix_form(MatrixForm<Scalar>* form, int order, Func<double>** base_fns, Func<double>** test_fns, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, AsmList<Scalar>** current_als, Traverse::State* current_state);

      virtual void assemble_matrix_form(MatrixForm<Scalar>* form, int order, Func<double>** base_fns, Func<double>** test_fns, Solution<Scalar>** current_u_ext, 
      AsmList<Scalar>* current_als_i, AsmList<Scalar>* current_als_j, Traverse::State* current_state, int n_quadrature_points, Geom<double>* geometry, double* jacobian_x_weights);

      /// Vector volumetric forms - calculate the integration order.
      int calc_order_vector_form(VectorForm<Scalar>* mfv, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, Traverse::State* current_state);

      /// Vector volumetric forms - assemble the form.
      void assemble_vector_form(VectorForm<Scalar>* form, int order, Func<double>** test_fns, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, AsmList<Scalar>** current_als, Traverse::State* current_state);

      void assemble_vector_form(VectorForm<Scalar>* form, int order, Func<double>** test_fns, Solution<Scalar>** current_u_ext, 
      AsmList<Scalar>* current_als, Traverse::State* current_state, int n_quadrature_points, Geom<double>* geometry, double* jacobian_x_weights);

      /// \ingroup Helper methods inside {calc_order_*, assemble_*}
      /// Init geometry, jacobian * weights, return the number of integration points.
      int init_geometry_points(RefMap* reference_mapping, int order, Geom<double>*& geometry, double*& jacobian_x_weights);
      int init_surface_geometry_points(RefMap* reference_mapping, int& order, Traverse::State* current_state, Geom<double>*& geometry, double*& jacobian_x_weights);

      /// \ingroup Helper methods inside {calc_order_*, assemble_*}
      /// Calculates orders for external functions.
      void init_ext_orders(Form<Scalar> *form, Func<Hermes::Ord>** oi, ExtData<Hermes::Ord>* oext, Solution<Scalar>** current_u_ext, Traverse::State* current_state);
      /// \ingroup Helper methods inside {calc_order_*, assemble_*}
      /// Cleans up after init_ext_orders.
      void deinit_ext_orders(Form<Scalar> *form, Func<Hermes::Ord>** oi, ExtData<Hermes::Ord>* oext);

      /// \ingroup Helper methods inside {calc_order_*, assemble_*}
      /// Calculates external functions.
      void init_ext(Form<Scalar> *form, Func<Scalar>** u_ext, ExtData<Scalar>* ext, int order, Solution<Scalar>** current_u_ext, Traverse::State* current_state);
      /// \ingroup Helper methods inside {calc_order_*, assemble_*}
      /// Cleans up after init_ext.
      void deinit_ext(Form<Scalar> *form, Func<Scalar>** u_ext, ExtData<Scalar>* ext);

      /// Init function. Common code for the constructors.
      void init();

      /// Matrix structure as well as spaces and weak formulation is up-to-date.
      bool is_up_to_date() const;

      /// Weak formulation.
      const WeakForm<Scalar>* wf;

      /// Seq number of the WeakForm.
      int wf_seq;

      /// Space instances for all equations in the system.
      Hermes::vector<const Space<Scalar>*> spaces;

      Hermes::vector<unsigned int> spaces_first_dofs;

      /// Seq numbers of Space instances in spaces.
      int* sp_seq;

      /// Number of DOFs of all Space instances in spaces.
      int ndof;

      /// Instance of the class Geom used in the calculation of integration order.
      Geom<Hermes::Ord> geom_ord;

      /// Fake weight used in the calculation of integration order.
      static double fake_wt;

      /// If the problem has only constant test functions, there is no need for order calculation,
      /// which saves time.
      bool is_fvm;

      /// Matrix structure can be reused.
      /// If other conditions apply.
      bool have_matrix;

      /// There is a matrix form set on DG_INNER_EDGE area or not.
      bool DG_matrix_forms_present;

      /// There is a vector form set on DG_INNER_EDGE area or not.
      bool DG_vector_forms_present;

      /// Turn on Runge-Kutta specific handling of external functions.
      bool RungeKutta;

      /// Number of spaces in the original problem in a Runge-Kutta method.
      int RK_original_spaces_count;

      /// Storing assembling info.
      SparseMatrix<Scalar>* current_mat;
      Vector<Scalar>* current_rhs;
      bool current_force_diagonal_blocks;
      Table* current_block_weights;

      /// Caching.
      class CacheRecordPerElement
      {
      public:
        void clear();
        int* asmlistIdx;
        int asmlistCnt;
      };

      class CacheRecordPerSubIdx
      {
      public:
        int order;
        int nvert;
        void clear();
        int asmlistCnt;
        Func<double>** fns;
        Func<double>*** fnsSurface;
        Geom<double>* geometry;
        Geom<double>** geometrySurface;
        double* jacobian_x_weights;
        double** jacobian_x_weightsSurface;
        int n_quadrature_points;
        int* n_quadrature_pointsSurface;
        int* orderSurface;
        int* asmlistSurfaceCnt;
      };

      std::map<uint64_t, CacheRecordPerSubIdx*>*** cache_records_sub_idx;
      CacheRecordPerElement*** cache_records_element;
      bool** cache_element_stored;
      int cache_size;
      bool doNotUseCache;

      /// Exception caught in a parallel region.
      Hermes::Exceptions::Exception* caughtException;
    
      
      ///* DG *///

      /// Assemble DG forms.
      void assemble_one_DG_state(PrecalcShapeset** current_pss, PrecalcShapeset** current_spss, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, AsmList<Scalar>** current_als,
        Traverse::State* current_state, MatrixFormSurf<Scalar>** current_mfsurf, VectorFormSurf<Scalar>** current_vfsurf, Transformable** fn);

      /// Assemble one DG neighbor.
      void assemble_DG_one_neighbor(bool edge_processed, unsigned int neighbor_i,
        PrecalcShapeset** current_pss, PrecalcShapeset** current_spss, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, AsmList<Scalar>** current_als,
        Traverse::State* current_state, MatrixFormSurf<Scalar>** current_mfsurf, VectorFormSurf<Scalar>** current_vfsurf, Transformable** fn,
        std::map<unsigned int, PrecalcShapeset *> npss, std::map<unsigned int, PrecalcShapeset *> nspss, std::map<unsigned int, RefMap *> nrefmap,
        LightArray<NeighborSearch<Scalar>*>& neighbor_searches, unsigned int min_dg_mesh_seq);

      /// Assemble DG matrix forms.
      void assemble_DG_matrix_forms(PrecalcShapeset** current_pss, PrecalcShapeset** current_spss, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, AsmList<Scalar>** current_als,
        Traverse::State* current_state, MatrixFormSurf<Scalar>** current_mfsurf, std::map<unsigned int, PrecalcShapeset*> npss,
        std::map<unsigned int, PrecalcShapeset*> nspss, std::map<unsigned int, RefMap*> nrefmap, LightArray<NeighborSearch<Scalar>*>& neighbor_searches);

      /// Assemble DG vector forms.
      void assemble_DG_vector_forms(PrecalcShapeset** current_spss, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, AsmList<Scalar>** current_als,
        Traverse::State* current_state, VectorFormSurf<Scalar>** current_vfsurf, std::map<unsigned int, PrecalcShapeset*> nspss,
        std::map<unsigned int, RefMap*> nrefmap, LightArray<NeighborSearch<Scalar>*>& neighbor_searches);

      DiscontinuousFunc<Hermes::Ord>* init_ext_fn_ord(NeighborSearch<Scalar>* ns, MeshFunction<Scalar>* fu);

      /// Calculates integration order for DG matrix forms.
      int calc_order_dg_matrix_form(MatrixFormSurf<Scalar>* mfs, Hermes::vector<Solution<Scalar>*> u_ext,
        PrecalcShapeset* fu, PrecalcShapeset* fv, RefMap* ru, SurfPos* surf_pos,
        bool neighbor_supp_u, bool neighbor_supp_v, LightArray<NeighborSearch<Scalar>*>& neighbor_searches, int neighbor_index_u);

      /// Evaluates DG matrix forms on an edge between elements identified by ru_actual, rv.
      Scalar eval_dg_form(MatrixFormSurf<Scalar>* mfs, Hermes::vector<Solution<Scalar>*> u_ext,
        PrecalcShapeset* fu, PrecalcShapeset* fv, RefMap* ru_central, RefMap* ru_actual, RefMap* rv,
        bool neighbor_supp_u, bool neighbor_supp_v,
        SurfPos* surf_pos, LightArray<NeighborSearch<Scalar>*>& neighbor_searches, int neighbor_index_u, int neighbor_index_v);

      /// Calculates integration order for DG vector forms.
      int calc_order_dg_vector_form(VectorFormSurf<Scalar>* vfs, Hermes::vector<Solution<Scalar>*> u_ext,
        PrecalcShapeset* fv, RefMap* ru, SurfPos* surf_pos,
        LightArray<NeighborSearch<Scalar>*>& neighbor_searches, int neighbor_index_v);

      /// Evaluates DG vector forms on an edge between elements identified by ru_actual, rv.
      Scalar eval_dg_form(VectorFormSurf<Scalar>* vfs, Hermes::vector<Solution<Scalar>*> u_ext,
        PrecalcShapeset* fv, RefMap* rv,
        SurfPos* surf_pos, LightArray<NeighborSearch<Scalar>*>& neighbor_searches, int neighbor_index_v);

      /// Initialize orders of external functions for DG forms.
      ExtData<Hermes::Ord>* init_ext_fns_ord(Hermes::vector<MeshFunction<Scalar>*> &ext,
        LightArray<NeighborSearch<Scalar>*>& neighbor_searches);

      /// Initialize external functions for DG forms.
      ExtData<Scalar>* init_ext_fns(Hermes::vector<MeshFunction<Scalar>*> &ext,
        LightArray<NeighborSearch<Scalar>*>& neighbor_searches,
        int order, unsigned int min_dg_mesh_seq);

      /// Initialize neighbors.
      bool init_neighbors(LightArray<NeighborSearch<Scalar>*>& neighbor_searches, Traverse::State* current_state, unsigned int min_dg_mesh_seq);

      /// Initialize the tree for traversing multimesh neighbors.
      void build_multimesh_tree(NeighborNode* root, LightArray<NeighborSearch<Scalar>*>& neighbor_searches);

      /// Recursive insertion function into the tree.
      void insert_into_multimesh_tree(NeighborNode* node, unsigned int* transformations, unsigned int transformation_count);

      /// Return a global (unified list of central element transformations representing the neighbors on the union mesh.
      Hermes::vector<Hermes::vector<unsigned int>*> get_multimesh_neighbors_transformations(NeighborNode* multimesh_tree);

      /// Traverse the multimesh tree. Used in the function get_multimesh_neighbors_transformations().
      void traverse_multimesh_tree(NeighborNode* node, Hermes::vector<Hermes::vector<unsigned int>*>& running_transformations);

      /// Update the NeighborSearch according to the multimesh tree.
      void update_neighbor_search(NeighborSearch<Scalar>* ns, NeighborNode* multimesh_tree);

      /// Finds a node in the multimesh tree that corresponds to the array transformations, with the length of transformation_count,
      /// starting to look for it in the NeighborNode node.
      NeighborNode* find_node(unsigned int* transformations, unsigned int transformation_count, NeighborNode* node);

      /// Updates the NeighborSearch ns according to the subtree of NeighborNode node.
      /// Returns 0 if no neighbor was deleted, -1 otherwise.
      int update_ns_subtree(NeighborSearch<Scalar>* ns, NeighborNode* node, unsigned int ith_neighbor);

      /// Traverse the multimesh subtree. Used in the function update_ns_subtree().
      void traverse_multimesh_subtree(NeighborNode* node, Hermes::vector<Hermes::vector<unsigned int>*>& running_central_transformations,
        Hermes::vector<Hermes::vector<unsigned int>*>& running_neighbor_transformations, const typename NeighborSearch<Scalar>::NeighborEdgeInfo& edge_info, const int& active_edge, const int& mode);

      template<typename T> friend class KellyTypeAdapt;
      template<typename T> friend class NewtonSolver;
      template<typename T> friend class PicardSolver;
      template<typename T> friend class RungeKutta;
    };
  }
}
#endif