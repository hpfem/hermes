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

namespace Hermes
{
  namespace Hermes2D
  {
    class PrecalcShapeset;

    /// Multimesh neighbors traversal class.
    class NeighborNode
    {
    public:
      NeighborNode(NeighborNode* parent, unsigned int transformation);
      ~NeighborNode();
      void set_left_son(NeighborNode* left_son);
      void set_right_son(NeighborNode* right_son);
      void set_transformation(unsigned int transformation);
      NeighborNode* get_left_son();
      NeighborNode* get_right_son();
      unsigned int get_transformation();
    private:
      NeighborNode* parent;
      NeighborNode* left_son;
      NeighborNode* right_son;
      unsigned int transformation;
    };

    /// Discrete problem class.
    ///
    /// This class does assembling into external matrix / vector structures.
    ///
    template<typename Scalar>
    class HERMES_API DiscreteProblem : public DiscreteProblemInterface<Scalar>
    {
    public:
      /// Constructor for multiple components / equations.
      DiscreteProblem(WeakForm<Scalar>* wf, Hermes::vector<Space<Scalar>*> spaces);

      /// Constructor for one equation.
      DiscreteProblem(WeakForm<Scalar>* wf, Space<Scalar>* space);

      /// Non-parameterized constructor (currently used only in KellyTypeAdapt to gain access to NeighborSearch methods).
      DiscreteProblem();

      /// Init function. Common code for the constructors.
      void init();

      /// Destuctor.
      virtual ~DiscreteProblem();
      void free();

      // GET functions.
      /// Get pointer to n-th space.
      Space<Scalar>* get_space(int n);

      /// Get the weak forms.
      WeakForm<Scalar>* get_weak_formulation();

      /// Get all spaces as a Hermes::vector.
      Hermes::vector<Space<Scalar>*> get_spaces();

      /// This is different from H3D.
      PrecalcShapeset* get_pss(int n);

      /// Get the number of unknowns.
      int get_num_dofs();

      /// Get info about presence of a matrix.
      bool is_matrix_free() ;

      /// Sets the storing of the previous matrix in adaptivity.
      void set_adaptivity_cache();
      void temp_enable_adaptivity_cache();
      void temp_disable_adaptivity_cache();

      /// Preassembling.
      /// Precalculate matrix sparse structure.
      /// If force_diagonal_block == true, then (zero) matrix
      /// antries are created in diagonal blocks even if corresponding matrix weak
      /// forms do not exist. This is useful if the matrix is later to be merged with
      /// a matrix that has nonzeros in these blocks. The Table serves for optional
      /// weighting of matrix blocks in systems.
      void create_sparse_structure(SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs = NULL,
        bool force_diagonal_blocks = false, Table* block_weights = NULL);

      /// Assembling utilities.
      /// Check whether it is sane to assemble.
      /// Throws errors if not.
      void assemble_sanity_checks(Table* block_weights);

      /// Converts coeff_vec to u_ext.
      /// If the supplied coeff_vec is NULL, it supplies the appropriate number of NULL pointers.
      void convert_coeff_vec(Scalar* coeff_vec, Hermes::vector<Solution<Scalar>*> & u_ext, bool add_dir_lift);

      /// Initializes psss.
      void initialize_psss(Hermes::vector<PrecalcShapeset*>& spss);

      /// Initializes refmaps.
      void initialize_refmaps(Hermes::vector<RefMap*>& refmap);

      /// Initialize a state, returns a non-NULL Element.
      Element* init_state(Stage<Scalar>& stage, Hermes::vector<PrecalcShapeset*>& spss, 
        Hermes::vector<RefMap*>& refmap, Element** e, Hermes::vector<AsmList<Scalar>*>& al);

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
        bool force_diagonal_blocks = false, bool add_dir_lift = true, Table* block_weights = NULL);

      /// Assembling.
      /// Without the matrix.
      void assemble(Scalar* coeff_vec, Vector<Scalar>* rhs = NULL,
        bool force_diagonal_blocks = false, bool add_dir_lift = true, Table* block_weights = NULL);

      /// Light version passing NULL for the coefficient vector. External solutions 
      /// are initialized with zeros.
      void assemble(SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs = NULL, bool force_diagonal_blocks = false, 
        Table* block_weights = NULL);

      /// Light version passing NULL for the coefficient vector. External solutions 
      /// are initialized with zeros.
      /// Without the matrix.
      void assemble(Vector<Scalar>* rhs = NULL, bool force_diagonal_blocks = false, 
        Table* block_weights = NULL);

      /// Assemble one stage.
      void assemble_one_stage(Stage<Scalar>& stage, 
        SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs, bool force_diagonal_blocks, Table* block_weights,
        Hermes::vector<PrecalcShapeset*>& spss, Hermes::vector<RefMap*>& refmap, Hermes::vector<Solution<Scalar>*>& u_ext);

      /// Assemble one state.
      void assemble_one_state(Stage<Scalar>& stage, 
        SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs, bool force_diagonal_blocks, Table* block_weights,
        Hermes::vector<PrecalcShapeset*>& spss, Hermes::vector<RefMap*>& refmap, Hermes::vector<Solution<Scalar>*>& u_ext, Element** e, 
        bool* bnd, SurfPos* surf_pos, Element* trav_base);

      /// Assemble volume matrix forms.
      void assemble_volume_matrix_forms(Stage<Scalar>& stage, 
        SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs, bool force_diagonal_blocks, Table* block_weights,
        Hermes::vector<PrecalcShapeset*>& spss, Hermes::vector<RefMap*>& refmap, Hermes::vector<Solution<Scalar>*>& u_ext, 
        int marker, Hermes::vector<AsmList<Scalar>*>& al);
      /// Assemble multicomponent volume matrix forms.
      void assemble_multicomponent_volume_matrix_forms(Stage<Scalar>& stage, 
        SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs, bool force_diagonal_blocks, Table* block_weights,
        Hermes::vector<PrecalcShapeset*>& spss, Hermes::vector<RefMap*>& refmap, Hermes::vector<Solution<Scalar>*>& u_ext, 
        int marker, Hermes::vector<AsmList<Scalar>*>& al);

      /// Assemble volume vector forms.
      void assemble_volume_vector_forms(Stage<Scalar>& stage, 
        SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs, bool force_diagonal_blocks, Table* block_weights,
        Hermes::vector<PrecalcShapeset*>& spss, Hermes::vector<RefMap*>& refmap, Hermes::vector<Solution<Scalar>*>& u_ext, 
        int marker, Hermes::vector<AsmList<Scalar>*>& al);
      /// Assemble multicomponent volume vector forms.
      void assemble_multicomponent_volume_vector_forms(Stage<Scalar>& stage, 
        SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs, bool force_diagonal_blocks, Table* block_weights,
        Hermes::vector<PrecalcShapeset*>& spss, Hermes::vector<RefMap*>& refmap, Hermes::vector<Solution<Scalar>*>& u_ext, 
        int marker, Hermes::vector<AsmList<Scalar>*>& al);

      /// Assemble surface and DG forms.
      void assemble_surface_integrals(Stage<Scalar>& stage, 
        SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs, bool force_diagonal_blocks, Table* block_weights,
        Hermes::vector<PrecalcShapeset*>& spss, Hermes::vector<RefMap*>& refmap, Hermes::vector<Solution<Scalar>*>& u_ext, 
        int marker, Hermes::vector<AsmList<Scalar>*>& al, bool bnd, SurfPos& surf_pos, Hermes::vector<bool>& nat, 
        int isurf, Element** e, Element* trav_base, Element* rep_element);

      /// Assemble surface matrix forms.
      void assemble_surface_matrix_forms(Stage<Scalar>& stage, 
        SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs, bool force_diagonal_blocks, Table* block_weights,
        Hermes::vector<PrecalcShapeset*>& spss, Hermes::vector<RefMap*>& refmap, Hermes::vector<Solution<Scalar>*>& u_ext, 
        int marker, Hermes::vector<AsmList<Scalar>*>& al, bool bnd, SurfPos& surf_pos, Hermes::vector<bool>& nat, 
        int isurf, Element** e, Element* trav_base);
      /// Assemble multicomponent surface matrix forms.
      void assemble_multicomponent_surface_matrix_forms(Stage<Scalar>& stage, 
        SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs, bool force_diagonal_blocks, Table* block_weights,
        Hermes::vector<PrecalcShapeset*>& spss, Hermes::vector<RefMap*>& refmap, Hermes::vector<Solution<Scalar>*>& u_ext, 
        int marker, Hermes::vector<AsmList<Scalar>*>& al, bool bnd, SurfPos& surf_pos, Hermes::vector<bool>& nat, 
        int isurf, Element** e, Element* trav_base);

      /// Assemble surface vector forms.
      void assemble_surface_vector_forms(Stage<Scalar>& stage, 
        SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs, bool force_diagonal_blocks, Table* block_weights,
        Hermes::vector<PrecalcShapeset*>& spss, Hermes::vector<RefMap*>& refmap, Hermes::vector<Solution<Scalar>*>& u_ext, 
        int marker, Hermes::vector<AsmList<Scalar>*>& al, bool bnd, SurfPos& surf_pos, Hermes::vector<bool>& nat, 
        int isurf, Element** e, Element* trav_base);
      /// Assemble multicomponent surface vector forms.
      void assemble_multicomponent_surface_vector_forms(Stage<Scalar>& stage, 
        SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs, bool force_diagonal_blocks, Table* block_weights,
        Hermes::vector<PrecalcShapeset*>& spss, Hermes::vector<RefMap*>& refmap, Hermes::vector<Solution<Scalar>*>& u_ext, 
        int marker, Hermes::vector<AsmList<Scalar>*>& al, bool bnd, SurfPos& surf_pos, Hermes::vector<bool>& nat, 
        int isurf, Element** e, Element* trav_base);

      /// Assemble DG forms.
      void assemble_DG_forms(Stage<Scalar>& stage, 
        SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs, bool force_diagonal_blocks, Table* block_weights,
        Hermes::vector<PrecalcShapeset*>& spss, Hermes::vector<RefMap*>& refmap, Hermes::vector<Solution<Scalar>*>& u_ext, 
        int marker, Hermes::vector<AsmList<Scalar>*>& al, bool bnd, SurfPos& surf_pos, Hermes::vector<bool>& nat, 
        int isurf, Element** e, Element* trav_base, Element* rep_element);

      /// Assemble one DG neighbor.
      void assemble_DG_one_neighbor(bool edge_processed, unsigned int neighbor_i, Stage<Scalar>& stage, 
        SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs, bool force_diagonal_blocks, Table* block_weights,
        Hermes::vector<PrecalcShapeset *>& spss, Hermes::vector<RefMap *>& refmap, std::map<unsigned int, PrecalcShapeset *> npss,
        std::map<unsigned int, PrecalcShapeset *> nspss, std::map<unsigned int, RefMap *> nrefmap, LightArray<NeighborSearch<Scalar>*>& neighbor_searches, Hermes::vector<Solution<Scalar>*>& u_ext, 
        int marker, Hermes::vector<AsmList<Scalar>*>& al, bool bnd, SurfPos& surf_pos, Hermes::vector<bool>& nat, 
        int isurf, Element** e, Element* trav_base, Element* rep_element);

      /// Assemble DG matrix forms.
      void assemble_DG_matrix_forms(Stage<Scalar>& stage, 
        SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs, bool force_diagonal_blocks, Table* block_weights,
        Hermes::vector<PrecalcShapeset*>& spss, Hermes::vector<RefMap*>& refmap, std::map<unsigned int, PrecalcShapeset*> npss,
        std::map<unsigned int, PrecalcShapeset*> nspss, std::map<unsigned int, RefMap*> nrefmap, LightArray<NeighborSearch<Scalar>*>& neighbor_searches, Hermes::vector<Solution<Scalar>*>& u_ext, 
        int marker, Hermes::vector<AsmList<Scalar>*>& al, bool bnd, SurfPos& surf_pos, Hermes::vector<bool>& nat, 
        int isurf, Element** e, Element* trav_base, Element* rep_element);
      /// Assemble multicomponent DG matrix forms.
      void assemble_multicomponent_DG_matrix_forms(Stage<Scalar>& stage, 
        SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs, bool force_diagonal_blocks, Table* block_weights,
        Hermes::vector<PrecalcShapeset*>& spss, Hermes::vector<RefMap*>& refmap, std::map<unsigned int, PrecalcShapeset*> npss,
        std::map<unsigned int, PrecalcShapeset*> nspss, std::map<unsigned int, RefMap*> nrefmap, LightArray<NeighborSearch<Scalar>*>& neighbor_searches, Hermes::vector<Solution<Scalar>*>& u_ext, 
        int marker, Hermes::vector<AsmList<Scalar>*>& al, bool bnd, SurfPos& surf_pos, Hermes::vector<bool>& nat, 
        int isurf, Element** e, Element* trav_base, Element* rep_element);

      /// Assemble DG vector forms.
      void assemble_DG_vector_forms(Stage<Scalar>& stage, 
        SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs, bool force_diagonal_blocks, Table* block_weights,
        Hermes::vector<PrecalcShapeset*>& spss, Hermes::vector<RefMap*>& refmap, LightArray<NeighborSearch<Scalar>*>& neighbor_searches, Hermes::vector<Solution<Scalar>*>& u_ext, 
        int marker, Hermes::vector<AsmList<Scalar>*>& al, bool bnd, SurfPos& surf_pos, Hermes::vector<bool>& nat, 
        int isurf, Element** e, Element* trav_base, Element* rep_element);
      /// Assemble multicomponent DG vector forms.
      void assemble_multicomponent_DG_vector_forms(Stage<Scalar>& stage, 
        SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs, bool force_diagonal_blocks, Table* block_weights,
        Hermes::vector<PrecalcShapeset*>& spss, Hermes::vector<RefMap*>& refmap, LightArray<NeighborSearch<Scalar>*>& neighbor_searches, Hermes::vector<Solution<Scalar>*>& u_ext, 
        int marker, Hermes::vector<AsmList<Scalar>*>& al, bool bnd, SurfPos& surf_pos, Hermes::vector<bool>& nat, 
        int isurf, Element** e, Element* trav_base, Element* rep_element);

      void invalidate_matrix();

      /// Set this problem to Finite Volume.
      void set_fvm();

      void set_spaces(Hermes::vector<Space<Scalar>*> spaces);

      void set_spaces(Space<Scalar>* space);

      /// Set the special handling of external functions of Runge-Kutta methods, including information how many spaces were there in the original problem.
      void set_RK(int original_spaces_count) { this->RungeKutta = true; RK_original_spaces_count = original_spaces_count; }

    protected:
      bool cache_for_adaptivity;
      bool temp_cache_for_adaptivity;

      DiscontinuousFunc<Hermes::Ord>* init_ext_fn_ord(NeighborSearch<Scalar>* ns, MeshFunction<Scalar>* fu);

      /// Main function for the evaluation of volumetric matrix forms. 
      /// Evaluates weak form on element given by the RefMap.
      /// Calls the function calc_order_matrix_form_vol to get the integration order.
      Scalar eval_form(MatrixFormVol<Scalar>* mfv, Hermes::vector<Solution<Scalar>*> u_ext,
        PrecalcShapeset* fu, PrecalcShapeset* fv, RefMap* ru, RefMap* rv);
      /// Main function for the evaluation of multicomponent volumetric matrix forms. 
      void eval_form(MultiComponentMatrixFormVol<Scalar>* mfv, Hermes::vector<Solution<Scalar>*> u_ext,
        PrecalcShapeset* fu, PrecalcShapeset* fv, RefMap* ru, RefMap* rv, Hermes::vector<Scalar>& result);

      /// Calculates the necessary integration order to use for a particular volumetric matrix form.
      int calc_order_matrix_form_vol(MatrixFormVol<Scalar>* mfv, Hermes::vector<Solution<Scalar>*> u_ext,
        PrecalcShapeset* fu, PrecalcShapeset* fv, RefMap* ru, RefMap* rv);
      /// Calculates the necessary integration order to use for a particular multicomponent volumetric matrix form.
      int calc_order_matrix_form_vol(MultiComponentMatrixFormVol<Scalar>* mfv, Hermes::vector<Solution<Scalar>*> u_ext,
        PrecalcShapeset* fu, PrecalcShapeset* fv, RefMap* ru, RefMap* rv);

      /// Elementary function used in eval_form() in adaptive mode for volumetric matrix forms.
      Scalar eval_form_subelement(int order, MatrixFormVol<Scalar>* mfv, 
        Hermes::vector<Solution<Scalar>*> u_ext,
        PrecalcShapeset* fu, PrecalcShapeset* fv, 
        RefMap* ru, RefMap* rv);

      /// Vector<Scalar> volume forms. The functions provide the same functionality as the
      /// parallel ones for matrix volume forms.
      Scalar eval_form(VectorFormVol<Scalar>* vfv, Hermes::vector<Solution<Scalar>*> u_ext,
        PrecalcShapeset* fv, RefMap* rv);
      void eval_form(MultiComponentVectorFormVol<Scalar>* vfv, Hermes::vector<Solution<Scalar>*> u_ext,
        PrecalcShapeset* fv, RefMap* rv, Hermes::vector<Scalar>& result);

      int calc_order_vector_form_vol(VectorFormVol<Scalar>* mfv, Hermes::vector<Solution<Scalar>*> u_ext,
        PrecalcShapeset* fv, RefMap* rv);
      int calc_order_vector_form_vol(MultiComponentVectorFormVol<Scalar>* mfv, Hermes::vector<Solution<Scalar>*> u_ext,
        PrecalcShapeset* fv, RefMap* rv);

      Scalar eval_form_subelement(int order, VectorFormVol<Scalar>* vfv, 
        Hermes::vector<Solution<Scalar>*> u_ext,
        PrecalcShapeset* fv, RefMap* rv);

      /// Matrix<Scalar> surface forms. The functions provide the same functionality as the
      /// parallel ones for matrix volume forms.
      Scalar eval_form(MatrixFormSurf<Scalar>* mfs, 
        Hermes::vector<Solution<Scalar>*> u_ext, 
        PrecalcShapeset* fu, PrecalcShapeset* fv, RefMap* ru, RefMap* rv, SurfPos* surf_pos);
      void eval_form(MultiComponentMatrixFormSurf<Scalar>* mfs, 
        Hermes::vector<Solution<Scalar>*> u_ext, 
        PrecalcShapeset* fu, PrecalcShapeset* fv, RefMap* ru, RefMap* rv, SurfPos* surf_pos, Hermes::vector<Scalar>& result);

      int calc_order_matrix_form_surf(MatrixFormSurf<Scalar>* mfs, 
        Hermes::vector<Solution<Scalar>*> u_ext,
        PrecalcShapeset* fu, PrecalcShapeset* fv, 
        RefMap* ru, RefMap* rv, SurfPos* surf_pos);
      int calc_order_matrix_form_surf(MultiComponentMatrixFormSurf<Scalar>* mfs, 
        Hermes::vector<Solution<Scalar>*> u_ext,
        PrecalcShapeset* fu, PrecalcShapeset* fv, 
        RefMap* ru, RefMap* rv, SurfPos* surf_pos);

      Scalar eval_form_subelement(int order, MatrixFormSurf<Scalar>* mfs, Hermes::vector<Solution<Scalar>*> u_ext,
        PrecalcShapeset* fu, PrecalcShapeset* fv, RefMap* ru, RefMap* rv, SurfPos* surf_pos);

      /// Vector<Scalar> surface forms. The functions provide the same functionality as the
      /// parallel ones for matrix volume forms.
      Scalar eval_form(VectorFormSurf<Scalar>* vfs, 
        Hermes::vector<Solution<Scalar>*> u_ext, 
        PrecalcShapeset* fv, RefMap* rv, SurfPos* surf_pos);
      void eval_form(MultiComponentVectorFormSurf<Scalar>* vfs, 
        Hermes::vector<Solution<Scalar>*> u_ext, 
        PrecalcShapeset* fv, RefMap* rv, SurfPos* surf_pos, Hermes::vector<Scalar>& result);

      int calc_order_vector_form_surf(VectorFormSurf<Scalar>* vfs, Hermes::vector<Solution<Scalar>*> u_ext,
        PrecalcShapeset* fv, RefMap* rv, SurfPos* surf_pos);
      int calc_order_vector_form_surf(MultiComponentVectorFormSurf<Scalar>* vfs, Hermes::vector<Solution<Scalar>*> u_ext,
        PrecalcShapeset* fv, RefMap* rv, SurfPos* surf_pos);

      Scalar eval_form_subelement(int order, VectorFormSurf<Scalar>* vfs, Hermes::vector<Solution<Scalar>*> u_ext,
        PrecalcShapeset* fv, RefMap* rv, SurfPos* surf_pos);

      /// Calculates integration order for DG matrix forms.
      int calc_order_dg_matrix_form(MatrixFormSurf<Scalar>* mfs, Hermes::vector<Solution<Scalar>*> u_ext,
        PrecalcShapeset* fu, PrecalcShapeset* fv, RefMap* ru, SurfPos* surf_pos,
        bool neighbor_supp_u, bool neighbor_supp_v, LightArray<NeighborSearch<Scalar>*>& neighbor_searches, int neighbor_index_u);
      /// Calculates integration order for multicomponent DG matrix forms.
      int calc_order_dg_matrix_form(MultiComponentMatrixFormSurf<Scalar>* mfs, Hermes::vector<Solution<Scalar>*> u_ext,
        PrecalcShapeset* fu, PrecalcShapeset* fv, RefMap* ru, SurfPos* surf_pos,
        bool neighbor_supp_u, bool neighbor_supp_v, LightArray<NeighborSearch<Scalar>*>& neighbor_searches, int neighbor_index_u);


      /// Evaluates DG matrix forms on an edge between elements identified by ru_actual, rv.
      Scalar eval_dg_form(MatrixFormSurf<Scalar>* mfs, Hermes::vector<Solution<Scalar>*> u_ext,
        PrecalcShapeset* fu, PrecalcShapeset* fv, RefMap* ru_central, RefMap* ru_actual, RefMap* rv, 
        bool neighbor_supp_u, bool neighbor_supp_v,
        SurfPos* surf_pos, LightArray<NeighborSearch<Scalar>*>& neighbor_searches, int neighbor_index_u, int neighbor_index_v);
      /// Evaluates multicomponent DG matrix forms on an edge between elements identified by ru_actual, rv.
      void eval_dg_form(MultiComponentMatrixFormSurf<Scalar>* mfs, Hermes::vector<Solution<Scalar>*> u_ext,
        PrecalcShapeset* fu, PrecalcShapeset* fv, RefMap* ru_central, RefMap* ru_actual, RefMap* rv, 
        bool neighbor_supp_u, bool neighbor_supp_v,
        SurfPos* surf_pos, LightArray<NeighborSearch<Scalar>*>& neighbor_searches, int neighbor_index_u, int neighbor_index_v, Hermes::vector<Scalar>& result);

      /// Calculates integration order for DG vector forms.
      int calc_order_dg_vector_form(VectorFormSurf<Scalar>* vfs, Hermes::vector<Solution<Scalar>*> u_ext,
        PrecalcShapeset* fv, RefMap* ru, SurfPos* surf_pos,
        LightArray<NeighborSearch<Scalar>*>& neighbor_searches, int neighbor_index_v);
      /// Calculates integration order for multicomponent DG vector forms.
      int calc_order_dg_vector_form(MultiComponentVectorFormSurf<Scalar>* vfs, Hermes::vector<Solution<Scalar>*> u_ext,
        PrecalcShapeset* fv, RefMap* ru, SurfPos* surf_pos,
        LightArray<NeighborSearch<Scalar>*>& neighbor_searches, int neighbor_index_v);

      /// Evaluates DG vector forms on an edge between elements identified by ru_actual, rv.
      Scalar eval_dg_form(VectorFormSurf<Scalar>* vfs, Hermes::vector<Solution<Scalar>*> u_ext,
        PrecalcShapeset* fv, RefMap* rv, 
        SurfPos* surf_pos, LightArray<NeighborSearch<Scalar>*>& neighbor_searches, int neighbor_index_v);
      /// Evaluates multicomponent DG vector forms on an edge between elements identified by ru_actual, rv.
      void eval_dg_form(MultiComponentVectorFormSurf<Scalar>* vfs, Hermes::vector<Solution<Scalar>*> u_ext,
        PrecalcShapeset* fv, RefMap* rv, 
        SurfPos* surf_pos, LightArray<NeighborSearch<Scalar>*>& neighbor_searches, int neighbor_index_v, Hermes::vector<Scalar>& result);

      /// Initialize orders of external functions for volumetric forms.
      ExtData<Hermes::Ord>* init_ext_fns_ord(Hermes::vector<MeshFunction<Scalar>*> &ext);

      /// Initialize orders of external functions for surface forms.
      ExtData<Hermes::Ord>* init_ext_fns_ord(Hermes::vector<MeshFunction<Scalar>*> &ext,
        int edge);

      /// Initialize orders of external functions for DG forms.
      ExtData<Hermes::Ord>* init_ext_fns_ord(Hermes::vector<MeshFunction<Scalar>*> &ext,
        LightArray<NeighborSearch<Scalar>*>& neighbor_searches);

      /// Initialize external functions for volumetric / surface forms.
      ExtData<Scalar>* init_ext_fns(Hermes::vector<MeshFunction<Scalar>*> &ext,  
        RefMap* rm, const int order);

      /// Initialize external functions for DG forms.
      ExtData<Scalar>* init_ext_fns(Hermes::vector<MeshFunction<Scalar>*> &ext, 
        LightArray<NeighborSearch<Scalar>*>& neighbor_searches,
        int order);

      /// Uses assembling_caches to get (possibly cached) precalculated shapeset function values.
      Func<double>* get_fn(PrecalcShapeset* fu, RefMap* rm, const int order);

      /// Uses assembling_caches to get (possibly chached) dummy function for calculation of the integration order.  
      Func<Hermes::Ord>* get_fn_ord(const int order);

      /// Initialize all caches.
      void init_cache();

      /// Deinitialize all caches.
      void delete_cache();

      /// Deinitialize a single geometry cache.
      void delete_single_geom_cache(int order);

      /// Initialize neighbors.
      void init_neighbors(LightArray<NeighborSearch<Scalar>*>& neighbor_searches, const Stage<Scalar>& stage, const int& isurf);

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
      unsigned int update_ns_subtree(NeighborSearch<Scalar>* ns, NeighborNode* node, unsigned int ith_neighbor);

      /// Traverse the multimesh subtree. Used in the function update_ns_subtree().
      void traverse_multimesh_subtree(NeighborNode* node, Hermes::vector<Hermes::vector<unsigned int>*>& running_central_transformations,
        Hermes::vector<Hermes::vector<unsigned int>*>& running_neighbor_transformations, const typename NeighborSearch<Scalar>::NeighborEdgeInfo& edge_info, const int& active_edge, const int& mode);

      /// Returns the matrix_buffer of the size n.
      Scalar** get_matrix_buffer(int n);

      /// Matrix structure as well as spaces and weak formulation is up-to-date.
      bool is_up_to_date();

      /// Minimum identifier of the meshes used in DG assembling in one stage.
      unsigned int min_dg_mesh_seq;

      /// Weak formulation.
      WeakForm<Scalar>* wf;

      /// Seq number of the WeakForm.
      int wf_seq;

      /// Space instances for all equations in the system.
      Hermes::vector<Space<Scalar>*> spaces;

      /// Seq numbers of Space instances in spaces.
      int* sp_seq;

      /// Number of DOFs of all Space instances in spaces.
      int ndof;

      /// Element usage flag: iempty[i] == true if the current state does not posses an active element in the i-th space.
      Hermes::vector<bool> isempty;

      /// Instance of the class Geom used in the calculation of integration order.
      Geom<Hermes::Ord> geom_ord;

      /// If the problem has only constant test functions, there is no need for order calculation,
      /// which saves time.
      bool is_fvm;

      Scalar** matrix_buffer;///< buffer for holding square matrix (during assembling)

      int matrix_buffer_dim;///< dimension of the matrix held by 'matrix_buffer'

      /// Matrix structure can be reused.
      /// If other conditions apply.
      bool have_matrix;

      /// PrecalcShapeset instances for the problem (as many as equations in the system).
      PrecalcShapeset** pss;

      /// Geometry cache.
      Geom<double>* cache_e[g_max_quad + 1 + 4* g_max_quad + 4];

      /// Jacobian * weights cache.
      double* cache_jwt[g_max_quad + 1 + 4* g_max_quad + 4];

      /// There is a matrix form set on DG_INNER_EDGE area or not.
      bool DG_matrix_forms_present;

      /// There is a vector form set on DG_INNER_EDGE area or not.
      bool DG_vector_forms_present;

      /// Turn on Runge-Kutta specific handling of external functions.
      bool RungeKutta;
      
      /// Number of spaces in the original problem in a Runge-Kutta method.
      int RK_original_spaces_count;

      /// Class handling various caches used in assembling.
      class AssemblingCaches 
      {
      public:
        /// Basic constructor.
        AssemblingCaches();

        /// Basic destructor.
        ~AssemblingCaches();

        /// Key for caching precalculated shapeset values on transformed elements with constant
        /// jacobians.
        struct KeyConst
        {
          int index;
          int order;
#ifdef _MSC_VER
          UINT64 sub_idx;
#else
          unsigned int sub_idx;
#endif
          int shapeset_type;
          double inv_ref_map[2][2];
#ifdef _MSC_VER
          KeyConst(int index, int order, UINT64 sub_idx, int shapeset_type, double2x2* inv_ref_map);
#else
          KeyConst(int index, int order, unsigned int sub_idx, int shapeset_type, double2x2* inv_ref_map);
#endif
        };

        /// Functor that compares two above keys (needed e.g. to create a std::map indexed by these keys);
        struct CompareConst 
        {
          bool operator()(KeyConst a, KeyConst b) const;
        };

        /// PrecalcShapeset stored values for Elements with constant jacobian of the reference mapping for triangles.
        std::map<KeyConst, Func<double>* , CompareConst> const_cache_fn_triangles;

        /// PrecalcShapeset stored values for Elements with constant jacobian of the reference mapping for quads.
        std::map<KeyConst, Func<double>* , CompareConst> const_cache_fn_quads;

        /// The same setup for elements with non-constant jacobians.
        /// This cache is deleted with every change of the state in assembling.
        struct KeyNonConst 
        {
          int index;
          int order;
#ifdef _MSC_VER
          UINT64 sub_idx;
#else
          unsigned int sub_idx;
#endif
          int shapeset_type;
#ifdef _MSC_VER
          KeyNonConst(int index, int order, UINT64 sub_idx, int shapeset_type);
#else
          KeyNonConst(int index, int order, unsigned int sub_idx, int shapeset_type);
#endif
        };

        /// Functor that compares two above keys (needed e.g. to create a std::map indexed by these keys);
        struct CompareNonConst 
        {
          bool operator()(KeyNonConst a, KeyNonConst b) const;
        };

        /// PrecalcShapeset stored values for Elements with non-constant jacobian of the reference mapping for triangles.
        std::map<KeyNonConst, Func<double>* , CompareNonConst> cache_fn_triangles;

        /// PrecalcShapeset stored values for Elements with non-constant jacobian of the reference mapping for quads.
        std::map<KeyNonConst, Func<double>* , CompareNonConst> cache_fn_quads;

        /// For caching of values calculated on the previous reference space during adaptivity.
        /// Contains true if the matrix entry for the Element has already been recalculated, and the value
        /// saved is therefore not possible to use.
        bool** element_reassembled_matrix;
        /// For caching of values calculated on the previous reference space during adaptivity.
        /// Contains pair<Element id, true> if the vector entry for the Element has already been recalculated, and the value
        /// saved is therefore not possible to use.
        bool** element_reassembled_vector;

        /// For caching of values calculated on the previous reference space during adaptivity.
        /// The matrix cache, coordinates are : [Space_1][Space_2][Element id on Space_1][DOF on Space_1][DOF on Space_2].
        Scalar***** previous_reference_dp_cache_matrix;

        /// For caching of values calculated on the previous reference space during adaptivity.
        /// The vector cache, coordinates are : [Space][Element id on Space][DOF on Space].
        Scalar*** previous_reference_dp_cache_vector;
        
        /// For caching of values calculated on the previous reference space during adaptivity.
        /// Current size of the [Element id on Space_1] dimension in previous_reference_dp_cache_matrix.
        unsigned int** cache_matrix_size;

        /// For caching of values calculated on the previous reference space during adaptivity.
        /// Current size of the [Element id on Space] dimension in previous_reference_dp_cache_vector.
        unsigned int* cache_vector_size;

        /// For caching of values calculated on the previous reference space during adaptivity.
        /// Previous spaces (incl. mesh) to be able to determine the previous element assembly lists etc.
        Hermes::vector<Space<Scalar>*> stored_spaces_for_adaptivity;

        LightArray<Func<Hermes::Ord>*> cache_fn_ord;
      };

      /// An AssemblingCaches instance for this instance of DiscreteProblem.
      AssemblingCaches assembling_caches;

      /// Class used for profiling of DiscreteProblem (assembling).
      class Profiling
      {
      private:
        Profiling();

        /// Utility time measurement.
        /// For measuring parts of assembling in methods "assemble_*".
        /// Also for initialization parts.
        Hermes::TimePeriod assemble_util_time;

        /// Utility time measurement.
        /// For measuring parts of assembling in methods "eval_*_form".
        /// Except for numerical integration.
        Hermes::TimePeriod eval_util_time;

        /// Total time measurement.
        /// For measuring the whole assembling.
        Hermes::TimePeriod total_time;

        /// Integration time measurement.
        /// For measuring the numerical integration.
        /// Also, the numerical integration when calculating the integration order is counted.
        Hermes::TimePeriod integration_time;

        /// One record stores information about one matrix assembling.
        class Record
        {
        public:
          Record();
          void reset();

          double total;
          double create_sparse_structure;
          double initialization;
          double state_init;
          /// In assemble_* methods.
          double form_preparation_assemble;
          /// In eval_*_form methods.
          double form_preparation_eval;
          double form_evaluation;
        };

        /// The Record instance that is currently being filled.
        Record current_record;

        /// Stores all Records.
        Hermes::vector<Record> profile;

        /// Returns the profiling info.
        void get_profiling_output(std::ostream & out, unsigned int order);

        friend class DiscreteProblem<Scalar>;
      };

      /// A Profiling instance for this instance of DiscreteProblem.
      Profiling profiling;

    public:
      /// Returns all profiling info.
      void get_all_profiling_output(std::ostream & out);

      /// Returns the profiling info from the last call to assemble().
      void get_last_profiling_output(std::ostream & out);

      friend class KellyTypeAdapt<Scalar>;
    };
  }
}
#endif
