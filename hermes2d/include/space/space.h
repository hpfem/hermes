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

#ifndef __H2D_SPACE_H
#define __H2D_SPACE_H

#include "../shapeset/shapeset.h"
#include "asmlist.h"
#include "../boundary_conditions/essential_boundary_conditions.h"
#include "../mesh/traverse.h"
#include "../quadrature/quad_all.h"
#include "algebra/dense_matrix_operations.h"

using namespace Hermes::Algebra::DenseMatrixOperations;

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar> class Space;
  }
}

template<typename Scalar>
#ifdef _WINDOWS
class HERMES_API SpaceSharedPtr : public std::shared_ptr<Hermes::Hermes2D::Space<Scalar> >
#else
class HERMES_API SpaceSharedPtr : public std::tr1::shared_ptr<Hermes::Hermes2D::Space<Scalar> >
#endif
{
public:
  SpaceSharedPtr(Hermes::Hermes2D::Space<Scalar>* ptr = NULL);

  SpaceSharedPtr(const SpaceSharedPtr<Scalar>& other);

  void operator=(const SpaceSharedPtr<Scalar>& other);

  ~SpaceSharedPtr();
};

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Mixins
    {
      /// \ingroup g_mixins2d
      /// Mixin for classes where Spaces can be (re-)set during their existence.
      template<typename Scalar>
      class HERMES_API SettableSpaces
      {
      public:
        /// Sets new spaces for the instance.
        virtual void set_spaces(Hermes::vector<SpaceSharedPtr<Scalar> >& spaces) = 0;
        virtual void set_space(SpaceSharedPtr<Scalar>& space);
        /// Get all spaces as a Hermes::vector.
        virtual Hermes::vector<SpaceSharedPtr<Scalar> >& get_spaces();
        virtual SpaceSharedPtr<Scalar>& get_space(int n);
      };
    }

    template<typename Scalar> class Adapt;
    template<typename Scalar> class DiscreteProblem;
    template<typename Scalar> class DiscreteProblemDGAssembler;
    template<typename Scalar> class DiscreteProblemThreadAssembler;
    namespace Views
    {
      template<typename Scalar> class BaseView;
      template<typename Scalar> class VectorBaseView;
      class Orderizer;
      class OrderView;
    };
    class Shapeset;

    template<typename Scalar> class L2Space;
    template<typename Scalar> class H1Space;
    template<typename Scalar> class HcurlSpace;
    template<typename Scalar> class HdivSpace;

    /** @defgroup spaces FEM Spaces
    * \brief Collection of classes that represent and specify FE spaces.
    */

    /// @ingroup spaces
    /// \brief Represents a finite element space over a domain.
    ///
    /// The Space class represents a finite element space over a domain defined by 'mesh', spanned
    /// by basis functions constructed using 'shapeset'. It serves as a base class for H1Space,
    /// HcurlSpace, HdivSpace and L2Space, since most of the functionality is common for all these spaces.
    ///
    /// There are four main functions the Space class provides:
    /// <ol>
    ///&nbsp; <li> It handles the Dirichlet (essential) boundary conditions. The user provides a pointer
    ///&nbsp;      to an instance of the EssentialBCs class that determines
    ///&nbsp;      which markers represent the Dirichlet part of the boundary, and which markers represent
    ///&nbsp;      the Neumann and Newton parts (all other than those with an Essential BC specified on them).
    ///&nbsp;      Handling of these conditions is done naturally - through weak formulation.
    ///
    ///&nbsp; <li> It stores element polynomial degrees, or 'orders'. All active elements need to have
    ///&nbsp;      an order set for the Space to be valid. Individual orders can be set by calling
    ///&nbsp;      set_element_order(). You can also set the same order for all elements using
    ///&nbsp;      set_uniform_order(). Quadrilateral elements can have different orders in the vertical
    ///&nbsp;      and horizontal directions. It is therefore necessary to form the order using the macro
    ///&nbsp;      H2D_MAKE_QUAD_ORDER() when calling the aforementioned functions.
    ///
    ///&nbsp; <li> It builds and enumerates the basis functions. After all element orders have been set,
    ///&nbsp;      you must call the function assign_dofs(). This function assigns the DOF (degree-of-
    ///&nbsp;      freedom) numbers to basis functions, starting with 'first_dof' (optional parameter).
    ///&nbsp;      It also determines constraining relationships in the mesh due to hanging nodes and
    ///&nbsp;      builds constrained basis functions. The total number of basis functions can then
    ///&nbsp;      be obtained by calling  (). Standard basis functions are assigned positive
    ///&nbsp;      numbers from 'first_dof' to ('first_dof' + (get_num_dofs() - 1) * 'stride'). All
    ///&nbsp;      shape functions belonging to the Dirichlet lift are assigned DOF number of -1. This
    ///&nbsp;      way the Dirichlet lift becomes a (virtual) basis function. This simplifies assembling.
    ///
    ///&nbsp; <li> Finally, and most importantly, the Space is able to deliver a list of shape functions
    ///&nbsp;      existing on each element. Such a list is called an "assembly list" and is represented
    ///&nbsp;      by the class AsmList<Scalar>. The assembly list contains the triplets (idx, dof, coef).
    ///&nbsp;      'idx' is the shape function index, as understood by the Shapeset. 'dof' is the number
    ///&nbsp;      of the basis function, whose part the shape function forms. 'coef' is a Real constant,
    ///&nbsp;      with which the shape function must be multiplied in order to fit into the basis
    ///&nbsp;      function. This is typically 1, but can be less than that for constrained functions.
    ///&nbsp;      Constrained vertex functions can belong to more than one basis functions. This
    ///&nbsp;      results in more triplets with the same 'idx'. However, the assembling procedure or
    ///&nbsp;      the Solution class do not have to worry about that. For instance, the Solution class
    ///&nbsp;      simply obtains the values of all shape functions contained in the list, multiplies them
    ///&nbsp;      by 'coef', and forms a linear combination of them by taking the solution vector values
    ///&nbsp;      at the 'dof' positions. This way the solution to the PDE is obtained.
    /// </ol>
    ///
    /// Space is an abstract class and cannot be instatiated. Use one of the specializations H1Space,
    /// HcurlSpace or L2Space instead.
    /// <br>
    /// The handling of irregular meshes is desribed in H1Space and HcurlSpace.<br>
    ///
    template<typename Scalar>
    class HERMES_API Space : public Hermes::Mixins::Loggable, public Hermes::Hermes2D::Mixins::StateQueryable, public Hermes::Hermes2D::Mixins::XMLParsing
    {
    public:
      Space();
      Space(MeshSharedPtr mesh, Shapeset* shapeset, EssentialBCs<Scalar>* essential_bcs);

      /// Common code for constructors.
      void init();

      /// State querying helpers.
      virtual bool isOkay() const;
      inline std::string getClassName() const { return "Space"; }

      /// Destructor.
      virtual ~Space();

      /// Sets element polynomial order. Can be called by the user. Should not be called
      /// for many elements at once, since assign_dofs() is called at the end of this function.
      virtual void set_element_order(int id, int order, int order_v = -1);

      /// Sets polynomial order to all elements.
      virtual void set_element_orders(int* elem_orders);

      /// Returns element polynomial order.
      int get_element_order(int id) const;

      /// Sets the same polynomial order for all elements in the mesh. Intended for
      /// the user and thus assign_dofs() is called at the end of this function.
      void set_uniform_order(int order, std::string marker = HERMES_ANY);

      /// Set the element order relative to the current order.
      /// The parameter min_order prevents decreasing polynomial order below this threshold.
      void adjust_element_order(int order_change, int min_order);

      /// Version for quads.
      void adjust_element_order(int horizontal_order_change, int vertical_order_change, unsigned int horizontal_min_order, unsigned int vertical_min_order);

      /// Recursively removes all son elements of the given element and
      /// makes it active. Also handles element orders.
      /// \param[in] keep_initial_refinements Refinements in Mesh can be marked as initial (to prevent taking them back), 
      /// this parameter serves to prevent taking them back with this method.
      void unrefine_all_mesh_elements(bool keep_initial_refinements = true);
      
      /// Recursively removes all son elements of the given element and
      /// Version for more spaces sharing the mesh
      static void unrefine_all_mesh_elements(Hermes::vector<SpaceSharedPtr<Scalar> > spaces, bool keep_initial_refinements = true);

      /// Updates element orders when the underlying mesh has been refined.
      void update_element_orders_after_refinement();

      /// Sets the shapeset.
      virtual void set_shapeset(Shapeset* shapeset) = 0;

      /// \brief Returns the number of basis functions contained in the space.
      int get_num_dofs() const;

      /// \brief Returns the number of basis functions contained in the spaces.
      static int get_num_dofs(Hermes::vector<SpaceSharedPtr<Scalar> > spaces);

      /// \brief Returns the number of basis functions contained in the space.
      static int get_num_dofs(SpaceSharedPtr<Scalar> space);

      MeshSharedPtr get_mesh() const;

      /// \brief Sets a (new) mesh and calls assign_dofs().
      void set_mesh(MeshSharedPtr mesh);

      /// \brief Sets a (new) mesh seq, and mesh_seq.
      void set_mesh_seq(int seq);

      /// Sets the boundary condition.
      void set_essential_bcs(EssentialBCs<Scalar>* essential_bcs);

      /// Obtains an boundary conditions
      EssentialBCs<Scalar>* get_essential_bcs() const;

      /// Updates essential BC values. Typically used for time-dependent
      /// essnetial boundary conditions.
      void update_essential_bc_values();

      Shapeset* get_shapeset() const;

      /// Saves this space into a file.
      void save(const char *filename) const;
#ifdef WITH_BSON
      void save_bson(const char* filename) const;
#endif

      /// Loads a space from a file in XML format.
      static SpaceSharedPtr<Scalar> load(const char *filename, MeshSharedPtr mesh, bool validate = false, EssentialBCs<Scalar>* essential_bcs = NULL, Shapeset* shapeset = NULL);
      /// This method is here for rapid re-loading.
      void load(const char *filename);

#ifdef WITH_BSON
      /// Loads a space from a file in BSON.
      static SpaceSharedPtr<Scalar> load_bson(const char *filename, MeshSharedPtr mesh, EssentialBCs<Scalar>* essential_bcs = NULL, Shapeset* shapeset = NULL);
      /// This method is here for rapid re-loading.
      void load_bson(const char *filename);
#endif

      /// Obtains an assembly list for the given element.
      virtual void get_element_assembly_list(Element* e, AsmList<Scalar>* al) const;

      /// Copy from Space instance 'space'
      /// \param[in] new_mesh Mesh where data will be copied to.
      virtual void copy(SpaceSharedPtr<Scalar> space, MeshSharedPtr new_mesh);

      /// Class for creating reference space.
      class HERMES_API ReferenceSpaceCreator
      {
      public:
        /// Constructor.
        /// \param[in] order_increase Increase of the polynomial order.
        ReferenceSpaceCreator(unsigned int order_increase = 1);
        /// Constructor.
        /// \param[in] coarse_space The coarse (original) space.
        /// \param[in] ref_mesh The refined mesh.
        /// \param[in] order_increase Increase of the polynomial order.
        ReferenceSpaceCreator(SpaceSharedPtr<Scalar> coarse_space, MeshSharedPtr ref_mesh, unsigned int order_increase = 1);

        /// Method that does the creation.
        /// THIS IS THE METHOD TO OVERLOAD FOR CUSTOM CREATING OF A REFERENCE SPACE.
        virtual void handle_orders(SpaceSharedPtr<Scalar> ref_space);

        /// Methods that user calls to get the reference space pointer (has to be properly casted if necessary).
        virtual SpaceSharedPtr<Scalar> create_ref_space(bool assign_dofs = true);
        virtual SpaceSharedPtr<Scalar> create_ref_space(SpaceSharedPtr<Scalar> coarse_space, MeshSharedPtr ref_mesh, bool assign_dofs = true);

      private:
        /// Construction initialization.
        SpaceSharedPtr<Scalar> init_construction_l2();
        SpaceSharedPtr<Scalar> init_construction_h1();
        SpaceSharedPtr<Scalar> init_construction_hcurl();
        SpaceSharedPtr<Scalar> init_construction_hdiv();

        /// Construction finalization.
        virtual void finish_construction(SpaceSharedPtr<Scalar> ref_space);

        /// Storage.
        SpaceSharedPtr<Scalar> coarse_space;
        MeshSharedPtr ref_mesh;
        unsigned int order_increase;
      };

      /// Sets element polynomial order. This version does not call assign_dofs() and is
      /// intended primarily for internal use.
      virtual void set_element_order_internal(int id, int order, int order_v = -1);

      /// \brief Builds basis functions and assigns DOF numbers to them.
      /// \details This functions must be called \b after assigning element orders, and \b before
      /// using the space in a computation, otherwise an error will occur.
      /// \param first_dof[in] The DOF number of the first basis function.
      /// \param stride[in] The difference between the DOF numbers of successive basis functions.
      /// \return The number of basis functions contained in the space.
      virtual int assign_dofs(int first_dof = 0, int stride = 1);

      /// \brief Assings the degrees of freedom to all Spaces in the Hermes::vector.
      static int assign_dofs(Hermes::vector<SpaceSharedPtr<Scalar> > spaces);

      virtual Scalar* get_bc_projection(SurfPos* surf_pos, int order, EssentialBoundaryCondition<Scalar> *bc) = 0;

      static void update_essential_bc_values(Hermes::vector<SpaceSharedPtr<Scalar> >& spaces, double time);

      static void update_essential_bc_values(SpaceSharedPtr<Scalar>& space, double time);

      /// Internal. Return type of this space (H1 = HERMES_H1_SPACE, Hcurl = HERMES_HCURL_SPACE,
      /// Hdiv = HERMES_HDIV_SPACE, L2 = HERMES_L2_SPACE)
      virtual SpaceType get_type() const = 0;

      static Node* get_mid_edge_vertex_node(Element* e, int i, int j);

      /// Sets polynomial orders to elements created by Mesh::regularize() using "parents".
      void distribute_orders(MeshSharedPtr mesh, int* parents);

      /// Internal. Obtains the order of an edge, according to the minimum rule.
      virtual int get_edge_order(Element* e, int edge) const;

      /// \brief Returns the DOF number of the last basis function.
      int get_max_dof() const;

      /// Returns true if the space is ready for computation, false otherwise.
      bool is_up_to_date() const;

      /// Obtains an edge assembly list (contains shape functions that are nonzero on the specified edge).
      void get_boundary_assembly_list(Element* e, int surf_num, AsmList<Scalar>* al) const;

      /// Sets the same polynomial order for all elements in the mesh. Does not
      /// call assign_dofs(). For internal use.
      void set_uniform_order_internal(int order, int marker);

      void free();

      /// Returns the total (global) number of vertex functions.
      /// The DOF ordering starts with vertex functions, so it it necessary to know how many of them there are.
      int get_vertex_functions_count();
      /// Returns the total (global) number of edge functions.
      int get_edge_functions_count();
      /// Returns the total (global) number of bubble functions.
      int get_bubble_functions_count();
      
      /// Internal. Used by DiscreteProblem to detect changes in the space.
      int get_seq() const;

    protected:
      /// Number of degrees of freedom (dimension of the space).
      int ndof;

      static const int H2D_UNASSIGNED_DOF = -2; ///< DOF which was not assigned yet.
      static const int H2D_CONSTRAINED_DOF = -1; ///< DOF which is constrained.

      Shapeset* shapeset;

      bool own_shapeset;  ///< true if default shapeset is created in the constructor, false if shapeset is supplied by user.

      /// Boundary conditions.
      EssentialBCs<Scalar>* essential_bcs;

      /// FE mesh.
      MeshSharedPtr mesh;

      int default_tri_order, default_quad_order;
      int vertex_functions_count, edge_functions_count, bubble_functions_count;
      int first_dof, next_dof;
      int stride;
      int seq, mesh_seq;
      int was_assigned;

      struct BaseComponent
      {
        int dof;
        Scalar coef;
      };

      union NodeData
      {
        struct // regular node
        {
          int dof;
          union {
            Scalar* edge_bc_proj;
            Scalar* vertex_bc_coef;
          };
          int n; ///< Number of dofs. Temporarily used during assignment
          ///< of DOFs to indicate nodes which were not processed yet.
        };
        struct // constrained vertex node
        {
          BaseComponent* baselist;
          int ncomponents;
        };
        struct // constrained edge node
        {
          Node* base;
          int part;
        };
        NodeData() : dof(0), edge_bc_proj(NULL) {}
      };

      class ElementData
      {
      public:
        ElementData() : changed_in_last_adaptation(true) {};
        int order;
        int bdof, n;
        bool changed_in_last_adaptation;
      };

      NodeData* ndata;    ///< node data table
      ElementData* edata; ///< element data table
      int nsize, ndata_allocated; ///< number of items in ndata, allocated space
      int esize;

      virtual int get_edge_order_internal(Node* en) const;

      /// Recursively removes all son elements of the given element and
      /// makes it active. Also handles element orders.
      /// \param[in] keep_initial_refinements Refinements in Mesh can be marked as initial (to prevent taking them back), 
      /// this parameter serves to prevent taking them back with this method.
      /// \param[in] only_unrefine_space_data Useful when more spaces share the mesh if one wants to unrefine the underlying
      /// Mesh only once, but wants other spaces know about the change.
      void unrefine_all_mesh_elements_internal(bool keep_initial_refinements, bool only_unrefine_space_data);
      
      /// \brief Updates internal node and element tables.
      /// \details Since meshes only contain geometric information, the Space class keeps two
      /// tables with FEM-related information. The first one, 'ndata', contains DOF numbers
      /// and other things for each node. The second table, 'edata', holds element orders
      /// and bubble DOF numbers. Both tables are directly indexed by the node and element
      /// IDs. The function resize_tables() is called to check whether the tables are large
      /// enough to contain all node and element id's, and to reallocate them if not.
      virtual void resize_tables();

      void update_orders_recurrent(Element* e, int order);

      virtual void reset_dof_assignment(); ///< Resets assignment of DOF to an unassigned state.
      virtual void assign_vertex_dofs() = 0;
      virtual void assign_edge_dofs() = 0;
      virtual void assign_bubble_dofs() = 0;

      virtual void get_vertex_assembly_list(Element* e, int iv, AsmList<Scalar>* al) const = 0;
      virtual void get_boundary_assembly_list_internal(Element* e, int surf_num, AsmList<Scalar>* al) const = 0;
      virtual void get_bubble_assembly_list(Element* e, AsmList<Scalar>* al) const;

      double** proj_mat;
      double*  chol_p;

      /// Used for bc projection.
      Hermes::vector<void*> bc_data;

      void precalculate_projection_matrix(int nv, double**& mat, double*& p);
      void update_edge_bc(Element* e, SurfPos* surf_pos);

      /// Called by Space to update constraining relationships between shape functions due
      /// to hanging nodes in the mesh. As this is space-specific, this function is reimplemented
      /// in H1Space and HcurlSpace.
      virtual void update_constraints();

      /// Auxiliary function the descendants may implement to perform additional tasks after
      /// the DOFs have been assigned.
      virtual void post_assign();

      void free_bc_data();

      /// Internal.
      /// Returns a new Space according to the type provided.
      /// Used in loading.
      static SpaceSharedPtr<Scalar> init_empty_space(const char* spaceType, MeshSharedPtr mesh, Shapeset* shapeset);

      template<typename T> friend class OGProjection;
      template<typename T> friend class NewtonSolver;
      template<typename T> friend class PicardSolver;
      template<typename T> friend class LinearSolver;
      template<typename T> friend class OGProjectionNOX;
      template<typename T> friend class Solution;
      template<typename T> friend class RungeKutta;
      template<typename T> friend class ExactSolution;
      template<typename T> friend class NeighborSearch;
      template<typename T> friend class ExactSolutionScalar;
      template<typename T> friend class ExactSolutionVector;
      template<typename T> friend class Views::BaseView;
      friend class Views::Orderizer;
      friend class Views::OrderView;
      template<typename T> friend class Views::VectorBaseView;
      friend class Adapt<Scalar>;
      friend class DiscreteProblem<Scalar>;
      friend class DiscreteProblemDGAssembler<Scalar>;
      friend class DiscreteProblemThreadAssembler<Scalar>;
      template<typename T> friend class CalculationContinuity;
    };
  }
}
#endif
