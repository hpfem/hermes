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

#include "global.h"
#include "integrals/h1.h"
#include "quadrature/limit_order.h"
#include "discrete_problem.h"
#include "mesh/traverse.h"
#include "space/space.h"
#include "shapeset/precalc.h"
#include "mesh/refmap.h"
#include "function/solution.h"
#include "neighbor.h"

using namespace Hermes::Algebra::DenseMatrixOperations;

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    double DiscreteProblem<Scalar>::fake_wt = 1.0;

    template<typename Scalar>
    DiscreteProblem<Scalar>::DiscreteProblem(const WeakForm<Scalar>* wf, Hermes::vector<const Space<Scalar> *> spaces) : wf(wf), wf_seq(-1)
    {
      _F_;
      if (spaces.empty()) throw Exceptions::NullException(2);
      unsigned int first_dof_running = 0;
      for(unsigned int i = 0; i < spaces.size(); i++)
      {
        this->spaces.push_back(spaces.at(i));
        this->spaces_first_dofs.push_back(first_dof_running);
        first_dof_running += spaces.at(i)->get_num_dofs();
      }
      init();
    }

    template<typename Scalar>
    DiscreteProblem<Scalar>::DiscreteProblem(const WeakForm<Scalar>* wf, const Space<Scalar>* space)
      : wf(wf), wf_seq(-1)
    {
      _F_;
      spaces.push_back(space);
      this->spaces_first_dofs.push_back(0);

      init();
    }

    template<typename Scalar>
    DiscreteProblem<Scalar>::DiscreteProblem() : wf(NULL), current_stage(NULL)
    {
      // Set all attributes for which we don't need to acces wf or spaces.
      // This is important for the destructor to properly detect what needs to be deallocated.
      sp_seq = NULL;
      is_fvm = false;
      RungeKutta = false;
      RK_original_spaces_count = 0;
      matrix_buffer = NULL;
      matrix_buffer_dim = 0;
      have_matrix = false;
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::init()
    {
      _F_

      // Initialize special variable for Runge-Kutta time integration.
      RungeKutta = false;
      RK_original_spaces_count = 0;

      ndof = Space<Scalar>::get_num_dofs(spaces);

      // Sanity checks.
      if(wf == NULL)
        error("WeakForm<Scalar>* wf can not be NULL in DiscreteProblem<Scalar>::DiscreteProblem.");

      if (spaces.size() != (unsigned) wf->get_neq())
        error("Bad number of spaces in DiscreteProblem.");
      if (spaces.size() == 0)
        error("Zero number of spaces in DiscreteProblem.");

      // Internal variables settings.
      sp_seq = new int[wf->get_neq()];
      memset(sp_seq, -1, sizeof(int) * wf->get_neq());

      // Matrix<Scalar> related settings.
      matrix_buffer = NULL;
      matrix_buffer_dim = 0;
      have_matrix = false;

      // There is a special function that sets a DiscreteProblem to be FVM.
      // Purpose is that this constructor looks cleaner and is simpler.
      this->is_fvm = false;

      Geom<Hermes::Ord> *tmp = init_geom_ord();
      geom_ord = *tmp;
      delete tmp;
      quad = &g_quad_2d_std;

      current_mat = NULL;
      current_stage = NULL;
      current_rhs = NULL;
      current_block_weights = NULL;
    }

    template<typename Scalar>
    DiscreteProblem<Scalar>::~DiscreteProblem()
    {
      _F_;
      if (wf != NULL)
        memset(sp_seq, -1, sizeof(int) * wf->get_neq());
      wf_seq = -1;
      if (sp_seq != NULL) delete [] sp_seq;
    }

    template<typename Scalar>
    int DiscreteProblem<Scalar>::get_num_dofs()
    {
      _F_;
      ndof = 0;
      for (unsigned int i = 0; i < wf->get_neq(); i++)
        ndof += spaces[i]->get_num_dofs();
      return ndof;
    }

    template<typename Scalar>
    Scalar** DiscreteProblem<Scalar>::get_matrix_buffer(int n)
    {
      _F_;
      if (n <= matrix_buffer_dim)
        return matrix_buffer;
      if (matrix_buffer != NULL)
        delete [] matrix_buffer;
      matrix_buffer_dim = n;
      return (matrix_buffer = new_matrix<Scalar>(n, n));
    }

    template<typename Scalar>
    const Space<Scalar>* DiscreteProblem<Scalar>::get_space(int n)
    {
      return this->spaces[n];
    }

    template<typename Scalar>
    const WeakForm<Scalar>* DiscreteProblem<Scalar>::get_weak_formulation()
    {
      return this->wf;
    }

    template<typename Scalar>
    Hermes::vector<const Space<Scalar>*> DiscreteProblem<Scalar>::get_spaces()
    {
      return this->spaces;
    }

    template<typename Scalar>
    bool DiscreteProblem<Scalar>::is_matrix_free()
    {
      return wf->is_matrix_free();
    }

    template<typename Scalar>
    bool DiscreteProblem<Scalar>::is_up_to_date()
    {
      _F_;
      // check if we can reuse the matrix structure
      bool up_to_date = true;
      if (!have_matrix)
        up_to_date = false;

      for (unsigned int i = 0; i < wf->get_neq(); i++)
      {
        if (spaces[i]->get_seq() != sp_seq[i])
        {
          up_to_date = false;
          break;
        }
      }

      if (wf->get_seq() != wf_seq)
        up_to_date = false;

      return up_to_date;
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::invalidate_matrix()
    {
      have_matrix = false;
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::set_fvm()
    {
      this->is_fvm = true;
    }
      
    template<typename Scalar>
    double DiscreteProblem<Scalar>::block_scaling_coeff(MatrixForm<Scalar>* form)
    {
      if(current_block_weights != NULL)
        return current_block_weights->get_A(form->i, form->j);
      return 1.0;
    }

    template<typename Scalar>
    bool DiscreteProblem<Scalar>::form_to_be_assembled(MatrixForm<Scalar>* form, Traverse::State* current_state)
    {
      if (current_state->e[form->i] == NULL || current_state->e[form->j] == NULL)
        return false;
      if (fabs(form->scaling_factor) < 1e-12)
        return false;

      // If a block scaling table is provided, and if the scaling coefficient
      // A_mn for this block is zero, then the form does not need to be assembled.
      if (current_block_weights != NULL)
        if (fabs(current_block_weights->get_A(form->i, form->j)) < 1e-12)
          return false;
      return true;
    }

    template<typename Scalar>
    bool DiscreteProblem<Scalar>::form_to_be_assembled(MatrixFormVol<Scalar>* form, Traverse::State* current_state)
    {
      if(!form_to_be_assembled((MatrixForm<Scalar>*)form, current_state))
        return false;
      
      // Assemble this form only if one of its areas is HERMES_ANY
      // of if the element marker coincides with one of the form's areas.
      bool assemble_this_form = false;
      for (unsigned int ss = 0; ss < form->areas.size(); ss++)
      {
        if(form->areas[ss] == HERMES_ANY)
        {
          assemble_this_form = true;
          break;
        }
        else
        {
          bool marker_on_space_m = this->spaces[form->i]->get_mesh()->get_element_markers_conversion().get_internal_marker(form->areas[ss]).valid;
          if(marker_on_space_m)
            marker_on_space_m = (this->spaces[form->i]->get_mesh()->get_element_markers_conversion().get_internal_marker(form->areas[ss]).marker == current_state->rep->marker);

          bool marker_on_space_n = this->spaces[form->j]->get_mesh()->get_element_markers_conversion().get_internal_marker(form->areas[ss]).valid;
          if(marker_on_space_n)
            marker_on_space_n = (this->spaces[form->j]->get_mesh()->get_element_markers_conversion().get_internal_marker(form->areas[ss]).marker == current_state->rep->marker);

          if (marker_on_space_m && marker_on_space_n)
          {
            assemble_this_form = true;
            break;
          }
        }
      }
      if (!assemble_this_form)
        return false;
      return true;
    }

    template<typename Scalar>
    bool DiscreteProblem<Scalar>::form_to_be_assembled(MatrixFormSurf<Scalar>* form, Traverse::State* current_state)
    {
      if(current_state->rep->en[current_state->isurf]->marker == 0)
        return false;

      if (form->areas[0] == H2D_DG_INNER_EDGE)
        return false;
      if(!form_to_be_assembled((MatrixForm<Scalar>*)form, current_state))
        return false;

      bool assemble_this_form = false;
      for (unsigned int ss = 0; ss < form->areas.size(); ss++)
      {
        if(form->areas[ss] == HERMES_ANY || form->areas[ss] == H2D_DG_BOUNDARY_EDGE)
        {
          assemble_this_form = true;
          break;
        }
        else
        {
          bool marker_on_space_m = this->spaces[form->i]->get_mesh()->get_boundary_markers_conversion().get_internal_marker(form->areas[ss]).valid;
          if(marker_on_space_m)
            marker_on_space_m = (this->spaces[form->i]->get_mesh()->get_boundary_markers_conversion().get_internal_marker(form->areas[ss]).marker == current_state->rep->en[current_state->isurf]->marker);

          bool marker_on_space_n = this->spaces[form->j]->get_mesh()->get_boundary_markers_conversion().get_internal_marker(form->areas[ss]).valid;
          if(marker_on_space_n)
            marker_on_space_n = (this->spaces[form->j]->get_mesh()->get_boundary_markers_conversion().get_internal_marker(form->areas[ss]).marker == current_state->rep->en[current_state->isurf]->marker);

          if (marker_on_space_m && marker_on_space_n)
          {
            assemble_this_form = true;
            break;
          }
        }
      }
      if (assemble_this_form == false)
        return false;
      return true;
    }
    
    template<typename Scalar>
    bool DiscreteProblem<Scalar>::form_to_be_assembled(VectorForm<Scalar>* form, Traverse::State* current_state)
    {
      if (current_state->e[form->i] == NULL)
        return false;
      if (fabs(form->scaling_factor) < 1e-12)
        return false;

      return true;
    }
    
    template<typename Scalar>
    bool DiscreteProblem<Scalar>::form_to_be_assembled(VectorFormVol<Scalar>* form, Traverse::State* current_state)
    {
      if(!form_to_be_assembled((VectorForm<Scalar>*)form, current_state))
        return false;
      
      // Assemble this form only if one of its areas is HERMES_ANY
      // of if the element marker coincides with one of the form's areas.
      bool assemble_this_form = false;
      for (unsigned int ss = 0; ss < form->areas.size(); ss++)
      {
        if(form->areas[ss] == HERMES_ANY)
        {
          assemble_this_form = true;
          break;
        }
        else
        {
          bool marker_on_space_m = this->spaces[form->i]->get_mesh()->get_element_markers_conversion().get_internal_marker(form->areas[ss]).valid;
          if(marker_on_space_m)
            marker_on_space_m = (this->spaces[form->i]->get_mesh()->get_element_markers_conversion().get_internal_marker(form->areas[ss]).marker == current_state->rep->marker);

          if (marker_on_space_m)
          {
            assemble_this_form = true;
            break;
          }
        }
      }
      if (!assemble_this_form)
        return false;
      return true;
    }
    
    template<typename Scalar>
    bool DiscreteProblem<Scalar>::form_to_be_assembled(VectorFormSurf<Scalar>* form, Traverse::State* current_state)
    {
      if(current_state->rep->en[current_state->isurf]->marker == 0)
        return false;

      if (form->areas[0] == H2D_DG_INNER_EDGE)
        return false;
      
      if(!form_to_be_assembled((VectorForm<Scalar>*)form, current_state))
        return false;
      
      bool assemble_this_form = false;
      for (unsigned int ss = 0; ss < form->areas.size(); ss++)
      {
        if(form->areas[ss] == HERMES_ANY || form->areas[ss] == H2D_DG_BOUNDARY_EDGE)
        {
          assemble_this_form = true;
          break;
        }
        else
        {
          bool marker_on_space_m = this->spaces[form->i]->get_mesh()->get_boundary_markers_conversion().get_internal_marker(form->areas[ss]).valid;
          if(marker_on_space_m)
            marker_on_space_m = (this->spaces[form->i]->get_mesh()->get_boundary_markers_conversion().get_internal_marker(form->areas[ss]).marker == current_state->rep->en[current_state->isurf]->marker);

          if (marker_on_space_m)
          {
            assemble_this_form = true;
            break;
          }
        }
      }
      if (assemble_this_form == false)
        return false;
      return true;
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::create_sparse_structure(SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs)
    {
      this->current_mat = mat;
      if(rhs != NULL)
        this->current_rhs = rhs;
      this->create_sparse_structure();
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::create_sparse_structure()
    {
      _F_;

      if (is_up_to_date())
      {
        if (current_mat != NULL)
        {
          verbose("Reusing matrix sparse structure.");
          current_mat->zero();
        }
        if (current_rhs != NULL)
        {
          // If we use e.g. a new NewtonSolver (providing a new Vector) for this instance of DiscreteProblem that already assembled a system,
          // we end up with everything up_to_date, but unallocated Vector.
          if(current_rhs->length() == 0)
            current_rhs->alloc(ndof);
          else
            current_rhs->zero();
        }
        return;
      }

      // For DG, the sparse structure is different as we have to
      // account for over-edge calculations.
      bool is_DG = false;
      for(unsigned int i = 0; i < this->wf->mfsurf.size(); i++)
      {
        if(this->wf->mfsurf[i]->areas[0] == H2D_DG_INNER_EDGE)
        {
          is_DG = true;
          break;
        }
      }
      for(unsigned int i = 0; i < this->wf->vfsurf.size() && is_DG == false; i++)
      {
        if(this->wf->vfsurf[i]->areas[0] == H2D_DG_INNER_EDGE)
        {
          is_DG = true;
          break;
        }
      }

      if (current_mat != NULL)
      {
        // Spaces have changed: create the matrix from scratch.
        have_matrix = true;
        current_mat->free();
        current_mat->prealloc(ndof);

        AsmList<Scalar>* al = new AsmList<Scalar>[wf->get_neq()];
        Mesh** meshes = new Mesh*[wf->get_neq()];
        bool **blocks = wf->get_blocks(current_force_diagonal_blocks);

        // Init multi-mesh traversal.
        for (unsigned int i = 0; i < wf->get_neq(); i++)
          meshes[i] = spaces[i]->get_mesh();

        Traverse trav;
        trav.begin(wf->get_neq(), meshes);

        if(is_DG)
        {
          Hermes::vector<Space<Scalar>*> mutable_spaces;
          for(unsigned int i = 0; i < this->spaces.size(); i++)
          {
            mutable_spaces.push_back(const_cast<Space<Scalar>*>(spaces.at(i)));
            spaces_first_dofs[i] = 0;
          }
          Space<Scalar>::assign_dofs(mutable_spaces);
        }

        Traverse::State* current_state;
        // Loop through all elements.
        while ((current_state = trav.get_next_state()) != NULL)
        {
          // Obtain assembly lists for the element at all spaces.
          /// \todo do not get the assembly list again if the element was not changed.
          for (unsigned int i = 0; i < wf->get_neq(); i++)
            if (current_state->e[i] != NULL)
              if(is_DG)
                spaces[i]->get_element_assembly_list(current_state->e[i], &(al[i]));
              else
                spaces[i]->get_element_assembly_list(current_state->e[i], &(al[i]), spaces_first_dofs[i]);

          if(is_DG)
          {
            // Number of edges ( =  number of vertices).
            int num_edges = current_state->e[0]->get_num_surf();

            // Allocation an array of arrays of neighboring elements for every mesh x edge.
            Element **** neighbor_elems_arrays = new Element *** [wf->get_neq()];
            for(unsigned int i = 0; i < wf->get_neq(); i++)
              neighbor_elems_arrays[i] = new Element ** [num_edges];

            // The same, only for number of elements
            int ** neighbor_elems_counts = new int * [wf->get_neq()];
            for(unsigned int i = 0; i < wf->get_neq(); i++)
              neighbor_elems_counts[i] = new int [num_edges];

            // Get the neighbors.
            for(unsigned int el = 0; el < wf->get_neq(); el++)
            {
              NeighborSearch<Scalar> ns(current_state->e[el], meshes[el]);

              // Ignoring errors (and doing nothing) in case the edge is a boundary one.
              ns.set_ignore_errors(true);

              for(int ed = 0; ed < num_edges; ed++)
              {
                ns.set_active_edge(ed);
                const Hermes::vector<Element *> *neighbors = ns.get_neighbors();

                neighbor_elems_counts[el][ed] = ns.get_num_neighbors();
                neighbor_elems_arrays[el][ed] = new Element * [neighbor_elems_counts[el][ed]];
                for(int neigh = 0; neigh < neighbor_elems_counts[el][ed]; neigh++)
                  neighbor_elems_arrays[el][ed][neigh] = (*neighbors)[neigh];
              }
            }

            // Pre-add into the stiffness matrix.
            for (unsigned int m = 0; m < wf->get_neq(); m++)
              for(unsigned int el = 0; el < wf->get_neq(); el++)
                for(int ed = 0; ed < num_edges; ed++)
                  for(int neigh = 0; neigh < neighbor_elems_counts[el][ed]; neigh++)
                    if ((blocks[m][el] || blocks[el][m]) && current_state->e[m] != NULL)
                    {
                      AsmList<Scalar>*am = &(al[m]);
                      AsmList<Scalar>*an = new AsmList<Scalar>;
                      spaces[el]->get_element_assembly_list(neighbor_elems_arrays[el][ed][neigh], an);

                      // pretend assembling of the element stiffness matrix
                      // register nonzero elements
                      for (unsigned int i = 0; i < am->cnt; i++)
                        if (am->dof[i] >= 0)
                          for (unsigned int j = 0; j < an->cnt; j++)
                            if (an->dof[j] >= 0)
                            {
                              if(blocks[m][el]) current_mat->pre_add_ij(am->dof[i], an->dof[j]);
                              if(blocks[el][m]) current_mat->pre_add_ij(an->dof[j], am->dof[i]);
                            }
                            delete an;
                    }

                    // Deallocation an array of arrays of neighboring elements
                    // for every mesh x edge.
                    for(unsigned int el = 0; el < wf->get_neq(); el++)
                    {
                      for(int ed = 0; ed < num_edges; ed++)
                        delete [] neighbor_elems_arrays[el][ed];
                      delete [] neighbor_elems_arrays[el];
                    }
                    delete [] neighbor_elems_arrays;

                    // The same, only for number of elements.
                    for(unsigned int el = 0; el < wf->get_neq(); el++)
                      delete [] neighbor_elems_counts[el];
                    delete [] neighbor_elems_counts;
          }

          // Go through all equation-blocks of the local stiffness matrix.
          for (unsigned int m = 0; m < wf->get_neq(); m++)
          {
            for (unsigned int n = 0; n < wf->get_neq(); n++)
            {
              if (blocks[m][n] && current_state->e[m] != NULL && current_state->e[n] != NULL)
              {
                AsmList<Scalar>*am = &(al[m]);
                AsmList<Scalar>*an = &(al[n]);

                // Pretend assembling of the element stiffness matrix.
                for (unsigned int i = 0; i < am->cnt; i++)
                  if (am->dof[i] >= 0)
                    for (unsigned int j = 0; j < an->cnt; j++)
                      if (an->dof[j] >= 0)
                        current_mat->pre_add_ij(am->dof[i], an->dof[j]);
              }
            }
          }
        }

        trav.finish();
        delete [] al;
        delete [] meshes;
        delete [] blocks;

        current_mat->alloc();
      }

      // WARNING: unlike Matrix<Scalar>::alloc(), Vector<Scalar>::alloc(ndof) frees the memory occupied
      // by previous vector before allocating
      if (current_rhs != NULL)
        current_rhs->alloc(ndof);

      // save space seq numbers and weakform seq number, so we can detect their changes
      for (unsigned int i = 0; i < wf->get_neq(); i++)
        sp_seq[i] = spaces[i]->get_seq();

      wf_seq = wf->get_seq();
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::assemble(SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs,
      bool force_diagonal_blocks, Table* block_weights)
    {
      _F_;
      Scalar* coeff_vec = NULL;
      assemble(coeff_vec, mat, rhs, force_diagonal_blocks, block_weights);
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::assemble(Vector<Scalar>* rhs,
      bool force_diagonal_blocks, Table* block_weights)
    {
      _F_;
      assemble(NULL, NULL, rhs, force_diagonal_blocks, block_weights);
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::assemble(Scalar* coeff_vec, SparseMatrix<Scalar>* mat,
      Vector<Scalar>* rhs,
      bool force_diagonal_blocks,
      Table* block_weights)
    {
      _F_;

      current_rhs = rhs;
      current_mat = mat;
      current_force_diagonal_blocks = force_diagonal_blocks;
      current_block_weights = block_weights;

      // Check that the block scaling table have proper dimension.
      if (block_weights != NULL)
        if (block_weights->get_size() != wf->get_neq())
          throw Exceptions::LengthException(6, block_weights->get_size(), wf->get_neq());

      // Creating matrix sparse structure.
      create_sparse_structure();

      // Reset the warnings about insufficiently high integration order.
      reset_warn_order();

      // Initialize matrix buffer.
      matrix_buffer = NULL;
      matrix_buffer_dim = 0;
      if (mat != NULL)
        get_matrix_buffer(9);

      // Create assembling stages.
      Hermes::vector<Stage<Scalar> > stages = Hermes::vector<Stage<Scalar> >();
      bool want_matrix = (mat != NULL);
      bool want_vector = (rhs != NULL);
      wf->get_stages(spaces, stages, want_matrix, want_vector);

      // Loop through all assembling stages -- the purpose of this is increased performance
      // in multi-mesh calculations, where, e.g., only the right hand side uses two meshes.
      // In such a case, the matrix forms are assembled over one mesh, and only the rhs
      // traverses through the union mesh. On the other hand, if you don't use multi-mesh
      // at all, there will always be only one stage in which all forms are assembled as usual.
      for (unsigned ss = 0; ss < stages.size(); ss++)
      {
        current_stage = &stages[ss];
        // Check that there is a DG form, so that the DG assembling procedure needs to be performed.
        is_DG_stage();

        // Assemble one stage. One stage is a collection of functions,
        // and meshes that can not be further minimized.
        // E.g. if a linear form uses two external solutions, each of
        // which is defined on a different mesh, and different to the
        // mesh of the current test function, then the stage would have
        // three meshes. By stage functions, all functions are meant: shape
        // functions (their precalculated values), and mesh functions.
        assemble_one_stage(coeff_vec);
      }
      // Deinitialize matrix buffer.
      if(matrix_buffer != NULL)
        delete [] matrix_buffer;
      matrix_buffer = NULL;
      matrix_buffer_dim = 0;
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::assemble(Scalar* coeff_vec, Vector<Scalar>* rhs,
      bool force_diagonal_blocks, Table* block_weights)
    {
      _F_;
      assemble(coeff_vec, NULL, rhs, force_diagonal_blocks, block_weights);
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::is_DG_stage()
    {
      DG_matrix_forms_present = false;
      DG_vector_forms_present = false;
      for(unsigned int i = 0; i < current_stage->mfsurf.size() && DG_matrix_forms_present == false; i++)
      {
        if (current_stage->mfsurf[i]->areas[0] == H2D_DG_INNER_EDGE)
        {
          DG_matrix_forms_present = true;
          break;
        }
      }
      for(unsigned int i = 0; i < current_stage->vfsurf.size() && DG_vector_forms_present == false; i++)
      {
        if (current_stage->vfsurf[i]->areas[0] == H2D_DG_INNER_EDGE)
        {
          DG_vector_forms_present = true;
          break;
        }
      }
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::assemble_one_stage(Scalar* coeff_vec)
    {
      _F_;

      PrecalcShapeset*** pss = new PrecalcShapeset**[omp_get_max_threads()];
      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
      {
        pss[i] = new PrecalcShapeset*[wf->get_neq()];
        for (unsigned int j = 0; j < wf->get_neq(); j++)
          pss[i][j] = new PrecalcShapeset(spaces[j]->shapeset);
      }
      PrecalcShapeset*** spss = new PrecalcShapeset**[omp_get_max_threads()];
      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
      {
        spss[i] = new PrecalcShapeset*[wf->get_neq()];
        for (unsigned int j = 0; j < wf->get_neq(); j++)
          spss[i][j] = new PrecalcShapeset(pss[i][j]);
      }
      RefMap*** refmaps = new RefMap**[omp_get_max_threads()];
      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
      {
        refmaps[i] = new RefMap*[wf->get_neq()];
        for (unsigned int j = 0; j < wf->get_neq(); j++)
        {
          refmaps[i][j] = new RefMap();
          refmaps[i][j]->set_quad_2d(&g_quad_2d_std);
        }
      }
      Solution<Scalar>*** u_ext = new Solution<Scalar>**[omp_get_max_threads()];
      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
      {
        if (coeff_vec != NULL)
        {
          u_ext[i] = new Solution<Scalar>*[wf->get_neq()];
          if(i == 0)
          {
            int first_dof = 0;
            for (int j = 0; j < wf->get_neq(); j++)
            {
              u_ext[i][j] = new Solution<Scalar>(spaces[j]->get_mesh());
              Solution<Scalar>::vector_to_solution(coeff_vec, spaces[j], u_ext[i][j], !RungeKutta, first_dof);
              first_dof += spaces[j]->get_num_dofs();
            }
          }
          else
          {
            for (int j = 0; j < wf->get_neq(); j++)
            {
              u_ext[i][j] = new Solution<Scalar>(spaces[j]->get_mesh());
              u_ext[i][j]->copy(u_ext[0][j]);
            }
          }
        }
      }
      AsmList<Scalar>*** als = new AsmList<Scalar>**[omp_get_max_threads()];
      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
      {
        als[i] = new AsmList<Scalar>*[wf->get_neq()];
        for (unsigned int j = 0; j < wf->get_neq(); j++)
          als[i][j] = new AsmList<Scalar>();
      }
      MeshFunction<Scalar>*** ext = new MeshFunction<Scalar>**[omp_get_max_threads()];
      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
      {
        ext[i] = new MeshFunction<Scalar>*[current_stage->ext.size()];
        for (int j = 0; j < current_stage->ext.size(); j++)
          ext[i][j] = current_stage->ext[j]->clone();
      }
      Hermes::vector<MatrixFormVol<Scalar>*>* mfvol = new Hermes::vector<MatrixFormVol<Scalar>*>[omp_get_max_threads()];
      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
      {
        for (int j = 0; j < current_stage->mfvol.size(); j++)
        {
          mfvol[i].push_back(current_stage->mfvol[j]->clone());
          // Inserting proper ext.
          for(int k = 0; k < current_stage->mfvol[j]->ext.size(); k++)
          {
            for (int l = 0; l < current_stage->ext.size(); l++)
            {
              if(current_stage->ext[l] == current_stage->mfvol[j]->ext[j])
              {
                mfvol[i][j]->ext[j] = ext[i][l];
                break;
              }
            }
          }
        }
      }
      Hermes::vector<MatrixFormSurf<Scalar>*>* mfsurf = new Hermes::vector<MatrixFormSurf<Scalar>*>[omp_get_max_threads()];
      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
      {
        for (int j = 0; j < current_stage->mfsurf.size(); j++)
        {
          mfsurf[i].push_back(current_stage->mfsurf[j]->clone());
          // Inserting proper ext.
          for(int k = 0; k < current_stage->mfsurf[j]->ext.size(); k++)
          {
            for (int l = 0; l < current_stage->ext.size(); l++)
            {
              if(current_stage->ext[l] == current_stage->mfsurf[j]->ext[j])
              {
                mfsurf[i][j]->ext[j] = ext[i][l];
                break;
              }
            }
          }
        }
      }
      Hermes::vector<VectorFormVol<Scalar>*>* vfvol = new Hermes::vector<VectorFormVol<Scalar>*>[omp_get_max_threads()];
      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
      {
        for (int j = 0; j < current_stage->vfvol.size(); j++)
        {
          vfvol[i].push_back(current_stage->vfvol[j]->clone());
          // Inserting proper ext.
          for(int k = 0; k < current_stage->vfvol[j]->ext.size(); k++)
          {
            for (int l = 0; l < current_stage->ext.size(); l++)
            {
              if(current_stage->ext[l] == current_stage->vfvol[j]->ext[j])
              {
                vfvol[i][j]->ext[j] = ext[i][l];
                break;
              }
            }
          }
        }
      }
      Hermes::vector<VectorFormSurf<Scalar>*>* vfsurf = new Hermes::vector<VectorFormSurf<Scalar>*>[omp_get_max_threads()];
      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
      {
        for (int j = 0; j < current_stage->vfsurf.size(); j++)
        {
          vfsurf[i].push_back(current_stage->vfsurf[j]->clone());
          // Inserting proper ext.
          for(int k = 0; k < current_stage->vfsurf[j]->ext.size(); k++)
          {
            for (int l = 0; l < current_stage->ext.size(); l++)
            {
              if(current_stage->ext[l] == current_stage->vfsurf[j]->ext[j])
              {
                vfsurf[i][j]->ext[j] = ext[i][l];
                break;
              }
            }
          }
        }
      }

      Traverse trav_master(true);
      unsigned int num_states = trav_master.get_num_states(current_stage->meshes);
      
      trav_master.begin(current_stage->meshes.size(), &(current_stage->meshes.front()));

      Traverse* trav = new Traverse[omp_get_max_threads()];
      Hermes::vector<Transformable *>* fns = new Hermes::vector<Transformable *>[omp_get_max_threads()];
      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
      {
        for (unsigned j = 0; j < current_stage->idx.size(); j++)
          fns[i].push_back(pss[i][current_stage->idx[j]]);
        for (unsigned j = 0; j < current_stage->ext.size(); j++)
        {
          fns[i].push_back(ext[i][j]);
          ext[i][j]->set_quad_2d(&g_quad_2d_std);
        }
        for (unsigned j = 0; j < current_stage->idx.size(); j++)
        {
          fns[i].push_back(u_ext[i][current_stage->idx[j]]);
          if(i == 0)
            current_stage->meshes.push_back(u_ext[i][current_stage->idx[j]]->get_mesh());
          u_ext[i][current_stage->idx[j]]->set_quad_2d(&g_quad_2d_std);
        }
        trav[i].begin(current_stage->meshes.size(), &(current_stage->meshes.front()), &(fns[i].front()));
        trav[i].stack = trav_master.stack;
      }
int state_i;
#define CHUNKSIZE 1
//#pragma omp parallel shared(trav_master) private(state_i)
      {
        //#pragma omp for schedule(dynamic, CHUNKSIZE)
        for(state_i = 0; state_i < num_states; state_i++)
        {
          Traverse::State* current_state;
          #pragma omp critical
          {
            current_state = trav[omp_get_thread_num()].get_next_state(&trav_master.top, &trav_master.id);
          }

          PrecalcShapeset** current_pss = pss[omp_get_thread_num()];
          PrecalcShapeset** current_spss = spss[omp_get_thread_num()];
          RefMap** current_refmaps = refmaps[omp_get_thread_num()];
          Solution<Scalar>** current_u_ext = u_ext[omp_get_thread_num()];
          AsmList<Scalar>** current_als = als[omp_get_thread_num()];

          MatrixFormVol<Scalar>** current_mfvol = mfvol[omp_get_thread_num()].size() == 0 ? NULL : &(mfvol[omp_get_thread_num()].front());
          MatrixFormSurf<Scalar>** current_mfsurf = mfsurf[omp_get_thread_num()].size() == 0 ? NULL : &(mfsurf[omp_get_thread_num()].front());
          VectorFormVol<Scalar>** current_vfvol = vfvol[omp_get_thread_num()].size() == 0 ? NULL : &(vfvol[omp_get_thread_num()].front());
          VectorFormSurf<Scalar>** current_vfsurf = vfsurf[omp_get_thread_num()].size() == 0 ? NULL : &(vfsurf[omp_get_thread_num()].front());

          // One state is a collection of (virtual) elements sharing
          // the same physical location on (possibly) different meshes.
          // This is then the same element of the virtual union mesh.
          // The proper sub-element mappings to all the functions of
          // this stage is supplied by the function Traverse::get_next_state()
          // called in the while loop.
          if(current_state->e != NULL)
            assemble_one_state(current_pss, current_spss, current_refmaps, current_u_ext, current_als, current_state, current_mfvol, current_mfsurf, current_vfvol, current_vfsurf);
        }
      }
      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
      {
        for (unsigned int j = 0; j < wf->get_neq(); j++)
          delete pss[i][j];
        delete [] pss[i];
      }
      delete [] pss;

      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
      {
        for (unsigned int j = 0; j < wf->get_neq(); j++)
          delete spss[i][j];
        delete [] spss[i];
      }
      delete [] spss;

      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
      {
        for (unsigned int j = 0; j < wf->get_neq(); j++)
          delete refmaps[i][j];
        delete [] refmaps[i];
      }
      delete [] refmaps;

      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
      {
        for (unsigned int j = 0; j < wf->get_neq(); j++)
          delete u_ext[i][j];
        delete [] u_ext[i];
      }
      delete [] u_ext;

      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
      {
        for (unsigned int j = 0; j < wf->get_neq(); j++)
          delete als[i][j];
        delete [] als[i];
      }
      delete [] als;

      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
      {
        for (unsigned int j = 0; j < current_stage->ext.size(); j++)
          delete ext[i][j];
        delete [] ext[i];
      }
      delete [] ext;

      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
      {
        for (int j = 0; j < current_stage->mfvol.size(); j++)
          delete mfvol[i][j];
        mfvol[i].clear();
      }
      delete [] mfvol;
      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
      {
        for (int j = 0; j < current_stage->mfsurf.size(); j++)
          delete mfsurf[i][j];
        mfsurf[i].clear();
      }
      delete [] mfsurf;
      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
      {
        for (int j = 0; j < current_stage->vfvol.size(); j++)
          delete vfvol[i][j];
        vfvol[i].clear();
      }
      delete [] vfvol;
      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
      {
        for (int j = 0; j < current_stage->vfsurf.size(); j++)
          delete vfsurf[i][j];
        vfsurf[i].clear();
      }
      delete [] vfsurf;

      /// \todo Should this be really here? Or in assemble()?
      if (current_mat != NULL)
        current_mat->finish();
      if (current_rhs != NULL)
        current_rhs->finish();

      trav_master.finish();
      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
        trav[i].finish();

      if(DG_matrix_forms_present || DG_vector_forms_present)
      {
        Element* element_to_set_nonvisited;
        for(unsigned int mesh_i = 0; mesh_i < current_stage->meshes.size(); mesh_i++)
          for_all_elements(element_to_set_nonvisited, current_stage->meshes[mesh_i])
          element_to_set_nonvisited->visited = false;
      }
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::init_state(PrecalcShapeset** current_pss, PrecalcShapeset** current_spss, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, AsmList<Scalar>** current_als, Traverse::State* current_state)
    {
      _F_;

      // Set maximum integration order for use in integrals, see limit_order()
      update_limit_table(current_state->rep->get_mode());

      // Obtain assembly lists for the element at all spaces of the stage, set appropriate mode for each pss.
      // NOTE: Active elements and transformations for external functions (including the solutions from previous
      // Newton's iteration) as well as basis functions (master PrecalcShapesets) have already been set in
      // trav.get_next_state(...).
      for (unsigned int i = 0; i < current_stage->idx.size(); i++)
      {
        int j = current_stage->idx[i];
        if (current_state->e[j] == NULL)
          continue;

        // \todo do not obtain again if the element was not changed.
        spaces[j]->get_element_assembly_list(current_state->e[j], current_als[j], spaces_first_dofs[j]);

        // Set active element to all test functions.
        current_spss[j]->set_active_element(current_state->e[j]);
        current_spss[j]->set_master_transform();

        // Set active element to reference mappings.
        current_refmaps[j]->set_active_element(current_state->e[j]);
        current_refmaps[j]->force_transform(current_pss[j]->get_transform(), current_pss[j]->get_ctm());

        // Mark the active element on each mesh in order to prevent assembling on its edges from the other side.
        if(DG_matrix_forms_present || DG_vector_forms_present)
          current_state->e[j]->visited = true;
      }
      return;
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::init_surface_state(AsmList<Scalar>** current_als, Traverse::State* current_state)
    {
      _F_;
      // Obtain the list of shape functions which are nonzero on this surface.
      for (unsigned int i = 0; i < current_stage->idx.size(); i++)
      {
        int j = current_stage->idx[i];
        if (current_state->e[j] == NULL)
          continue;
        
        spaces[j]->get_boundary_assembly_list(current_state->e[j], current_state->isurf, current_als[j], spaces_first_dofs[j]);
      }
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::assemble_one_state(PrecalcShapeset** current_pss, PrecalcShapeset** current_spss, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, AsmList<Scalar>** current_als, Traverse::State* current_state,
        MatrixFormVol<Scalar>** current_mfvol, MatrixFormSurf<Scalar>** current_mfsurf, VectorFormVol<Scalar>** current_vfvol, VectorFormSurf<Scalar>** current_vfsurf)
    {
      _F_;
        
      // Initialize the state, return a non-NULL element; if no such element found, return.
      init_state(current_pss, current_spss, current_refmaps, current_u_ext, current_als, current_state);

      init_cache();

      if (current_mat != NULL)
      {
        for(int current_mfvol_i = 0; current_mfvol_i < current_stage->mfvol.size(); current_mfvol_i++)
        {
          if(!form_to_be_assembled(current_mfvol[current_mfvol_i], current_state))
            continue;
        
          Func<double>** base_fns = new Func<double>*[current_als[current_mfvol[current_mfvol_i]->j]->cnt];
          Func<double>** test_fns = new Func<double>*[current_als[current_mfvol[current_mfvol_i]->i]->cnt];
          
          int order = calc_order_matrix_form(current_mfvol[current_mfvol_i], current_refmaps, current_u_ext, current_state);
        
          for (unsigned int i = 0; i < current_als[current_mfvol[current_mfvol_i]->i]->cnt; i++)
          {
            if (std::abs(current_als[current_mfvol[current_mfvol_i]->i]->coef[i]) < 1e-12)
              continue;
            if (current_als[current_mfvol[current_mfvol_i]->i]->dof[i] >= 0)
            {
              current_spss[current_mfvol[current_mfvol_i]->i]->set_active_shape(current_als[current_mfvol[current_mfvol_i]->i]->idx[i]);
              test_fns[i] = get_fn(current_spss[current_mfvol[current_mfvol_i]->i], current_refmaps[current_mfvol[current_mfvol_i]->i], order);
            }
          }

          for (unsigned int j = 0; j < current_als[current_mfvol[current_mfvol_i]->j]->cnt; j++)
          {
            if (std::abs(current_als[current_mfvol[current_mfvol_i]->j]->coef[j]) < 1e-12)
                continue;
            if (current_als[current_mfvol[current_mfvol_i]->j]->dof[j] >= 0)
            {
              current_pss[current_mfvol[current_mfvol_i]->j]->set_active_shape(current_als[current_mfvol[current_mfvol_i]->j]->idx[j]);
              base_fns[j] = get_fn(current_pss[current_mfvol[current_mfvol_i]->j], current_refmaps[current_mfvol[current_mfvol_i]->j], order);
            }
          }

          assemble_matrix_form(current_mfvol[current_mfvol_i], order, base_fns, test_fns, current_refmaps, current_u_ext, current_als, current_state);

          delete [] base_fns;
          delete [] test_fns;
        }
      }
    
      if (current_rhs != NULL)
      {
        for(int current_vfvol_i = 0; current_vfvol_i < current_stage->vfvol.size(); current_vfvol_i++)
        {
          if(!form_to_be_assembled(current_vfvol[current_vfvol_i], current_state))
              continue;
        
          Func<double>** test_fns = new Func<double>*[current_als[current_vfvol[current_vfvol_i]->i]->cnt];
          
          int order = calc_order_vector_form(current_vfvol[current_vfvol_i], current_refmaps, current_u_ext, current_state);
        
          for (unsigned int i = 0; i < current_als[current_vfvol[current_vfvol_i]->i]->cnt; i++)
          {
            if (std::abs(current_als[current_vfvol[current_vfvol_i]->i]->coef[i]) < 1e-12)
              continue;
            if (current_als[current_vfvol[current_vfvol_i]->i]->dof[i] >= 0)
            {
              current_spss[current_vfvol[current_vfvol_i]->i]->set_active_shape(current_als[current_vfvol[current_vfvol_i]->i]->idx[i]);
              test_fns[i] = get_fn(current_spss[current_vfvol[current_vfvol_i]->i], current_refmaps[current_vfvol[current_vfvol_i]->i], order);
            }
          }

          assemble_vector_form(current_vfvol[current_vfvol_i], order, test_fns, current_refmaps, current_u_ext, current_als, current_state);

          delete [] test_fns;
        }
      }

      // Assemble surface integrals now: loop through surfaces of the element.
      for (current_state->isurf = 0; current_state->isurf < current_state->rep->get_num_surf(); current_state->isurf++)
      {
        // \todo DG.
        if(!current_state->bnd[current_state->isurf])
          continue;

        init_surface_state(current_als, current_state);

        if (current_mat != NULL)
        {
          for(int current_mfsurf_i = 0; current_mfsurf_i < current_stage->mfsurf.size(); current_mfsurf_i++)
          {
            if(!form_to_be_assembled(current_mfsurf[current_mfsurf_i], current_state))
              continue;
        
            Func<double>** base_fns = new Func<double>*[current_als[current_mfsurf[current_mfsurf_i]->j]->cnt];
            Func<double>** test_fns = new Func<double>*[current_als[current_mfsurf[current_mfsurf_i]->i]->cnt];
          
            int order = calc_order_matrix_form(current_mfsurf[current_mfsurf_i], current_refmaps, current_u_ext, current_state);
        
            for (unsigned int i = 0; i < current_als[current_mfsurf[current_mfsurf_i]->i]->cnt; i++)
            {
              if (std::abs(current_als[current_mfsurf[current_mfsurf_i]->i]->coef[i]) < 1e-12)
                continue;
              if (current_als[current_mfsurf[current_mfsurf_i]->i]->dof[i] >= 0)
              {
                current_spss[current_mfsurf[current_mfsurf_i]->i]->set_active_shape(current_als[current_mfsurf[current_mfsurf_i]->i]->idx[i]);
                test_fns[i] = get_fn(current_spss[current_mfsurf[current_mfsurf_i]->i], current_refmaps[current_mfsurf[current_mfsurf_i]->i], quad->get_edge_points(current_state->isurf, order));
              }
            }

            for (unsigned int j = 0; j < current_als[current_mfsurf[current_mfsurf_i]->j]->cnt; j++)
            {
              if (std::abs(current_als[current_mfsurf[current_mfsurf_i]->j]->coef[j]) < 1e-12)
                continue;
              if (current_als[current_mfsurf[current_mfsurf_i]->j]->dof[j] >= 0)
              {
                current_pss[current_mfsurf[current_mfsurf_i]->j]->set_active_shape(current_als[current_mfsurf[current_mfsurf_i]->j]->idx[j]);
                base_fns[j] = get_fn(current_pss[current_mfsurf[current_mfsurf_i]->j], current_refmaps[current_mfsurf[current_mfsurf_i]->j], quad->get_edge_points(current_state->isurf, order));
              }
            }

            assemble_matrix_form(current_mfsurf[current_mfsurf_i], order, base_fns, test_fns, current_refmaps, current_u_ext, current_als, current_state);

            delete [] base_fns;
            delete [] test_fns;
          }
        }
    
        if (current_rhs != NULL)
        {
          for(int current_vfsurf_i = 0; current_vfsurf_i < current_stage->vfsurf.size(); current_vfsurf_i++)
          {
            if(!form_to_be_assembled(current_vfsurf[current_vfsurf_i], current_state))
                continue;
        
            Func<double>** test_fns = new Func<double>*[current_als[current_vfsurf[current_vfsurf_i]->i]->cnt];
          
            int order = calc_order_vector_form(current_vfsurf[current_vfsurf_i], current_refmaps, current_u_ext, current_state);
        
            for (unsigned int i = 0; i < current_als[current_vfsurf[current_vfsurf_i]->i]->cnt; i++)
            {
              if (std::abs(current_als[current_vfsurf[current_vfsurf_i]->i]->coef[i]) < 1e-12)
                continue;
              if (current_als[current_vfsurf[current_vfsurf_i]->i]->dof[i] >= 0)
              {
                current_spss[current_vfsurf[current_vfsurf_i]->i]->set_active_shape(current_als[current_vfsurf[current_vfsurf_i]->i]->idx[i]);
                test_fns[i] = get_fn(current_spss[current_vfsurf[current_vfsurf_i]->i], current_refmaps[current_vfsurf[current_vfsurf_i]->i], quad->get_edge_points(current_state->isurf, order));
              }
            }

            assemble_vector_form(current_vfsurf[current_vfsurf_i], order, test_fns, current_refmaps, current_u_ext, current_als, current_state);

            delete [] test_fns;
          }
        }
      }

      delete_cache();
    }

    template<typename Scalar>
    int DiscreteProblem<Scalar>::calc_order_matrix_form(MatrixForm<Scalar> *form, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, Traverse::State* current_state)
    {
      _F_;

      int order;

      if(is_fvm)
        order = current_refmaps[form->i]->get_inv_ref_order();
      else
      {
        // order of solutions from the previous Newton iteration etc..
        Func<Hermes::Ord>** u_ext_ord = new Func<Hermes::Ord>*[RungeKutta ? RK_original_spaces_count : this->wf->get_neq() - form->u_ext_offset];
        ExtData<Hermes::Ord> ext_ord;
        init_ext_orders(form, u_ext_ord, &ext_ord, current_u_ext);

        // Order of shape functions.
        Func<Hermes::Ord>* ou = get_fn_ord(this->spaces[form->j]->get_element_order(current_state->e[form->j]->id));
        Func<Hermes::Ord>* ov = get_fn_ord(this->spaces[form->i]->get_element_order(current_state->e[form->i]->id));

        // Total order of the vector form.
        Hermes::Ord o = form->ord(1, &fake_wt, u_ext_ord, ou, ov, &geom_ord, &ext_ord);

        adjust_order_to_refmaps(form, order, &o, current_refmaps);

        // Cleanup.
        deinit_ext_orders(form, u_ext_ord, &ext_ord);
      }
      return order;
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::assemble_matrix_form(MatrixForm<Scalar>* form, int order, Func<double>** base_fns, Func<double>** test_fns, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, AsmList<Scalar>** current_als, Traverse::State* current_state)
    {
      _F_;
      bool surface_form = (dynamic_cast<MatrixFormVol<Scalar>*>(form) == NULL);

      double block_scaling_coef = this->block_scaling_coeff(form);

      bool tra = (form->i != form->j) && (form->sym != 0);
      bool sym = (form->i == form->j) && (form->sym == 1);

      // Assemble the local stiffness matrix for the form form.
      Scalar **local_stiffness_matrix = get_matrix_buffer(std::max(current_als[form->i]->cnt, current_als[form->j]->cnt));

      // Init external functions.
      Func<Scalar>** u_ext = new Func<Scalar>*[RungeKutta ? RK_original_spaces_count : this->wf->get_neq() - form->u_ext_offset];
      ExtData<Scalar> ext;
      init_ext(form, u_ext, &ext, order, current_u_ext, current_state);
      
      // Add the previous time level solution previously inserted at the back of ext.
      if(RungeKutta)
        for(int ext_i = 0; ext_i < this->RK_original_spaces_count; ext_i++)
          u_ext[ext_i]->add(*ext.fn[form->ext.size() - this->RK_original_spaces_count + ext_i]);

      // Init geometry.
      int n_quadrature_points;
      if(surface_form)
        n_quadrature_points = init_surface_geometry_points(current_refmaps[form->i], order, current_state);
      else
        n_quadrature_points = init_geometry_points(current_refmaps[form->i], order);

      // Actual form-specific calculation.
      for (unsigned int i = 0; i < current_als[form->i]->cnt; i++)
      {
        if (current_als[form->i]->dof[i] < 0)
          continue;

        if ((!tra || surface_form) && current_als[form->i]->dof[i] < 0) 
          continue;
        if(std::abs(current_als[form->i]->coef[i]) < 1e-12)
          continue;
        if (!sym)
        {
          for (unsigned int j = 0; j < current_als[form->j]->cnt; j++)
          {
            if (current_als[form->j]->dof[j] >= 0)
            {
              // Is this necessary, i.e. is there a coefficient smaller than 1e-12?
              if (std::abs(current_als[form->j]->coef[j]) < 1e-12)
                continue;
              
              Func<double>* u = base_fns[j];
              Func<double>* v = test_fns[i];

              if(surface_form)
                local_stiffness_matrix[i][j] = 0.5 * block_scaling_coeff(form) * form->value(n_quadrature_points, jacobian_x_weights_cache[order], u_ext, u, v, geometry_cache[order], &ext) * form->scaling_factor * current_als[form->j]->coef[j] * current_als[form->i]->coef[i];
              else
                local_stiffness_matrix[i][j] = block_scaling_coeff(form) * form->value(n_quadrature_points, jacobian_x_weights_cache[order], u_ext, u, v, geometry_cache[order], &ext) * form->scaling_factor * current_als[form->j]->coef[j] * current_als[form->i]->coef[i];
            }
          }
        }
        // Symmetric block.
        else
        {
          for (unsigned int j = 0; j < current_als[form->j]->cnt; j++)
          {
            if (j < i && current_als[form->j]->dof[j] >= 0)
              continue;
            if (current_als[form->j]->dof[j] >= 0)
            {
              // Is this necessary, i.e. is there a coefficient smaller than 1e-12?
              if (std::abs(current_als[form->j]->coef[j]) < 1e-12)
                continue;

              Func<double>* u = base_fns[j];
              Func<double>* v = test_fns[i];

              Scalar val = block_scaling_coeff(form) * form->value(n_quadrature_points, jacobian_x_weights_cache[order], u_ext, u, v, geometry_cache[order], &ext) * form->scaling_factor * current_als[form->j]->coef[j] * current_als[form->i]->coef[i];

              local_stiffness_matrix[i][j] = local_stiffness_matrix[j][i] = val;
            }
          }
        }
      }

      // Insert the local stiffness matrix into the global one.
      current_mat->add(current_als[form->i]->cnt, current_als[form->j]->cnt, local_stiffness_matrix, current_als[form->i]->dof, current_als[form->j]->dof);

      // Insert also the off-diagonal (anti-)symmetric block, if required.
      if (tra)
      {
        if (form->sym < 0)
          chsgn(local_stiffness_matrix, current_als[form->i]->cnt, current_als[form->j]->cnt);
        transpose(local_stiffness_matrix, current_als[form->i]->cnt, current_als[form->j]->cnt);
        current_mat->add(current_als[form->j]->cnt, current_als[form->i]->cnt, local_stiffness_matrix, current_als[form->j]->dof, current_als[form->i]->dof);
      }

      // Cleanup.
      deinit_ext(form, u_ext, &ext);
    }
    
    template<typename Scalar>
    int DiscreteProblem<Scalar>::calc_order_vector_form(VectorForm<Scalar> *form, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, Traverse::State* current_state)
    {
      _F_;

      int order;

      if(is_fvm)
        order = current_refmaps[form->i]->get_inv_ref_order();
      else
      {
        // order of solutions from the previous Newton iteration etc..
        Func<Hermes::Ord>** u_ext_ord = new Func<Hermes::Ord>*[RungeKutta ? RK_original_spaces_count : this->wf->get_neq() - form->u_ext_offset];
        ExtData<Hermes::Ord> ext_ord;
        init_ext_orders(form, u_ext_ord, &ext_ord, current_u_ext);

        // Order of shape functions.
        Func<Hermes::Ord>* ov = get_fn_ord(this->spaces[form->i]->get_element_order(current_state->e[form->i]->id));

        // Total order of the vector form.
        Hermes::Ord o = form->ord(1, &fake_wt, u_ext_ord, ov, &geom_ord, &ext_ord);

        adjust_order_to_refmaps(form, order, &o, current_refmaps);

        // Cleanup.
        deinit_ext_orders(form, u_ext_ord, &ext_ord);
      }
      return order;
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::assemble_vector_form(VectorForm<Scalar>* form, int order, Func<double>** test_fns, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, AsmList<Scalar>** current_als, Traverse::State* current_state)
    {
      _F_;
      bool surface_form = (dynamic_cast<VectorFormVol<Scalar>*>(form) == NULL);

      // Init geometry.
      int n_quadrature_points;
      if(surface_form)
        n_quadrature_points = init_surface_geometry_points(current_refmaps[form->i], order, current_state);
      else
        n_quadrature_points = init_geometry_points(current_refmaps[form->i], order);

      // Init external functions.
      Func<Scalar>** u_ext = new Func<Scalar>*[RungeKutta ? RK_original_spaces_count : this->wf->get_neq() - form->u_ext_offset];
      ExtData<Scalar> ext;
      init_ext(form, u_ext, &ext, order, current_u_ext, current_state);
      
      // Add the previous time level solution previously inserted at the back of ext.
      if(RungeKutta)
        for(int ext_i = 0; ext_i < this->RK_original_spaces_count; ext_i++)
          u_ext[ext_i]->add(*ext.fn[form->ext.size() - this->RK_original_spaces_count + ext_i]);

      // Actual form-specific calculation.
      for (unsigned int i = 0; i < current_als[form->i]->cnt; i++)
      {
        if (current_als[form->i]->dof[i] < 0)
          continue;

        // Is this necessary, i.e. is there a coefficient smaller than 1e-12?
        if (std::abs(current_als[form->i]->coef[i]) < 1e-12)
          continue;

        Func<double>* v = test_fns[i];
        
        if(surface_form)
          current_rhs->add(current_als[form->i]->dof[i], 0.5 * form->value(n_quadrature_points, jacobian_x_weights_cache[order], u_ext, v, geometry_cache[order], &ext) * form->scaling_factor * current_als[form->i]->coef[i]);
        else
          current_rhs->add(current_als[form->i]->dof[i], form->value(n_quadrature_points, jacobian_x_weights_cache[order], u_ext, v, geometry_cache[order], &ext) * form->scaling_factor * current_als[form->i]->coef[i]);
      }

      // Cleanup.
      deinit_ext(form, u_ext, &ext);
    }

    template<typename Scalar>
    int DiscreteProblem<Scalar>::init_geometry_points(RefMap* reference_mapping, int order)
    {
      _F_;
      double3* pt = quad->get_points(order);
      int np = quad->get_num_points(order);

      // Init geometry and jacobian*weights.
      if (geometry_cache[order] == NULL)
      {
        geometry_cache[order] = init_geom_vol(reference_mapping, order);
        double* jac = NULL;
        if(!reference_mapping->is_jacobian_const())
          jac = reference_mapping->get_jacobian(order);
        jacobian_x_weights_cache[order] = new double[np];
        for(int i = 0; i < np; i++)
        {
          if(reference_mapping->is_jacobian_const())
            jacobian_x_weights_cache[order][i] = pt[i][2] * reference_mapping->get_const_jacobian();
          else
            jacobian_x_weights_cache[order][i] = pt[i][2] * jac[i];
        }
      }
      Geom<double>* e = geometry_cache[order];
      double* jwt = jacobian_x_weights_cache[order];
      return np;
    }

    template<typename Scalar>
    int DiscreteProblem<Scalar>::init_surface_geometry_points(RefMap* reference_mapping, int& order, Traverse::State* current_state)
    {
      _F_;
      int eo = quad->get_edge_points(current_state->isurf, order);
      double3* pt = quad->get_points(eo);
      int np = quad->get_num_points(eo);

      // Init geometry and jacobian*weights.
      if (geometry_cache[eo] == NULL)
      {
        geometry_cache[eo] = init_geom_surf(reference_mapping, current_state->isurf, current_state->rep->marker, eo);
        double3* tan = reference_mapping->get_tangent(current_state->isurf, eo);
        jacobian_x_weights_cache[eo] = new double[np];
        for(int i = 0; i < np; i++)
          jacobian_x_weights_cache[eo][i] = pt[i][2] * tan[i][2];
      }
      Geom<double>* e = geometry_cache[eo];
      double* jwt = jacobian_x_weights_cache[eo];

      order = eo;
      return np;
    }
    
    template<typename Scalar>
    void DiscreteProblem<Scalar>::init_ext_orders(Form<Scalar> *form, Func<Hermes::Ord>** oi, ExtData<Hermes::Ord>* oext, Solution<Scalar>** current_u_ext)
    {
      _F_;
      unsigned int prev_size = RungeKutta ? RK_original_spaces_count : this->wf->get_neq() - form->u_ext_offset;

      if (current_u_ext != NULL)
        for(int i = 0; i < prev_size; i++)
          if (current_u_ext[i + form->u_ext_offset] != NULL)
            oi[i] = get_fn_ord(current_u_ext[i + form->u_ext_offset]->get_fn_order());
          else
            oi[i] = get_fn_ord(0);
      else
        for(int i = 0; i < prev_size; i++)
          oi[i] = get_fn_ord(0);
      
      oext->nf = form->ext.size();
      oext->fn = new Func<Hermes::Ord>*[oext->nf];
      for (int i = 0; i < oext->nf; i++)
        oext->fn[i] = get_fn_ord(form->ext[i]->get_fn_order());
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::deinit_ext_orders(Form<Scalar> *form, Func<Hermes::Ord>** oi, ExtData<Hermes::Ord>* oext)
    {
      _F_;
      delete [] oi;

      if (oext != NULL)
        oext->free_ord();
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::init_ext(Form<Scalar> *form, Func<Scalar>** u_ext, ExtData<Scalar>* ext, int order, Solution<Scalar>** current_u_ext, Traverse::State* current_state)
    {
      _F_;
      unsigned int prev_size = RungeKutta ? RK_original_spaces_count : this->wf->get_neq() - form->u_ext_offset;
      
      if (current_u_ext != NULL)
        for(int i = 0; i < prev_size; i++)
          if (current_u_ext[i + form->u_ext_offset] != NULL)
            u_ext[i] = current_state->e[i] == NULL ? NULL : init_fn(current_u_ext[i + form->u_ext_offset], order);
          else
            u_ext[i] = NULL;
      else
        for(int i = 0; i < prev_size; i++)
          u_ext[i] = NULL;
      
      ext->nf = form->ext.size();
      ext->fn = new Func<Scalar>*[ext->nf];
      for (unsigned i = 0; i < ext->nf; i++)
      {
        if (form->ext[i] != NULL) 
          ext->fn[i] = init_fn(form->ext[i], order);
        else ext->fn[i] = NULL;
      }
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::deinit_ext(Form<Scalar> *form, Func<Scalar>** u_ext, ExtData<Scalar>* ext)
    {
      _F_;
      // Values of the previous Newton iteration, shape functions
      // and external functions in quadrature points.
      int prev_size = this->wf->get_neq() - form->u_ext_offset;
      // In case of Runge-Kutta, this is time-saving, as it is known how many functions are there for the user.
      if(this->RungeKutta)
        prev_size = this->RK_original_spaces_count;
      
      for(int i = 0; i < prev_size; i++)
      if (u_ext[i] != NULL)
      {
        u_ext[i]->free_fn();
        delete u_ext[i];
      }

      delete [] u_ext;

      if (ext != NULL)
        ext->free();
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::adjust_order_to_refmaps(Form<Scalar> *form, int& order, Hermes::Ord* o, RefMap** current_refmaps)
    {
      _F_;
      // Increase due to reference map.
      order = current_refmaps[form->i]->get_inv_ref_order();
      order += o->get_order();
      limit_order(order);
    }

    template<typename Scalar>
    Func<double>* DiscreteProblem<Scalar>::get_fn(PrecalcShapeset *fu, RefMap *rm, const int order)
    {
      _F_;

      if(rm->is_jacobian_const())
      {
        typename AssemblingCaches::KeyConst key(256 - fu->get_active_shape(), order, fu->get_transform(), fu->get_shapeset()->get_id(), rm->get_const_inv_ref_map());
        if(rm->get_active_element()->get_mode() == HERMES_MODE_TRIANGLE)
        {
          if(assembling_caches.const_cache_fn_triangles.find(key) == assembling_caches.const_cache_fn_triangles.end())
            assembling_caches.const_cache_fn_triangles[key] = init_fn(fu, rm, order);
          return assembling_caches.const_cache_fn_triangles[key];
        }
        else
        {
          if(assembling_caches.const_cache_fn_quads.find(key) == assembling_caches.const_cache_fn_quads.end())
            assembling_caches.const_cache_fn_quads[key] = init_fn(fu, rm, order);
          return assembling_caches.const_cache_fn_quads[key];
        }
      }
      else
      {
        typename AssemblingCaches::KeyNonConst key(256 - fu->get_active_shape(), order,
          fu->get_transform(), fu->get_shapeset()->get_id());
        if(rm->get_active_element()->get_mode() == HERMES_MODE_TRIANGLE)
        {
          if(assembling_caches.cache_fn_triangles.find(key) == assembling_caches.cache_fn_triangles.end())
            assembling_caches.cache_fn_triangles[key] = init_fn(fu, rm, order);
          return assembling_caches.cache_fn_triangles[key];
        }
        else
        {
          if(assembling_caches.cache_fn_quads.find(key) == assembling_caches.cache_fn_quads.end())
            assembling_caches.cache_fn_quads[key] = init_fn(fu, rm, order);
          return assembling_caches.cache_fn_quads[key];
        }
      }
    }

    template<typename Scalar>
    Func<Hermes::Ord>* DiscreteProblem<Scalar>::get_fn_ord(const int order)
    {
      _F_;
      assert(order >= 0);
      unsigned int cached_order = (unsigned int) order;
      if(!assembling_caches.cache_fn_ord.present(cached_order))
        assembling_caches.cache_fn_ord.add(init_fn_ord(cached_order), cached_order);
      return assembling_caches.cache_fn_ord.get(cached_order);
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::init_cache()
    {
      _F_;
      for (int i = 0; i < g_max_quad + 1 + 4 * g_max_quad + 4; i++)
      {
        geometry_cache[i] = NULL;
        jacobian_x_weights_cache[i] = NULL;
      }
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::delete_single_geom_cache(int order)
    {
      if (geometry_cache[order] != NULL)
      {
        geometry_cache[order]->free();
        delete geometry_cache[order];
        geometry_cache[order] = NULL;
        delete [] jacobian_x_weights_cache[order];
      }
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::delete_cache()
    {
      _F_;
      for (int i = 0; i < g_max_quad + 1 + 4 * g_max_quad + 4; i++)
      {
        if (geometry_cache[i] != NULL)
        {
          geometry_cache[i]->free(); delete geometry_cache[i];
          delete [] jacobian_x_weights_cache[i];
        }
      }

      for (typename std::map<typename AssemblingCaches::KeyNonConst, Func<double>*, typename AssemblingCaches::CompareNonConst>::const_iterator it = assembling_caches.cache_fn_quads.begin();
        it != assembling_caches.cache_fn_quads.end(); it++)
      {
        (it->second)->free_fn(); delete (it->second);
      }
      assembling_caches.cache_fn_quads.clear();

      for (typename std::map<typename AssemblingCaches::KeyNonConst, Func<double>*, typename AssemblingCaches::CompareNonConst>::const_iterator it = assembling_caches.cache_fn_triangles.begin();
        it != assembling_caches.cache_fn_triangles.end(); it++)
      {
        (it->second)->free_fn();
        delete (it->second);
      }
      assembling_caches.cache_fn_triangles.clear();
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::set_quad_2d(Quad2D* quad)
    {
      this->quad = quad;
    }

    template<typename Scalar>
    DiscreteProblem<Scalar>::AssemblingCaches::AssemblingCaches()
    {
    };

    template<typename Scalar>
    DiscreteProblem<Scalar>::AssemblingCaches::~AssemblingCaches()
    {
      _F_;
      for (typename std::map<KeyConst, Func<double>*, CompareConst>::const_iterator it = const_cache_fn_triangles.begin();
        it != const_cache_fn_triangles.end(); it++)
      {
        (it->second)->free_fn(); delete (it->second);
      }
      const_cache_fn_triangles.clear();

      for (typename std::map<KeyConst, Func<double>*, CompareConst>::const_iterator it = const_cache_fn_quads.begin();
        it != const_cache_fn_quads.end(); it++)
      {
        (it->second)->free_fn(); delete (it->second);
      }
      const_cache_fn_quads.clear();

      for(unsigned int i = 0; i < cache_fn_ord.get_size(); i++)
        if(cache_fn_ord.present(i))
        {
          cache_fn_ord.get(i)->free_ord();
          delete cache_fn_ord.get(i);
        }
    };

#ifdef _MSC_VER
    template<typename Scalar>
    DiscreteProblem<Scalar>::AssemblingCaches::KeyConst::KeyConst(int index, int order, UINT64 sub_idx, int shapeset_type, double2x2* inv_ref_map)
    {
      this->index = index;
      this->order = order;
      this->sub_idx = sub_idx;
      this->shapeset_type = shapeset_type;
      this->inv_ref_map[0][0] = (* inv_ref_map)[0][0];
      this->inv_ref_map[0][1] = (* inv_ref_map)[0][1];
      this->inv_ref_map[1][0] = (* inv_ref_map)[1][0];
      this->inv_ref_map[1][1] = (* inv_ref_map)[1][1];
    }
#else
    template<typename Scalar>
    DiscreteProblem<Scalar>::AssemblingCaches::KeyConst::KeyConst(int index, int order, unsigned int sub_idx, int shapeset_type, double2x2* inv_ref_map)
    {
      this->index = index;
      this->order = order;
      this->sub_idx = sub_idx;
      this->shapeset_type = shapeset_type;
      this->inv_ref_map[0][0] = (* inv_ref_map)[0][0];
      this->inv_ref_map[0][1] = (* inv_ref_map)[0][1];
      this->inv_ref_map[1][0] = (* inv_ref_map)[1][0];
      this->inv_ref_map[1][1] = (* inv_ref_map)[1][1];
    }
#endif

    template<typename Scalar>
    bool DiscreteProblem<Scalar>::AssemblingCaches::CompareConst::operator()(KeyConst a, KeyConst b) const
    {
      if(a.inv_ref_map[0][0] < b.inv_ref_map[0][0]) return true;
      else if(a.inv_ref_map[0][0] > b.inv_ref_map[0][0]) return false;
      else
        if(a.inv_ref_map[0][1] < b.inv_ref_map[0][1]) return true;
        else if(a.inv_ref_map[0][1] > b.inv_ref_map[0][1]) return false;
        else
          if(a.inv_ref_map[1][0] < b.inv_ref_map[1][0]) return true;
          else if(a.inv_ref_map[1][0] > b.inv_ref_map[1][0]) return false;
          else
            if(a.inv_ref_map[1][1] < b.inv_ref_map[1][1]) return true;
            else if(a.inv_ref_map[1][1] > b.inv_ref_map[1][1]) return false;
            else
              if (a.index < b.index) return true;
              else if (a.index > b.index) return false;
              else
                if (a.order < b.order) return true;
                else if (a.order > b.order) return false;
                else
                  if (a.sub_idx < b.sub_idx) return true;
                  else if (a.sub_idx > b.sub_idx) return false;
                  else
                    if (a.shapeset_type < b.shapeset_type) return true;
                    else return false;
    }

#ifdef _MSC_VER
    template<typename Scalar>
    DiscreteProblem<Scalar>::AssemblingCaches::KeyNonConst::KeyNonConst(int index, int order, UINT64 sub_idx, int shapeset_type)
    {
      this->index = index;
      this->order = order;
      this->sub_idx = sub_idx;
      this->shapeset_type = shapeset_type;
    }
#else
    template<typename Scalar>
    DiscreteProblem<Scalar>::AssemblingCaches::KeyNonConst::KeyNonConst(int index, int order, unsigned int sub_idx, int shapeset_type)
    {
      this->index = index;
      this->order = order;
      this->sub_idx = sub_idx;
      this->shapeset_type = shapeset_type;
    }
#endif

    template<typename Scalar> bool DiscreteProblem<Scalar>::AssemblingCaches::CompareNonConst::operator()(KeyNonConst a, KeyNonConst b) const
    {
      if (a.index < b.index) return true;
      else if (a.index > b.index) return false;
      else
      {
        if (a.order < b.order) return true;
        else if (a.order > b.order) return false;
        else
        {
          if (a.sub_idx < b.sub_idx) return true;
          else if (a.sub_idx > b.sub_idx) return false;
          else
          {
            if (a.shapeset_type < b.shapeset_type) return true;
            else return false;
          }
        }
      }
    }

    template class HERMES_API DiscreteProblem<double>;
    template class HERMES_API DiscreteProblem<std::complex<double> >;
  }
}
