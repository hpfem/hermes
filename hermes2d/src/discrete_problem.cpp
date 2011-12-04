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
    DiscreteProblem<Scalar>::DiscreteProblem() : wf(NULL), pss(NULL), current_stage(NULL)
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

      // Initialize precalc shapesets according to spaces provided.
      pss = new PrecalcShapeset*[wf->get_neq()];

      for (unsigned int i = 0; i < wf->get_neq(); i++)
      {
        pss[i] = NULL;
        Shapeset *shapeset = spaces[i]->shapeset;
        if (shapeset == NULL) error("Internal in DiscreteProblem<Scalar>::init_spaces().");
        PrecalcShapeset *p = new PrecalcShapeset(shapeset);
        if (p == NULL) error("New PrecalcShapeset could not be allocated in DiscreteProblem<Scalar>::init_spaces().");
        pss[i] = p;
      }

      // There is a special function that sets a DiscreteProblem to be FVM.
      // Purpose is that this constructor looks cleaner and is simpler.
      this->is_fvm = false;

      Geom<Hermes::Ord> *tmp = init_geom_ord();
      geom_ord = *tmp;
      delete tmp;
      quad = &g_quad_2d_std;

      current_al.resize(wf->get_neq());

      current_mat = NULL;
      current_stage = NULL;
      current_rhs = NULL;
      current_block_weights = NULL;
      current_isurf = -1;
    }

    template<typename Scalar>
    DiscreteProblem<Scalar>::~DiscreteProblem()
    {
      _F_;
      if (wf != NULL)
        memset(sp_seq, -1, sizeof(int) * wf->get_neq());
      wf_seq = -1;
      if (sp_seq != NULL) delete [] sp_seq;
      if (pss != NULL)
      {
        for(unsigned int i = 0; i < wf->get_neq(); i++)
          delete pss[i];
        delete [] pss;
      }
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
    PrecalcShapeset* DiscreteProblem<Scalar>::get_pss(int n)
    {
      return this->pss[n];
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
    bool DiscreteProblem<Scalar>::form_to_be_assembled(MatrixForm<Scalar>* form)
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
    bool DiscreteProblem<Scalar>::form_to_be_assembled(MatrixFormVol<Scalar>* form)
    {
      if(!form_to_be_assembled((MatrixForm<Scalar>*)form))
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
    bool DiscreteProblem<Scalar>::form_to_be_assembled(MatrixFormSurf<Scalar>* form)
    {
      if(current_state->rep->en[this->current_isurf]->marker == 0)
        return false;

      if (form->areas[0] == H2D_DG_INNER_EDGE)
        return false;
      if(!form_to_be_assembled((MatrixForm<Scalar>*)form))
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
            marker_on_space_m = (this->spaces[form->i]->get_mesh()->get_boundary_markers_conversion().get_internal_marker(form->areas[ss]).marker == current_state->rep->en[this->current_isurf]->marker);

          bool marker_on_space_n = this->spaces[form->j]->get_mesh()->get_boundary_markers_conversion().get_internal_marker(form->areas[ss]).valid;
          if(marker_on_space_n)
            marker_on_space_n = (this->spaces[form->j]->get_mesh()->get_boundary_markers_conversion().get_internal_marker(form->areas[ss]).marker == current_state->rep->en[this->current_isurf]->marker);

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
    bool DiscreteProblem<Scalar>::form_to_be_assembled(VectorForm<Scalar>* form)
    {
      if (current_state->e[form->i] == NULL)
        return false;
      if (fabs(form->scaling_factor) < 1e-12)
        return false;

      return true;
    }
    
    template<typename Scalar>
    bool DiscreteProblem<Scalar>::form_to_be_assembled(VectorFormVol<Scalar>* form)
    {
      if(!form_to_be_assembled((VectorForm<Scalar>*)form))
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
    bool DiscreteProblem<Scalar>::form_to_be_assembled(VectorFormSurf<Scalar>* form)
    {
      if(current_state->rep->en[this->current_isurf]->marker == 0)
        return false;

      if (form->areas[0] == H2D_DG_INNER_EDGE)
        return false;
      
      if(!form_to_be_assembled((VectorForm<Scalar>*)form))
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
            marker_on_space_m = (this->spaces[form->i]->get_mesh()->get_boundary_markers_conversion().get_internal_marker(form->areas[ss]).marker == current_state->rep->en[this->current_isurf]->marker);

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
    void DiscreteProblem<Scalar>::init_psss()
    {
      _F_;
      for (unsigned int i = 0; i < wf->get_neq(); i++)
      {
        current_spss.push_back(new PrecalcShapeset(pss[i]));
        current_spss[i]->set_quad_2d(&g_quad_2d_std);
      }
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::init_refmaps()
    {
      _F_;
      for (unsigned int i = 0; i < wf->get_neq(); i++)
      {
        current_refmap.push_back(new RefMap());
        current_refmap[i]->set_quad_2d(&g_quad_2d_std);
      }
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::init_u_ext(Scalar* coeff_vec)
    {
      _F_;
      int first_dof = 0;
      if (coeff_vec != NULL)
        for (int i = 0; i < wf->get_neq(); i++)
        {
          Solution<Scalar>* external_solution_i = new Solution<Scalar>(spaces[i]->get_mesh());
          Solution<Scalar>::vector_to_solution(coeff_vec, spaces[i], external_solution_i, true, first_dof);
          current_u_ext.push_back(external_solution_i);
          first_dof += spaces[i]->get_num_dofs();
        }
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::deinit_psss()
    {
      _F_;
      for (unsigned int i = 0; i < wf->get_neq(); i++)
        delete current_spss[i];
      current_spss.clear();
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::deinit_refmaps()
    {
      _F_;
      for (unsigned int i = 0; i < wf->get_neq(); i++)
        delete current_refmap[i];
      current_refmap.clear();
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::deinit_u_ext()
    {
      _F_;
      for(typename Hermes::vector<Solution<Scalar>*>::iterator it = current_u_ext.begin(); it != current_u_ext.end(); it++)
        delete *it;
      current_u_ext.clear();
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

      // Convert the coefficient vector into vector of external solutions.

      // Reset the warnings about insufficiently high integration order.
      reset_warn_order();

      // Initialize slave pss's, refmaps, u_ext.
      init_psss();
      init_refmaps();
      init_u_ext(coeff_vec);

      // Initialize matrix buffer.
      matrix_buffer = NULL;
      matrix_buffer_dim = 0;
      if (mat != NULL)
        get_matrix_buffer(9);

      // Create assembling stages.
      Hermes::vector<Stage<Scalar> > stages = Hermes::vector<Stage<Scalar> >();
      bool want_matrix = (mat != NULL);
      bool want_vector = (rhs != NULL);
      wf->get_stages(spaces, current_u_ext, stages, want_matrix, want_vector);

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

        // Assemble one current_stage-> One stage is a collection of functions,
        // and meshes that can not be further minimized.
        // E.g. if a linear form uses two external solutions, each of
        // which is defined on a different mesh, and different to the
        // mesh of the current test function, then the stage would have
        // three meshes. By stage functions, all functions are meant: shape
        // functions (their precalculated values), and mesh functions.
        assemble_one_stage();
      }
      // Deinitialize matrix buffer.
      if(matrix_buffer != NULL)
        delete [] matrix_buffer;
      matrix_buffer = NULL;
      matrix_buffer_dim = 0;

      deinit_psss();
      deinit_refmaps();
      deinit_u_ext();
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
    void DiscreteProblem<Scalar>::assemble_one_stage()
    {
      _F_;

      // Create the assembling states.
      Traverse trav;
      for (unsigned i = 0; i < current_stage->idx.size(); i++)
        current_stage->fns[i] = pss[current_stage->idx[i]];
      for (unsigned i = 0; i < current_stage->ext.size(); i++)
        current_stage->ext[i]->set_quad_2d(&g_quad_2d_std);
      trav.begin(current_stage->meshes.size(), &(current_stage->meshes.front()), &(current_stage->fns.front()));

      // Loop through all assembling states.
      // Assemble each one.
      while ((current_state = trav.get_next_state()) != NULL)
      {
        // One state is a collection of (virtual) elements sharing
        // the same physical location on (possibly) different meshes.
        // This is then the same element of the virtual union mesh.
        // The proper sub-element mappings to all the functions of
        // this stage is supplied by the function Traverse::get_next_state()
        // called in the while loop.
        if(current_state->e != NULL)
          assemble_one_state();
      }

      if (current_mat != NULL)
        current_mat->finish();
      if (current_rhs != NULL)
        current_rhs->finish();
      trav.finish();

      if(DG_matrix_forms_present || DG_vector_forms_present)
      {
        Element* element_to_set_nonvisited;
        for(unsigned int mesh_i = 0; mesh_i < current_stage->meshes.size(); mesh_i++)
          for_all_elements(element_to_set_nonvisited, current_stage->meshes[mesh_i])
          element_to_set_nonvisited->visited = false;
      }
    }

    template<typename Scalar>
    Element* DiscreteProblem<Scalar>::init_state()
    {
      _F_;
      // Find a non-NULL e[i].
      Element* e0 = NULL;
      for (unsigned int i = 0; i < current_stage->idx.size(); i++)
        if ((e0 = current_state->e[i]) != NULL)
          break;
      if(e0 == NULL)
        return NULL;

      // Set maximum integration order for use in integrals, see limit_order()
      update_limit_table(e0->get_mode());

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
        spaces[j]->get_element_assembly_list(current_state->e[j], current_al[j], spaces_first_dofs[j]);

        // Set active element to all test functions.
        current_spss[j]->set_active_element(current_state->e[j]);
        current_spss[j]->set_master_transform();

        // Set active element to reference mappings.
        current_refmap[j]->set_active_element(current_state->e[j]);
        current_refmap[j]->force_transform(pss[j]->get_transform(), pss[j]->get_ctm());

        // Mark the active element on each mesh in order to prevent assembling on its edges from the other side.
        if(DG_matrix_forms_present || DG_vector_forms_present)
          current_state->e[j]->visited = true;
      }
      return e0;
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::init_surface_state()
    {
      _F_;
      // Obtain the list of shape functions which are nonzero on this surface.
      for (unsigned int i = 0; i < current_stage->idx.size(); i++)
      {
        int j = current_stage->idx[i];
        if (current_state->e[j] == NULL)
          continue;
        
        spaces[j]->get_boundary_assembly_list(current_state->e[j], current_isurf, current_al[j], spaces_first_dofs[j]);
      }
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::assemble_one_state()
    {
      _F_;

      // Assembly list vector.
      for(unsigned int i = 0; i < wf->get_neq(); i++)
        current_al[i] = new AsmList<Scalar>;
        
      // Initialize the state, return a non-NULL element; if no such element found, return.
      Element* rep_element = init_state();

      init_cache();

      if (current_mat != NULL)
      {
        for(typename Hermes::vector<MatrixFormVol<Scalar> *>::iterator it = current_stage->mfvol.begin(); it != current_stage->mfvol.end(); it++)
        {
          if(!form_to_be_assembled(*it))
            continue;
        
          int order = calc_order_matrix_form(*it);
        
          assemble_matrix_form(*it, order);
        }
      }
    
      if (current_rhs != NULL)
      {
        for(typename Hermes::vector<VectorFormVol<Scalar> *>::iterator it = current_stage->vfvol.begin(); it != current_stage->vfvol.end(); it++)
        {
          if(!form_to_be_assembled(*it))
              continue;
        
          int order = calc_order_vector_form(*it);
        
          assemble_vector_form(*it, order);
        }
      }

      // Assemble surface integrals now: loop through surfaces of the element.
      for (this->current_isurf = 0; this->current_isurf < rep_element->get_num_surf(); this->current_isurf++)
      {
        // \todo DG.
        if(!this->current_state->bnd[current_isurf])
          continue;
        
        init_surface_state();

        if (current_mat != NULL)
        {
          for(typename Hermes::vector<MatrixFormSurf<Scalar> *>::iterator it = current_stage->mfsurf.begin(); it != current_stage->mfsurf.end(); it++)
          {
            if(!form_to_be_assembled(*it))
              continue;
        
            int order = calc_order_matrix_form(*it);
        
            assemble_matrix_form(*it, order);
          }
        }
    
        if (current_rhs != NULL)
        {
          for(typename Hermes::vector<VectorFormSurf<Scalar> *>::iterator it = current_stage->vfsurf.begin(); it != current_stage->vfsurf.end(); it++)
          {
            if(!form_to_be_assembled(*it))
                continue;
        
            int order = calc_order_vector_form(*it);
        
            assemble_vector_form(*it, order);
          }
        }
      }

      // Delete assembly lists.
      for(unsigned int i = 0; i < wf->get_neq(); i++)
        delete current_al[i];

      delete_cache();
    }
    



    template<typename Scalar>
    int DiscreteProblem<Scalar>::calc_order_matrix_form(MatrixForm<Scalar> *form)
    {
      _F_;

      int order;

      if(is_fvm)
        order = current_refmap[form->i]->get_inv_ref_order();
      else
      {
        // order of solutions from the previous Newton iteration etc..
        Func<Hermes::Ord>** u_ext_ord = new Func<Hermes::Ord>*[RungeKutta ? RK_original_spaces_count : current_u_ext.size() - form->u_ext_offset];
        ExtData<Hermes::Ord> ext_ord;
        init_ext_orders(form, u_ext_ord, &ext_ord);

        // Order of shape functions.
        Func<Hermes::Ord>* ou = get_fn_ord(this->spaces[form->j]->get_element_order(current_state->e[form->j]->id));
        Func<Hermes::Ord>* ov = get_fn_ord(this->spaces[form->i]->get_element_order(current_state->e[form->i]->id));

        // Total order of the vector form.
        Hermes::Ord o = form->ord(1, &fake_wt, u_ext_ord, ou, ov, &geom_ord, &ext_ord);

        adjust_order_to_refmaps(form, order, &o);

        // Cleanup.
        deinit_ext_orders(form, u_ext_ord, &ext_ord);
      }
      return order;
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::assemble_matrix_form(MatrixForm<Scalar>* form, int order)
    {
      _F_;
      bool surface_form = (dynamic_cast<MatrixFormVol<Scalar>*>(form) == NULL);

      double block_scaling_coef = this->block_scaling_coeff(form);

      bool tra = (form->i != form->j) && (form->sym != 0);
      bool sym = (form->i == form->j) && (form->sym == 1);

      // Assemble the local stiffness matrix for the form form.
      Scalar **local_stiffness_matrix = get_matrix_buffer(std::max(current_al[form->i]->cnt, current_al[form->j]->cnt));

      // Init external functions.
      Func<Scalar>** u_ext = new Func<Scalar>*[RungeKutta ? RK_original_spaces_count : current_u_ext.size() - form->u_ext_offset];
      ExtData<Scalar> ext;
      init_ext(form, u_ext, &ext, order);

      // Init geometry.
      int n_quadrature_points;
      if(surface_form)
        n_quadrature_points = init_surface_geometry_points(current_refmap[form->i], order);
      else
        n_quadrature_points = init_geometry_points(current_refmap[form->i], order);

      // Actual form-specific calculation.
      for (unsigned int i = 0; i < current_al[form->i]->cnt; i++)
      {
        if ((!tra || surface_form) && current_al[form->i]->dof[i] < 0) 
          continue;
        current_spss[form->i]->set_active_shape(current_al[form->i]->idx[i]);
        if (!sym)
        {
          for (unsigned int j = 0; j < current_al[form->j]->cnt; j++)
          {
            this->pss[form->j]->set_active_shape(current_al[form->j]->idx[j]);
            if (current_al[form->j]->dof[j] >= 0)
            {
              // Is this necessary, i.e. is there a coefficient smaller than 1e-12?
              if (std::abs(current_al[form->i]->coef[i]) < 1e-12 || std::abs(current_al[form->j]->coef[j]) < 1e-12)
                continue;
              Func<double>* u = get_fn(this->pss[form->j], this->current_refmap[form->j], order);
              Func<double>* v = get_fn(this->current_spss[form->i], this->current_refmap[form->i], order);

              // Add the previous time level solution previously inserted at the back of ext.
              if(RungeKutta)
                for(int ext_i = 0; ext_i < this->RK_original_spaces_count; ext_i++)
                  u_ext[ext_i]->add(*ext.fn[form->ext.size() - this->RK_original_spaces_count + ext_i]);

              if(surface_form)
                local_stiffness_matrix[i][j] = 0.5 * block_scaling_coeff(form) * form->value(n_quadrature_points, jacobian_x_weights_cache[order], u_ext, u, v, geometry_cache[order], &ext) * form->scaling_factor * current_al[form->j]->coef[j] * current_al[form->i]->coef[i];
              else
                local_stiffness_matrix[i][j] = block_scaling_coeff(form) * form->value(n_quadrature_points, jacobian_x_weights_cache[order], u_ext, u, v, geometry_cache[order], &ext) * form->scaling_factor * current_al[form->j]->coef[j] * current_al[form->i]->coef[i];
            }
          }
        }
        // Symmetric block.
        else
        {
          for (unsigned int j = 0; j < current_al[form->j]->cnt; j++)
          {
            if (j < i && current_al[form->j]->dof[j] >= 0)
              continue;
            this->pss[form->j]->set_active_shape(current_al[form->j]->idx[j]);
            if (current_al[form->j]->dof[j] >= 0)
            {
              Scalar val = 0;
              // Is this necessary, i.e. is there a coefficient smaller than 1e-12?
              if (std::abs(current_al[form->i]->coef[i]) < 1e-12 || std::abs(current_al[form->j]->coef[j]) < 1e-12)
                continue;

              Func<double>* u = get_fn(this->pss[form->j], this->current_refmap[form->j], order);
              Func<double>* v = get_fn(this->current_spss[form->i], this->current_refmap[form->i], order);

              // Add the previous time level solution previously inserted at the back of ext.
              if(RungeKutta)
                for(int ext_i = 0; ext_i < this->RK_original_spaces_count; ext_i++)
                  u_ext[ext_i]->add(*ext.fn[form->ext.size() - this->RK_original_spaces_count + ext_i]);

              val = block_scaling_coeff(form) * form->value(n_quadrature_points, jacobian_x_weights_cache[order], u_ext, u, v, geometry_cache[order], &ext) * form->scaling_factor * current_al[form->j]->coef[j] * current_al[form->i]->coef[i];

              local_stiffness_matrix[i][j] = local_stiffness_matrix[j][i] = val;
            }
          }
        }
      }

      // Insert the local stiffness matrix into the global one.
      current_mat->add(current_al[form->i]->cnt, current_al[form->j]->cnt, local_stiffness_matrix, current_al[form->i]->dof, current_al[form->j]->dof);

      // Insert also the off-diagonal (anti-)symmetric block, if required.
      if (tra)
      {
        if (form->sym < 0)
          chsgn(local_stiffness_matrix, current_al[form->i]->cnt, current_al[form->j]->cnt);
        transpose(local_stiffness_matrix, current_al[form->i]->cnt, current_al[form->j]->cnt);
        current_mat->add(current_al[form->j]->cnt, current_al[form->i]->cnt, local_stiffness_matrix, current_al[form->j]->dof, current_al[form->i]->dof);
      }

      // Cleanup.
      deinit_ext(form, u_ext, &ext);
    }
    
    
    template<typename Scalar>
    int DiscreteProblem<Scalar>::calc_order_vector_form(VectorForm<Scalar> *form)
    {
      _F_;

      int order;

      if(is_fvm)
        order = current_refmap[form->i]->get_inv_ref_order();
      else
      {
        // order of solutions from the previous Newton iteration etc..
        Func<Hermes::Ord>** u_ext_ord = new Func<Hermes::Ord>*[RungeKutta ? RK_original_spaces_count : current_u_ext.size() - form->u_ext_offset];
        ExtData<Hermes::Ord> ext_ord;
        init_ext_orders(form, u_ext_ord, &ext_ord);

        // Order of shape functions.
        Func<Hermes::Ord>* ov = get_fn_ord(this->spaces[form->i]->get_element_order(current_state->e[form->i]->id));

        // Total order of the vector form.
        Hermes::Ord o = form->ord(1, &fake_wt, u_ext_ord, ov, &geom_ord, &ext_ord);

        adjust_order_to_refmaps(form, order, &o);

        // Cleanup.
        deinit_ext_orders(form, u_ext_ord, &ext_ord);
      }
      return order;
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::assemble_vector_form(VectorForm<Scalar>* form, int order)
    {
      _F_;
      bool surface_form = (dynamic_cast<VectorFormVol<Scalar>*>(form) == NULL);

      // Init geometry.
      int n_quadrature_points;
      if(surface_form)
        n_quadrature_points = init_surface_geometry_points(current_refmap[form->i], order);
      else
        n_quadrature_points = init_geometry_points(current_refmap[form->i], order);

      // Init external functions.
      Func<Scalar>** u_ext = new Func<Scalar>*[RungeKutta ? RK_original_spaces_count : current_u_ext.size() - form->u_ext_offset];
      ExtData<Scalar> ext;
      init_ext(form, u_ext, &ext, order);

      // Actual form-specific calculation.
      for (unsigned int i = 0; i < current_al[form->i]->cnt; i++)
      {
        if (current_al[form->i]->dof[i] < 0)
          continue;
        current_spss[form->i]->set_active_shape(current_al[form->i]->idx[i]);
        
        // Is this necessary, i.e. is there a coefficient smaller than 1e-12?
        if (std::abs(current_al[form->i]->coef[i]) < 1e-12)
          continue;

        Func<double>* v = get_fn(this->current_spss[form->i], this->current_refmap[form->i], order);

        // Add the previous time level solution previously inserted at the back of ext.
        if(RungeKutta)
          for(int ext_i = 0; ext_i < this->RK_original_spaces_count; ext_i++)
            u_ext[ext_i]->add(*ext.fn[form->ext.size() - this->RK_original_spaces_count + ext_i]);
        
        if(surface_form)
          current_rhs->add(current_al[form->i]->dof[i], 0.5 * form->value(n_quadrature_points, jacobian_x_weights_cache[order], u_ext, v, geometry_cache[order], &ext) * form->scaling_factor * current_al[form->i]->coef[i]);
        else
          current_rhs->add(current_al[form->i]->dof[i], form->value(n_quadrature_points, jacobian_x_weights_cache[order], u_ext, v, geometry_cache[order], &ext) * form->scaling_factor * current_al[form->i]->coef[i]);
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
    int DiscreteProblem<Scalar>::init_surface_geometry_points(RefMap* reference_mapping, int& order)
    {
      _F_;
      int eo = quad->get_edge_points(current_isurf, order);
      double3* pt = quad->get_points(eo);
      int np = quad->get_num_points(eo);

      // Init geometry and jacobian*weights.
      if (geometry_cache[eo] == NULL)
      {
        geometry_cache[eo] = init_geom_surf(reference_mapping, current_isurf, current_state->rep->marker, eo);
        double3* tan = reference_mapping->get_tangent(current_isurf, eo);
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
    void DiscreteProblem<Scalar>::init_ext_orders(Form<Scalar> *form, Func<Hermes::Ord>** oi, ExtData<Hermes::Ord>* oext)
    {
      _F_;
      unsigned int prev_size = RungeKutta ? RK_original_spaces_count : current_u_ext.size() - form->u_ext_offset;

      if (current_u_ext != Hermes::vector<Solution<Scalar>*>())
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
    void DiscreteProblem<Scalar>::init_ext(Form<Scalar> *form, Func<Scalar>** u_ext, ExtData<Scalar>* ext, int order)
    {
      _F_;
      unsigned int prev_size = RungeKutta ? RK_original_spaces_count : current_u_ext.size() - form->u_ext_offset;
      
      if (current_u_ext != Hermes::vector<Solution<Scalar>*>())
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
      int prev_size = current_u_ext.size() - form->u_ext_offset;
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
    void DiscreteProblem<Scalar>::adjust_order_to_refmaps(Form<Scalar> *form, int& order, Hermes::Ord* o)
    {
      _F_;
      // Increase due to reference map.
      order = current_refmap[form->i]->get_inv_ref_order();
      order += o->get_order();
      limit_order(order);
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::assemble_DG_forms(Stage<Scalar>& stage,
      SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs, bool force_diagonal_blocks,
      Table* block_weights, Hermes::vector<PrecalcShapeset *>& spss,
      Hermes::vector<RefMap *>& refmap, Hermes::vector<Solution<Scalar>*>& u_ext,
      int marker, Hermes::vector<AsmList<Scalar>*>& al,
      bool bnd, SurfPos& surf_pos, Hermes::vector<bool>& nat,
      int isurf, Element** e, Element* trav_base, Element* rep_element)
    {
      _F_;
      // Determine the minimum mesh seq in this current_stage->
      min_dg_mesh_seq = 0;
      for(unsigned int i = 0; i < current_stage->meshes.size(); i++)
        if(current_stage->meshes[i]->get_seq() < min_dg_mesh_seq || i == 0)
          min_dg_mesh_seq = current_stage->meshes[i]->get_seq();

      // Initialize the NeighborSearches.
      // 5 is for bits per page in the array.
      LightArray<NeighborSearch<Scalar>*> neighbor_searches(5);
      init_neighbors(neighbor_searches, stage, isurf);

      // Create a multimesh tree;
      NeighborNode* root = new NeighborNode(NULL, 0);
      build_multimesh_tree(root, neighbor_searches);

      // Update all NeighborSearches according to the multimesh tree.
      // After this, all NeighborSearches in neighbor_searches should have the same count
      // of neighbors and proper set of transformations
      // for the central and the neighbor element(s) alike.
      // Also check that every NeighborSearch has the same number of neighbor elements.
      unsigned int num_neighbors = 0;
      for(unsigned int i = 0; i < neighbor_searches.get_size(); i++)
        if(neighbor_searches.present(i))
        {
          NeighborSearch<Scalar>* ns = neighbor_searches.get(i);
          update_neighbor_search(ns, root);
          if(num_neighbors == 0)
            num_neighbors = ns->n_neighbors;
          if(ns->n_neighbors != num_neighbors)
            error("Num_neighbors of different NeighborSearches not matching in DiscreteProblem<Scalar>::assemble_surface_integrals().");
        }

        // Create neighbor psss, refmaps.
        std::map<unsigned int, PrecalcShapeset *> npss;
        std::map<unsigned int, PrecalcShapeset *> nspss;
        std::map<unsigned int, RefMap *> nrefmap;

        // Initialize neighbor precalc shapesets and refmaps.
        // This is only needed when there are matrix DG forms present.
        if(DG_matrix_forms_present)
          for (unsigned int i = 0; i < current_stage->idx.size(); i++)
          {
            PrecalcShapeset* new_ps = new PrecalcShapeset(pss[i]->get_shapeset());
            new_ps->set_quad_2d(&g_quad_2d_std);
            npss.insert(std::pair<unsigned int, PrecalcShapeset*>(current_stage->idx[i], new_ps));

            PrecalcShapeset* new_pss = new PrecalcShapeset(new_ps);
            new_pss->set_quad_2d(&g_quad_2d_std);
            nspss.insert(std::pair<unsigned int, PrecalcShapeset*>(current_stage->idx[i], new_pss));

            RefMap* new_rm = new RefMap();
            new_rm->set_quad_2d(&g_quad_2d_std);
            nrefmap.insert(std::pair<unsigned int, RefMap*>(current_stage->idx[i], new_rm));
          }

          for(unsigned int neighbor_i = 0; neighbor_i < num_neighbors; neighbor_i++)
          {
            // If the active segment has already been processed (when the neighbor element was assembled), it is skipped.
            // We test all neighbor searches, because in the case of intra-element edge, the neighboring (the same as central) element
            // will be marked as visited, even though the edge was not calculated.
            bool processed = true;
            for(unsigned int i = 0; i < neighbor_searches.get_size(); i++)
              if(neighbor_searches.present(i))
                if(!neighbor_searches.get(i)->neighbors.at(neighbor_i)->visited)
                {
                  processed = false;
                  break;
                }

                if(!DG_vector_forms_present && processed)
                  continue;

                // For every neighbor we want to delete the geometry caches and create new ones.
                for (int i = 0; i < g_max_quad + 1 + 4 * g_max_quad + 4; i++)
                  if (geometry_cache[i] != NULL)
                  {
                    geometry_cache[i]->free();
                    delete geometry_cache[i];
                    geometry_cache[i] = NULL;
                    delete [] jacobian_x_weights_cache[i];
                  }

                  assemble_DG_one_neighbor(processed, neighbor_i, stage, mat, rhs,
                    force_diagonal_blocks, block_weights, spss, refmap,
                    npss, nspss, nrefmap, neighbor_searches, u_ext,
                    marker, al, bnd, surf_pos, nat, isurf, e, trav_base, rep_element);
          }

          // Delete the multimesh tree;
          delete root;

          // Deinitialize neighbor pss's, refmaps.
          if(DG_matrix_forms_present)
          {
            for(std::map<unsigned int, PrecalcShapeset *>::iterator it = nspss.begin(); it != nspss.end(); it++)
              delete it->second;
            for(std::map<unsigned int, PrecalcShapeset *>::iterator it = npss.begin(); it != npss.end(); it++)
              delete it->second;
            for(std::map<unsigned int, RefMap *>::iterator it = nrefmap.begin(); it != nrefmap.end(); it++)
              delete it->second;
          }

          // Delete the neighbor_searches array.
          for(unsigned int i = 0; i < neighbor_searches.get_size(); i++)
            if(neighbor_searches.present(i))
              delete neighbor_searches.get(i);
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::assemble_DG_one_neighbor(bool edge_processed, unsigned int neighbor_i, Stage<Scalar>& stage,
      SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs, bool force_diagonal_blocks, Table* block_weights,
      Hermes::vector<PrecalcShapeset *>& spss, Hermes::vector<RefMap *>& refmap, std::map<unsigned int, PrecalcShapeset *> npss,
      std::map<unsigned int, PrecalcShapeset *> nspss, std::map<unsigned int, RefMap *> nrefmap, LightArray<NeighborSearch<Scalar>*>& neighbor_searches, Hermes::vector<Solution<Scalar>*>& u_ext,
      int marker, Hermes::vector<AsmList<Scalar>*>& al, bool bnd, SurfPos& surf_pos, Hermes::vector<bool>& nat,
      int isurf, Element** e, Element* trav_base, Element* rep_element)
    {
      _F_;
      // Set the active segment in all NeighborSearches
      for(unsigned int i = 0; i < neighbor_searches.get_size(); i++)
        if(neighbor_searches.present(i))
        {
          neighbor_searches.get(i)->active_segment = neighbor_i;
          neighbor_searches.get(i)->neighb_el = neighbor_searches.get(i)->neighbors[neighbor_i];
          neighbor_searches.get(i)->neighbor_edge = neighbor_searches.get(i)->neighbor_edges[neighbor_i];
        }

        // Push all the necessary transformations to all functions of this current_stage->
        // The important thing is that the transformations to the current subelement are already there.
        for(unsigned int fns_i = 0; fns_i < current_stage->fns.size(); fns_i++)
        {
          NeighborSearch<Scalar>* ns = neighbor_searches.get(current_stage->meshes[fns_i]->get_seq() - min_dg_mesh_seq);
          if (ns->central_transformations.present(neighbor_i))
            ns->central_transformations.get(neighbor_i)->apply_on(current_stage->fns[fns_i]);
        }

        // For neighbor psss.
        if(DG_matrix_forms_present && !edge_processed)
          for(unsigned int idx_i = 0; idx_i < current_stage->idx.size(); idx_i++)
          {
            NeighborSearch<Scalar>* ns = neighbor_searches.get(current_stage->meshes[idx_i]->get_seq() - min_dg_mesh_seq);
            npss[current_stage->idx[idx_i]]->set_active_element((*ns->get_neighbors())[neighbor_i]);
            if (ns->neighbor_transformations.present(neighbor_i))
              ns->neighbor_transformations.get(neighbor_i)->apply_on(npss[current_stage->idx[idx_i]]);
          }

          // Also push the transformations to the slave psss and refmaps.
          for (unsigned int i = 0; i < current_stage->idx.size(); i++)
          {
            if(current_state->e[current_stage->idx[i]] == NULL)
              continue;
            spss[current_stage->idx[i]]->set_master_transform();
            refmap[current_stage->idx[i]]->force_transform(pss[current_stage->idx[i]]->get_transform(), pss[current_stage->idx[i]]->get_ctm());

            // Neighbor.
            if(DG_matrix_forms_present && !edge_processed)
            {
              nspss[current_stage->idx[i]]->set_active_element(npss[current_stage->idx[i]]->get_active_element());
              nspss[current_stage->idx[i]]->set_master_transform();
              nrefmap[current_stage->idx[i]]->set_active_element(npss[current_stage->idx[i]]->get_active_element());
              nrefmap[current_stage->idx[i]]->force_transform(npss[current_stage->idx[i]]->get_transform(), npss[current_stage->idx[i]]->get_ctm());
            }
          }

          /***/

          // The computation takes place here.
          if(DG_matrix_forms_present && !edge_processed)
          {
            assemble_DG_matrix_forms(stage, mat, rhs, force_diagonal_blocks, block_weights, spss, refmap, npss, nspss, nrefmap, neighbor_searches, u_ext,
              marker, al, bnd, surf_pos, nat, isurf, e, trav_base, rep_element);
          }
          if (DG_vector_forms_present && rhs != NULL)
          {
            assemble_DG_vector_forms(stage, mat, rhs, force_diagonal_blocks, block_weights, spss, refmap, neighbor_searches, u_ext,
              marker, al, bnd, surf_pos, nat, isurf, e, trav_base, rep_element);
          }
          /***/

          // This is just cleaning after ourselves.
          // Clear the transformations from the RefMaps and all functions.
          for(unsigned int fns_i = 0; fns_i < current_stage->fns.size(); fns_i++)
            current_stage->fns[fns_i]->set_transform(neighbor_searches.get(current_stage->meshes[fns_i]->get_seq() - min_dg_mesh_seq)->original_central_el_transform);

          // Also clear the transformations from the slave psss and refmaps.
          for (unsigned int i = 0; i < current_stage->idx.size(); i++)
          {
            if(current_state->e[current_stage->idx[i]] == NULL)
              continue;
            spss[current_stage->idx[i]]->set_master_transform();
            refmap[current_stage->idx[i]]->force_transform(pss[current_stage->idx[i]]->get_transform(), pss[current_stage->idx[i]]->get_ctm());
          }
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::init_neighbors(LightArray<NeighborSearch<Scalar>*>& neighbor_searches,
      const Stage<Scalar>& stage, const int& isurf)
    {
      _F_;
      // Initialize the NeighborSearches.
      for(unsigned int i = 0; i < current_stage->meshes.size(); i++)
      {
        if(i > 0 && current_stage->meshes.at(i - 1)->get_seq() == current_stage->meshes.at(i)->get_seq())
          continue;
        else
          if(!neighbor_searches.present(current_stage->meshes[i]->get_seq() - min_dg_mesh_seq))
          {
            NeighborSearch<Scalar>* ns = new NeighborSearch<Scalar>(current_stage->fns[i]->get_active_element(), current_stage->meshes[i]);
            ns->original_central_el_transform = current_stage->fns[i]->get_transform();
            neighbor_searches.add(ns, current_stage->meshes[i]->get_seq() - min_dg_mesh_seq);
          }
      }

      // Calculate respective neighbors.
      // Also clear the initial_sub_idxs from the central element transformations
      // of NeighborSearches with multiple neighbors.
      for(unsigned int i = 0; i < neighbor_searches.get_size(); i++)
      {
        if(i > 0 && current_stage->meshes.at(i - 1)->get_seq() == current_stage->meshes.at(i)->get_seq())
          continue;
        if(neighbor_searches.present(i))
        {
          neighbor_searches.get(i)->set_active_edge_multimesh(isurf);
          neighbor_searches.get(i)->clear_initial_sub_idx();
        }
      }
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::build_multimesh_tree(NeighborNode* root,
      LightArray<NeighborSearch<Scalar>*>& neighbor_searches)
    {
      _F_;
        for(unsigned int i = 0; i < neighbor_searches.get_size(); i++)
        if(neighbor_searches.present(i))
        {
          NeighborSearch<Scalar>* ns = neighbor_searches.get(i);
          if (ns->n_neighbors == 1 &&
              (ns->central_transformations.get_size() == 0 || ns->central_transformations.get(0)->num_levels == 0))
            continue;
          for(unsigned int j = 0; j < ns->n_neighbors; j++)
            if (ns->central_transformations.present(j))
              insert_into_multimesh_tree(root, ns->central_transformations.get(j)->transf, ns->central_transformations.get(j)->num_levels);
        }
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::insert_into_multimesh_tree(NeighborNode* node,
      unsigned int* transformations,
      unsigned int transformation_count)
    {
      _F_;
      // If we are already in the leaf.
      if(transformation_count == 0)
        return;
      // Both sons are null. We have to add a new Node. Let us do it for the left sone of node.
      if(node->get_left_son() == NULL && node->get_right_son() == NULL)
      {
        node->set_left_son(new NeighborNode(node, transformations[0]));
        insert_into_multimesh_tree(node->get_left_son(), transformations + 1, transformation_count - 1);
      }
      // At least the left son is not null (it is impossible only for the right one to be not null, because
      // the left one always gets into the tree first, as seen above).
      else
      {
        // The existing left son is the right one to continue through.
        if(node->get_left_son()->get_transformation() == transformations[0])
          insert_into_multimesh_tree(node->get_left_son(), transformations + 1, transformation_count - 1);
        // The right one also exists, check that it is the right one, or return an error.
        else if(node->get_right_son() != NULL)
        {
          if(node->get_right_son()->get_transformation() == transformations[0])
            insert_into_multimesh_tree(node->get_right_son(), transformations + 1, transformation_count - 1);
          else error("More than two possible sons in insert_into_multimesh_tree().");
        }
        // If the right one does not exist and the left one was not correct, create a right son and continue this way.
        else
        {
          node->set_right_son(new NeighborNode(node, transformations[0]));
          insert_into_multimesh_tree(node->get_right_son(), transformations + 1, transformation_count - 1);
        }
      }
    }

    template<typename Scalar>
    Hermes::vector<Hermes::vector<unsigned int>*> DiscreteProblem<Scalar>::get_multimesh_neighbors_transformations(NeighborNode* multimesh_tree)
    {
      _F_;
      // Initialize the vector.
      Hermes::vector<Hermes::vector<unsigned int>*> running_transformations;
      // Prepare the first neighbor's vector.
      running_transformations.push_back(new Hermes::vector<unsigned int>);
      // Fill the vector.
      traverse_multimesh_tree(multimesh_tree, running_transformations);
      return running_transformations;
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::traverse_multimesh_tree(NeighborNode* node,
      Hermes::vector<Hermes::vector<unsigned int>*>& running_transformations)
    {
      _F_;
      // If we are in the root.
      if(node->get_transformation() == 0)
      {
        if(node->get_left_son() != NULL)
          traverse_multimesh_tree(node->get_left_son(), running_transformations);
        if(node->get_right_son() != NULL)
          traverse_multimesh_tree(node->get_right_son(), running_transformations);
        // Delete the vector prepared by the last accessed leaf.
        delete running_transformations.back();
        running_transformations.pop_back();
        return;
      }
      // If we are in a leaf.
      if(node->get_left_son() == NULL && node->get_right_son() == NULL)
      {
        // Create a vector for the new neighbor.
        Hermes::vector<unsigned int>* new_neighbor_transformations = new Hermes::vector<unsigned int>;
        // Copy there the whole path except for this leaf.
        for(unsigned int i = 0; i < running_transformations.back()->size(); i++)
          new_neighbor_transformations->push_back((*running_transformations.back())[i]);
        // Insert this leaf into the current running transformation, thus complete it.
        running_transformations.back()->push_back(node->get_transformation());
        // Make the new_neighbor_transformations the current running transformation.
        running_transformations.push_back(new_neighbor_transformations);
        return;
      }
      else
      {
        running_transformations.back()->push_back(node->get_transformation());
        if(node->get_left_son() != NULL)
          traverse_multimesh_tree(node->get_left_son(), running_transformations);
        if(node->get_right_son() != NULL)
          traverse_multimesh_tree(node->get_right_son(), running_transformations);
        running_transformations.back()->pop_back();
        return;
      }
      return;
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::update_neighbor_search(NeighborSearch<Scalar>* ns, NeighborNode* multimesh_tree)
    {
      _F_;
      // This has to be done, because we pass ns by reference and the number of neighbors is changing.
      unsigned int num_neighbors = ns->get_num_neighbors();

      for(unsigned int i = 0; i < num_neighbors; i++)
      {
        // Find the node corresponding to this neighbor in the tree.
        NeighborNode* node;
        if (ns->central_transformations.present(i))
          node = find_node(ns->central_transformations.get(i)->transf, ns->central_transformations.get(i)->num_levels, multimesh_tree);
        else
          node = multimesh_tree;

        // Update the NeighborSearch.
        unsigned int added = update_ns_subtree(ns, node, i);
        i += added;
        num_neighbors += added;
      }
    }

    template<typename Scalar>
    NeighborNode* DiscreteProblem<Scalar>::find_node(unsigned int* transformations,
      unsigned int transformation_count,
      NeighborNode* node)
    {
      _F_;
      // If there are no transformations left.
      if(transformation_count == 0)
        return node;
      else
      {
        if(node->get_left_son() != NULL)
        {
          if(node->get_left_son()->get_transformation() == transformations[0])
            return find_node(transformations + 1, transformation_count - 1, node->get_left_son());
        }
        if(node->get_right_son() != NULL)
        {
          if(node->get_right_son()->get_transformation() == transformations[0])
            return find_node(transformations + 1, transformation_count - 1, node->get_right_son());
        }
      }
      // We always should be able to empty the transformations array.
      error("Transformation of a central element not found in the multimesh tree.");
      return NULL;
    }

    template<typename Scalar>
    unsigned int DiscreteProblem<Scalar>::update_ns_subtree(NeighborSearch<Scalar>* ns,
      NeighborNode* node, unsigned int ith_neighbor)
    {
      _F_;
      // No subtree => no work.
      // Also check the assertion that if one son is null, then the other too.
      if(node->get_left_son() == NULL)
      {
        if(node->get_right_son() != NULL)
          error("Only one son (right) not null in DiscreteProblem<Scalar>::update_ns_subtree.");
        return 0;
      }

      // Key part.
      // Begin with storing the info about the current neighbor.
      Element* neighbor = ns->neighbors[ith_neighbor];
      typename NeighborSearch<Scalar>::NeighborEdgeInfo edge_info = ns->neighbor_edges[ith_neighbor];

      // Initialize the vector for central transformations->
      Hermes::vector<Hermes::vector<unsigned int>*> running_central_transformations;
      // Prepare the first new neighbor's vector. Push back the current transformations (in case of GO_DOWN neighborhood).
      running_central_transformations.push_back(new Hermes::vector<unsigned int>);
      if (ns->central_transformations.present(ith_neighbor))
        ns->central_transformations.get(ith_neighbor)->copy_to(running_central_transformations.back());

      // Initialize the vector for neighbor transformations->
      Hermes::vector<Hermes::vector<unsigned int>*> running_neighbor_transformations;
      // Prepare the first new neighbor's vector. Push back the current transformations (in case of GO_UP/NO_TRF neighborhood).
      running_neighbor_transformations.push_back(new Hermes::vector<unsigned int>);
      if (ns->neighbor_transformations.present(ith_neighbor))
        ns->neighbor_transformations.get(ith_neighbor)->copy_to(running_neighbor_transformations.back());

      // Delete the current neighbor.
      ns->delete_neighbor(ith_neighbor);

      // Move down the subtree.
      if(node->get_left_son() != NULL)
        traverse_multimesh_subtree(node->get_left_son(), running_central_transformations,
        running_neighbor_transformations, edge_info, ns->active_edge,
        ns->central_el->get_mode());
      if(node->get_right_son() != NULL)
        traverse_multimesh_subtree(node->get_right_son(), running_central_transformations,
        running_neighbor_transformations, edge_info, ns->active_edge,
        ns->central_el->get_mode());

      // Delete the last neighbors' info (this is a dead end, caused by the function traverse_multimesh_subtree.
      delete running_central_transformations.back();
      running_central_transformations.pop_back();
      delete running_neighbor_transformations.back();
      running_neighbor_transformations.pop_back();

      // Insert new neighbors.
      for(unsigned int i = 0; i < running_central_transformations.size(); i++)
      {
        ns->neighbors.push_back(neighbor);
        ns->neighbor_edges.push_back(edge_info);

        if (!ns->central_transformations.present(ns->n_neighbors))
          ns->central_transformations.add(new typename NeighborSearch<Scalar>::Transformations, ns->n_neighbors);
        if (!ns->neighbor_transformations.present(ns->n_neighbors))
          ns->neighbor_transformations.add(new typename NeighborSearch<Scalar>::Transformations, ns->n_neighbors);
        ns->central_transformations.get(ns->n_neighbors)->copy_from(*running_central_transformations[i]);
        ns->neighbor_transformations.get(ns->n_neighbors)->copy_from(*running_neighbor_transformations[i]);

        ns->n_neighbors++;
      }

      for(unsigned int i = 0; i < running_central_transformations.size(); i++)
        delete running_central_transformations[i];
      for(unsigned int i = 0; i < running_neighbor_transformations.size(); i++)
        delete running_neighbor_transformations[i];

      // Return the number of neighbors deleted.
      return -1;
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::traverse_multimesh_subtree(NeighborNode* node,
      Hermes::vector<Hermes::vector<unsigned int>*>& running_central_transformations,
      Hermes::vector<Hermes::vector<unsigned int>*>& running_neighbor_transformations,
      const typename NeighborSearch<Scalar>::NeighborEdgeInfo& edge_info, const int& active_edge, const int& mode)
    {
      _F_;
      // If we are in a leaf.
      if(node->get_left_son() == NULL && node->get_right_son() == NULL)
      {
        // Create vectors for the new neighbor.
        Hermes::vector<unsigned int>* new_neighbor_central_transformations = new Hermes::vector<unsigned int>;
        Hermes::vector<unsigned int>* new_neighbor_neighbor_transformations = new Hermes::vector<unsigned int>;

        // Copy there the whole path except for this leaf.
        for(unsigned int i = 0; i < running_central_transformations.back()->size(); i++)
          new_neighbor_central_transformations->push_back((*running_central_transformations.back())[i]);
        for(unsigned int i = 0; i < running_neighbor_transformations.back()->size(); i++)
          new_neighbor_neighbor_transformations->push_back((*running_neighbor_transformations.back())[i]);

        // Insert this leaf into the current running central transformation, thus complete it.
        running_central_transformations.back()->push_back(node->get_transformation());

        // Make the new_neighbor_central_transformations the current running central transformation.
        running_central_transformations.push_back(new_neighbor_central_transformations);

        // Take care of the neighbor transformation.
        // Insert appropriate info from this leaf into the current running neighbor transformation, thus complete it.
        if(mode == HERMES_MODE_TRIANGLE)
          if ((active_edge == 0 && node->get_transformation() == 0) ||
            (active_edge == 1 && node->get_transformation() == 1) ||
            (active_edge == 2 && node->get_transformation() == 2))
            running_neighbor_transformations.back()->push_back((!edge_info.orientation ? edge_info.local_num_of_edge : (edge_info.local_num_of_edge + 1) % 3));
          else
            running_neighbor_transformations.back()->push_back((edge_info.orientation ? edge_info.local_num_of_edge : (edge_info.local_num_of_edge + 1) % 3));
        // Quads.
        else
          if ((active_edge == 0 && (node->get_transformation() == 0 || node->get_transformation() == 6)) ||
            (active_edge == 1 && (node->get_transformation() == 1 || node->get_transformation() == 4)) ||
            (active_edge == 2 && (node->get_transformation() == 2 || node->get_transformation() == 7)) ||
            (active_edge == 3 && (node->get_transformation() == 3 || node->get_transformation() == 5)))
            running_neighbor_transformations.back()->push_back((!edge_info.orientation ? edge_info.local_num_of_edge : (edge_info.local_num_of_edge + 1) % 4));
          else
            running_neighbor_transformations.back()->push_back((edge_info.orientation ? edge_info.local_num_of_edge : (edge_info.local_num_of_edge + 1) % 4));

        // Make the new_neighbor_neighbor_transformations the current running neighbor transformation.
        running_neighbor_transformations.push_back(new_neighbor_neighbor_transformations);

        return;
      }
      else
      {
        // Insert this leaf into the current running central transformation, thus complete it.
        running_central_transformations.back()->push_back(node->get_transformation());

        // Insert appropriate info from this leaf into the current running neighbor transformation, thus complete it.
        // Triangles.
        if(mode == HERMES_MODE_TRIANGLE)
          if ((active_edge == 0 && node->get_transformation() == 0) ||
            (active_edge == 1 && node->get_transformation() == 1) ||
            (active_edge == 2 && node->get_transformation() == 2))
            running_neighbor_transformations.back()->push_back((!edge_info.orientation ? edge_info.local_num_of_edge : (edge_info.local_num_of_edge + 1) % 3));
          else
            running_neighbor_transformations.back()->push_back((edge_info.orientation ? edge_info.local_num_of_edge : (edge_info.local_num_of_edge + 1) % 3));
        // Quads.
        else
          if ((active_edge == 0 && (node->get_transformation() == 0 || node->get_transformation() == 6)) ||
            (active_edge == 1 && (node->get_transformation() == 1 || node->get_transformation() == 4)) ||
            (active_edge == 2 && (node->get_transformation() == 2 || node->get_transformation() == 7)) ||
            (active_edge == 3 && (node->get_transformation() == 3 || node->get_transformation() == 5)))
            running_neighbor_transformations.back()->push_back((!edge_info.orientation ? edge_info.local_num_of_edge : (edge_info.local_num_of_edge + 1) % 4));
          else
            running_neighbor_transformations.back()->push_back((edge_info.orientation ? edge_info.local_num_of_edge : (edge_info.local_num_of_edge + 1) % 4));

        // Move down.
        if(node->get_left_son() != NULL)
          traverse_multimesh_subtree(node->get_left_son(), running_central_transformations, running_neighbor_transformations,
          edge_info, active_edge, mode);
        if(node->get_right_son() != NULL)
          traverse_multimesh_subtree(node->get_right_son(), running_central_transformations, running_neighbor_transformations,
          edge_info, active_edge, mode);

        // Take this transformation out.
        running_central_transformations.back()->pop_back();
        running_neighbor_transformations.back()->pop_back();
        return;
      }
      return;
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::assemble_DG_matrix_forms(Stage<Scalar>& stage,
      SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs, bool force_diagonal_blocks,
      Table* block_weights, Hermes::vector<PrecalcShapeset *>& spss,
      Hermes::vector<RefMap *>& refmap, std::map<unsigned int, PrecalcShapeset *> npss,
      std::map<unsigned int, PrecalcShapeset *> nspss, std::map<unsigned int, RefMap *> nrefmap,
      LightArray<NeighborSearch<Scalar>*>& neighbor_searches, Hermes::vector<Solution<Scalar>*>& u_ext,
      int marker, Hermes::vector<AsmList<Scalar>*>& al, bool bnd,
      SurfPos& surf_pos, Hermes::vector<bool>& nat, int isurf, Element** e,
      Element* trav_base, Element* rep_element)
    {
      _F_;
      for (unsigned int ww = 0; ww < current_stage->mfsurf.size(); ww++)
      {
        MatrixFormSurf<Scalar>* mfs = current_stage->mfsurf[ww];
        if (mfs->areas[0] != H2D_DG_INNER_EDGE)
          continue;
        int m = mfs->i;
        int n = mfs->j;

        if (current_state->e[m] == NULL || current_state->e[n] == NULL)
          continue;
        if (fabs(mfs->scaling_factor) < 1e-12)
          continue;

        surf_pos.base = trav_base;

        // Create the extended shapeset on the union of the central element and its current neighbor.
        typename NeighborSearch<Scalar>::ExtendedShapeset* ext_asmlist_u = neighbor_searches.get(spaces[n]->get_mesh()->get_seq() - min_dg_mesh_seq)->create_extended_asmlist(spaces[n], al[n]);
        typename NeighborSearch<Scalar>::ExtendedShapeset* ext_asmlist_v = neighbor_searches.get(spaces[m]->get_mesh()->get_seq() - min_dg_mesh_seq)->create_extended_asmlist(spaces[m], al[m]);

        // If a block scaling table is provided, and if the scaling coefficient
        // A_mn for this block is zero, then the form does not need to be assembled.
        double block_scaling_coeff = 1.;
        if (block_weights != NULL)
        {
          block_scaling_coeff = block_weights->get_A(m, n);
          if (fabs(block_scaling_coeff) < 1e-12)
            continue;
        }

        // Precalc shapeset and refmaps used for the evaluation.
        PrecalcShapeset* fu;
        PrecalcShapeset* fv;
        RefMap* ru;
        RefMap* rv;
        bool support_neigh_u, support_neigh_v;

        Scalar **local_stiffness_matrix = get_matrix_buffer(std::max(ext_asmlist_u->cnt, ext_asmlist_v->cnt));
        for (int i = 0; i < ext_asmlist_v->cnt; i++)
        {
          if (ext_asmlist_v->dof[i] < 0)
            continue;
          // Choose the correct shapeset for the test function.
          if (!ext_asmlist_v->has_support_on_neighbor(i))
          {
            spss[m]->set_active_shape(ext_asmlist_v->central_al->idx[i]);
            fv = spss[m];
            rv = refmap[m];
            support_neigh_v = false;
          }
          else
          {
            nspss[m]->set_active_shape(ext_asmlist_v->neighbor_al->idx[i - ext_asmlist_v->central_al->cnt]);
            fv = nspss[m];
            rv = nrefmap[m];
            support_neigh_v = true;
          }
          for (int j = 0; j < ext_asmlist_u->cnt; j++)
          {
            // Choose the correct shapeset for the solution function.
            if (!ext_asmlist_u->has_support_on_neighbor(j))
            {
              pss[n]->set_active_shape(ext_asmlist_u->central_al->idx[j]);
              fu = pss[n];
              ru = refmap[n];
              support_neigh_u = false;
            }
            else
            {
              npss[n]->set_active_shape(ext_asmlist_u->neighbor_al->idx[j - ext_asmlist_u->central_al->cnt]);
              fu = npss[n];
              ru = nrefmap[n];
              support_neigh_u = true;
            }

            if (ext_asmlist_u->dof[j] >= 0)
            {
              if (mat != NULL)
              {
                Scalar val = block_scaling_coeff * eval_dg_form(mfs, u_ext, fu, fv, refmap[n], ru, rv, support_neigh_u, support_neigh_v, &surf_pos, neighbor_searches, current_stage->meshes[n]->get_seq() - min_dg_mesh_seq, current_stage->meshes[m]->get_seq() - min_dg_mesh_seq)
                  * (support_neigh_u ? ext_asmlist_u->neighbor_al->coef[j - ext_asmlist_u->central_al->cnt]: ext_asmlist_u->central_al->coef[j])
                  * (support_neigh_v ? ext_asmlist_v->neighbor_al->coef[i - ext_asmlist_v->central_al->cnt]: ext_asmlist_v->central_al->coef[i]);
                local_stiffness_matrix[i][j] = val;
              }
            }
          }
        }
        if (mat != NULL)
          current_mat->add(ext_asmlist_v->cnt, ext_asmlist_u->cnt, local_stiffness_matrix, ext_asmlist_v->dof, ext_asmlist_u->dof);
      }
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::assemble_DG_vector_forms(Stage<Scalar>& stage,
      SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs, bool force_diagonal_blocks, Table* block_weights,
      Hermes::vector<PrecalcShapeset *>& spss, Hermes::vector<RefMap *>& refmap, LightArray<NeighborSearch<Scalar>*>& neighbor_searches, Hermes::vector<Solution<Scalar>*>& u_ext,
      int marker, Hermes::vector<AsmList<Scalar>*>& al, bool bnd, SurfPos& surf_pos, Hermes::vector<bool>& nat,
      int isurf, Element** e, Element* trav_base, Element* rep_element)
    {
      _F_;
      for (unsigned int ww = 0; ww < current_stage->vfsurf.size(); ww++)
      {
        VectorFormSurf<Scalar>* vfs = current_stage->vfsurf[ww];
        if (vfs->areas[0] != H2D_DG_INNER_EDGE)
          continue;
        int m = vfs->i;
        if (current_state->e[m] == NULL)
          continue;
        if (fabs(vfs->scaling_factor) < 1e-12)
          continue;

        // Here we use the standard pss, possibly just transformed by NeighborSearch.
        for (unsigned int i = 0; i < al[m]->cnt; i++)
        {
          if (al[m]->dof[i] < 0)
            continue;
          spss[m]->set_active_shape(al[m]->idx[i]);
          rhs->add(al[m]->dof[i], eval_dg_form(vfs, u_ext, spss[m], refmap[m], &surf_pos, neighbor_searches, current_stage->meshes[m]->get_seq() - min_dg_mesh_seq) * al[m]->coef[i]);
        }
      }
    }

    template<typename Scalar>
    ExtData<Scalar>* DiscreteProblem<Scalar>::init_ext_fns(Hermes::vector<MeshFunction<Scalar>*> &ext,
      LightArray<NeighborSearch<Scalar>*>& neighbor_searches, int order)
    {
      _F_;
      Func<Scalar>** ext_fns = new Func<Scalar>*[ext.size()];
      for(unsigned int j = 0; j < ext.size(); j++)
      {
        neighbor_searches.get(ext[j]->get_mesh()->get_seq() - min_dg_mesh_seq)->set_quad_order(order);
        ext_fns[j] = neighbor_searches.get(ext[j]->get_mesh()->get_seq() - min_dg_mesh_seq)->init_ext_fn(ext[j]);
      }

      ExtData<Scalar>* ext_data = new ExtData<Scalar>;
      ext_data->fn = ext_fns;
      ext_data->nf = ext.size();

      return ext_data;
    }

    template<typename Scalar>
    ExtData<Hermes::Ord>* DiscreteProblem<Scalar>::init_ext_fns_ord(Hermes::vector<MeshFunction<Scalar>*> &ext,
      LightArray<NeighborSearch<Scalar>*>& neighbor_searches)
    {
      _F_;
      Func<Hermes::Ord>** fake_ext_fns = new Func<Hermes::Ord>*[ext.size()];
      for (unsigned int j = 0; j < ext.size(); j++)
        fake_ext_fns[j] = init_ext_fn_ord(neighbor_searches.get(ext[j]->get_mesh()->get_seq() - min_dg_mesh_seq), ext[j]);

      ExtData<Hermes::Ord>* fake_ext = new ExtData<Hermes::Ord>;
      fake_ext->fn = fake_ext_fns;
      fake_ext->nf = ext.size();

      return fake_ext;
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
    DiscontinuousFunc<Hermes::Ord>* DiscreteProblem<Scalar>::init_ext_fn_ord(NeighborSearch<Scalar>* ns, MeshFunction<Scalar>* fu)
    {
      _F_;
      int inc = (fu->get_num_components() == 2) ? 1 : 0;
      int central_order = fu->get_edge_fn_order(ns->active_edge) + inc;
      int neighbor_order = fu->get_edge_fn_order(ns->neighbor_edge.local_num_of_edge) + inc;
      return new DiscontinuousFunc<Hermes::Ord>(get_fn_ord(central_order), get_fn_ord(neighbor_order));
    }
    
    template<typename Scalar>
    void DiscreteProblem<Scalar>::set_quad_2d(Quad2D* quad)
    {
      this->quad = quad;
    }

    template<typename Scalar>
    int DiscreteProblem<Scalar>::calc_order_dg_matrix_form(MatrixFormSurf<Scalar> *mfs, Hermes::vector<Solution<Scalar>*> u_ext,
      PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, SurfPos* surf_pos,
      bool neighbor_supp_u, bool neighbor_supp_v, LightArray<NeighborSearch<Scalar>*>& neighbor_searches, int neighbor_index_u)
    {
      NeighborSearch<Scalar>* nbs_u = neighbor_searches.get(neighbor_index_u);
      // Hermes::Order that will be returned.
      int order;

      if(this->is_fvm)
        order = ru->get_inv_ref_order();
      else
      {
        // Hermes::Order of solutions from the previous Newton iteration.
        int prev_size = u_ext.size() - mfs->u_ext_offset;
        Func<Hermes::Ord>** oi = new Func<Hermes::Ord>*[prev_size];
        if (u_ext != Hermes::vector<Solution<Scalar>*>())
          for (int i = 0; i < prev_size; i++)
            if (u_ext[i + mfs->u_ext_offset] != NULL)
              oi[i] = init_ext_fn_ord(neighbor_searches.get(u_ext[i]->get_mesh()->get_seq() - min_dg_mesh_seq), u_ext[i]);
            else
              oi[i] = get_fn_ord(0);
        else
          for (int i = 0; i < prev_size; i++)
            oi[i] = get_fn_ord(0);

        // Order of shape functions.
        int inc = (fv->get_num_components() == 2) ? 1 : 0;
        DiscontinuousFunc<Hermes::Ord>* ou = new DiscontinuousFunc<Hermes::Ord>(get_fn_ord(fu->get_edge_fn_order(surf_pos->surf_num) + inc), neighbor_supp_u);
        DiscontinuousFunc<Hermes::Ord>* ov = new DiscontinuousFunc<Hermes::Ord>(get_fn_ord(fv->get_edge_fn_order(surf_pos->surf_num) + inc), neighbor_supp_v);

        // Hermes::Order of additional external functions.
        ExtData<Hermes::Ord>* fake_ext = init_ext_fns_ord(mfs->ext, neighbor_searches);

        // Hermes::Order of geometric attributes (eg. for multiplication of a solution with coordinates, normals, etc.).
        Geom<Hermes::Ord>* fake_e = new InterfaceGeom<Hermes::Ord>(&geom_ord, nbs_u->neighb_el->marker,
          nbs_u->neighb_el->id, Hermes::Ord(nbs_u->neighb_el->get_diameter()));
        double fake_wt = 1.0;

        // Total order of the vector form.
        Hermes::Ord o = mfs->ord(1, &fake_wt, oi, ou, ov, fake_e, fake_ext);

        // Increase due to reference maps.
        order = ru->get_inv_ref_order();

        order += o.get_order();
        limit_order(order);

        // Clean up.
        if (u_ext != Hermes::vector<Solution<Scalar>*>())
          for (int i = 0; i < prev_size; i++)
            if (u_ext[i + mfs->u_ext_offset] != NULL)
              delete oi[i];
        delete [] oi;
        delete fake_e;
        delete ou;
        delete ov;
        if (fake_ext != NULL)
        {
          for (int i = 0; i < fake_ext->nf; i++)
          {
            delete fake_ext->fn[i];
          }
          fake_ext->free_ord();
          delete fake_ext;
        }
      }

      return order;
    }
    
    template<typename Scalar>
    Scalar DiscreteProblem<Scalar>::eval_dg_form(MatrixFormSurf<Scalar>* mfs, Hermes::vector<Solution<Scalar>*> u_ext,
      PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru_central, RefMap *ru_actual, RefMap *rv,
      bool neighbor_supp_u, bool neighbor_supp_v,
      SurfPos* surf_pos, LightArray<NeighborSearch<Scalar>*>& neighbor_searches, int neighbor_index_u, int neighbor_index_v)
    {
      _F_;

      NeighborSearch<Scalar>* nbs_u = neighbor_searches.get(neighbor_index_u);
      NeighborSearch<Scalar>* nbs_v = neighbor_searches.get(neighbor_index_v);

      // Determine the integration order.
      int order = calc_order_dg_matrix_form(mfs, u_ext, fu, fv, ru_actual, surf_pos, neighbor_supp_u, neighbor_supp_v, neighbor_searches, neighbor_index_u);

      // Evaluate the form using just calculated order.
      Quad2D* quad = fu->get_quad_2d();
      int eo = quad->get_edge_points(surf_pos->surf_num, order);
      int np = quad->get_num_points(eo);
      double3* pt = quad->get_points(eo);

      // A (debug) check.
      assert(surf_pos->surf_num == neighbor_searches.get(neighbor_index_u)->active_edge);

      // Init geometry and jacobian*weights.
      if (geometry_cache[eo] == NULL)
      {
        geometry_cache[eo] = init_geom_surf(ru_central, surf_pos->surf_num, surf_pos->marker, eo);
        double3* tan = ru_central->get_tangent(surf_pos->surf_num, eo);
        jacobian_x_weights_cache[eo] = new double[np];
        for(int i = 0; i < np; i++)
          jacobian_x_weights_cache[eo][i] = pt[i][2] * tan[i][2];
      }

      Geom<double>* e = new InterfaceGeom<double>(geometry_cache[eo], nbs_u->neighb_el->marker,
        nbs_u->neighb_el->id, nbs_u->neighb_el->get_diameter());
      double* jwt = jacobian_x_weights_cache[eo];

      // Values of the previous Newton iteration, shape functions and external functions in quadrature points.
      int prev_size = u_ext.size() - mfs->u_ext_offset;
      Func<Scalar>** prev = new Func<Scalar>*[prev_size];
      if (u_ext != Hermes::vector<Solution<Scalar>*>())
        for (int i = 0; i < prev_size; i++)
          if (u_ext[i + mfs->u_ext_offset] != NULL)
          {
            neighbor_searches.get(u_ext[i]->get_mesh()->get_seq() - min_dg_mesh_seq)->set_quad_order(order);
            prev[i]  = neighbor_searches.get(u_ext[i]->get_mesh()->get_seq() - min_dg_mesh_seq)->init_ext_fn(u_ext[i]);
          }
          else prev[i] = NULL;
      else
        for (int i = 0; i < prev_size; i++)
          prev[i] = NULL;

      // Values of the previous Newton iteration, shape functions and external functions in quadrature points.
      nbs_u->set_quad_order(order);
      DiscontinuousFunc<double>* u = new DiscontinuousFunc<double>(get_fn(fu, ru_actual, nbs_u->get_quad_eo(neighbor_supp_u)),
        neighbor_supp_u, nbs_u->neighbor_edge.orientation);
      nbs_v->set_quad_order(order);
      DiscontinuousFunc<double>* v = new DiscontinuousFunc<double>(get_fn(fv, rv, nbs_v->get_quad_eo(neighbor_supp_v)),
        neighbor_supp_v, nbs_v->neighbor_edge.orientation);

      ExtData<Scalar>* ext = init_ext_fns(mfs->ext, neighbor_searches, order);

      Scalar res = mfs->value(np, jwt, prev, u, v, e, ext);

      // Clean up.
      for (int i = 0; i < prev_size; i++)
      {
        if (prev[i] != NULL)
        {
          prev[i]->free_fn();
          delete prev[i];
        }
      }

      delete [] prev;


      if (ext != NULL)
      {
        ext->free();
        delete ext;
      }

      delete u;
      delete v;
      delete e;

      // Scaling.
      res *= mfs->scaling_factor;

      return 0.5 * res; // Edges are parametrized from 0 to 1 while integration weights
      // are defined in (-1, 1). Thus multiplying with 0.5 to correct
      // the weights.
    }

    template<typename Scalar>
    int DiscreteProblem<Scalar>::calc_order_dg_vector_form(VectorFormSurf<Scalar> *vfs, Hermes::vector<Solution<Scalar>*> u_ext,
      PrecalcShapeset *fv, RefMap *rv, SurfPos* surf_pos,
      LightArray<NeighborSearch<Scalar>*>& neighbor_searches, int neighbor_index_v)
    {
      NeighborSearch<Scalar>* nbs_v = neighbor_searches.get(neighbor_index_v);
      // Hermes::Order that will be returned.
      int order;
      int u_ext_length = u_ext.size();      // Number of external solutions.
      int u_ext_offset = vfs->u_ext_offset; // External solutions will start with u_ext[u_ext_offset]
      // and there will be only u_ext_length - u_ext_offset of them.
      if(this->is_fvm)
        order = rv->get_inv_ref_order();
      else
      {
        // Hermes::Order of solutions from the previous Newton iteration.
        int prev_size = u_ext.size() - vfs->u_ext_offset;
        Func<Hermes::Ord>** oi = new Func<Hermes::Ord>*[prev_size];
        if (u_ext != Hermes::vector<Solution<Scalar>*>())
          for (int i = 0; i < prev_size; i++)
            if (u_ext[i + vfs->u_ext_offset] != NULL)
              oi[i] = init_ext_fn_ord(neighbor_searches.get(u_ext[i]->get_mesh()->get_seq() - min_dg_mesh_seq), u_ext[i]);
            else
              oi[i] = get_fn_ord(0);
        else
          for (int i = 0; i < prev_size; i++)
            oi[i] = get_fn_ord(0);

        // Hermes::Order of the shape function.
        // Determine the integration order.
        int inc = (fv->get_num_components() == 2) ? 1 : 0;
        Func<Hermes::Ord>* ov = get_fn_ord(fv->get_edge_fn_order(surf_pos->surf_num) + inc);

        // Hermes::Order of additional external functions.
        ExtData<Hermes::Ord>* fake_ext = init_ext_fns_ord(vfs->ext, neighbor_searches);

        // Hermes::Order of geometric attributes (eg. for multiplication of a solution with coordinates, normals, etc.).
        Geom<Hermes::Ord>* fake_e = new InterfaceGeom<Hermes::Ord>(&geom_ord,
          nbs_v->neighb_el->marker, nbs_v->neighb_el->id, Hermes::Ord(nbs_v->neighb_el->get_diameter()));
        double fake_wt = 1.0;

        // Total order of the vector form.
        Hermes::Ord o = vfs->ord(1, &fake_wt, oi, ov, fake_e, fake_ext);

        // Increase due to reference map.
        order = rv->get_inv_ref_order();
        order += o.get_order();
        limit_order(order);

        // Clean up.
        if (u_ext != Hermes::vector<Solution<Scalar>*>())
          for (int i = 0; i < prev_size; i++)
            if (u_ext[i + vfs->u_ext_offset] != NULL)
              delete oi[i];
        delete [] oi;
        if (fake_ext != NULL)
        {
          for (int i = 0; i < fake_ext->nf; i++)
          {
            delete fake_ext->fn[i];
          }
          fake_ext->free_ord();
          delete fake_ext;
        }

        delete fake_e;
      }

      return order;
    }

    template<typename Scalar>
    Scalar DiscreteProblem<Scalar>::eval_dg_form(VectorFormSurf<Scalar>* vfs, Hermes::vector<Solution<Scalar>*> u_ext,
      PrecalcShapeset *fv, RefMap *rv,
      SurfPos* surf_pos, LightArray<NeighborSearch<Scalar>*>& neighbor_searches, int neighbor_index_v)
    {
      _F_;

      NeighborSearch<Scalar>* nbs_v = (neighbor_searches.get(neighbor_index_v));
      int order = calc_order_dg_vector_form(vfs, u_ext, fv, rv, surf_pos, neighbor_searches, neighbor_index_v);

      // Evaluate the form using just calculated order.
      Quad2D* quad = fv->get_quad_2d();
      int eo = quad->get_edge_points(surf_pos->surf_num, order);
      int np = quad->get_num_points(eo);
      double3* pt = quad->get_points(eo);

      // A (debug) check.
      assert(surf_pos->surf_num == nbs_v->active_edge);

      // Init geometry and jacobian*weights.
      if (geometry_cache[eo] == NULL)
      {
        geometry_cache[eo] = init_geom_surf(rv, surf_pos->surf_num, surf_pos->marker, eo);
        double3* tan = rv->get_tangent(surf_pos->surf_num, eo);
        jacobian_x_weights_cache[eo] = new double[np];
        for(int i = 0; i < np; i++)
          jacobian_x_weights_cache[eo][i] = pt[i][2] * tan[i][2];
      }

      Geom<double>* e = new InterfaceGeom<double>(geometry_cache[eo], nbs_v->neighb_el->marker,
        nbs_v->neighb_el->id, nbs_v->neighb_el->get_diameter());
      double* jwt = jacobian_x_weights_cache[eo];

      // Values of the previous Newton iteration, shape functions and external functions in quadrature points.
      int prev_size = u_ext.size() - vfs->u_ext_offset;
      Func<Scalar>** prev = new Func<Scalar>*[prev_size];
      if (u_ext != Hermes::vector<Solution<Scalar>*>())
        for (int i = 0; i < prev_size; i++)
          if (u_ext[i + vfs->u_ext_offset] != NULL)
          {
            neighbor_searches.get(u_ext[i]->get_mesh()->get_seq() - min_dg_mesh_seq)->set_quad_order(order);
            prev[i]  = neighbor_searches.get(u_ext[i]->get_mesh()->get_seq() - min_dg_mesh_seq)->init_ext_fn(u_ext[i]);
          }
          else prev[i] = NULL;
      else
        for (int i = 0; i < prev_size; i++)
          prev[i] = NULL;

      Func<double>* v = get_fn(fv, rv, eo);
      ExtData<Scalar>* ext = init_ext_fns(vfs->ext, neighbor_searches, order);

      Scalar res = vfs->value(np, jwt, prev, v, e, ext) * vfs->scaling_factor;

      // Clean up.
      for (int i = 0; i < prev_size; i++)
      {
        if (prev[i] != NULL)
        {
          prev[i]->free_fn();
          delete prev[i];
        }
      }

      delete [] prev;

      if (ext != NULL)
      {
        ext->free();
        delete ext;
      }

      delete e;

      // Scaling.
      res *= vfs->scaling_factor;

      return 0.5 * res; // Edges are parametrized from 0 to 1 while integration weights
      // are defined in (-1, 1). Thus multiplying with 0.5 to correct
      // the weights.
    }

    NeighborNode::NeighborNode(NeighborNode* parent, unsigned int transformation) : parent(parent), transformation(transformation)
    {
      left_son = right_son = NULL;
    }

    NeighborNode::~NeighborNode()
    {
      if(left_son != NULL)
      {
        delete left_son;
        left_son = NULL;
      }
      if(right_son != NULL)
      {
        delete right_son;
        right_son = NULL;
      }
    }

    void NeighborNode::set_left_son(NeighborNode* left_son)
    {
      this->left_son = left_son;
    }
    void NeighborNode::set_right_son(NeighborNode* right_son)
    {
      this->right_son = right_son;
    }
    void NeighborNode::set_transformation(unsigned int transformation)
    {
      this->transformation = transformation;
    }
    NeighborNode* NeighborNode::get_left_son()
    {
      return left_son;
    }
    NeighborNode* NeighborNode::get_right_son()
    {
      return right_son;
    }
    unsigned int NeighborNode::get_transformation()
    {
      return this->transformation;
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
