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
    DiscreteProblem<Scalar>::DiscreteProblem() : wf(NULL)
    {
      // Set all attributes for which we don't need to acces wf or spaces.
      // This is important for the destructor to properly detect what needs to be deallocated.
      sp_seq = NULL;
      is_fvm = false;
      RungeKutta = false;
      RK_original_spaces_count = 0;
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
      have_matrix = false;

      // There is a special function that sets a DiscreteProblem to be FVM.
      // Purpose is that this constructor looks cleaner and is simpler.
      this->is_fvm = false;

      Geom<Hermes::Ord> *tmp = init_geom_ord();
      geom_ord = *tmp;
      delete tmp;

      current_mat = NULL;
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
    void DiscreteProblem<Scalar>::init_assembling(Scalar* coeff_vec, PrecalcShapeset*** pss , PrecalcShapeset*** spss, RefMap*** refmaps, Solution<Scalar>*** u_ext, AsmList<Scalar>*** als, Hermes::vector<MeshFunction<Scalar>*>& ext_functions, MeshFunction<Scalar>*** ext, 
          Hermes::vector<MatrixFormVol<Scalar>*>* mfvol, Hermes::vector<MatrixFormSurf<Scalar>*>* mfsurf, Hermes::vector<VectorFormVol<Scalar>*>* vfvol, Hermes::vector<VectorFormSurf<Scalar>*>* vfsurf)
    {
      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
      {
        pss[i] = new PrecalcShapeset*[wf->get_neq()];
        for (unsigned int j = 0; j < wf->get_neq(); j++)
          pss[i][j] = new PrecalcShapeset(spaces[j]->shapeset);
      }
      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
      {
        spss[i] = new PrecalcShapeset*[wf->get_neq()];
        for (unsigned int j = 0; j < wf->get_neq(); j++)
          spss[i][j] = new PrecalcShapeset(pss[i][j]);
      }
      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
      {
        refmaps[i] = new RefMap*[wf->get_neq()];
        for (unsigned int j = 0; j < wf->get_neq(); j++)
        {
          refmaps[i][j] = new RefMap();
          refmaps[i][j]->set_quad_2d(&g_quad_2d_std);
        }
      }
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
        else
          u_ext[i] = NULL;
      }
      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
      {
        als[i] = new AsmList<Scalar>*[wf->get_neq()];
        for (unsigned int j = 0; j < wf->get_neq(); j++)
          als[i][j] = new AsmList<Scalar>();
      }
      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
      {
        ext[i] = new MeshFunction<Scalar>*[ext_functions.size()];
        for (int j = 0; j < ext_functions.size(); j++)
          ext[i][j] = ext_functions[j]->clone();
      }
      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
      {
        for (int j = 0; j < wf->mfvol.size(); j++)
        {
          mfvol[i].push_back(wf->mfvol[j]->clone());
          // Inserting proper ext.
          for(int k = 0; k < wf->mfvol[j]->ext.size(); k++)
          {
            for (int l = 0; l < ext_functions.size(); l++)
            {
              if(ext_functions[l] == wf->mfvol[j]->ext[k])
              {
                while(k >= mfvol[i][j]->ext.size())
                  mfvol[i][j]->ext.push_back(NULL);
                mfvol[i][j]->ext[k] = ext[i][l];
                break;
              }
            }
          }
        }
      }
      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
      {
        for (int j = 0; j < wf->mfsurf.size(); j++)
        {
          mfsurf[i].push_back(wf->mfsurf[j]->clone());
          // Inserting proper ext.
          for(int k = 0; k < wf->mfsurf[j]->ext.size(); k++)
          {
            for (int l = 0; l < ext_functions.size(); l++)
            {
              if(ext_functions[l] == wf->mfsurf[j]->ext[k])
              {
                while(k >= mfsurf[i][j]->ext.size())
                  mfsurf[i][j]->ext.push_back(NULL);
                mfsurf[i][j]->ext[k] = ext[i][l];
                break;
              }
            }
          }
        }
      }
      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
      {
        for (int j = 0; j < wf->vfvol.size(); j++)
        {
          vfvol[i].push_back(wf->vfvol[j]->clone());
          // Inserting proper ext.
          for(int k = 0; k < wf->vfvol[j]->ext.size(); k++)
          {
            for (int l = 0; l < ext_functions.size(); l++)
            {
              if(ext_functions[l] == wf->vfvol[j]->ext[k])
              {
                while(k >= vfvol[i][j]->ext.size())
                  vfvol[i][j]->ext.push_back(NULL);

                vfvol[i][j]->ext[k] = ext[i][l];
                break;
              }
            }
          }
        }
      }
      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
      {
        for (int j = 0; j < wf->vfsurf.size(); j++)
        {
          vfsurf[i].push_back(wf->vfsurf[j]->clone());
          // Inserting proper ext.
          for(int k = 0; k < wf->vfsurf[j]->ext.size(); k++)
          {
            for (int l = 0; l < ext_functions.size(); l++)
            {
              if(ext_functions[l] == wf->vfsurf[j]->ext[k])
              {
                while(k >= vfsurf[i][j]->ext.size())
                  vfsurf[i][j]->ext.push_back(NULL);
                vfsurf[i][j]->ext[k] = ext[i][l];
                break;
              }
            }
          }
        }
      }
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::deinit_assembling(PrecalcShapeset*** pss , PrecalcShapeset*** spss, RefMap*** refmaps, Solution<Scalar>*** u_ext, AsmList<Scalar>*** als, Hermes::vector<MeshFunction<Scalar>*>& ext_functions, MeshFunction<Scalar>*** ext, 
          Hermes::vector<MatrixFormVol<Scalar>*>* mfvol, Hermes::vector<MatrixFormSurf<Scalar>*>* mfsurf, Hermes::vector<VectorFormVol<Scalar>*>* vfvol, Hermes::vector<VectorFormSurf<Scalar>*>* vfsurf)
    {
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
        if(u_ext[i] != NULL)
        {
          for (unsigned int j = 0; j < wf->get_neq(); j++)
            delete u_ext[i][j];
          delete [] u_ext[i];
        }
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
        for (unsigned int j = 0; j < ext_functions.size(); j++)
          delete ext[i][j];
        delete [] ext[i];
      }
      delete [] ext;

      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
      {
        for (int j = 0; j < wf->mfvol.size(); j++)
          delete mfvol[i][j];
        mfvol[i].clear();
      }
      delete [] mfvol;
      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
      {
        for (int j = 0; j < wf->mfsurf.size(); j++)
          delete mfsurf[i][j];
        mfsurf[i].clear();
      }
      delete [] mfsurf;
      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
      {
        for (int j = 0; j < wf->vfvol.size(); j++)
          delete vfvol[i][j];
        vfvol[i].clear();
      }
      delete [] vfvol;
      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
      {
        for (int j = 0; j < wf->vfsurf.size(); j++)
          delete vfsurf[i][j];
        vfsurf[i].clear();
      }
      delete [] vfsurf;
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::assemble(Scalar* coeff_vec, SparseMatrix<Scalar>* mat,
      Vector<Scalar>* rhs,
      bool force_diagonal_blocks,
      Table* block_weights)
    {
      _F_;

      current_mat = mat;
      current_rhs = rhs;
      current_force_diagonal_blocks = force_diagonal_blocks;
      current_block_weights = block_weights;

      // Check that the block scaling table have proper dimension.
      if (block_weights != NULL)
        if (block_weights->get_size() != wf->get_neq())
          throw Exceptions::LengthException(6, block_weights->get_size(), wf->get_neq());

      // Creating matrix sparse structure.
      create_sparse_structure();

      Hermes::vector<MeshFunction<Scalar>*> ext_functions;
      for(unsigned int form_i = 0; form_i < wf->mfvol.size(); form_i++)
        for(unsigned int ext_i = 0; ext_i < wf->mfvol.at(form_i)->ext.size(); ext_i++)
          ext_functions.push_back(wf->mfvol.at(form_i)->ext[ext_i]);
      for(unsigned int form_i = 0; form_i < wf->mfsurf.size(); form_i++)
        for(unsigned int ext_i = 0; ext_i < wf->mfsurf.at(form_i)->ext.size(); ext_i++)
          ext_functions.push_back(wf->mfsurf.at(form_i)->ext[ext_i]);
      for(unsigned int form_i = 0; form_i < wf->vfvol.size(); form_i++)
        for(unsigned int ext_i = 0; ext_i < wf->vfvol.at(form_i)->ext.size(); ext_i++)
          ext_functions.push_back(wf->vfvol.at(form_i)->ext[ext_i]);
      for(unsigned int form_i = 0; form_i < wf->vfsurf.size(); form_i++)
        for(unsigned int ext_i = 0; ext_i < wf->vfsurf.at(form_i)->ext.size(); ext_i++)
          ext_functions.push_back(wf->vfsurf.at(form_i)->ext[ext_i]);

      // Structures that cloning will be done into.
      PrecalcShapeset*** pss = new PrecalcShapeset**[omp_get_max_threads()];
      PrecalcShapeset*** spss = new PrecalcShapeset**[omp_get_max_threads()];
      RefMap*** refmaps = new RefMap**[omp_get_max_threads()];
      Solution<Scalar>*** u_ext = new Solution<Scalar>**[omp_get_max_threads()];
      AsmList<Scalar>*** als = new AsmList<Scalar>**[omp_get_max_threads()];
      MeshFunction<Scalar>*** ext = new MeshFunction<Scalar>**[omp_get_max_threads()];
      Hermes::vector<MatrixFormVol<Scalar>*>* mfvol = new Hermes::vector<MatrixFormVol<Scalar>*>[omp_get_max_threads()];
      Hermes::vector<MatrixFormSurf<Scalar>*>* mfsurf = new Hermes::vector<MatrixFormSurf<Scalar>*>[omp_get_max_threads()];
      Hermes::vector<VectorFormVol<Scalar>*>* vfvol = new Hermes::vector<VectorFormVol<Scalar>*>[omp_get_max_threads()];
      Hermes::vector<VectorFormSurf<Scalar>*>* vfsurf = new Hermes::vector<VectorFormSurf<Scalar>*>[omp_get_max_threads()];

      // Fill these structures.
      init_assembling(coeff_vec, pss, spss, refmaps, u_ext, als, ext_functions, ext, mfvol, mfsurf, vfvol, vfsurf);

      // Vector of meshes.
      Hermes::vector<Mesh*> meshes;
      for(unsigned int space_i = 0; space_i < spaces.size(); space_i++)
        meshes.push_back(spaces[space_i]->get_mesh());
      for (unsigned j = 0; j < ext_functions.size(); j++)
        meshes.push_back(ext_functions[j]->get_mesh());
      if (coeff_vec != NULL)
        for(unsigned int space_i = 0; space_i < spaces.size(); space_i++)
          meshes.push_back(spaces[space_i]->get_mesh());

      Traverse trav_master(true);
      unsigned int num_states = trav_master.get_num_states(meshes);
      
      trav_master.begin(meshes.size(), &(meshes.front()));

      Traverse* trav = new Traverse[omp_get_max_threads()];
      Hermes::vector<Transformable *>* fns = new Hermes::vector<Transformable *>[omp_get_max_threads()];
      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
      {
        for (unsigned j = 0; j < spaces.size(); j++)
          fns[i].push_back(pss[i][j]);
        for (unsigned j = 0; j < ext_functions.size(); j++)
        {
          fns[i].push_back(ext[i][j]);
          ext[i][j]->set_quad_2d(&g_quad_2d_std);
        }
        if (coeff_vec != NULL)
          for (unsigned j = 0; j < wf->get_neq(); j++)
          {
            fns[i].push_back(u_ext[i][j]);
            u_ext[i][j]->set_quad_2d(&g_quad_2d_std);
          }
        trav[i].begin(meshes.size(), &(meshes.front()), &(fns[i].front()));
        trav[i].stack = trav_master.stack;
      }

int state_i;

PrecalcShapeset** current_pss;
PrecalcShapeset** current_spss;
RefMap** current_refmaps;
Solution<Scalar>** current_u_ext;
AsmList<Scalar>** current_als;

MatrixFormVol<Scalar>** current_mfvol;
MatrixFormSurf<Scalar>** current_mfsurf;
VectorFormVol<Scalar>** current_vfvol;
VectorFormSurf<Scalar>** current_vfsurf;

#define CHUNKSIZE 1
#pragma omp parallel shared(trav_master, mat, rhs) private(state_i, current_pss, current_spss, current_refmaps, current_u_ext, current_als, current_mfvol, current_mfsurf, current_vfvol, current_vfsurf)
      {
        #pragma omp for schedule(dynamic, CHUNKSIZE)
        for(state_i = 0; state_i < num_states; state_i++)
        {
          Traverse::State current_state;
          #pragma omp critical (get_next_state)
            current_state = trav[omp_get_thread_num()].get_next_state(&trav_master.top, &trav_master.id);

          current_pss = pss[omp_get_thread_num()];
          current_spss = spss[omp_get_thread_num()];
          current_refmaps = refmaps[omp_get_thread_num()];
          current_u_ext = u_ext[omp_get_thread_num()];
          current_als = als[omp_get_thread_num()];

          current_mfvol = mfvol[omp_get_thread_num()].size() == 0 ? NULL : &(mfvol[omp_get_thread_num()].front());
          current_mfsurf = mfsurf[omp_get_thread_num()].size() == 0 ? NULL : &(mfsurf[omp_get_thread_num()].front());
          current_vfvol = vfvol[omp_get_thread_num()].size() == 0 ? NULL : &(vfvol[omp_get_thread_num()].front());
          current_vfsurf = vfsurf[omp_get_thread_num()].size() == 0 ? NULL : &(vfsurf[omp_get_thread_num()].front());

          // One state is a collection of (virtual) elements sharing
          // the same physical location on (possibly) different meshes.
          // This is then the same element of the virtual union mesh.
          // The proper sub-element mappings to all the functions of
          // this stage is supplied by the function Traverse::get_next_state()
          // called in the while loop.
          assemble_one_state(current_pss, current_spss, current_refmaps, current_u_ext, current_als, &current_state, current_mfvol, current_mfsurf, current_vfvol, current_vfsurf);
        }
      }

      deinit_assembling(pss, spss, refmaps, u_ext, als, ext_functions, ext, mfvol, mfsurf, vfvol, vfsurf);
      
      trav_master.finish();
      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
        trav[i].finish();
      
      for(unsigned int i = 0; i < omp_get_max_threads(); i++)
      {
        fns[i].clear();
      }
      delete [] fns;
      delete [] trav;

      /// \todo Should this be really here? Or in assemble()?
      if (current_mat != NULL)
        current_mat->finish();
      if (current_rhs != NULL)
        current_rhs->finish();

      if(DG_matrix_forms_present || DG_vector_forms_present)
      {
        Element* element_to_set_nonvisited;
        for(unsigned int mesh_i = 0; mesh_i < meshes.size(); mesh_i++)
          for_all_elements(element_to_set_nonvisited, meshes[mesh_i])
          element_to_set_nonvisited->visited = false;
      }
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
      for(unsigned int i = 0; i < wf->mfsurf.size() && DG_matrix_forms_present == false; i++)
      {
        if (wf->mfsurf[i]->areas[0] == H2D_DG_INNER_EDGE)
        {
          DG_matrix_forms_present = true;
          break;
        }
      }
      for(unsigned int i = 0; i < wf->vfsurf.size() && DG_vector_forms_present == false; i++)
      {
        if (wf->vfsurf[i]->areas[0] == H2D_DG_INNER_EDGE)
        {
          DG_vector_forms_present = true;
          break;
        }
      }
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::init_state(PrecalcShapeset** current_pss, PrecalcShapeset** current_spss, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, AsmList<Scalar>** current_als, Traverse::State* current_state)
    {
      _F_;

      // Obtain assembly lists for the element at all spaces of the stage, set appropriate mode for each pss.
      // NOTE: Active elements and transformations for external functions (including the solutions from previous
      // Newton's iteration) as well as basis functions (master PrecalcShapesets) have already been set in
      // trav.get_next_state(...).
      for (unsigned int i = 0; i < spaces.size(); i++)
      {
        if (current_state->e[i] == NULL)
          continue;

        // \todo do not obtain again if the element was not changed.
        spaces[i]->get_element_assembly_list(current_state->e[i], current_als[i], spaces_first_dofs[i]);

        // Set active element to all test functions.
        current_spss[i]->set_active_element(current_state->e[i]);
        current_spss[i]->set_master_transform();

        // Set active element to reference mappings.
        current_refmaps[i]->set_active_element(current_state->e[i]);
        current_refmaps[i]->force_transform(current_pss[i]->get_transform(), current_pss[i]->get_ctm());

        // Mark the active element on each mesh in order to prevent assembling on its edges from the other side.
        if(DG_matrix_forms_present || DG_vector_forms_present)
          current_state->e[i]->visited = true;
      }
      return;
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::init_surface_state(AsmList<Scalar>** current_als, Traverse::State* current_state)
    {
      _F_;
      // Obtain the list of shape functions which are nonzero on this surface.
      for (unsigned int i = 0; i < spaces.size(); i++)
      {
        if (current_state->e[i] == NULL)
          continue;
        
        spaces[i]->get_boundary_assembly_list(current_state->e[i], current_state->isurf, current_als[i], spaces_first_dofs[i]);
      }
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::assemble_one_state(PrecalcShapeset** current_pss, PrecalcShapeset** current_spss, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, AsmList<Scalar>** current_als, Traverse::State* current_state,
        MatrixFormVol<Scalar>** current_mfvol, MatrixFormSurf<Scalar>** current_mfsurf, VectorFormVol<Scalar>** current_vfvol, VectorFormSurf<Scalar>** current_vfsurf)
    {
      _F_;

      // Initialize the state, return a non-NULL element; if no such element found, return.
      init_state(current_pss, current_spss, current_refmaps, current_u_ext, current_als, current_state);
      
      if (current_mat != NULL)
      {
        for(int current_mfvol_i = 0; current_mfvol_i < wf->mfvol.size(); current_mfvol_i++)
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
              test_fns[i] = init_fn(current_spss[current_mfvol[current_mfvol_i]->i], current_refmaps[current_mfvol[current_mfvol_i]->i], order);
            }
          }

          for (unsigned int j = 0; j < current_als[current_mfvol[current_mfvol_i]->j]->cnt; j++)
          {
            if (std::abs(current_als[current_mfvol[current_mfvol_i]->j]->coef[j]) < 1e-12)
                continue;
            if (current_als[current_mfvol[current_mfvol_i]->j]->dof[j] >= 0)
            {
              current_pss[current_mfvol[current_mfvol_i]->j]->set_active_shape(current_als[current_mfvol[current_mfvol_i]->j]->idx[j]);

              base_fns[j] = init_fn(current_pss[current_mfvol[current_mfvol_i]->j], current_refmaps[current_mfvol[current_mfvol_i]->j], order);

            }
          }
          
          assemble_matrix_form(current_mfvol[current_mfvol_i], order, base_fns, test_fns, current_refmaps, current_u_ext, current_als, current_state);

          for (unsigned int j = 0; j < current_als[current_mfvol[current_mfvol_i]->j]->cnt; j++)
            if (std::abs(current_als[current_mfvol[current_mfvol_i]->j]->coef[j]) >= 1e-12)
              if (current_als[current_mfvol[current_mfvol_i]->j]->dof[j] >= 0)
              {
                base_fns[j]->free_fn();
                delete base_fns[j];
              }
          delete [] base_fns;
          for (unsigned int i = 0; i < current_als[current_mfvol[current_mfvol_i]->i]->cnt; i++)
            if (std::abs(current_als[current_mfvol[current_mfvol_i]->i]->coef[i]) >= 1e-12)
              if (current_als[current_mfvol[current_mfvol_i]->i]->dof[i] >= 0)
              {
                test_fns[i]->free_fn();
                delete test_fns[i];
              }
          delete [] test_fns;
        }
      }
      if (current_rhs != NULL)
      {
        for(int current_vfvol_i = 0; current_vfvol_i < wf->vfvol.size(); current_vfvol_i++)
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

              test_fns[i] = init_fn(current_spss[current_vfvol[current_vfvol_i]->i], current_refmaps[current_vfvol[current_vfvol_i]->i], order);
            }
          }

          assemble_vector_form(current_vfvol[current_vfvol_i], order, test_fns, current_refmaps, current_u_ext, current_als, current_state);

          for (unsigned int i = 0; i < current_als[current_vfvol[current_vfvol_i]->i]->cnt; i++)
            if (std::abs(current_als[current_vfvol[current_vfvol_i]->i]->coef[i]) >= 1e-12)
              if (current_als[current_vfvol[current_vfvol_i]->i]->dof[i] >= 0)
              {
                test_fns[i]->free_fn();
                delete test_fns[i];
              }
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
          for(int current_mfsurf_i = 0; current_mfsurf_i < wf->mfsurf.size(); current_mfsurf_i++)
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

                test_fns[i] = init_fn(current_spss[current_mfsurf[current_mfsurf_i]->i], current_refmaps[current_mfsurf[current_mfsurf_i]->i], current_refmaps[current_mfsurf[current_mfsurf_i]->i]->get_quad_2d()->get_edge_points(current_state->isurf, order, current_state->e[0]->get_mode()));
              }
            }

            for (unsigned int j = 0; j < current_als[current_mfsurf[current_mfsurf_i]->j]->cnt; j++)
            {
              if (std::abs(current_als[current_mfsurf[current_mfsurf_i]->j]->coef[j]) < 1e-12)
                continue;
              if (current_als[current_mfsurf[current_mfsurf_i]->j]->dof[j] >= 0)
              {
                current_pss[current_mfsurf[current_mfsurf_i]->j]->set_active_shape(current_als[current_mfsurf[current_mfsurf_i]->j]->idx[j]);

                base_fns[j] = init_fn(current_pss[current_mfsurf[current_mfsurf_i]->j], current_refmaps[current_mfsurf[current_mfsurf_i]->j], current_refmaps[current_mfsurf[current_mfsurf_i]->j]->get_quad_2d()->get_edge_points(current_state->isurf, order,current_state->e[0]->get_mode()));
              }
            }
            
            assemble_matrix_form(current_mfsurf[current_mfsurf_i], order, base_fns, test_fns, current_refmaps, current_u_ext, current_als, current_state);

            for (unsigned int j = 0; j < current_als[current_mfsurf[current_mfsurf_i]->j]->cnt; j++)
            if (std::abs(current_als[current_mfsurf[current_mfsurf_i]->j]->coef[j]) >= 1e-12)
              if (current_als[current_mfsurf[current_mfsurf_i]->j]->dof[j] >= 0)
              {
                base_fns[j]->free_fn();
                delete base_fns[j];
              }
          delete [] base_fns;
          for (unsigned int i = 0; i < current_als[current_mfsurf[current_mfsurf_i]->i]->cnt; i++)
            if (std::abs(current_als[current_mfsurf[current_mfsurf_i]->i]->coef[i]) >= 1e-12)
              if (current_als[current_mfsurf[current_mfsurf_i]->i]->dof[i] >= 0)
              {
                test_fns[i]->free_fn();
                delete test_fns[i];
              }
          delete [] test_fns;
          }
        }
    
        if (current_rhs != NULL)
        {
          for(int current_vfsurf_i = 0; current_vfsurf_i < wf->vfsurf.size(); current_vfsurf_i++)
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

                test_fns[i] = init_fn(current_spss[current_vfsurf[current_vfsurf_i]->i], current_refmaps[current_vfsurf[current_vfsurf_i]->i], current_refmaps[current_vfsurf[current_vfsurf_i]->i]->get_quad_2d()->get_edge_points(current_state->isurf, order, current_state->e[0]->get_mode()));
              }
            }
            
            assemble_vector_form(current_vfsurf[current_vfsurf_i], order, test_fns, current_refmaps, current_u_ext, current_als, current_state);

          for (unsigned int i = 0; i < current_als[current_vfsurf[current_vfsurf_i]->i]->cnt; i++)
            if (std::abs(current_als[current_vfsurf[current_vfsurf_i]->i]->coef[i]) >= 1e-12)
              if (current_als[current_vfsurf[current_vfsurf_i]->i]->dof[i] >= 0)
              {
                test_fns[i]->free_fn();
                delete test_fns[i];
              }
          delete [] test_fns;
          }
        }
      }
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
        init_ext_orders(form, u_ext_ord, &ext_ord, current_u_ext, current_state);

        // Order of shape functions.
        int max_order_j = this->spaces[form->j]->get_element_order(current_state->e[form->j]->id);
        int max_order_i = this->spaces[form->i]->get_element_order(current_state->e[form->i]->id);
        if(H2D_GET_V_ORDER(max_order_i) > H2D_GET_H_ORDER(max_order_i))
          max_order_i = H2D_GET_V_ORDER(max_order_i);
        else
          max_order_i = H2D_GET_H_ORDER(max_order_i);
        if(H2D_GET_V_ORDER(max_order_j) > H2D_GET_H_ORDER(max_order_j))
          max_order_j = H2D_GET_V_ORDER(max_order_j);
        else
          max_order_j = H2D_GET_H_ORDER(max_order_j);

        for (unsigned int k = 0; k < current_state->rep->get_num_surf(); k++)
        {
          int eo = this->spaces[form->i]->get_edge_order(current_state->e[form->i], k);
          if (eo > max_order_i) 
            max_order_i = eo;
          eo = this->spaces[form->j]->get_edge_order(current_state->e[form->j], k);
          if (eo > max_order_j) 
            max_order_j = eo;
        }

        Func<Hermes::Ord>* ou = init_fn_ord(max_order_j + (spaces[form->j]->get_shapeset()->get_num_components() > 1 ? 1 : 0));
        Func<Hermes::Ord>* ov = init_fn_ord(max_order_i + (spaces[form->i]->get_shapeset()->get_num_components() > 1 ? 1 : 0));

        // Total order of the vector form.
        Hermes::Ord o = form->ord(1, &fake_wt, u_ext_ord, ou, ov, &geom_ord, &ext_ord);

        adjust_order_to_refmaps(form, order, &o, current_refmaps);

        // Cleanup.
        deinit_ext_orders(form, u_ext_ord, &ext_ord);
        ou->free_ord();
        delete ou;
        ov->free_ord();
        delete ov;
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
      Scalar **local_stiffness_matrix = new_matrix<Scalar>(std::max(current_als[form->i]->cnt, current_als[form->j]->cnt));

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
      Geom<double>* geometry = NULL;
      double* jacobian_x_weights = NULL;
      if(surface_form)
        n_quadrature_points = init_surface_geometry_points(current_refmaps[form->i], order, current_state, geometry, jacobian_x_weights);
      else
        n_quadrature_points = init_geometry_points(current_refmaps[form->i], order, geometry, jacobian_x_weights);

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
                local_stiffness_matrix[i][j] = 0.5 * block_scaling_coeff(form) * form->value(n_quadrature_points, jacobian_x_weights, u_ext, u, v, geometry, &ext) * form->scaling_factor * current_als[form->j]->coef[j] * current_als[form->i]->coef[i];
              else
                local_stiffness_matrix[i][j] = block_scaling_coeff(form) * form->value(n_quadrature_points, jacobian_x_weights, u_ext, u, v, geometry, &ext) * form->scaling_factor * current_als[form->j]->coef[j] * current_als[form->i]->coef[i];
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

              Scalar val = block_scaling_coeff(form) * form->value(n_quadrature_points, jacobian_x_weights, u_ext, u, v, geometry, &ext) * form->scaling_factor * current_als[form->j]->coef[j] * current_als[form->i]->coef[i];

              local_stiffness_matrix[i][j] = local_stiffness_matrix[j][i] = val;
            }
          }
        }
      }

      // Insert the local stiffness matrix into the global one.
#pragma omp critical (mat)
      current_mat->add(current_als[form->i]->cnt, current_als[form->j]->cnt, local_stiffness_matrix, current_als[form->i]->dof, current_als[form->j]->dof);

      // Insert also the off-diagonal (anti-)symmetric block, if required.
      if (tra)
      {
        if (form->sym < 0)
          chsgn(local_stiffness_matrix, current_als[form->i]->cnt, current_als[form->j]->cnt);
        transpose(local_stiffness_matrix, current_als[form->i]->cnt, current_als[form->j]->cnt);
#pragma omp critical (mat)
        current_mat->add(current_als[form->j]->cnt, current_als[form->i]->cnt, local_stiffness_matrix, current_als[form->j]->dof, current_als[form->i]->dof);
      }

      // Cleanup.
      deinit_ext(form, u_ext, &ext);
      delete [] local_stiffness_matrix;
      delete [] jacobian_x_weights;
      geometry->free();
      delete geometry;
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
        init_ext_orders(form, u_ext_ord, &ext_ord, current_u_ext, current_state);

        // Order of shape functions.
        int max_order_i = this->spaces[form->i]->get_element_order(current_state->e[form->i]->id);
        if(H2D_GET_V_ORDER(max_order_i) > H2D_GET_H_ORDER(max_order_i))
          max_order_i = H2D_GET_V_ORDER(max_order_i);
        else
          max_order_i = H2D_GET_H_ORDER(max_order_i);

        for (unsigned int k = 0; k < current_state->rep->get_num_surf(); k++)
        {
          int eo = this->spaces[form->i]->get_edge_order(current_state->e[form->i], k);
          if (eo > max_order_i) 
            max_order_i = eo;
        }
        Func<Hermes::Ord>* ov = init_fn_ord(max_order_i + (spaces[form->i]->get_shapeset()->get_num_components() > 1 ? 1 : 0));

        // Total order of the vector form.
        Hermes::Ord o = form->ord(1, &fake_wt, u_ext_ord, ov, &geom_ord, &ext_ord);

        adjust_order_to_refmaps(form, order, &o, current_refmaps);
        
        // Cleanup.
        deinit_ext_orders(form, u_ext_ord, &ext_ord);
        ov->free_ord();
        delete ov;
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
      Geom<double>* geometry = NULL;
      double* jacobian_x_weights = NULL;

      if(surface_form)
        n_quadrature_points = init_surface_geometry_points(current_refmaps[form->i], order, current_state, geometry, jacobian_x_weights);
      else
        n_quadrature_points = init_geometry_points(current_refmaps[form->i], order, geometry, jacobian_x_weights);

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
        
        Scalar val;
        if(surface_form)
          val = 0.5 * form->value(n_quadrature_points, jacobian_x_weights, u_ext, v, geometry, &ext) * form->scaling_factor * current_als[form->i]->coef[i];
        else
          val = form->value(n_quadrature_points, jacobian_x_weights, u_ext, v, geometry, &ext) * form->scaling_factor * current_als[form->i]->coef[i];
#pragma omp critical (rhs)
        current_rhs->add(current_als[form->i]->dof[i], val);
      }

      // Cleanup.
      deinit_ext(form, u_ext, &ext);
      delete [] jacobian_x_weights;
      geometry->free();
      delete geometry;
    }

    template<typename Scalar>
    int DiscreteProblem<Scalar>::init_geometry_points(RefMap* reference_mapping, int order, Geom<double>*& geometry, double*& jacobian_x_weights)
    {
      _F_;
      double3* pt = reference_mapping->get_quad_2d()->get_points(order, reference_mapping->get_active_element()->get_mode());
      int np = reference_mapping->get_quad_2d()->get_num_points(order, reference_mapping->get_active_element()->get_mode());

      // Init geometry and jacobian*weights.
      geometry = init_geom_vol(reference_mapping, order);
      double* jac = NULL;
      if(!reference_mapping->is_jacobian_const())
        jac = reference_mapping->get_jacobian(order);
      jacobian_x_weights = new double[np];
      for(int i = 0; i < np; i++)
      {
        if(reference_mapping->is_jacobian_const())
          jacobian_x_weights[i] = pt[i][2] * reference_mapping->get_const_jacobian();
        else
          jacobian_x_weights[i] = pt[i][2] * jac[i];
      }
      return np;
    }

    template<typename Scalar>
    int DiscreteProblem<Scalar>::init_surface_geometry_points(RefMap* reference_mapping, int& order, Traverse::State* current_state, Geom<double>*& geometry, double*& jacobian_x_weights)
    {
      _F_;
      int eo = reference_mapping->get_quad_2d()->get_edge_points(current_state->isurf, order, reference_mapping->get_active_element()->get_mode());
      double3* pt = reference_mapping->get_quad_2d()->get_points(eo, reference_mapping->get_active_element()->get_mode());
      int np = reference_mapping->get_quad_2d()->get_num_points(eo, reference_mapping->get_active_element()->get_mode());

      // Init geometry and jacobian*weights.
      double3* tan;
      geometry = init_geom_surf(reference_mapping, current_state->isurf, current_state->rep->marker, eo, tan);
      jacobian_x_weights = new double[np];
      for(int i = 0; i < np; i++)
        jacobian_x_weights[i] = pt[i][2] * tan[i][2];
      
      order = eo;
      return np;
    }
    
    template<typename Scalar>
    void DiscreteProblem<Scalar>::init_ext_orders(Form<Scalar> *form, Func<Hermes::Ord>** oi, ExtData<Hermes::Ord>* oext, Solution<Scalar>** current_u_ext, Traverse::State* current_state)
    {
      _F_;
      unsigned int prev_size = RungeKutta ? RK_original_spaces_count : this->wf->get_neq() - form->u_ext_offset;
      bool surface_form = (current_state->isurf > -1);

      if (current_u_ext != NULL)
        for(int i = 0; i < prev_size; i++)
          if (current_u_ext[i + form->u_ext_offset] != NULL)
            if(surface_form)
              oi[i] = init_fn_ord(current_u_ext[i + form->u_ext_offset]->get_edge_fn_order(current_state->isurf) + (current_u_ext[i + form->u_ext_offset]->get_num_components() > 1 ? 1 : 0));
            else
              oi[i] = init_fn_ord(current_u_ext[i + form->u_ext_offset]->get_fn_order() + (current_u_ext[i + form->u_ext_offset]->get_num_components() > 1 ? 1 : 0));
          else

          oi[i] = init_fn_ord(0);
      else
        for(int i = 0; i < prev_size; i++)
          oi[i] = init_fn_ord(0);
      
      oext->nf = form->ext.size();
      oext->fn = new Func<Hermes::Ord>*[oext->nf];
      for (int i = 0; i < oext->nf; i++)
        if(surface_form)
          oext->fn[i] = init_fn_ord(form->ext[i]->get_edge_fn_order(current_state->isurf) + (form->ext[i]->get_num_components() > 1 ? 1 : 0));
        else
          oext->fn[i] = init_fn_ord(form->ext[i]->get_fn_order() + (form->ext[i]->get_num_components() > 1 ? 1 : 0));
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
      limit_order(order, current_refmaps[form->i]->get_active_element()->get_mode());
    }

    template class HERMES_API DiscreteProblem<double>;
    template class HERMES_API DiscreteProblem<std::complex<double> >;
  }
}
