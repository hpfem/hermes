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

#include "h2d_common.h"
#include "integrals/integrals_h1.h"
#include "quadrature/limit_order.h"
#include "discrete_problem.h"
#include "mesh/traverse.h"
#include "space/space.h"
#include "shapeset/precalc.h"
#include "../../hermes_common/matrix.h"
#include "../../hermes_common/solver/umfpack_solver.h"
#include "mesh/refmap.h"
#include "function/solution.h"
#include "config.h"
#include "neighbor.h"
#include "views/scalar_view.h"
#include "views/base_view.h"
#include "boundaryconditions/essential_bcs.h"

DiscreteProblem::DiscreteProblem(WeakForm* wf, Hermes::vector<Space *> spaces, 
         bool is_linear) : wf(wf), is_linear(is_linear), wf_seq(-1), spaces(spaces)
{
  _F_
  init();
}

DiscreteProblem::DiscreteProblem(WeakForm* wf, Space* space, bool is_linear)
   : wf(wf), is_linear(is_linear), wf_seq(-1)
{
  _F_
  spaces.push_back(space);
  init();
}

void DiscreteProblem::init()
{
  _F_
  // Sanity checks.
  if(wf == NULL)
    error("WeakForm* wf can not be NULL in DiscreteProblem::DiscreteProblem.");

  if (spaces.size() != (unsigned) wf->get_neq()) error("Bad number of spaces in DiscreteProblem.");
  if (spaces.size() > 0) have_spaces = true;
  else error("Zero number of spaces in DiscreteProblem.");

  // Internal variables settings.
  sp_seq = new int[wf->get_neq()];
  memset(sp_seq, -1, sizeof(int) * wf->get_neq());

  // Matrix related settings.
  matrix_buffer = NULL;
  matrix_buffer_dim = 0;
  have_matrix = false;
  values_changed = true;
  struct_changed = true;

  // Initialize precalc shapesets according to spaces provided.
  pss = new PrecalcShapeset*[wf->get_neq()];
  for (unsigned int i = 0; i < wf->get_neq(); i++) pss[i] = NULL;
  num_user_pss = 0;
  for (unsigned int i = 0; i < wf->get_neq(); i++){
    Shapeset *shapeset = spaces[i]->get_shapeset();
    if (shapeset == NULL) error("Internal in DiscreteProblem::init_spaces().");
    PrecalcShapeset *p = new PrecalcShapeset(shapeset);
    if (p == NULL) error("New PrecalcShapeset could not be allocated in DiscreteProblem::init_spaces().");
    pss[i] = p;
    num_user_pss++;
  }

  // Create global enumeration of dof and fill the ndof variable.
  ndof = Space::assign_dofs(spaces);

  // Update the weak formulation with the user-supplied string markers
  // according to the conversion table contained in the mesh.
  element_markers_conversion = &spaces[0]->get_mesh()->element_markers_conversion;
  boundary_markers_conversion = &spaces[0]->get_mesh()->boundary_markers_conversion;
  wf->set_markers_conversion(&spaces[0]->get_mesh()->element_markers_conversion, 
                             &spaces[0]->get_mesh()->boundary_markers_conversion);

  // There is a special function that sets a DiscreteProblem to be FVM.
  // Purpose is that this constructor looks cleaner and is simpler.
  this->is_fvm = false;

  vector_valued_forms = false;

  Geom<Ord> *tmp = init_geom_ord();
  geom_ord = *tmp;
  delete tmp;
}

DiscreteProblem::~DiscreteProblem()
{
  _F_
  free();
  if (sp_seq != NULL) delete [] sp_seq;
  if (pss != NULL) {
    for(int i = 0; i < num_user_pss; i++)
      delete pss[i];
    delete [] pss;
  }
}

void DiscreteProblem::free()
{
  _F_
  struct_changed = values_changed = true;
  if (wf != NULL)
    memset(sp_seq, -1, sizeof(int) * wf->get_neq());
  wf_seq = -1;
}

int DiscreteProblem::get_num_dofs()
{
  _F_
  ndof = 0;
  for (unsigned int i = 0; i < wf->get_neq(); i++)
    ndof += spaces[i]->get_num_dofs();
  return ndof;
}

scalar** DiscreteProblem::get_matrix_buffer(int n)
{
  _F_
  if (n <= matrix_buffer_dim) 
    return matrix_buffer;
  if (matrix_buffer != NULL) 
    delete [] matrix_buffer;
  matrix_buffer_dim = n;
  return (matrix_buffer = new_matrix<scalar>(n, n));
}

//// matrix structure precalculation /////////////////////////////////////////

// This functions is identical in H2D and H3D.
bool DiscreteProblem::is_up_to_date()
{
  _F_
  // check if we can reuse the matrix structure
  bool up_to_date = true;
  if (!have_matrix) up_to_date = false;

  for (unsigned int i = 0; i < wf->get_neq(); i++) {
    if (spaces[i]->get_seq() != sp_seq[i]) {
      up_to_date = false;
      break;
    }
  }

  if (wf->get_seq() != wf_seq)
    up_to_date = false;

  return up_to_date;
}

//// matrix creation /////////////////////////////////////////////////////////
void DiscreteProblem::create_sparse_structure(SparseMatrix* mat, Vector* rhs, 
                      bool force_diagonal_blocks, Table* block_weights)
{
  _F_

  if (is_up_to_date())
  {
    if (mat != NULL)
    {
      verbose("Reusing matrix sparse structure.");
      mat->zero();
    }
    if (rhs != NULL) rhs->zero();
    return;
  }

  // For DG, the sparse structure is different as we have to 
  // account for over-edge calculations.
  bool is_DG = false;
  for(unsigned int i = 0; i < this->wf->mfsurf.size(); i++) {
    if(this->wf->mfsurf[i]->area == H2D_DG_INNER_EDGE) {
      is_DG = true;
      break;
    }
  }
  for(unsigned int i = 0; i < this->wf->vfsurf.size(); i++) {
    if(this->wf->vfsurf[i]->area == H2D_DG_INNER_EDGE) {
      is_DG = true;
      break;
    }
  }

  for(unsigned int i = 0; i < this->wf->mfsurf_mc.size(); i++) {
    if(this->wf->mfsurf_mc[i]->area == H2D_DG_INNER_EDGE) {
      is_DG = true;
      break;
    }
  }
  for(unsigned int i = 0; i < this->wf->vfsurf_mc.size(); i++) {
    if(this->wf->vfsurf_mc[i]->area == H2D_DG_INNER_EDGE) {
      is_DG = true;
      break;
    }
  }

  int ndof = get_num_dofs();

  if (mat != NULL)  
  {
    // Spaces have changed: create the matrix from scratch.
    have_matrix = true;
    mat->free();
    mat->prealloc(ndof);

    AsmList* al = new AsmList[wf->get_neq()];
    Mesh** meshes = new Mesh*[wf->get_neq()];
    bool **blocks = wf->get_blocks(force_diagonal_blocks);

    // Init multi-mesh traversal.
    for (unsigned int i = 0; i < wf->get_neq(); i++) meshes[i] = spaces[i]->get_mesh();

    Traverse trav;
    trav.begin(wf->get_neq(), meshes);

    // Loop through all elements.
    Element **e;
    while ((e = trav.get_next_state(NULL, NULL)) != NULL) {
      // Obtain assembly lists for the element at all spaces.
      for (unsigned int i = 0; i < wf->get_neq(); i++) {
        // TODO: do not get the assembly list again if the element was not changed.
        if (e[i] != NULL) spaces[i]->get_element_assembly_list(e[i], &(al[i]));
      }

      if(is_DG) {
        // Number of edges (= number of vertices).
        int num_edges = e[0]->get_num_surf();

        // Allocation an array of arrays of neighboring elements for every mesh x edge.
        Element **** neighbor_elems_arrays = new Element *** [wf->get_neq()];
        for(unsigned int i = 0; i < wf->get_neq(); i++)
          neighbor_elems_arrays[i] = new Element ** [num_edges];

        // The same, only for number of elements
        int ** neighbor_elems_counts = new int * [wf->get_neq()];
        for(unsigned int i = 0; i < wf->get_neq(); i++)
          neighbor_elems_counts[i] = new int [num_edges];

        // Get the neighbors.
        for(unsigned int el = 0; el < wf->get_neq(); el++) {
          NeighborSearch ns(e[el], meshes[el]);

          // Ignoring errors (and doing nothing) in case the edge is a boundary one.
          ns.set_ignore_errors(true);

          for(int ed = 0; ed < num_edges; ed++) {
            ns.set_active_edge(ed);
            std::vector<Element *> *neighbors = ns.get_neighbors();

            neighbor_elems_counts[el][ed] = ns.get_num_neighbors();
            neighbor_elems_arrays[el][ed] = new Element * [neighbor_elems_counts[el][ed]];
            for(int neigh = 0; neigh < neighbor_elems_counts[el][ed]; neigh++)
              neighbor_elems_arrays[el][ed][neigh] = (*neighbors)[neigh];
          }
        }

        // Pre-add into the stiffness matrix.
        for (unsigned int m = 0; m < wf->get_neq(); m++) {
          for(unsigned int el = 0; el < wf->get_neq(); el++) {

            // Do not include blocks with zero weight except if
            // (force_diagonal_blocks == true && this is a diagonal block).
            bool is_diagonal_block = (m == el);
            if (is_diagonal_block == false || force_diagonal_blocks == false) {
              if (block_weights != NULL) {
                if (fabs(block_weights->get_A(m, el)) < 1e-12) continue;
              }
            }

            for(int ed = 0; ed < num_edges; ed++) {
              for(int neigh = 0; neigh < neighbor_elems_counts[el][ed]; neigh++) {
                if ((blocks[m][el] || blocks[el][m]) && e[m] != NULL)  {
                  AsmList *am = &(al[m]);
                  AsmList *an = new AsmList;
                  spaces[el]->get_element_assembly_list(neighbor_elems_arrays[el][ed][neigh], an);

                  // pretend assembling of the element stiffness matrix
                  // register nonzero elements
                  for (unsigned int i = 0; i < am->cnt; i++) {
                    if (am->dof[i] >= 0) {
                      for (unsigned int j = 0; j < an->cnt; j++) {
                        if (an->dof[j] >= 0) {
                          if(blocks[m][el]) mat->pre_add_ij(am->dof[i], an->dof[j]);
                          if(blocks[el][m]) mat->pre_add_ij(an->dof[j], am->dof[i]);
                        }
                      }
                    }
                  }
                  delete an;
                }
              }
            }
          }
        }

        // Deallocation an array of arrays of neighboring elements 
        // for every mesh x edge.
        for(unsigned int el = 0; el < wf->get_neq(); el++) {
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
      for (unsigned int m = 0; m < wf->get_neq(); m++) {
        for (unsigned int n = 0; n < wf->get_neq(); n++) {

          // Do not include blocks with zero weight except if
          // (force_diagonal_blocks == true && this is a diagonal block).
          bool is_diagonal_block = (m == n);
          if (is_diagonal_block == false || force_diagonal_blocks == false) {
            if (block_weights != NULL) {
              if (fabs(block_weights->get_A(m, n)) < 1e-12) continue;
            }
          }

          if (blocks[m][n] && e[m] != NULL && e[n] != NULL) {
            AsmList *am = &(al[m]);
            AsmList *an = &(al[n]);

            // Pretend assembling of the element stiffness matrix.
            for (unsigned int i = 0; i < am->cnt; i++) {
              if (am->dof[i] >= 0) {
                for (unsigned int j = 0; j < an->cnt; j++) {
                  if (an->dof[j] >= 0) {
                    mat->pre_add_ij(am->dof[i], an->dof[j]);
                  }
                }
              }
            }
          }
        }
      }
    }

    trav.finish();
    delete [] blocks;

    mat->alloc();
  }

  // WARNING: unlike Matrix::alloc(), Vector::alloc(ndof) frees the memory occupied
  // by previous vector before allocating
  if (rhs != NULL) rhs->alloc(ndof);

  // save space seq numbers and weakform seq number, so we can detect their changes
  for (unsigned int i = 0; i < wf->get_neq(); i++) sp_seq[i] = spaces[i]->get_seq();

  wf_seq = wf->get_seq();

  struct_changed = true;
}

//// assembly ////////////////////////////////////////////////////////////////////

// Light version for linear problems.
// The Table is here for optional weighting of matrix blocks in systems.
void DiscreteProblem::assemble(SparseMatrix* mat, Vector* rhs,
                               bool force_diagonal_blocks, Table* block_weights)
{
  _F_
  scalar* coeff_vec = NULL;
  bool add_dir_lift = false;
  assemble(coeff_vec, mat, rhs, force_diagonal_blocks, 
           add_dir_lift, block_weights);
}


void DiscreteProblem::assemble_sanity_checks(Table* block_weights) 
{
  _F_
  // Check that spaces have been set either when constructing or by a call to set_spaces().
  if (!have_spaces) 
    error("You have to call DiscreteProblem::set_spaces() before calling assemble().");
  for (unsigned int i = 0; i < wf->get_neq(); i++)
    if (this->spaces[i] == NULL) error("A space is NULL in assemble().");
  
  // Check that the block scaling table have proper dimension.
  if (block_weights != NULL)
    if (block_weights->get_size() != wf->get_neq())
      error ("Bad dimension of block scaling table in DiscreteProblem::assemble().");
}

void DiscreteProblem::convert_coeff_vec(scalar* coeff_vec, Hermes::vector<Solution *> & u_ext,
                                        bool add_dir_lift) 
{
  _F_
  if (coeff_vec != NULL) {
    for (unsigned int i = 0; i < wf->get_neq(); i++) {
      Solution* external_solution_i = new Solution(spaces[i]->get_mesh());
      Solution::vector_to_solution(coeff_vec, spaces[i], external_solution_i, add_dir_lift);
      u_ext.push_back(external_solution_i);
    }
  }
  else for (unsigned int i = 0; i < wf->get_neq(); i++) u_ext.push_back(NULL);
}

void DiscreteProblem::initialize_psss(Hermes::vector<PrecalcShapeset *>& spss) 
{
  _F_
  for (unsigned int i = 0; i < wf->get_neq(); i++) {
    spss.push_back(new PrecalcShapeset(pss[i]));
    spss[i]->set_quad_2d(&g_quad_2d_std);
  }
}

void DiscreteProblem::initialize_refmaps(Hermes::vector<RefMap *>& refmap) 
{
  _F_
  for (unsigned int i = 0; i < wf->get_neq(); i++) {
    refmap.push_back(new RefMap());
    refmap[i]->set_quad_2d(&g_quad_2d_std);
  }
}

// General assembling function for nonlinear problems. For linear problems use the
// light version above.
// The Table is here for optional weighting of matrix blocks in systems.

void DiscreteProblem::assemble(scalar* coeff_vec, SparseMatrix* mat, 
                               Vector* rhs, bool force_diagonal_blocks, 
                               bool add_dir_lift, Table* block_weights)
{
  _F_
  // Sanity checks.
  assemble_sanity_checks(block_weights);

  // Creating matrix sparse structure.
  create_sparse_structure(mat, rhs, force_diagonal_blocks, block_weights);
 
  // Convert the coefficient vector 'coeff_vec' into solutions Hermes::vector 'u_ext'.
  Hermes::vector<Solution *> u_ext = Hermes::vector<Solution *>();
  convert_coeff_vec(coeff_vec, u_ext, add_dir_lift);
  
  // Reset the warnings about insufficiently high integration order.
  reset_warn_order();

  // Create slave pss's, refmaps.
  Hermes::vector<PrecalcShapeset *> spss;
  Hermes::vector<RefMap *> refmap;

  // Initialize slave pss's, refmaps.
  initialize_psss(spss);
  initialize_refmaps(refmap);

  // Initialize matrix buffer.
  matrix_buffer = NULL;
  matrix_buffer_dim = 0;
  if (mat != NULL) get_matrix_buffer(9);
  
  // Create assembling stages.
  std::vector<WeakForm::Stage> stages = std::vector<WeakForm::Stage>();
  bool want_matrix = (mat != NULL);
  bool want_vector = (rhs != NULL);
  wf->get_stages(spaces, u_ext, stages, want_matrix, want_vector);

  // Loop through all assembling stages -- the purpose of this is increased performance
  // in multi-mesh calculations, where, e.g., only the right hand side uses two meshes.
  // In such a case, the matrix forms are assembled over one mesh, and only the rhs
  // traverses through the union mesh. On the other hand, if you don't use multi-mesh
  // at all, there will always be only one stage in which all forms are assembled as usual.
  for (unsigned ss = 0; ss < stages.size(); ss++) {
    // Assemble one stage. One stage is a collection of functions, 
    // and meshes that can not be further minimized.
    // E.g. if a linear form uses two external solutions, each of 
    // which is defined on a different mesh, and different to the
    // mesh of the current test function, then the stage would have 
    // three meshes. By stage functions, all functions are meant: shape 
    // functions (their precalculated values), and mesh functions.
    assemble_one_stage(stages[ss], mat, rhs, force_diagonal_blocks, 
                       block_weights, spss, refmap, u_ext);
  }

  // Deinitialize matrix buffer.
  if(matrix_buffer != NULL)
    delete [] matrix_buffer;
  matrix_buffer = NULL;
  matrix_buffer_dim = 0;

  // Deinitialize slave pss's, refmaps.
  for(std::vector<PrecalcShapeset *>::iterator it = spss.begin(); it != spss.end(); it++)
    delete *it;
  for(std::vector<RefMap *>::iterator it = refmap.begin(); it != refmap.end(); it++)
    delete *it;

  // Delete the vector u_ext.
  for(std::vector<Solution *>::iterator it = u_ext.begin(); it != u_ext.end(); it++)
    delete *it;
}

void DiscreteProblem::assemble_one_stage(WeakForm::Stage& stage, 
					 SparseMatrix* matrix, Vector* rhs,
                                         bool force_diagonal_blocks, Table* block_weights,
                                         Hermes::vector<PrecalcShapeset *>& spss, 
                                         Hermes::vector<RefMap *>& refmap, 
                                         Hermes::vector<Solution *>& u_ext)
{
  _F_
  // Boundary flags. bnd[i] == true if i-th edge of the current 
  // element is a boundary edge.
  bool bnd[4];

  // Info about the boundary edge.
  SurfPos surf_pos[4];

  // Create the assembling states.
  Traverse trav;
  for (unsigned i = 0; i < stage.idx.size(); i++)
    stage.fns[i] = pss[stage.idx[i]];
  for (unsigned i = 0; i < stage.ext.size(); i++)
    stage.ext[i]->set_quad_2d(&g_quad_2d_std);
  trav.begin(stage.meshes.size(), &(stage.meshes.front()), &(stage.fns.front()));

  // Check that there is a DG form, so that the DG assembling procedure needs to be performed.
  DG_matrix_forms_present = false;
  DG_vector_forms_present = false;
  for(unsigned int i = 0; i < stage.mfsurf.size(); i++) {
    if (stage.mfsurf[i]->area == H2D_DG_INNER_EDGE) {
      DG_matrix_forms_present = true;
      break;
    }
  }
  for(unsigned int i = 0; i < stage.vfsurf.size(); i++) {
    if (stage.vfsurf[i]->area == H2D_DG_INNER_EDGE) {
      DG_vector_forms_present = true;
      break;
    }
  }

  for(unsigned int i = 0; i < stage.mfsurf_mc.size(); i++) {
    if (stage.mfsurf_mc[i]->area == H2D_DG_INNER_EDGE) {
      DG_matrix_forms_present = true;
      break;
    }
  }
  for(unsigned int i = 0; i < stage.vfsurf_mc.size(); i++) {
    if (stage.vfsurf_mc[i]->area == H2D_DG_INNER_EDGE) {
      DG_vector_forms_present = true;
      break;
    }
  }

  // Loop through all assembling states.
  // Assemble each one.
  Element** e;
  while ((e = trav.get_next_state(bnd, surf_pos)) != NULL) {
    // One state is a collection of (virtual) elements sharing 
    // the same physical location on (possibly) different meshes.
    // This is then the same element of the virtual union mesh. 
    // The proper sub-element mappings to all the functions of
    // this stage is supplied by the function Traverse::get_next_state() 
    // called in the while loop.
    assemble_one_state(stage, matrix, rhs, force_diagonal_blocks, 
                       block_weights, spss, refmap, 
                       u_ext, e, bnd, surf_pos, trav.get_base());
  }

  if (matrix != NULL) matrix->finish();
  if (rhs != NULL) rhs->finish();
  trav.finish();

  if(DG_matrix_forms_present || DG_vector_forms_present) {
    Element* element_to_set_nonvisited;
    for(unsigned int mesh_i = 0; mesh_i < stage.meshes.size(); mesh_i++)
      for_all_elements(element_to_set_nonvisited, stage.meshes[mesh_i])
        element_to_set_nonvisited->visited = false;
  }
}

Element* DiscreteProblem::init_state(WeakForm::Stage& stage, Hermes::vector<PrecalcShapeset *>& spss, 
  Hermes::vector<RefMap *>& refmap, Element** e, Hermes::vector<bool>& isempty, Hermes::vector<AsmList *>& al)
{
  _F_
  // Find a non-NULL e[i].
  Element* e0 = NULL;
  for (unsigned int i = 0; i < stage.idx.size(); i++)
    if ((e0 = e[i]) != NULL) 
      break;
  if(e0 == NULL) 
    return NULL;

  // Set maximum integration order for use in integrals, see limit_order()
  update_limit_table(e0->get_mode());

  // Obtain assembly lists for the element at all spaces of the stage, set appropriate mode for each pss.
  // NOTE: Active elements and transformations for external functions (including the solutions from previous
  // Newton's iteration) as well as basis functions (master PrecalcShapesets) have already been set in
  // trav.get_next_state(...).
  for (unsigned int i = 0; i < stage.idx.size(); i++) {
    int j = stage.idx[i];
    if (e[i] == NULL) {
      isempty[j] = true;
      continue;
    }

    // TODO: do not obtain again if the element was not changed.
    spaces[j]->get_element_assembly_list(e[i], al[j]);

    // Set active element to all test functions.
    spss[j]->set_active_element(e[i]);
    spss[j]->set_master_transform();

    // Set active element to reference mappings.
    refmap[j]->set_active_element(e[i]);
    refmap[j]->force_transform(pss[j]->get_transform(), pss[j]->get_ctm());

    // Mark the active element on each mesh in order to prevent assembling on its edges from the other side.
    if(DG_matrix_forms_present || DG_vector_forms_present)
      e[i]->visited = true;
  }
  return e0;
}

void DiscreteProblem::assemble_one_state(WeakForm::Stage& stage,
      SparseMatrix* matrix, Vector* rhs,
      bool force_diagonal_blocks, Table* block_weights,
      Hermes::vector<PrecalcShapeset *>& spss, Hermes::vector<RefMap *>& refmap,
      Hermes::vector<Solution *>& u_ext, Element** e,
      bool* bnd, SurfPos* surf_pos, Element* trav_base)
{
  _F_
  // Assembly list vector.
  Hermes::vector<AsmList *> al;
  for(unsigned int i = 0; i < wf->get_neq(); i++) 
    al.push_back(new AsmList);
  
  // Natural boundary condition flag.
  Hermes::vector<bool> nat;
  for(unsigned int i = 0; i < wf->get_neq(); i++) 
    nat.push_back(false);
  
  // Element usage flag: iempty[i] == true if the current state does not posses an active element in the i-th space.
  Hermes::vector<bool> isempty;
  for(unsigned int i = 0; i < wf->get_neq(); i++) 
    isempty.push_back(false);
    
  // Initialize the state, return a non-NULL element; if no such element found, return.
  Element* rep_element = init_state(stage, spss, refmap, e, isempty, al);
  if(rep_element == NULL)
    return;

  init_cache();

  /// Assemble volume matrix forms.
   assemble_volume_matrix_forms(stage, matrix, rhs, force_diagonal_blocks,
                               block_weights, spss, refmap, u_ext, isempty,
                               rep_element->marker, al);
  if(!stage.mfvol_mc.empty())
    assemble_multicomponent_volume_matrix_forms(stage, matrix, rhs, force_diagonal_blocks,
                               block_weights, spss, refmap, u_ext, isempty,
                               rep_element->marker, al);

  /// Assemble volume vector forms.
  if (rhs != NULL) {
    assemble_volume_vector_forms(stage, matrix, rhs, force_diagonal_blocks,
                                 block_weights, spss, refmap, u_ext, isempty,
                                 rep_element->marker, al);
    if(!stage.vfvol_mc.empty()) {
      assemble_multicomponent_volume_vector_forms(stage, matrix, rhs, force_diagonal_blocks,
                                 block_weights, spss, refmap, u_ext, isempty,
                                 rep_element->marker, al);
    }
  }

  // Assemble surface integrals now: loop through surfaces of the element.
  for (int isurf = 0; isurf < e[0]->get_num_surf(); isurf++) {
    assemble_surface_integrals(stage, matrix, rhs, force_diagonal_blocks,
                              block_weights, spss, refmap, u_ext, isempty,
                              surf_pos[isurf].marker, al, bnd[isurf], surf_pos[isurf],
                              nat, isurf, e, trav_base, rep_element);
  }

  // Delete assembly lists.
  for(unsigned int i = 0; i < wf->get_neq(); i++) 
    delete al[i];

  delete_cache();
}

void DiscreteProblem::assemble_volume_matrix_forms(WeakForm::Stage& stage, 
                      SparseMatrix* mat, Vector* rhs, bool force_diagonal_blocks, 
                      Table* block_weights,
                      Hermes::vector<PrecalcShapeset *>& spss, Hermes::vector<RefMap *>& refmap, 
                      Hermes::vector<Solution *>& u_ext, Hermes::vector<bool>& isempty, 
                      int marker, Hermes::vector<AsmList *>& al)
{
  _F_
  for (unsigned ww = 0; ww < stage.mfvol.size(); ww++) {
    WeakForm::MatrixFormVol* mfv = stage.mfvol[ww];
    int m = mfv->i;
    int n = mfv->j;
    if (isempty[m] || isempty[n])
      continue;
    if (fabs(mfv->scaling_factor) < 1e-12)
      continue;
    if (mfv->area != HERMES_ANY && !(marker == element_markers_conversion->get_internal_marker(mfv->area)))
      continue;

    // If a block scaling table is provided, and if the scaling coefficient
    // A_mn for this block is zero, then the form does not need to be assembled.
    double block_scaling_coeff = 1.;
    if (block_weights != NULL) {
      block_scaling_coeff = block_weights->get_A(m, n);
      if (fabs(block_scaling_coeff) < 1e-12)
        continue;
    }
    bool tra = (m != n) && (mfv->sym != 0);
    bool sym = (m == n) && (mfv->sym == 1);

    // Assemble the local stiffness matrix for the form mfv.
    scalar **local_stiffness_matrix = NULL;
    local_stiffness_matrix = get_matrix_buffer(std::max(al[m]->cnt, al[n]->cnt));

    for (unsigned int i = 0; i < al[m]->cnt; i++) {
      if (!tra && al[m]->dof[i] < 0) 
        continue;
      spss[m]->set_active_shape(al[m]->idx[i]);
        
      // Unsymmetric block.
      if (!sym) {
        for (unsigned int j = 0; j < al[n]->cnt; j++) {
          pss[n]->set_active_shape(al[n]->idx[j]);
          
          if (al[n]->dof[j] < 0) {
          // Linear problems only: Subtracting Dirichlet lift contribution from the RHS:
          if (rhs != NULL && this->is_linear) {
              // Numerical integration performed only if all 
              // coefficients multiplying the form are nonzero
              // and if the basis function is active.
              if (std::abs(al[m]->coef[i]) > 1e-12 && std::abs(al[n]->coef[j]) > 1e-12
                  && al[m]->dof[i] >= 0) {
                scalar val = eval_form(mfv, u_ext, pss[n], spss[m], refmap[n],
                                       refmap[m]) * al[n]->coef[j] * al[m]->coef[i];
                rhs->add(al[m]->dof[i], -val);
              }
            }
          }
          else if (mat != NULL) {
            scalar val = 0;
            // Numerical integration performed only if all 
            // coefficients multiplying the form are nonzero.
            if (std::abs(al[m]->coef[i]) > 1e-12 && std::abs(al[n]->coef[j]) > 1e-12) {
              val = block_scaling_coeff * eval_form(mfv, u_ext, pss[n], spss[m], refmap[n],
                                                    refmap[m]) * al[n]->coef[j] * al[m]->coef[i];
            }
            local_stiffness_matrix[i][j] = val;
          }
        }
      }
      // Symmetric block.
      else {
        for (unsigned int j = 0; j < al[n]->cnt; j++) {
          if (j < i && al[n]->dof[j] >= 0) 
            continue;
            
          pss[n]->set_active_shape(al[n]->idx[j]);
          
          if (al[n]->dof[j] < 0) {
            // Linear problems only: Subtracting Dirichlet lift contribution from the RHS:
            if (rhs != NULL && this->is_linear) {
              // Numerical integration performed only if all 
              // coefficients multiplying the form are nonzero
              // and if basis function is active.
              if (std::abs(al[m]->coef[i]) > 1e-12 && std::abs(al[n]->coef[j]) > 1e-12
                  && al[m]->dof[i] >= 0) {

                scalar val = eval_form(mfv, u_ext, pss[n], spss[m], refmap[n], refmap[m]) * 
                                       al[n]->coef[j] * al[m]->coef[i];
                rhs->add(al[m]->dof[i], -val);
	      }
            }
          }
          else if (mat != NULL) {
            scalar val = 0;
            // Numerical integration performed only if all coefficients 
            // multiplying the form are nonzero.
            if (std::abs(al[m]->coef[i]) > 1e-12 && std::abs(al[n]->coef[j]) > 1e-12) {
              val = block_scaling_coeff * eval_form(mfv, u_ext, pss[n], spss[m], refmap[n],
                                                    refmap[m]) * al[n]->coef[j] * al[m]->coef[i];
            }
            local_stiffness_matrix[i][j] = local_stiffness_matrix[j][i] = val;
          }
        }
      }
    }

    // Insert the local stiffness matrix into the global one.
    if (mat != NULL) {
      mat->add(al[m]->cnt, al[n]->cnt, local_stiffness_matrix, al[m]->dof, al[n]->dof);
    }

    // Insert also the off-diagonal (anti-)symmetric block, if required.
    if (tra) {
      if (mfv->sym < 0)
        chsgn(local_stiffness_matrix, al[m]->cnt, al[n]->cnt);

      transpose(local_stiffness_matrix, al[m]->cnt, al[n]->cnt);

      if (mat != NULL) {
        mat->add(al[n]->cnt, al[m]->cnt, local_stiffness_matrix, al[n]->dof, al[m]->dof);
      }

      // Linear problems only: Subtracting Dirichlet lift contribution from the RHS:
      if (rhs != NULL && this->is_linear)
        for (unsigned int j = 0; j < al[m]->cnt; j++)
          if (al[m]->dof[j] < 0)
            for (unsigned int i = 0; i < al[n]->cnt; i++)
              if (al[n]->dof[i] >= 0)
                rhs->add(al[n]->dof[i], -local_stiffness_matrix[i][j]);
    }
  }
}

void DiscreteProblem::assemble_multicomponent_volume_matrix_forms(WeakForm::Stage& stage, 
                      SparseMatrix* mat, Vector* rhs, bool force_diagonal_blocks, 
                      Table* block_weights,
                      Hermes::vector<PrecalcShapeset *>& spss, Hermes::vector<RefMap *>& refmap, 
                      Hermes::vector<Solution *>& u_ext, Hermes::vector<bool>& isempty, 
                      int marker, Hermes::vector<AsmList *>& al)
{
  _F_
  for (unsigned ww = 0; ww < stage.mfvol_mc.size(); ww++) {
    WeakForm::MultiComponentMatrixFormVol* mfv = stage.mfvol_mc[ww];
    if(fabs(mfv->scaling_factor) < 1e-12)
      continue;
    if (mfv->area != HERMES_ANY && !(marker == element_markers_conversion->get_internal_marker(mfv->area)))
      continue;

    // If a block scaling table is provided, and if the scaling coefficient
    // A_mn for this block is zero, then the form does not need to be assembled.
    Hermes::vector<double> block_scaling_coeffs;
    for(unsigned int coordinate_i = 0; coordinate_i < mfv->coordinates.size(); coordinate_i++)
      if (block_weights != NULL)
        block_scaling_coeffs.push_back(block_weights->get_A(mfv->coordinates[coordinate_i].first, mfv->coordinates[coordinate_i].second));
      else
        block_scaling_coeffs.push_back(1);

    // Assemble the local stiffness matrix for the form mfv.
    unsigned int m = mfv->coordinates[0].first;
    unsigned int n = mfv->coordinates[0].second;

    if(mfv->sym)
      for(unsigned int coordinate_i = 0; coordinate_i < mfv->coordinates.size(); coordinate_i++)
        if(mfv->coordinates[coordinate_i].first != mfv->coordinates[coordinate_i].second)
          error("Symmetrical multicomponent forms need to take both basis function from one space.");

    for (unsigned int i = 0; i < al[m]->cnt; i++) {
      spss[m]->set_active_shape(al[m]->idx[i]);
      if(al[m]->dof[i] < 0 && mfv->sym == HERMES_NONSYM)
        continue;
      if(mfv->sym == HERMES_NONSYM) {
        for (unsigned int j = 0; j < al[n]->cnt; j++) {
          pss[n]->set_active_shape(al[n]->idx[j]);
          if (al[n]->dof[j] < 0) {
            // Linear problems only: Subtracting Dirichlet lift contribution from the RHS:
            if (rhs != NULL && this->is_linear)
              if (std::abs(al[m]->coef[i]) > 1e-12 && std::abs(al[n]->coef[j]) > 1e-12 && al[m]->dof[i] >= 0) {
                Hermes::vector<scalar> result;
                eval_form(mfv, u_ext, pss[n], spss[m], refmap[n], refmap[m], result);
                for(unsigned int coordinate_i = 0; coordinate_i < mfv->coordinates.size(); coordinate_i++)
                  rhs->add(al[mfv->coordinates[coordinate_i].first]->dof[i], -result[coordinate_i] * al[n]->coef[j] * al[m]->coef[i]);
              }
          }
          else if (mat != NULL) {
            if (std::abs(al[m]->coef[i]) > 1e-12 && std::abs(al[n]->coef[j]) > 1e-12) {
              Hermes::vector<scalar> result;
              eval_form(mfv, u_ext, pss[n], spss[m], refmap[n], refmap[m], result);
              for(unsigned int coordinate_i = 0; coordinate_i < mfv->coordinates.size(); coordinate_i++) {
                mat->add(al[mfv->coordinates[coordinate_i].first]->dof[i], al[mfv->coordinates[coordinate_i].second]->dof[j],
                result[coordinate_i] * block_scaling_coeffs[coordinate_i] * al[n]->coef[j] * al[m]->coef[i]);
              }
            }
          }
        }
      }
      else {
        for (unsigned int j = 0; j < al[n]->cnt; j++) {
          if(j < i && al[n]->dof[j] >= 0)
            continue;
          pss[n]->set_active_shape(al[n]->idx[j]);

          if (al[n]->dof[j] < 0) {
            // Linear problems only: Subtracting Dirichlet lift contribution from the RHS:
            if (rhs != NULL && this->is_linear)
              if (std::abs(al[m]->coef[i]) > 1e-12 && std::abs(al[n]->coef[j]) > 1e-12 && al[m]->dof[i] >= 0) {
                Hermes::vector<scalar> result;
                eval_form(mfv, u_ext, pss[n], spss[m], refmap[n], refmap[m], result);
                for(unsigned int coordinate_i = 0; coordinate_i < mfv->coordinates.size(); coordinate_i++)
                  rhs->add(al[mfv->coordinates[coordinate_i].first]->dof[i], -result[coordinate_i] * al[n]->coef[j] * al[m]->coef[i]);
              }
          }
          else {
            if (al[m]->dof[i] < 0) {
              if (rhs != NULL && this->is_linear)
              if (std::abs(al[m]->coef[i]) > 1e-12 && std::abs(al[n]->coef[j]) > 1e-12) {
                Hermes::vector<scalar> result;
                eval_form(mfv, u_ext, pss[n], spss[m], refmap[n], refmap[m], result);
                for(unsigned int coordinate_i = 0; coordinate_i < mfv->coordinates.size(); coordinate_i++)
                  if(mfv->sym == HERMES_SYM)
                    rhs->add(al[mfv->coordinates[coordinate_i].second]->dof[j], - result[coordinate_i] * al[n]->coef[j] * al[m]->coef[i]);
                  else
                    rhs->add(al[mfv->coordinates[coordinate_i].second]->dof[j], result[coordinate_i] * al[n]->coef[j] * al[m]->coef[i]);
              }
            }
            else if (mat != NULL) {
              if (std::abs(al[m]->coef[i]) > 1e-12 && std::abs(al[n]->coef[j]) > 1e-12) {
                Hermes::vector<scalar> result;
                eval_form(mfv, u_ext, pss[n], spss[m], refmap[n], refmap[m], result);
                for(unsigned int coordinate_i = 0; coordinate_i < mfv->coordinates.size(); coordinate_i++) {
                  mat->add(al[mfv->coordinates[coordinate_i].first]->dof[i], al[mfv->coordinates[coordinate_i].second]->dof[j],
                  result[coordinate_i] * block_scaling_coeffs[coordinate_i] * al[n]->coef[j] * al[m]->coef[i]);
                  if(i != j)
                    mat->add(al[mfv->coordinates[coordinate_i].first]->dof[j], al[mfv->coordinates[coordinate_i].second]->dof[i],
                    result[coordinate_i] * block_scaling_coeffs[coordinate_i] * al[n]->coef[j] * al[m]->coef[i]);
                }
              }
            }
          }
        }
      }
    }
  }
}

void DiscreteProblem::assemble_volume_vector_forms(WeakForm::Stage& stage, 
      SparseMatrix* mat, Vector* rhs, bool force_diagonal_blocks, 
      Table* block_weights, Hermes::vector<PrecalcShapeset *>& spss, 
      Hermes::vector<RefMap *>& refmap, Hermes::vector<Solution *>& u_ext, 
      Hermes::vector<bool>& isempty, int marker, Hermes::vector<AsmList *>& al)
{
  _F_

  if (rhs == NULL) return;

  for (unsigned int ww = 0; ww < stage.vfvol.size(); ww++) {
    WeakForm::VectorFormVol* vfv = stage.vfvol[ww];
    int m = vfv->i;
    if (isempty[vfv->i]) 
      continue;
    if (fabs(vfv->scaling_factor) < 1e-12) 
      continue;
    if (vfv->area != HERMES_ANY && !(marker == element_markers_conversion->get_internal_marker(vfv->area))) 
      continue;

    for (unsigned int i = 0; i < al[m]->cnt; i++) {
      if (al[m]->dof[i] < 0) continue;
      
      spss[m]->set_active_shape(al[m]->idx[i]);

      // Numerical integration performed only if the coefficient 
      // multiplying the form is nonzero.
      if (std::abs(al[m]->coef[i]) > 1e-12) {   
        rhs->add(al[m]->dof[i], eval_form(vfv, u_ext, spss[m], refmap[m]) * al[m]->coef[i]);
      }
    }
  }
}
void DiscreteProblem::assemble_multicomponent_volume_vector_forms(WeakForm::Stage& stage, 
      SparseMatrix* mat, Vector* rhs, bool force_diagonal_blocks, 
      Table* block_weights, Hermes::vector<PrecalcShapeset *>& spss, 
      Hermes::vector<RefMap *>& refmap, Hermes::vector<Solution *>& u_ext, 
      Hermes::vector<bool>& isempty, int marker, Hermes::vector<AsmList *>& al)
{
  _F_

  if (rhs == NULL) return;

  for (unsigned int ww = 0; ww < stage.vfvol_mc.size(); ww++) {
    WeakForm::MultiComponentVectorFormVol* vfv = stage.vfvol_mc[ww];
    if (fabs(vfv->scaling_factor) < 1e-12) 
      continue;
    if (vfv->area != HERMES_ANY && !(marker == element_markers_conversion->get_internal_marker(vfv->area))) 
      continue;

    unsigned int m = vfv->coordinates[0];

    for (unsigned int i = 0; i < al[m]->cnt; i++) {
      if (al[m]->dof[i] < 0) continue;
      
      spss[m]->set_active_shape(al[m]->idx[i]);

      // Numerical integration performed only if the coefficient 
      // multiplying the form is nonzero.
      if (std::abs(al[m]->coef[i]) > 1e-12) {   
        Hermes::vector<scalar> result;
        eval_form(vfv, u_ext, spss[m], refmap[m], result);
        for(unsigned int coordinate_i = 0; coordinate_i < vfv->coordinates.size(); coordinate_i++)
          rhs->add(al[vfv->coordinates[coordinate_i]]->dof[i], result[coordinate_i] * al[vfv->coordinates[coordinate_i]]->coef[i]);
      }
    }
  }
}

void DiscreteProblem::assemble_surface_integrals(WeakForm::Stage& stage, 
       SparseMatrix* mat, Vector* rhs, bool force_diagonal_blocks, 
       Table* block_weights, Hermes::vector<PrecalcShapeset *>& spss, 
       Hermes::vector<RefMap *>& refmap, Hermes::vector<Solution *>& u_ext, 
       Hermes::vector<bool>& isempty, int marker, Hermes::vector<AsmList *>& al, 
       bool bnd, SurfPos& surf_pos, Hermes::vector<bool>& nat, int isurf, 
       Element** e, Element* trav_base, Element* rep_element)
{
  _F_
  // Obtain the list of shape functions which are nonzero on this surface.
  for (unsigned int i = 0; i < stage.idx.size(); i++) {
    int j = stage.idx[i];
    if (isempty[j])
      continue;
    if(marker > 0) {
      nat[j] = true;
      if(spaces[j]->get_essential_bcs() != NULL)
        if(spaces[j]->get_essential_bcs()->get_boundary_condition(boundary_markers_conversion->get_user_marker(marker)) != NULL)
          nat[j] = false;
    }
    spaces[j]->get_boundary_assembly_list(e[i], isurf, al[j]);
  }

  // Assemble boundary edges: 
  if(bnd == 1) {
    assemble_surface_matrix_forms(stage, mat, rhs, force_diagonal_blocks, block_weights, 
                                  spss, refmap, u_ext, isempty, 
                                  marker, al, bnd, surf_pos, nat, isurf, e, trav_base);
    if(!stage.mfsurf_mc.empty())
      assemble_multicomponent_surface_matrix_forms(stage, mat, rhs, force_diagonal_blocks, block_weights, 
                                  spss, refmap, u_ext, isempty, 
                                  marker, al, bnd, surf_pos, nat, isurf, e, trav_base);
    if (rhs != NULL) {
      assemble_surface_vector_forms(stage, mat, rhs, force_diagonal_blocks, block_weights, 
                                    spss, refmap, u_ext, isempty, 
                                    marker, al, bnd, surf_pos, nat, isurf, e, trav_base);
      if(!stage.vfsurf_mc.empty())
        assemble_multicomponent_surface_vector_forms(stage, mat, rhs, force_diagonal_blocks, block_weights, 
                                    spss, refmap, u_ext, isempty, 
                                    marker, al, bnd, surf_pos, nat, isurf, e, trav_base);
    }
  }
  // Assemble inner edges (in discontinuous Galerkin discretization): 
  else
    if(DG_vector_forms_present || DG_matrix_forms_present)
      assemble_DG_forms(stage, mat, rhs, force_diagonal_blocks, block_weights, 
                        spss, refmap, u_ext, isempty, 
                        marker, al, bnd, surf_pos, nat, isurf, e, trav_base, rep_element);
  }

void DiscreteProblem::assemble_DG_forms(WeakForm::Stage& stage, 
       SparseMatrix* mat, Vector* rhs, bool force_diagonal_blocks, 
       Table* block_weights, Hermes::vector<PrecalcShapeset *>& spss, 
       Hermes::vector<RefMap *>& refmap, Hermes::vector<Solution *>& u_ext, 
       Hermes::vector<bool>& isempty, int marker, Hermes::vector<AsmList *>& al, 
       bool bnd, SurfPos& surf_pos, Hermes::vector<bool>& nat, 
       int isurf, Element** e, Element* trav_base, Element* rep_element)
{
  _F_
  // Determine the minimum mesh seq in this stage.
  min_dg_mesh_seq = 0;
  for(unsigned int i = 0; i < stage.meshes.size(); i++)
    if(stage.meshes[i]->get_seq() < min_dg_mesh_seq || i == 0)
      min_dg_mesh_seq = stage.meshes[i]->get_seq();
  
  // Initialize the NeighborSearches.
  // 5 is for bits per page in the array.
  LightArray<NeighborSearch*> neighbor_searches(5);
  init_neighbors(neighbor_searches, stage, isurf);

  // Create a multimesh tree;
  DiscreteProblem::NeighborNode* root = new DiscreteProblem::NeighborNode(NULL, 0);
  build_multimesh_tree(root, neighbor_searches);
    
  // Update all NeighborSearches according to the multimesh tree.
  // After this, all NeighborSearches in neighbor_searches should have the same count 
  // of neighbors and proper set of transformations
  // for the central and the neighbor element(s) alike.
  // Also check that every NeighborSearch has the same number of neighbor elements.
  unsigned int num_neighbors = 0;
  for(unsigned int i = 0; i < neighbor_searches.get_size(); i++)
    if(neighbor_searches.present(i)) {
      NeighborSearch* ns = neighbor_searches.get(i);
      update_neighbor_search(ns, root);
      if(num_neighbors == 0)
        num_neighbors = ns->n_neighbors;
      if(ns->n_neighbors != num_neighbors)
        error("Num_neighbors of different NeighborSearches not matching in DiscreteProblem::assemble_surface_integrals().");
    }

  // Create neighbor psss, refmaps.
  Hermes::vector<PrecalcShapeset *> npss;
  Hermes::vector<PrecalcShapeset *> nspss;
  Hermes::vector<RefMap *> nrefmap;

  // Initialize neighbor precalc shapesets and refmaps.      
  // This is only needed when there are matrix DG forms present.
  if(DG_matrix_forms_present)
    for (unsigned int i = 0; i < stage.idx.size(); i++) {
      npss.push_back(new PrecalcShapeset(pss[i]->get_shapeset()));
      npss[i]->set_quad_2d(&g_quad_2d_std);
      nspss.push_back(new PrecalcShapeset(npss[i]));
      nspss[i]->set_quad_2d(&g_quad_2d_std);
      nrefmap.push_back(new RefMap());
      nrefmap[i]->set_quad_2d(&g_quad_2d_std);
    }

  for(unsigned int neighbor_i = 0; neighbor_i < num_neighbors; neighbor_i++) {
    // If the active segment has already been processed (when the neighbor element was assembled), it is skipped.
    // We test all neighbor searches, because in the case of intra-element edge, the neighboring (the same as central) element
    // will be marked as visited, even though the edge was not calculated.
    bool processed = true;
    for(unsigned int i = 0; i < neighbor_searches.get_size(); i++)
      if(neighbor_searches.present(i))
        if(!neighbor_searches.get(i)->neighbors.at(neighbor_i)->visited) {
            processed = false;
            break;
        }

    if(!DG_vector_forms_present && processed)
      continue;

    // For every neighbor we want to delete the geometry caches and create new ones.
    for (int i = 0; i < g_max_quad + 1 + 4 * g_max_quad + 4; i++)
      if (cache_e[i] != NULL) {
        cache_e[i]->free(); 
        delete cache_e[i];
        cache_e[i] = NULL;
        delete [] cache_jwt[i];
      }

    assemble_DG_one_neighbor(processed, neighbor_i, stage, mat, rhs, 
                             force_diagonal_blocks, block_weights, spss, refmap, 
                             npss, nspss, nrefmap, neighbor_searches, u_ext, isempty, 
                             marker, al, bnd, surf_pos, nat, isurf, e, trav_base, rep_element);
  }

  // Delete the multimesh tree;
  delete root;

  // Deinitialize neighbor pss's, refmaps.
  if(DG_matrix_forms_present) {
    for(std::vector<PrecalcShapeset *>::iterator it = nspss.begin(); it != nspss.end(); it++)
      delete *it;
    for(std::vector<PrecalcShapeset *>::iterator it = npss.begin(); it != npss.end(); it++)
      delete *it;
    for(std::vector<RefMap *>::iterator it = nrefmap.begin(); it != nrefmap.end(); it++)
      delete *it;
  }

  // Delete the neighbor_searches array.
  for(unsigned int i = 0; i < neighbor_searches.get_size(); i++) 
    if(neighbor_searches.present(i))
      delete neighbor_searches.get(i);
}


void DiscreteProblem::assemble_DG_one_neighbor(bool edge_processed, unsigned int neighbor_i, WeakForm::Stage& stage, 
      SparseMatrix* mat, Vector* rhs, bool force_diagonal_blocks, Table* block_weights,
       Hermes::vector<PrecalcShapeset *>& spss, Hermes::vector<RefMap *>& refmap, Hermes::vector<PrecalcShapeset *>& npss, 
       Hermes::vector<PrecalcShapeset *>& nspss, Hermes::vector<RefMap *>& nrefmap, LightArray<NeighborSearch*>& neighbor_searches, Hermes::vector<Solution *>& u_ext, 
       Hermes::vector<bool>& isempty, int marker, Hermes::vector<AsmList *>& al, bool bnd, SurfPos& surf_pos, Hermes::vector<bool>& nat, 
       int isurf, Element** e, Element* trav_base, Element* rep_element)
{
  _F_
  // Set the active segment in all NeighborSearches
  for(unsigned int i = 0; i < neighbor_searches.get_size(); i++)
    if(neighbor_searches.present(i)) {
      neighbor_searches.get(i)->active_segment = neighbor_i;
      neighbor_searches.get(i)->neighb_el = neighbor_searches.get(i)->neighbors[neighbor_i];
      neighbor_searches.get(i)->neighbor_edge = neighbor_searches.get(i)->neighbor_edges[neighbor_i];
    }

  // Push all the necessary transformations to all functions of this stage.
  // The important thing is that the transformations to the current subelement are already there.
  for(unsigned int fns_i = 0; fns_i < stage.fns.size(); fns_i++)
    for(unsigned int trf_i = 0; trf_i < neighbor_searches.get(stage.meshes[fns_i]->get_seq() - min_dg_mesh_seq)->central_n_trans[neighbor_i]; trf_i++)
      stage.fns[fns_i]->push_transform(neighbor_searches.get(stage.meshes[fns_i]->get_seq() - min_dg_mesh_seq)->central_transformations[neighbor_i][trf_i]);
  
  // For neighbor psss.
  if(DG_matrix_forms_present && !edge_processed)
    for(unsigned int idx_i = 0; idx_i < stage.idx.size(); idx_i++) {
      npss[idx_i]->set_active_element((*neighbor_searches.get(stage.meshes[idx_i]->get_seq() - min_dg_mesh_seq)->get_neighbors())[neighbor_i]);
      for(unsigned int trf_i = 0; trf_i < neighbor_searches.get(stage.meshes[idx_i]->get_seq() - min_dg_mesh_seq)->neighbor_n_trans[neighbor_i]; trf_i++)
        npss[idx_i]->push_transform(neighbor_searches.get(stage.meshes[idx_i]->get_seq() - min_dg_mesh_seq)->neighbor_transformations[neighbor_i][trf_i]);
    }

  // Also push the transformations to the slave psss and refmaps.
  for (unsigned int i = 0; i < stage.idx.size(); i++) {
    if(isempty[stage.idx[i]])
      continue;
    spss[stage.idx[i]]->set_master_transform();
    refmap[stage.idx[i]]->force_transform(pss[stage.idx[i]]->get_transform(), pss[stage.idx[i]]->get_ctm());

    // Neighbor.
    if(DG_matrix_forms_present && !edge_processed) {
      nspss[stage.idx[i]]->set_active_element(npss[stage.idx[i]]->get_active_element());
      nspss[stage.idx[i]]->set_master_transform();
      nrefmap[stage.idx[i]]->set_active_element(npss[stage.idx[i]]->get_active_element());
      nrefmap[stage.idx[i]]->force_transform(npss[stage.idx[i]]->get_transform(), npss[stage.idx[i]]->get_ctm());
    }
  }

  /***/

  // The computation takes place here.
  if(DG_matrix_forms_present && !edge_processed) {
    assemble_DG_matrix_forms(stage, mat, rhs, force_diagonal_blocks, block_weights, spss, refmap, npss, nspss, nrefmap, neighbor_searches, u_ext, isempty, 
      marker, al, bnd, surf_pos, nat, isurf, e, trav_base, rep_element);
    if(!stage.mfsurf_mc.empty())
      assemble_multicomponent_DG_matrix_forms(stage, mat, rhs, force_diagonal_blocks, block_weights, spss, refmap, npss, nspss, nrefmap, neighbor_searches, u_ext, isempty, 
      marker, al, bnd, surf_pos, nat, isurf, e, trav_base, rep_element);
  }
  if (DG_vector_forms_present && rhs != NULL) {
    assemble_DG_vector_forms(stage, mat, rhs, force_diagonal_blocks, block_weights, spss, refmap, neighbor_searches, u_ext, isempty, 
      marker, al, bnd, surf_pos, nat, isurf, e, trav_base, rep_element);
      if(!stage.vfsurf_mc.empty())
        assemble_multicomponent_DG_vector_forms(stage, mat, rhs, force_diagonal_blocks, block_weights, spss, refmap, neighbor_searches, u_ext, isempty, 
        marker, al, bnd, surf_pos, nat, isurf, e, trav_base, rep_element);
  }

  /***/

  // This is just cleaning after ourselves.
  // Clear the transformations from the RefMaps and all functions.
  for(unsigned int fns_i = 0; fns_i < stage.fns.size(); fns_i++)
      stage.fns[fns_i]->set_transform(neighbor_searches.get(stage.meshes[fns_i]->get_seq() 
                                      - min_dg_mesh_seq)->original_central_el_transform);

  // Also clear the transformations from the slave psss and refmaps.
  for (unsigned int i = 0; i < stage.idx.size(); i++) {
    if(isempty[stage.idx[i]])
      continue;
    spss[stage.idx[i]]->set_master_transform();
    refmap[stage.idx[i]]->force_transform(pss[stage.idx[i]]->get_transform(), pss[stage.idx[i]]->get_ctm());
  }
}

void DiscreteProblem::init_neighbors(LightArray<NeighborSearch*>& neighbor_searches, 
                                     const WeakForm::Stage& stage, const int& isurf)
{
  _F_
  // Initialize the NeighborSearches.
  for(unsigned int i = 0; i < stage.meshes.size(); i++) {
    if(!neighbor_searches.present(stage.meshes[i]->get_seq() - min_dg_mesh_seq)) {
      NeighborSearch* ns = new NeighborSearch(stage.fns[i]->get_active_element(), stage.meshes[i]);
      ns->original_central_el_transform = stage.fns[i]->get_transform();
      neighbor_searches.add(ns, stage.meshes[i]->get_seq() - min_dg_mesh_seq);
    }
  }

  // Calculate respective neighbors.
  // Also clear the initial_sub_idxs from the central element transformations 
  // of NeighborSearches with multiple neighbors.
  for(unsigned int i = 0; i < neighbor_searches.get_size(); i++) 
    if(neighbor_searches.present(i)) {
      neighbor_searches.get(i)->set_active_edge_multimesh(isurf);
      neighbor_searches.get(i)->clear_initial_sub_idx();
  }
  return;
}

void DiscreteProblem::build_multimesh_tree(DiscreteProblem::NeighborNode* root, 
                                           LightArray<NeighborSearch*>& neighbor_searches)
{
  _F_
    for(unsigned int i = 0; i < neighbor_searches.get_size(); i++) 
      if(neighbor_searches.present(i)) {
        NeighborSearch* ns = neighbor_searches.get(i);
        if(ns->n_neighbors == 1 && ns->central_n_trans[0] == 0)
          continue;
        for(unsigned int j = 0; j < ns->n_neighbors; j++)
          insert_into_multimesh_tree(root, ns->central_transformations[j], ns->central_n_trans[j]);
      }
}

void DiscreteProblem::insert_into_multimesh_tree(NeighborNode* node, 
                      unsigned int transformations [NeighborSearch::max_n_trans],  
                      unsigned int transformation_count)
{
  _F_
  // If we are already in the leaf.
  if(transformation_count == 0)
    return;
  // Both sons are null. We have to add a new Node. Let us do it for the left sone of node.
  if(node->get_left_son() == NULL && node->get_right_son() == NULL) {
    node->set_left_son(new NeighborNode(node, transformations[0]));
    insert_into_multimesh_tree(node->get_left_son(), transformations + 1, transformation_count - 1);
  }
  // At least the left son is not null (it is impossible only for the right one to be not null, because
  // the left one always gets into the tree first, as seen above).
  else {
    // The existing left son is the right one to continue through.
    if(node->get_left_son()->get_transformation() == transformations[0])
      insert_into_multimesh_tree(node->get_left_son(), transformations + 1, transformation_count - 1);
    // The right one also exists, check that it is the right one, or return an error.
    else if(node->get_right_son() != NULL) {
      if(node->get_right_son()->get_transformation() == transformations[0])
        insert_into_multimesh_tree(node->get_right_son(), transformations + 1, transformation_count - 1);
      else error("More than two possible sons in insert_into_multimesh_tree().");
    }
    // If the right one does not exist and the left one was not correct, create a right son and continue this way.
    else {
      node->set_right_son(new NeighborNode(node, transformations[0]));
      insert_into_multimesh_tree(node->get_right_son(), transformations + 1, transformation_count - 1);
    }
  }
}

/* This function is not used. It may be used when the implementation changes. 
   Will be deleted when found not necessary.*/
Hermes::vector<Hermes::vector<unsigned int>*> DiscreteProblem::get_multimesh_neighbors_transformations(DiscreteProblem::NeighborNode* multimesh_tree)
{
  _F_
  // Initialize the vector.
  Hermes::vector<Hermes::vector<unsigned int>*> running_transformations;
  // Prepare the first neighbor's vector.
  running_transformations.push_back(new Hermes::vector<unsigned int>);
  // Fill the vector.
  traverse_multimesh_tree(multimesh_tree, running_transformations);
  return running_transformations;
}

/* This function is not used. It may be used when the implementation changes. 
   Will be deleted when found not necessary.*/
void DiscreteProblem::traverse_multimesh_tree(DiscreteProblem::NeighborNode* node, 
                      Hermes::vector<Hermes::vector<unsigned int>*>& running_transformations)
{
  _F_
  // If we are in the root.
  if(node->get_transformation() == 0) {
    if(node->get_left_son() != NULL)
      traverse_multimesh_tree(node->get_left_son(), running_transformations);
    if(node->get_right_son() != NULL)
      traverse_multimesh_tree(node->get_right_son(), running_transformations);
    // Delete the vector prepared by the last accessed leaf.
    running_transformations.pop_back();
    return;
  }
  // If we are in a leaf.
  if(node->get_left_son() == NULL && node->get_right_son() == NULL) {
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
  else {
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

void DiscreteProblem::update_neighbor_search(NeighborSearch* ns, NeighborNode* multimesh_tree)
{
  _F_
  // This has to be done, because we pass ns by reference and the number of neighbors is changing.
  unsigned int num_neighbors = ns->get_num_neighbors();

  for(unsigned int i = 0; i < num_neighbors; i++) {
    // Find the node corresponding to this neighbor in the tree.
    NeighborNode* node = find_node(ns->central_transformations[i], ns->central_n_trans[i], multimesh_tree);
    // Update the NeighborSearch.
    unsigned int added = update_ns_subtree(ns, node, i);
    i += added;
    num_neighbors += added;
  }
}

DiscreteProblem::NeighborNode* DiscreteProblem::find_node(unsigned int* transformations, 
                                                          unsigned int transformation_count, 
                                                          DiscreteProblem::NeighborNode* node)
{
  _F_
  // If there are no transformations left.
  if(transformation_count == 0)
    return node;
  else {
    if(node->get_left_son() != NULL) {
      if(node->get_left_son()->get_transformation() == transformations[0])
        return find_node(transformations + 1, transformation_count - 1, node->get_left_son());
    }
    if(node->get_right_son() != NULL) {
      if(node->get_right_son()->get_transformation() == transformations[0])
        return find_node(transformations + 1, transformation_count - 1, node->get_right_son());
    }
  }
  // We always should be able to empty the transformations array.
  error("Transformation of a central element not found in the multimesh tree.");
  return NULL;
}
    
unsigned int DiscreteProblem::update_ns_subtree(NeighborSearch* ns, 
             DiscreteProblem::NeighborNode* node, unsigned int ith_neighbor)
{
  _F_
  // No subtree => no work.
  // Also check the assertion that if one son is null, then the other too.
  if(node->get_left_son() == NULL) {
    if(node->get_right_son() != NULL)
      error("Only one son (right) not null in DiscreteProblem::update_ns_subtree.");
    return 0;
  }

  // Key part.
  // Begin with storing the info about the current neighbor.
  Element* neighbor = ns->neighbors[ith_neighbor];
  NeighborSearch::NeighborEdgeInfo edge_info = ns->neighbor_edges[ith_neighbor];

  // Initialize the vector for central transformations->
  Hermes::vector<Hermes::vector<unsigned int>*> running_central_transformations;
  // Prepare the first new neighbor's vector. Push back the current transformations (in case of GO_DOWN neighborhood).
  running_central_transformations.push_back(new Hermes::vector<unsigned int>);
  for(unsigned int i = 0; i < ns->central_n_trans[ith_neighbor]; i++)
    running_central_transformations.back()->push_back(ns->central_transformations[ith_neighbor][i]);

  // Initialize the vector for neighbor transformations->
  Hermes::vector<Hermes::vector<unsigned int>*> running_neighbor_transformations;
  // Prepare the first new neighbor's vector. Push back the current transformations (in case of GO_UP/NO_TRF neighborhood).
  running_neighbor_transformations.push_back(new Hermes::vector<unsigned int>);
  for(unsigned int i = 0; i < ns->neighbor_n_trans[ith_neighbor]; i++)
    running_neighbor_transformations.back()->push_back(ns->neighbor_transformations[ith_neighbor][i]);

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
  running_central_transformations.pop_back();
  running_neighbor_transformations.pop_back();

  // Insert new neighbors.
  for(unsigned int i = 0; i < running_central_transformations.size(); i++) {
    ns->neighbors.push_back(neighbor);
    ns->neighbor_edges.push_back(edge_info);
    ns->central_n_trans[ns->n_neighbors] = running_central_transformations[i]->size();
    ns->neighbor_n_trans[ns->n_neighbors] = running_neighbor_transformations[i]->size();
    for(unsigned int ii = 0; ii < ns->central_n_trans[ns->n_neighbors]; ii++)
      ns->central_transformations[ns->n_neighbors][ii] = (*running_central_transformations[i])[ii];
    for(unsigned int ii = 0; ii < ns->neighbor_n_trans[ns->n_neighbors]; ii++)
      ns->neighbor_transformations[ns->n_neighbors][ii] = (*running_neighbor_transformations[i])[ii];
    ns->n_neighbors++;
  }

  // Return the number of neighbors deleted.
  return - 1;
}

void DiscreteProblem::traverse_multimesh_subtree(DiscreteProblem::NeighborNode* node, 
      Hermes::vector<Hermes::vector<unsigned int>*>& running_central_transformations,
      Hermes::vector<Hermes::vector<unsigned int>*>& running_neighbor_transformations, 
      const NeighborSearch::NeighborEdgeInfo& edge_info, const int& active_edge, const int& mode)
{
  _F_
  // If we are in a leaf.
  if(node->get_left_son() == NULL && node->get_right_son() == NULL) {
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
  else {
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


void DiscreteProblem::assemble_surface_matrix_forms(WeakForm::Stage& stage, 
       SparseMatrix* mat, Vector* rhs, bool force_diagonal_blocks, 
       Table* block_weights, Hermes::vector<PrecalcShapeset *>& spss, 
       Hermes::vector<RefMap *>& refmap, Hermes::vector<Solution *>& u_ext, 
       Hermes::vector<bool>& isempty, int marker, Hermes::vector<AsmList *>& al, 
       bool bnd, SurfPos& surf_pos, Hermes::vector<bool>& nat, 
       int isurf, Element** e, Element* trav_base)
{
  _F_
  for (unsigned int ww = 0; ww < stage.mfsurf.size(); ww++) {
    WeakForm::MatrixFormSurf* mfs = stage.mfsurf[ww];
    int m = mfs->i;
    int n = mfs->j;
    if (isempty[m] || isempty[n]) continue;
    if (!nat[m] || !nat[n]) continue;
    if (fabs(mfs->scaling_factor) < 1e-12) continue;
    if (mfs->area == H2D_DG_INNER_EDGE) continue;
    if (mfs->area != HERMES_ANY && mfs->area != H2D_DG_BOUNDARY_EDGE 
      && !(marker == boundary_markers_conversion->get_internal_marker(mfs->area))) continue;

    // If a block scaling table is provided, and if the scaling coefficient
    // A_mn for this block is zero, then the form does not need to be assembled.
    double block_scaling_coeff = 1.;
    if (block_weights != NULL) {
      block_scaling_coeff = block_weights->get_A(m, n);
      if (fabs(block_scaling_coeff) < 1e-12) 
        continue;
    }

    surf_pos.base = trav_base;
    surf_pos.space_v = spaces[m];
    surf_pos.space_u = spaces[n];

    scalar **local_stiffness_matrix = get_matrix_buffer(std::max(al[m]->cnt, al[n]->cnt));
    for (unsigned int i = 0; i < al[m]->cnt; i++) {
      if (al[m]->dof[i] < 0) continue;
      spss[m]->set_active_shape(al[m]->idx[i]);
      for (unsigned int j = 0; j < al[n]->cnt; j++) {
        pss[n]->set_active_shape(al[n]->idx[j]);
        if (al[n]->dof[j] < 0) {
          // Linear problems only: Subtracting Dirichlet lift contribution from the RHS:
          if (rhs != NULL && this->is_linear) {
            // Numerical integration performed only if all coefficients multiplying the form are nonzero
            // and if the basis function is active.
            if (std::abs(al[m]->coef[i]) > 1e-12 && std::abs(al[n]->coef[j]) > 1e-12
                && al[m]->dof[i] >= 0) {
              scalar val = eval_form(mfs, u_ext, pss[n], spss[m], refmap[n],
                                     refmap[m], &surf_pos) * al[n]->coef[j] * al[m]->coef[i];
              rhs->add(al[m]->dof[i], -val);
            }
          }
        }
        else if (mat != NULL) {
          scalar val = 0;
          // Numerical integration performed only if all coefficients multiplying the form are nonzero.
          if (std::abs(al[m]->coef[i]) > 1e-12 && std::abs(al[n]->coef[j]) > 1e-12) {
            val = block_scaling_coeff * eval_form(mfs, u_ext, pss[n], spss[m], refmap[n],
                                                  refmap[m], &surf_pos) * al[n]->coef[j] * al[m]->coef[i];
          }
          local_stiffness_matrix[i][j] = val;
        }
      }
    }
    if (mat != NULL)
      mat->add(al[m]->cnt, al[n]->cnt, local_stiffness_matrix, al[m]->dof, al[n]->dof);
  }
}
void DiscreteProblem::assemble_multicomponent_surface_matrix_forms(WeakForm::Stage& stage, 
       SparseMatrix* mat, Vector* rhs, bool force_diagonal_blocks, 
       Table* block_weights, Hermes::vector<PrecalcShapeset *>& spss, 
       Hermes::vector<RefMap *>& refmap, Hermes::vector<Solution *>& u_ext, 
       Hermes::vector<bool>& isempty, int marker, Hermes::vector<AsmList *>& al, 
       bool bnd, SurfPos& surf_pos, Hermes::vector<bool>& nat, 
       int isurf, Element** e, Element* trav_base)
{
  _F_
  for (unsigned int ww = 0; ww < stage.mfsurf_mc.size(); ww++) {
    WeakForm::MultiComponentMatrixFormSurf* mfs = stage.mfsurf_mc[ww];
    unsigned int m = mfs->coordinates[0].first;
    unsigned int n = mfs->coordinates[0].second;

    if (!nat[m] || !nat[n]) continue;
    if (fabs(mfs->scaling_factor) < 1e-12) continue;
    if (mfs->area == H2D_DG_INNER_EDGE) continue;
    if (mfs->area != HERMES_ANY && mfs->area != H2D_DG_BOUNDARY_EDGE 
      && !(marker == boundary_markers_conversion->get_internal_marker(mfs->area))) continue;

    // If a block scaling table is provided, and if the scaling coefficient
    // A_mn for this block is zero, then the form does not need to be assembled.
    // If a block scaling table is provided, and if the scaling coefficient
    // A_mn for this block is zero, then the form does not need to be assembled.
    Hermes::vector<double> block_scaling_coeffs;
    for(unsigned int coordinate_i = 0; coordinate_i < mfs->coordinates.size(); coordinate_i++)
      if (block_weights != NULL)
        block_scaling_coeffs.push_back(block_weights->get_A(mfs->coordinates[coordinate_i].first, mfs->coordinates[coordinate_i].second));
      else
        block_scaling_coeffs.push_back(1);

    surf_pos.base = trav_base;
    surf_pos.space_v = spaces[m];
    surf_pos.space_u = spaces[n];

    for (unsigned int i = 0; i < al[m]->cnt; i++) {
      if (al[m]->dof[i] < 0) continue;
      spss[m]->set_active_shape(al[m]->idx[i]);
      for (unsigned int j = 0; j < al[n]->cnt; j++) {
        pss[n]->set_active_shape(al[n]->idx[j]);
        if (al[n]->dof[j] < 0) {
          // Linear problems only: Subtracting Dirichlet lift contribution from the RHS:
          if (rhs != NULL && this->is_linear) {
            // Numerical integration performed only if all coefficients multiplying the form are nonzero
            // and if the basis function is active.
            if (std::abs(al[m]->coef[i]) > 1e-12 && std::abs(al[n]->coef[j]) > 1e-12 && al[m]->dof[i] >= 0) {
              Hermes::vector<scalar> result;
              eval_form(mfs, u_ext, pss[n], spss[m], refmap[n], refmap[m], &surf_pos, result);
              for(unsigned int coordinate_i = 0; coordinate_i < mfs->coordinates.size(); coordinate_i++)
                rhs->add(al[mfs->coordinates[coordinate_i].first]->dof[i], -result[coordinate_i] * al[n]->coef[j] * al[m]->coef[i]);
            }
          }
        }
        else 
          if (mat != NULL) {
            // Numerical integration performed only if all coefficients multiplying the form are nonzero.
            if (std::abs(al[m]->coef[i]) > 1e-12 && std::abs(al[n]->coef[j]) > 1e-12) {
              Hermes::vector<scalar> result;
              eval_form(mfs, u_ext, pss[n], spss[m], refmap[n], refmap[m], &surf_pos, result);
              for(unsigned int coordinate_i = 0; coordinate_i < mfs->coordinates.size(); coordinate_i++)
                mat->add(al[mfs->coordinates[coordinate_i].first]->dof[i], al[mfs->coordinates[coordinate_i].second]->dof[j],
                result[coordinate_i] * block_scaling_coeffs[coordinate_i] * al[n]->coef[j] * al[m]->coef[i]);
            }
        }
      }
    }
  }
}

void DiscreteProblem::assemble_surface_vector_forms(WeakForm::Stage& stage, 
       SparseMatrix* mat, Vector* rhs, bool force_diagonal_blocks, 
       Table* block_weights, Hermes::vector<PrecalcShapeset *>& spss, 
       Hermes::vector<RefMap *>& refmap, Hermes::vector<Solution *>& u_ext, 
       Hermes::vector<bool>& isempty, int marker, Hermes::vector<AsmList *>& al, 
       bool bnd, SurfPos& surf_pos, Hermes::vector<bool>& nat, 
       int isurf, Element** e, Element* trav_base)
{
  _F_

  if (rhs == NULL) return;

  for (unsigned int ww = 0; ww < stage.vfsurf.size(); ww++) {
    WeakForm::VectorFormSurf* vfs = stage.vfsurf[ww];
    int m = vfs->i;
    if (isempty[m]) continue;
    if (fabs(vfs->scaling_factor) < 1e-12) continue;
    if (vfs->area == H2D_DG_INNER_EDGE) continue;
    if (vfs->area != HERMES_ANY && vfs->area != H2D_DG_BOUNDARY_EDGE 
        && !(marker == boundary_markers_conversion->get_internal_marker(vfs->area))) continue;

    if (vfs->area == HERMES_ANY && !nat[m]) continue;

    surf_pos.base = trav_base;
    surf_pos.space_v = spaces[m];

    for (unsigned int i = 0; i < al[m]->cnt; i++) {
      if (al[m]->dof[i] < 0) continue;

      spss[m]->set_active_shape(al[m]->idx[i]);

      // Numerical integration performed only if the coefficient multiplying the form is nonzero.
      if (std::abs(al[m]->coef[i]) > 1e-12) {   
        rhs->add(al[m]->dof[i], eval_form(vfs, u_ext, spss[m], refmap[m], &surf_pos) * al[m]->coef[i]);
      }
    }
  }
}
void DiscreteProblem::assemble_multicomponent_surface_vector_forms(WeakForm::Stage& stage, 
       SparseMatrix* mat, Vector* rhs, bool force_diagonal_blocks, 
       Table* block_weights, Hermes::vector<PrecalcShapeset *>& spss, 
       Hermes::vector<RefMap *>& refmap, Hermes::vector<Solution *>& u_ext, 
       Hermes::vector<bool>& isempty, int marker, Hermes::vector<AsmList *>& al, 
       bool bnd, SurfPos& surf_pos, Hermes::vector<bool>& nat, 
       int isurf, Element** e, Element* trav_base)
{
  _F_

  if (rhs == NULL) return;

  for (unsigned int ww = 0; ww < stage.vfsurf_mc.size(); ww++) {
    WeakForm::MultiComponentVectorFormSurf* vfs = stage.vfsurf_mc[ww];
    unsigned int m = vfs->coordinates[0];
    if (fabs(vfs->scaling_factor) < 1e-12) continue;
    if (vfs->area == H2D_DG_INNER_EDGE) continue;
    if (vfs->area != HERMES_ANY && vfs->area != H2D_DG_BOUNDARY_EDGE 
        && !(marker == boundary_markers_conversion->get_internal_marker(vfs->area))) continue;

    if (vfs->area == HERMES_ANY && !nat[m]) continue;

    surf_pos.base = trav_base;
    surf_pos.space_v = spaces[m];

    for (unsigned int i = 0; i < al[m]->cnt; i++) {
      if (al[m]->dof[i] < 0) continue;

      spss[m]->set_active_shape(al[m]->idx[i]);

      // Numerical integration performed only if the coefficient multiplying the form is nonzero.
      if (std::abs(al[m]->coef[i]) > 1e-12) {   
        Hermes::vector<scalar> result;
        eval_form(vfs, u_ext, spss[m], refmap[m], &surf_pos, result);
        for(unsigned int coordinate_i = 0; coordinate_i < vfs->coordinates.size(); coordinate_i++)
          rhs->add(al[vfs->coordinates[coordinate_i]]->dof[i], result[coordinate_i] * al[vfs->coordinates[coordinate_i]]->coef[i]);
      }
    }
  }
}

void DiscreteProblem::assemble_DG_matrix_forms(WeakForm::Stage& stage, 
       SparseMatrix* mat, Vector* rhs, bool force_diagonal_blocks, 
       Table* block_weights, Hermes::vector<PrecalcShapeset *>& spss, 
       Hermes::vector<RefMap *>& refmap, Hermes::vector<PrecalcShapeset *>& npss, 
       Hermes::vector<PrecalcShapeset *>& nspss, Hermes::vector<RefMap *>& nrefmap, 
       LightArray<NeighborSearch*>& neighbor_searches, Hermes::vector<Solution *>& u_ext, 
       Hermes::vector<bool>& isempty, int marker, Hermes::vector<AsmList *>& al, bool bnd, 
       SurfPos& surf_pos, Hermes::vector<bool>& nat, int isurf, Element** e, 
       Element* trav_base, Element* rep_element)
{ 
  _F_
  for (unsigned int ww = 0; ww < stage.mfsurf.size(); ww++) {
    WeakForm::MatrixFormSurf* mfs = stage.mfsurf[ww];
    if (mfs->area != H2D_DG_INNER_EDGE) continue;
    int m = mfs->i;
    int n = mfs->j;

    if (isempty[m] || isempty[n]) continue;
    if (fabs(mfs->scaling_factor) < 1e-12) continue;

    surf_pos.base = trav_base;
    surf_pos.space_v = spaces[m];
    surf_pos.space_u = spaces[n];
            
    // Create the extended shapeset on the union of the central element and its current neighbor.
    NeighborSearch::ExtendedShapeset* ext_asmlist_u = neighbor_searches.get(stage.meshes[n]->get_seq() - min_dg_mesh_seq)->create_extended_asmlist(spaces[n], al[n]);
    NeighborSearch::ExtendedShapeset* ext_asmlist_v = neighbor_searches.get(stage.meshes[m]->get_seq() - min_dg_mesh_seq)->create_extended_asmlist(spaces[m], al[m]);

    // If a block scaling table is provided, and if the scaling coefficient
    // A_mn for this block is zero, then the form does not need to be assembled.
    double block_scaling_coeff = 1.;
    if (block_weights != NULL) {
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

    scalar **local_stiffness_matrix = get_matrix_buffer(std::max(ext_asmlist_u->cnt, ext_asmlist_v->cnt));
    for (int i = 0; i < ext_asmlist_v->cnt; i++) {
      if (ext_asmlist_v->dof[i] < 0) 
        continue;
      // Choose the correct shapeset for the test function.
      if (!ext_asmlist_v->has_support_on_neighbor(i)) {
        spss[m]->set_active_shape(ext_asmlist_v->central_al->idx[i]);
        fv = spss[m];
        rv = refmap[m];
        support_neigh_v = false;
      }
      else {
        nspss[m]->set_active_shape(ext_asmlist_v->neighbor_al->idx[i - ext_asmlist_v->central_al->cnt]);
        fv = nspss[m];
        rv = nrefmap[m];
        support_neigh_v = true;
      }
      for (int j = 0; j < ext_asmlist_u->cnt; j++) {
        // Choose the correct shapeset for the solution function.
        if (!ext_asmlist_u->has_support_on_neighbor(j)) {
          pss[n]->set_active_shape(ext_asmlist_u->central_al->idx[j]);
          fu = pss[n];
          ru = refmap[n];
          support_neigh_u = false;
        }
        else {
          npss[n]->set_active_shape(ext_asmlist_u->neighbor_al->idx[j - ext_asmlist_u->central_al->cnt]);
          fu = npss[n];
          ru = nrefmap[n];
          support_neigh_u = true;
        }

        if (ext_asmlist_u->dof[j] < 0) {
          if (rhs != NULL && this->is_linear) {
            // Evaluate the form with the activated discontinuous shape functions.
            scalar val = eval_dg_form(mfs, u_ext, fu, fv, refmap[n], ru, rv, support_neigh_u, support_neigh_v, 
                                      &surf_pos, neighbor_searches, stage.meshes[n]->get_seq() - min_dg_mesh_seq, 
                                      stage.meshes[m]->get_seq() - min_dg_mesh_seq)
              * (support_neigh_u ? ext_asmlist_u->neighbor_al->coef[j - ext_asmlist_u->central_al->cnt]: ext_asmlist_u->central_al->coef[j])
              * (support_neigh_v ? ext_asmlist_v->neighbor_al->coef[i - ext_asmlist_v->central_al->cnt]: ext_asmlist_v->central_al->coef[i]);
            // Add the contribution to the global dof index.
            rhs->add(ext_asmlist_v->dof[i], -val);
          }
        }
        else if (mat != NULL) {
          scalar val = block_scaling_coeff * eval_dg_form(mfs, u_ext, fu, fv, refmap[n], ru, rv, support_neigh_u, support_neigh_v, &surf_pos, neighbor_searches, stage.meshes[n]->get_seq() - min_dg_mesh_seq, stage.meshes[m]->get_seq() - min_dg_mesh_seq)
            * (support_neigh_u ? ext_asmlist_u->neighbor_al->coef[j - ext_asmlist_u->central_al->cnt]: ext_asmlist_u->central_al->coef[j])
            * (support_neigh_v ? ext_asmlist_v->neighbor_al->coef[i - ext_asmlist_v->central_al->cnt]: ext_asmlist_v->central_al->coef[i]);
          local_stiffness_matrix[i][j] = val;
        }
      }
    }
    if (mat != NULL)
      mat->add(ext_asmlist_v->cnt, ext_asmlist_u->cnt, local_stiffness_matrix, ext_asmlist_v->dof, ext_asmlist_u->dof);
  }
}
void DiscreteProblem::assemble_multicomponent_DG_matrix_forms(WeakForm::Stage& stage, 
       SparseMatrix* mat, Vector* rhs, bool force_diagonal_blocks, 
       Table* block_weights, Hermes::vector<PrecalcShapeset *>& spss, 
       Hermes::vector<RefMap *>& refmap, Hermes::vector<PrecalcShapeset *>& npss, 
       Hermes::vector<PrecalcShapeset *>& nspss, Hermes::vector<RefMap *>& nrefmap, 
       LightArray<NeighborSearch*>& neighbor_searches, Hermes::vector<Solution *>& u_ext, 
       Hermes::vector<bool>& isempty, int marker, Hermes::vector<AsmList *>& al, bool bnd, 
       SurfPos& surf_pos, Hermes::vector<bool>& nat, int isurf, Element** e, 
       Element* trav_base, Element* rep_element)
{ 
  _F_
  for (unsigned int ww = 0; ww < stage.mfsurf_mc.size(); ww++) {
    WeakForm::MultiComponentMatrixFormSurf* mfs = stage.mfsurf_mc[ww];
    if (mfs->area != H2D_DG_INNER_EDGE) continue;

    if (fabs(mfs->scaling_factor) < 1e-12) continue;

    unsigned int m = mfs->coordinates[0].first;
    unsigned int n = mfs->coordinates[0].second;

    surf_pos.base = trav_base;
    surf_pos.space_v = spaces[m];
    surf_pos.space_u = spaces[n];
            
    // Create the extended shapeset on the union of the central element and its current neighbor.
    Hermes::vector<NeighborSearch::ExtendedShapeset*> ext_asmlists;
    Hermes::vector<unsigned int> coordinates_processed;
    for(unsigned int coordinate_i = 0; coordinate_i < mfs->coordinates.size(); coordinate_i++) {
      bool new_coordinate = true;
      for(unsigned int i = 0; i < coordinates_processed.size(); i++)
        if(coordinates_processed[i] == mfs->coordinates[coordinate_i].first)
          new_coordinate = false;
      if(new_coordinate) {
        coordinates_processed.push_back(mfs->coordinates[coordinate_i].first);
        ext_asmlists.push_back(neighbor_searches.get(stage.meshes[mfs->coordinates[coordinate_i].first]->get_seq() - min_dg_mesh_seq)->create_extended_asmlist(spaces[mfs->coordinates[coordinate_i].first], al[mfs->coordinates[coordinate_i].first]));
      }

      new_coordinate = true;
      for(unsigned int i = 0; i < coordinates_processed.size(); i++)
        if(coordinates_processed[i] == mfs->coordinates[coordinate_i].second)
          new_coordinate = false;
      if(new_coordinate) {
        coordinates_processed.push_back(mfs->coordinates[coordinate_i].second);
        ext_asmlists.push_back(neighbor_searches.get(stage.meshes[mfs->coordinates[coordinate_i].second]->get_seq() - min_dg_mesh_seq)->create_extended_asmlist(spaces[mfs->coordinates[coordinate_i].second], al[mfs->coordinates[coordinate_i].second]));
      }
    }
      
    NeighborSearch::ExtendedShapeset* ext_asmlist_u = neighbor_searches.get(stage.meshes[n]->get_seq() - min_dg_mesh_seq)->create_extended_asmlist(spaces[n], al[n]);
    NeighborSearch::ExtendedShapeset* ext_asmlist_v = neighbor_searches.get(stage.meshes[m]->get_seq() - min_dg_mesh_seq)->create_extended_asmlist(spaces[m], al[m]);

    // If a block scaling table is provided, and if the scaling coefficient
    // A_mn for this block is zero, then the form does not need to be assembled.
    Hermes::vector<double> block_scaling_coeffs;
    for(unsigned int coordinate_i = 0; coordinate_i < mfs->coordinates.size(); coordinate_i++)
      if (block_weights != NULL)
        block_scaling_coeffs.push_back(block_weights->get_A(mfs->coordinates[coordinate_i].first, mfs->coordinates[coordinate_i].second));
      else
        block_scaling_coeffs.push_back(1);

    // Precalc shapeset and refmaps used for the evaluation.
    PrecalcShapeset* fu;
    PrecalcShapeset* fv;
    RefMap* ru;
    RefMap* rv;
    bool support_neigh_u, support_neigh_v;

    for (int i = 0; i < ext_asmlist_v->cnt; i++) {
      if (ext_asmlist_v->dof[i] < 0) 
        continue;
      // Choose the correct shapeset for the test function.
      if (!ext_asmlist_v->has_support_on_neighbor(i)) {
        spss[m]->set_active_shape(ext_asmlist_v->central_al->idx[i]);
        fv = spss[m];
        rv = refmap[m];
        support_neigh_v = false;
      }
      else {
        nspss[m]->set_active_shape(ext_asmlist_v->neighbor_al->idx[i - ext_asmlist_v->central_al->cnt]);
        fv = nspss[m];
        rv = nrefmap[m];
        support_neigh_v = true;
      }
      for (int j = 0; j < ext_asmlist_u->cnt; j++) {
        // Choose the correct shapeset for the solution function.
        if (!ext_asmlist_u->has_support_on_neighbor(j)) {
          pss[n]->set_active_shape(ext_asmlist_u->central_al->idx[j]);
          fu = pss[n];
          ru = refmap[n];
          support_neigh_u = false;
        }
        else {
          npss[n]->set_active_shape(ext_asmlist_u->neighbor_al->idx[j - ext_asmlist_u->central_al->cnt]);
          fu = npss[n];
          ru = nrefmap[n];
          support_neigh_u = true;
        }

        if (ext_asmlist_u->dof[j] < 0) {
          if (rhs != NULL && this->is_linear) {
            // Evaluate the form with the activated discontinuous shape functions.
            if (std::abs(al[m]->coef[i]) > 1e-12 && std::abs(al[n]->coef[j]) > 1e-12 && al[m]->dof[i] >= 0) {
              Hermes::vector<scalar> result;
              eval_dg_form(mfs, u_ext, fu, fv, refmap[n], ru, rv, support_neigh_u, support_neigh_v, 
                                      &surf_pos, neighbor_searches, stage.meshes[n]->get_seq() - min_dg_mesh_seq, 
                                      stage.meshes[m]->get_seq() - min_dg_mesh_seq, result);
              for(unsigned int coordinate_i = 0; coordinate_i < mfs->coordinates.size(); coordinate_i++)
                rhs->add(ext_asmlists[mfs->coordinates[coordinate_i].first]->dof[i], -result[coordinate_i] 
                * (support_neigh_u ? ext_asmlist_u->neighbor_al->coef[j - ext_asmlist_u->central_al->cnt]: ext_asmlist_u->central_al->coef[j])
                * (support_neigh_v ? ext_asmlist_v->neighbor_al->coef[i - ext_asmlist_v->central_al->cnt]: ext_asmlist_v->central_al->coef[i]));
            }
          }
        }
        else if (mat != NULL) {
          if (std::abs(al[m]->coef[i]) > 1e-12 && std::abs(al[n]->coef[j]) > 1e-12) {
              Hermes::vector<scalar> result;
              eval_dg_form(mfs, u_ext, fu, fv, refmap[n], ru, rv, support_neigh_u, support_neigh_v, &surf_pos, neighbor_searches, stage.meshes[n]->get_seq() - min_dg_mesh_seq, stage.meshes[m]->get_seq() - min_dg_mesh_seq, result);
              for(unsigned int coordinate_i = 0; coordinate_i < mfs->coordinates.size(); coordinate_i++)
                mat->add(ext_asmlists[mfs->coordinates[coordinate_i].first]->dof[i], ext_asmlists[mfs->coordinates[coordinate_i].second]->dof[j],
                result[coordinate_i] * block_scaling_coeffs[coordinate_i] 
                * (support_neigh_u ? ext_asmlist_u->neighbor_al->coef[j - ext_asmlist_u->central_al->cnt]: ext_asmlist_u->central_al->coef[j])
                * (support_neigh_v ? ext_asmlist_v->neighbor_al->coef[i - ext_asmlist_v->central_al->cnt]: ext_asmlist_v->central_al->coef[i]));
          }
        }
      }
    }
  }
}

void DiscreteProblem::assemble_DG_vector_forms(WeakForm::Stage& stage, 
      SparseMatrix* mat, Vector* rhs, bool force_diagonal_blocks, Table* block_weights,
       Hermes::vector<PrecalcShapeset *>& spss, Hermes::vector<RefMap *>& refmap, LightArray<NeighborSearch*>& neighbor_searches, Hermes::vector<Solution *>& u_ext, 
       Hermes::vector<bool>& isempty, int marker, Hermes::vector<AsmList *>& al, bool bnd, SurfPos& surf_pos, Hermes::vector<bool>& nat, 
       int isurf, Element** e, Element* trav_base, Element* rep_element)
{
  _F_

  if (rhs == NULL) return;

  for (unsigned int ww = 0; ww < stage.vfsurf.size(); ww++) {
    WeakForm::VectorFormSurf* vfs = stage.vfsurf[ww];
    if (vfs->area != H2D_DG_INNER_EDGE) continue;
    int m = vfs->i;
    if (isempty[m]) continue;
    if (fabs(vfs->scaling_factor) < 1e-12) continue;
            
    // Here we use the standard pss, possibly just transformed by NeighborSearch.
    for (unsigned int i = 0; i < al[m]->cnt; i++) {
      if (al[m]->dof[i] < 0) continue;
      spss[m]->set_active_shape(al[m]->idx[i]);
      rhs->add(al[m]->dof[i], eval_dg_form(vfs, u_ext, spss[m], refmap[m], &surf_pos, neighbor_searches, stage.meshes[m]->get_seq() - min_dg_mesh_seq) * al[m]->coef[i]);
    }
  }
}
void DiscreteProblem::assemble_multicomponent_DG_vector_forms(WeakForm::Stage& stage, 
      SparseMatrix* mat, Vector* rhs, bool force_diagonal_blocks, Table* block_weights,
       Hermes::vector<PrecalcShapeset *>& spss, Hermes::vector<RefMap *>& refmap, LightArray<NeighborSearch*>& neighbor_searches, Hermes::vector<Solution *>& u_ext, 
       Hermes::vector<bool>& isempty, int marker, Hermes::vector<AsmList *>& al, bool bnd, SurfPos& surf_pos, Hermes::vector<bool>& nat, 
       int isurf, Element** e, Element* trav_base, Element* rep_element)
{
  _F_

  if (rhs == NULL) return;

  for (unsigned int ww = 0; ww < stage.vfsurf_mc.size(); ww++) {
    WeakForm::MultiComponentVectorFormSurf* vfs = stage.vfsurf_mc[ww];
    if (vfs->area != H2D_DG_INNER_EDGE) continue;
    if (fabs(vfs->scaling_factor) < 1e-12) continue;
            
    // Here we use the standard pss, possibly just transformed by NeighborSearch.
    unsigned int m = vfs->coordinates[0];
    for (unsigned int i = 0; i < al[m]->cnt; i++) {
      if (al[m]->dof[i] < 0) continue;
      Hermes::vector<scalar> result;
      spss[m]->set_active_shape(al[m]->idx[i]);
      eval_dg_form(vfs, u_ext, spss[m], refmap[m], &surf_pos, neighbor_searches, stage.meshes[m]->get_seq() - min_dg_mesh_seq, result);
      for(unsigned int coordinate_i = 0; coordinate_i < vfs->coordinates.size(); coordinate_i++)
        rhs->add(al[vfs->coordinates[coordinate_i]]->dof[i], result[coordinate_i] * al[vfs->coordinates[coordinate_i]]->coef[i]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Initialize integration order for external functions
ExtData<Ord>* DiscreteProblem::init_ext_fns_ord(Hermes::vector<MeshFunction *> &ext)
{
  _F_
  ExtData<Ord>* fake_ext = new ExtData<Ord>;

  // External functions.
  fake_ext->nf = ext.size();
  Func<Ord>** fake_ext_fn = new Func<Ord>*[fake_ext->nf];
  for (int i = 0; i < fake_ext->nf; i++)
    fake_ext_fn[i] = get_fn_ord(ext[i]->get_fn_order());
  fake_ext->fn = fake_ext_fn;
  
  return fake_ext;
}

// Initialize external functions (obtain values, derivatives,...)
ExtData<scalar>* DiscreteProblem::init_ext_fns(Hermes::vector<MeshFunction *> &ext, 
                                               RefMap *rm, const int order)
{
  _F_
  ExtData<scalar>* ext_data = new ExtData<scalar>;

  // Copy external functions.
  Func<scalar>** ext_fn = new Func<scalar>*[ext.size()];
  for (unsigned i = 0; i < ext.size(); i++) {
    if (ext[i] != NULL) ext_fn[i] = init_fn(ext[i], order);
    else ext_fn[i] = NULL;
  }
  ext_data->nf = ext.size();
  ext_data->fn = ext_fn;

  return ext_data;
}

// Initialize integration order on a given edge for external functions
ExtData<Ord>* DiscreteProblem::init_ext_fns_ord(Hermes::vector<MeshFunction *> &ext, 
                                                int edge)
{
  _F_
  ExtData<Ord>* fake_ext = new ExtData<Ord>;

  // External functions.
  fake_ext->nf = ext.size();
  Func<Ord>** fake_ext_fn = new Func<Ord>*[fake_ext->nf];
  for (int i = 0; i < fake_ext->nf; i++)
    fake_ext_fn[i] = get_fn_ord(ext[i]->get_edge_fn_order(edge));
  fake_ext->fn = fake_ext_fn;

  return fake_ext;
}

// Initialize discontinuous external functions (obtain values, derivatives,... on both sides of the
// supplied NeighborSearch's active edge).
ExtData<scalar>* DiscreteProblem::init_ext_fns(Hermes::vector<MeshFunction *> &ext, 
                                               LightArray<NeighborSearch*>& neighbor_searches, int order)
{
  _F_
  Func<scalar>** ext_fns = new Func<scalar>*[ext.size()];
  for(unsigned int j = 0; j < ext.size(); j++) {
    neighbor_searches.get(ext[j]->get_mesh()->get_seq() - min_dg_mesh_seq)->set_quad_order(order);
    ext_fns[j] = neighbor_searches.get(ext[j]->get_mesh()->get_seq() - min_dg_mesh_seq)->init_ext_fn(ext[j]);
  }

  ExtData<scalar>* ext_data = new ExtData<scalar>;
  ext_data->fn = ext_fns;
  ext_data->nf = ext.size();

  return ext_data;
}

// Initialize integration order for discontinuous external functions.
ExtData<Ord>* DiscreteProblem::init_ext_fns_ord(Hermes::vector<MeshFunction *> &ext, 
                                                LightArray<NeighborSearch*>& neighbor_searches)
{
  _F_
  Func<Ord>** fake_ext_fns = new Func<Ord>*[ext.size()];
  for (unsigned int j = 0; j < ext.size(); j++)
    fake_ext_fns[j] = init_ext_fn_ord(neighbor_searches.get(ext[j]->get_mesh()->get_seq() - min_dg_mesh_seq), ext[j]);

  ExtData<Ord>* fake_ext = new ExtData<Ord>;
  fake_ext->fn = fake_ext_fns;
  fake_ext->nf = ext.size();

  return fake_ext;
}

// Initialize shape function values and derivatives (fill in the cache)
Func<double>* DiscreteProblem::get_fn(PrecalcShapeset *fu, RefMap *rm, const int order)
{
  _F_
  if(rm->is_jacobian_const()) {
    AssemblingCaches::KeyConst key(256 - fu->get_active_shape(), order, fu->get_transform(), fu->get_shapeset()->get_id(), rm->get_const_inv_ref_map());
    if(rm->get_active_element()->get_mode() == HERMES_MODE_TRIANGLE) {
      if(assembling_caches.const_cache_fn_triangles.find(key) == assembling_caches.const_cache_fn_triangles.end())
        assembling_caches.const_cache_fn_triangles[key] = init_fn(fu, rm, order);
      return assembling_caches.const_cache_fn_triangles[key];
    }
    else {
      if(assembling_caches.const_cache_fn_quads.find(key) == assembling_caches.const_cache_fn_quads.end())
        assembling_caches.const_cache_fn_quads[key] = init_fn(fu, rm, order);
      return assembling_caches.const_cache_fn_quads[key];
    }
  }
  else {
    AssemblingCaches::KeyNonConst key(256 - fu->get_active_shape(), order, 
                                      fu->get_transform(), fu->get_shapeset()->get_id());
    if(rm->get_active_element()->get_mode() == HERMES_MODE_TRIANGLE) {
      if(assembling_caches.cache_fn_triangles.find(key) == assembling_caches.cache_fn_triangles.end())
        assembling_caches.cache_fn_triangles[key] = init_fn(fu, rm, order);
      return assembling_caches.cache_fn_triangles[key];
    }
    else {
      if(assembling_caches.cache_fn_quads.find(key) == assembling_caches.cache_fn_quads.end())
        assembling_caches.cache_fn_quads[key] = init_fn(fu, rm, order);
      return assembling_caches.cache_fn_quads[key];
    }
  }
}

// Initialize shape function values and derivatives (fill in the cache)
Func<Ord>* DiscreteProblem::get_fn_ord(const int order)
{
  _F_
  assert(order >= 0);
  unsigned int cached_order = (unsigned int) order;
  if(!assembling_caches.cache_fn_ord.present(cached_order))
    assembling_caches.cache_fn_ord.add(init_fn_ord(cached_order), cached_order);
  return assembling_caches.cache_fn_ord.get(cached_order);
}

// Caching transformed values
void DiscreteProblem::init_cache()
{
  _F_
  for (int i = 0; i < g_max_quad + 1 + 4 * g_max_quad + 4; i++)
  {
    cache_e[i] = NULL;
    cache_jwt[i] = NULL;
  }
}

void DiscreteProblem::delete_single_geom_cache(int order)
{
  if (cache_e[order] != NULL) {
    cache_e[order]->free();
    delete cache_e[order];
    cache_e[order] = NULL;
    delete [] cache_jwt[order];
  }
}

void DiscreteProblem::delete_cache()
{
  _F_
  for (int i = 0; i < g_max_quad + 1 + 4 * g_max_quad + 4; i++)
  {
    if (cache_e[i] != NULL)
    {
      cache_e[i]->free(); delete cache_e[i];
      delete [] cache_jwt[i];
    }
  }
  
  for (std::map<AssemblingCaches::KeyNonConst, Func<double>*, AssemblingCaches::CompareNonConst>::const_iterator it = assembling_caches.cache_fn_quads.begin();
       it != assembling_caches.cache_fn_quads.end(); it++)
  {
    (it->second)->free_fn(); delete (it->second);
  }
  assembling_caches.cache_fn_quads.clear();

  for (std::map<AssemblingCaches::KeyNonConst, Func<double>*, AssemblingCaches::CompareNonConst>::const_iterator it = assembling_caches.cache_fn_triangles.begin();
       it != assembling_caches.cache_fn_triangles.end(); it++)
  {
    (it->second)->free_fn(); delete (it->second);
  }
  assembling_caches.cache_fn_triangles.clear();
}

DiscontinuousFunc<Ord>* DiscreteProblem::init_ext_fn_ord(NeighborSearch* ns, MeshFunction* fu)
{
  _F_
  int inc = (fu->get_num_components() == 2) ? 1 : 0;
  int central_order = fu->get_edge_fn_order(ns->active_edge) + inc;
  int neighbor_order = fu->get_edge_fn_order(ns->neighbor_edge.local_num_of_edge) + inc;
  return new DiscontinuousFunc<Ord>(get_fn_ord(central_order), get_fn_ord(neighbor_order));
}

//  Evaluation of forms  ///////////////////////////////////////////////////////////////////////

// Volume matrix forms.

scalar DiscreteProblem::eval_form(WeakForm::MatrixFormVol *mfv, 
                                  Hermes::vector<Solution *> u_ext,
                                  PrecalcShapeset *fu, PrecalcShapeset *fv, 
                                  RefMap *ru, RefMap *rv)
{
  _F_
  scalar result = 0;

  if (mfv->adapt_eval == false) {
    // Determine the integration order by parsing the form.
    int order = calc_order_matrix_form_vol(mfv, u_ext, fu, fv, ru, rv);
    // Perform non-adaptive numerical quadrature of order "order".
    result = eval_form_subelement(order, mfv, u_ext, fu, fv, ru, rv);
  }
  else {
    // Perform adaptive numerical quadrature starting with order = 2.
    // TODO: The choice of initial order matters a lot for efficiency,
    //       this needs more research.
    Shapeset* fu_shapeset = fu->get_shapeset();
    Shapeset* fv_shapeset = fv->get_shapeset();
    int fu_index = fu->get_active_shape();
    int fv_index = fv->get_active_shape();
    int fu_order = fu_shapeset->get_order(fu_index);
    int fv_order = fv_shapeset->get_order(fv_index);
    int fu_order_h = H2D_GET_H_ORDER(fu_order);
    int fu_order_v = H2D_GET_V_ORDER(fu_order);
    int fv_order_h = H2D_GET_H_ORDER(fv_order);
    int fv_order_v = H2D_GET_V_ORDER(fv_order);

    // FIXME - this needs more research.
    int order_init = (fu_order_h + fu_order_v) / 2 + (fv_order_h + fv_order_v) / 2;

    // Calculate initial value of the form on coarse element.
    scalar result_init = eval_form_subelement(order_init, mfv, u_ext, fu, fv, ru, rv);

    // Calculate the value of the form using adaptive quadrature.    
    result = eval_form_adaptive(order_init, result_init,
                                mfv, u_ext, fu, fv, ru, rv);
  }

  return result;
}
void DiscreteProblem::eval_form(WeakForm::MultiComponentMatrixFormVol *mfv, 
                                  Hermes::vector<Solution *> u_ext,
                                  PrecalcShapeset *fu, PrecalcShapeset *fv, 
                                  RefMap *ru, RefMap *rv, Hermes::vector<scalar>& result)
{
  _F_
  // Determine the integration order by parsing the form.
  int order = calc_order_matrix_form_vol(mfv, u_ext, fu, fv, ru, rv);
  
  // Evaluate the form using numerical quadrature of order "order".
  Quad2D* quad = fu->get_quad_2d();
  double3* pt = quad->get_points(order);
  int np = quad->get_num_points(order);

  // Init geometry and jacobian*weights.
  if (cache_e[order] == NULL)
  {
    cache_e[order] = init_geom_vol(ru, order);
    double* jac = NULL;
    if(!ru->is_jacobian_const()) 
      jac = ru->get_jacobian(order);
    cache_jwt[order] = new double[np];
    for(int i = 0; i < np; i++) {
      if(ru->is_jacobian_const())
        cache_jwt[order][i] = pt[i][2] * ru->get_const_jacobian();
      else
        cache_jwt[order][i] = pt[i][2] * jac[i];
    }
  }
  Geom<double>* e = cache_e[order];
  double* jwt = cache_jwt[order];

  // Values of the previous Newton iteration, shape functions 
  // and external functions in quadrature points.
  int prev_size = u_ext.size() - mfv->u_ext_offset;
  Func<scalar>** prev = new Func<scalar>*[prev_size];
  if (u_ext != Hermes::vector<Solution *>())
    for (int i = 0; i < prev_size; i++)
      if (u_ext[i + mfv->u_ext_offset] != NULL) 
        prev[i] = init_fn(u_ext[i + mfv->u_ext_offset], order);
      else 
        prev[i] = NULL;
  else
    for (int i = 0; i < prev_size; i++) 
      prev[i] = NULL;

  Func<double>* u = get_fn(fu, ru, order);
  Func<double>* v = get_fn(fv, rv, order);

  ExtData<scalar>* ext = init_ext_fns(mfv->ext, rv, order);

  // The actual calculation takes place here.
  mfv->value(np, jwt, prev, u, v, e, ext, result);

  for(unsigned int i = 0; i < result.size(); i++)
    result[i] *= mfv->scaling_factor;

  // Clean up.
  for(int i = 0; i < prev_size; i++)
    if (prev[i] != NULL) { 
      prev[i]->free_fn(); 
      delete prev[i]; 
    }
  delete [] prev;

  if (ext != NULL) {
    ext->free(); 
    delete ext;
  }
}

int DiscreteProblem::calc_order_matrix_form_vol(WeakForm::MatrixFormVol *mfv, Hermes::vector<Solution *> u_ext,
                                                PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv)
{
  _F_
  // Order that will be returned.
  int order;

  if(is_fvm) 
    order = ru->get_inv_ref_order();
  else {
    int u_ext_length = u_ext.size();      // Number of external solutions.
    int u_ext_offset = mfv->u_ext_offset; // External solutions will start with u_ext[u_ext_offset]
                                          // and there will be only u_ext_length - u_ext_offset of them.
  
    // Increase for multi-valued shape functions.
    int inc = (fu->get_num_components() == 2) ? 1 : 0;

    // Order of solutions from the previous Newton iteration.
    Func<Ord>** oi = new Func<Ord>*[u_ext_length - u_ext_offset];
    if (u_ext != Hermes::vector<Solution *>())
      for(int i = 0; i < u_ext_length - u_ext_offset; i++)
        if (u_ext[i + u_ext_offset] != NULL)
          oi[i] = get_fn_ord(u_ext[i + u_ext_offset]->get_fn_order() + inc);
        else
          oi[i] = get_fn_ord(0);
    else
      for(int i = 0; i < u_ext_length - u_ext_offset; i++)
        oi[i] = get_fn_ord(0);

    // Order of shape functions.
    Func<Ord>* ou = get_fn_ord(fu->get_fn_order() + inc);
    Func<Ord>* ov = get_fn_ord(fv->get_fn_order() + inc);

    // Order of additional external functions.
    ExtData<Ord>* fake_ext = init_ext_fns_ord(mfv->ext);

    // Order of geometric attributes (eg. for multiplication of a solution with coordinates, normals, etc.).
    double fake_wt = 1.0;

    // Total order of the matrix form.
    Ord o = mfv->ord(1, &fake_wt, oi, ou, ov, &geom_ord, fake_ext);

    // Increase due to reference map.
    order = ru->get_inv_ref_order();
    order += o.get_order();
    limit_order(order);
    
    // Cleanup.
    delete [] oi;

    if (fake_ext != NULL) {
      fake_ext->free_ord(); 
      delete fake_ext;
    }
  }
  return order;
}
int DiscreteProblem::calc_order_matrix_form_vol(WeakForm::MultiComponentMatrixFormVol *mfv, Hermes::vector<Solution *> u_ext,
                                                PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv)
{
  _F_
  // Order that will be returned.
  int order;

  if(is_fvm) 
    order = ru->get_inv_ref_order();
  else {
    int u_ext_length = u_ext.size();      // Number of external solutions.
    int u_ext_offset = mfv->u_ext_offset; // External solutions will start with u_ext[u_ext_offset]
                                          // and there will be only u_ext_length - u_ext_offset of them.
  
    // Increase for multi-valued shape functions.
    int inc = (fu->get_num_components() == 2) ? 1 : 0;

    // Order of solutions from the previous Newton iteration.
    Func<Ord>** oi = new Func<Ord>*[u_ext_length - u_ext_offset];
    if (u_ext != Hermes::vector<Solution *>())
      for(int i = 0; i < u_ext_length - u_ext_offset; i++)
        if (u_ext[i + u_ext_offset] != NULL)
          oi[i] = get_fn_ord(u_ext[i + u_ext_offset]->get_fn_order() + inc);
        else
          oi[i] = get_fn_ord(0);
    else
      for(int i = 0; i < u_ext_length - u_ext_offset; i++)
        oi[i] = get_fn_ord(0);

    // Order of shape functions.
    Func<Ord>* ou = get_fn_ord(fu->get_fn_order() + inc);
    Func<Ord>* ov = get_fn_ord(fv->get_fn_order() + inc);

    // Order of additional external functions.
    ExtData<Ord>* fake_ext = init_ext_fns_ord(mfv->ext);

    // Order of geometric attributes (eg. for multiplication of a solution with coordinates, normals, etc.).
    double fake_wt = 1.0;

    // Total order of the matrix form.
    Ord o = mfv->ord(1, &fake_wt, oi, ou, ov, &geom_ord, fake_ext);

    // Increase due to reference map.
    order = ru->get_inv_ref_order();
    order += o.get_order();
    limit_order(order);
    
    // Cleanup.
    delete [] oi;

    if (fake_ext != NULL) {
      fake_ext->free_ord(); 
      delete fake_ext;
    }
  }
  return order;
}

scalar DiscreteProblem::eval_form_subelement(int order, WeakForm::MatrixFormVol *mfv, 
                                             Hermes::vector<Solution *> u_ext,
                                             PrecalcShapeset *fu, PrecalcShapeset *fv, 
                                             RefMap *ru, RefMap *rv)
{
  // Evaluate the form using numerical quadrature of order "order".
  Quad2D* quad = fu->get_quad_2d();
  double3* pt = quad->get_points(order);
  int np = quad->get_num_points(order);

  // Init geometry and jacobian*weights.
  if (cache_e[order] == NULL)
  {
    cache_e[order] = init_geom_vol(ru, order);
    double* jac = NULL;
    if(!ru->is_jacobian_const()) 
      jac = ru->get_jacobian(order);
    cache_jwt[order] = new double[np];
    for(int i = 0; i < np; i++) {
      if(ru->is_jacobian_const())
        cache_jwt[order][i] = pt[i][2] * ru->get_const_jacobian();
      else
        cache_jwt[order][i] = pt[i][2] * jac[i];
    }
  }
  Geom<double>* e = cache_e[order];
  double* jwt = cache_jwt[order];

  // Values of the previous Newton iteration, shape functions 
  // and external functions in quadrature points.
  int prev_size = u_ext.size() - mfv->u_ext_offset;
  Func<scalar>** prev = new Func<scalar>*[prev_size];
  if (u_ext != Hermes::vector<Solution *>())
    for (int i = 0; i < prev_size; i++)
      if (u_ext[i + mfv->u_ext_offset] != NULL) 
        prev[i] = init_fn(u_ext[i + mfv->u_ext_offset], order);
      else 
        prev[i] = NULL;
  else
    for (int i = 0; i < prev_size; i++) 
      prev[i] = NULL;

  Func<double>* u = get_fn(fu, ru, order);
  Func<double>* v = get_fn(fv, rv, order);

  ExtData<scalar>* ext = init_ext_fns(mfv->ext, rv, order);

  // The actual calculation takes place here.
  scalar res = mfv->value(np, jwt, prev, u, v, e, ext) * mfv->scaling_factor;

  // Clean up.
  for(int i = 0; i < prev_size; i++)
    if (prev[i] != NULL) { 
      prev[i]->free_fn(); 
      delete prev[i]; 
    }
  delete [] prev;

  if (ext != NULL) {
    ext->free(); 
    delete ext;
  }

  return res;
}

scalar DiscreteProblem::eval_form_adaptive(int order_init, scalar result_init,
                                           WeakForm::MatrixFormVol *mfv, 
                                           Hermes::vector<Solution *> u_ext,
                                           PrecalcShapeset *fu, PrecalcShapeset *fv, 
                                           RefMap *ru, RefMap *rv)
{
  // Initialize set of all transformable entities.
  std::set<Transformable *> transformable_entities;
  transformable_entities.insert(fu);
  transformable_entities.insert(fv);
  transformable_entities.insert(ru);
  transformable_entities.insert(rv);
  transformable_entities.insert(mfv->ext.begin(), mfv->ext.end());
  transformable_entities.insert(u_ext.begin(), u_ext.end());

  scalar result = 0;
  
  int order_increase = mfv->adapt_order_increase;

  scalar subs_value[4];

#if !defined(H2D_COMPLEX)
  scalar result_current_subelements = 0;
#else
  scalar result_current_subelements = std::complex<double>(0, 0);
#endif
  
  // Clear of the geometry cache for the current active subelement. This needs to be done before AND after
  // as the further division to subelements must be only local.
  this->delete_single_geom_cache(order_init + order_increase);
  for(unsigned int sons_i = 0; sons_i < 4; sons_i++) {
    // Push the transformation to the current son to all functions and reference mappings involved.
    Transformable::push_transforms(transformable_entities, sons_i);

    // The actual calculation.
    subs_value[sons_i] = eval_form_subelement(order_init + order_increase, mfv, 
                                              u_ext, fu, fv, ru, rv);
    
    // Clear of the geometry cache for the current active subelement.
    this->delete_single_geom_cache(order_init + order_increase);

    result_current_subelements += subs_value[sons_i];

    // Reset the transformation.
    Transformable::pop_transforms(transformable_entities);
  }

  // Tolerance checking.
  // First, if the result is negligible.
  if (std::abs(result_current_subelements) < 1e-6)
    return result_current_subelements;

  // Relative error.
  double rel_error = std::abs(result_current_subelements - result_init) / std::abs(result_current_subelements);

  if (rel_error < mfv->adapt_rel_error_tol)
    // Relative error in tolerance.
    return result_current_subelements;
  else {
#if !defined(H2D_COMPLEX)
    scalar result_recursion = 0;
#else
    scalar result_recursion = std::complex<double>(0, 0);
#endif
    // Relative error exceeds allowed tolerance: Call the function 
    // eval_form_adaptive() in each subelement, with initial values
    // subs_value[sons_i].
    
    // Call the function eval_form_adaptive() recursively in all sons.
    for(unsigned int sons_i = 0; sons_i < 4; sons_i++) {
      // Push the transformation to the current son to all functions and reference mappings involved.
      Transformable::push_transforms(transformable_entities, sons_i);
    
      // Recursion.
      subs_value[sons_i] = eval_form_adaptive(order_init + order_increase, subs_value[sons_i],
                                          mfv, u_ext, fu, fv, ru, rv);

      result_recursion += subs_value[sons_i];

      // Reset the transformation.
      Transformable::pop_transforms(transformable_entities);
    }

    return result_recursion;
  }
}


// Volume vector forms.

scalar DiscreteProblem::eval_form(WeakForm::VectorFormVol *vfv, 
                                  Hermes::vector<Solution *> u_ext, 
                                  PrecalcShapeset *fv, RefMap *rv)
{
  _F_
  scalar result = 0;

  if (vfv->adapt_eval == false) {
    // Determine the integration order by parsing the form.
    int order = calc_order_vector_form_vol(vfv, u_ext, fv, rv);
    // Perform non-adaptive numerical quadrature of order "order".
    result = eval_form_subelement(order, vfv, u_ext, fv, rv);
  }
  else {
    // Perform adaptive numerical quadrature starting with order = 2.
    // TODO: The choice of initial order matters a lot for efficiency,
    //       this needs more research.
    Shapeset* fv_shapeset = fv->get_shapeset();
    int fv_index = fv->get_active_shape();
    int fv_order = fv_shapeset->get_order(fv_index);
    int fv_order_h = H2D_GET_H_ORDER(fv_order);
    int fv_order_v = H2D_GET_V_ORDER(fv_order);

    // FIXME - this needs more research.
    int order_init = (fv_order_h + fv_order_v) / 2;

    // Calculate initial value of the form on coarse element.
    scalar result_init = eval_form_subelement(order_init, vfv, u_ext, fv, rv);

    // Calculate the value of the form using adaptive quadrature.    
    result = eval_form_adaptive(order_init, result_init,
                                vfv, u_ext, fv, rv);
  }

  return result;
}
void DiscreteProblem::eval_form(WeakForm::MultiComponentVectorFormVol *vfv, 
                                  Hermes::vector<Solution *> u_ext, 
                                  PrecalcShapeset *fv, RefMap *rv, Hermes::vector<scalar>& result)
{
  _F_
  // Determine the integration order by parsing the form.
  int order = calc_order_vector_form_vol(vfv, u_ext, fv, rv);
  // Evaluate the form using numerical quadrature of order "order".
  Quad2D* quad = fv->get_quad_2d();
  double3* pt = quad->get_points(order);
  int np = quad->get_num_points(order);

  // Init geometry and jacobian*weights.
  if (cache_e[order] == NULL)
  {
    cache_e[order] = init_geom_vol(rv, order);
    double* jac = NULL;
    if(!rv->is_jacobian_const()) 
      jac = rv->get_jacobian(order);
    cache_jwt[order] = new double[np];
    for(int i = 0; i < np; i++) {
      if(rv->is_jacobian_const())
        cache_jwt[order][i] = pt[i][2] * rv->get_const_jacobian();
      else
        cache_jwt[order][i] = pt[i][2] * jac[i];
    }
  }
  Geom<double>* e = cache_e[order];
  double* jwt = cache_jwt[order];

  // Values of the previous Newton iteration, shape functions and external functions in quadrature points.
  int prev_size = u_ext.size() - vfv->u_ext_offset;
  Func<scalar>** prev = new Func<scalar>*[prev_size];
  if (u_ext != Hermes::vector<Solution *>())
    for (int i = 0; i < prev_size; i++)
      if (u_ext[i + vfv->u_ext_offset] != NULL) 
        prev[i] = init_fn(u_ext[i + vfv->u_ext_offset], order);
      else 
        prev[i] = NULL;
  else
    for (int i = 0; i < prev_size; i++) 
      prev[i] = NULL;

  Func<double>* v = get_fn(fv, rv, order);
  ExtData<scalar>* ext = init_ext_fns(vfv->ext, rv, order);

  // The actual calculation takes place here.
  vfv->value(np, jwt, prev, v, e, ext, result);

  for(unsigned int i = 0; i < result.size(); i++)
    result[i] *= vfv->scaling_factor;

  // Clean up.
  for(int i = 0; i < prev_size; i++)
    if (prev[i] != NULL) { 
      prev[i]->free_fn(); 
      delete prev[i]; 
    }
  delete [] prev;

  if (ext != NULL) {
    ext->free(); 
    delete ext;
  }
}

int DiscreteProblem::calc_order_vector_form_vol(WeakForm::VectorFormVol *vfv, 
                                                Hermes::vector<Solution *> u_ext,
                                                PrecalcShapeset *fv, RefMap *rv)
{
  _F_
  // Order that will be returned.
  int order;

  if(is_fvm) 
    order = rv->get_inv_ref_order();
  else {
    int u_ext_length = u_ext.size();      // Number of external solutions.
    int u_ext_offset = vfv->u_ext_offset; // External solutions will start with u_ext[u_ext_offset]
                                          // and there will be only u_ext_length - u_ext_offset of them.
  
    // Increase for multi-valued shape functions.
    int inc = (fv->get_num_components() == 2) ? 1 : 0;

    // Order of solutions from the previous Newton iteration.
    Func<Ord>** oi = new Func<Ord>*[u_ext_length - u_ext_offset];
    if (u_ext != Hermes::vector<Solution *>())
      for(int i = 0; i < u_ext_length - u_ext_offset; i++)
        if (u_ext[i + u_ext_offset] != NULL)
          oi[i] = get_fn_ord(u_ext[i + u_ext_offset]->get_fn_order() + inc);
        else
          oi[i] = get_fn_ord(0);
    else
      for(int i = 0; i < u_ext_length - u_ext_offset; i++)
        oi[i] = get_fn_ord(0);

    // Order of the shape function.
    Func<Ord>* ov = get_fn_ord(fv->get_fn_order() + inc);

    // Order of additional external functions.
    ExtData<Ord>* fake_ext = init_ext_fns_ord(vfv->ext);

    // Order of geometric attributes (eg. for multiplication of 
    // a solution with coordinates, normals, etc.).
    double fake_wt = 1.0;

    // Total order of the vector form.
    Ord o = vfv->ord(1, &fake_wt, oi, ov, &geom_ord, fake_ext);

    // Increase due to reference map.
    order = rv->get_inv_ref_order();
    order += o.get_order();
    limit_order(order);
    
    // Cleanup.
    delete [] oi;

    if (fake_ext != NULL) {
      fake_ext->free_ord(); 
      delete fake_ext;
    }
  }
  return order;
}
int DiscreteProblem::calc_order_vector_form_vol(WeakForm::MultiComponentVectorFormVol *vfv, 
                                                Hermes::vector<Solution *> u_ext,
                                                PrecalcShapeset *fv, RefMap *rv)
{
  _F_
  // Order that will be returned.
  int order;

  if(is_fvm) 
    order = rv->get_inv_ref_order();
  else {
    int u_ext_length = u_ext.size();      // Number of external solutions.
    int u_ext_offset = vfv->u_ext_offset; // External solutions will start with u_ext[u_ext_offset]
                                          // and there will be only u_ext_length - u_ext_offset of them.
  
    // Increase for multi-valued shape functions.
    int inc = (fv->get_num_components() == 2) ? 1 : 0;

    // Order of solutions from the previous Newton iteration.
    Func<Ord>** oi = new Func<Ord>*[u_ext_length - u_ext_offset];
    if (u_ext != Hermes::vector<Solution *>())
      for(int i = 0; i < u_ext_length - u_ext_offset; i++)
        if (u_ext[i + u_ext_offset] != NULL)
          oi[i] = get_fn_ord(u_ext[i + u_ext_offset]->get_fn_order() + inc);
        else
          oi[i] = get_fn_ord(0);
    else
      for(int i = 0; i < u_ext_length - u_ext_offset; i++)
        oi[i] = get_fn_ord(0);

    // Order of the shape function.
    Func<Ord>* ov = get_fn_ord(fv->get_fn_order() + inc);

    // Order of additional external functions.
    ExtData<Ord>* fake_ext = init_ext_fns_ord(vfv->ext);

    // Order of geometric attributes (eg. for multiplication of 
    // a solution with coordinates, normals, etc.).
    double fake_wt = 1.0;

    // Total order of the vector form.
    Ord o = vfv->ord(1, &fake_wt, oi, ov, &geom_ord, fake_ext);

    // Increase due to reference map.
    order = rv->get_inv_ref_order();
    order += o.get_order();
    limit_order(order);
    
    // Cleanup.
    delete [] oi;

    if (fake_ext != NULL) {
      fake_ext->free_ord(); 
      delete fake_ext;
    }
  }
  return order;
}

scalar DiscreteProblem::eval_form_subelement(int order, WeakForm::VectorFormVol *vfv, 
                                  Hermes::vector<Solution *> u_ext, 
                                  PrecalcShapeset *fv, RefMap *rv)
{
  _F_
  // Evaluate the form using numerical quadrature of order "order".
  Quad2D* quad = fv->get_quad_2d();
  double3* pt = quad->get_points(order);
  int np = quad->get_num_points(order);

  // Init geometry and jacobian*weights.
  if (cache_e[order] == NULL)
  {
    cache_e[order] = init_geom_vol(rv, order);
    double* jac = NULL;
    if(!rv->is_jacobian_const()) 
      jac = rv->get_jacobian(order);
    cache_jwt[order] = new double[np];
    for(int i = 0; i < np; i++) {
      if(rv->is_jacobian_const())
        cache_jwt[order][i] = pt[i][2] * rv->get_const_jacobian();
      else
        cache_jwt[order][i] = pt[i][2] * jac[i];
    }
  }
  Geom<double>* e = cache_e[order];
  double* jwt = cache_jwt[order];

  // Values of the previous Newton iteration, shape functions and external functions in quadrature points.
  int prev_size = u_ext.size() - vfv->u_ext_offset;
  Func<scalar>** prev = new Func<scalar>*[prev_size];
  if (u_ext != Hermes::vector<Solution *>())
    for (int i = 0; i < prev_size; i++)
      if (u_ext[i + vfv->u_ext_offset] != NULL) 
        prev[i] = init_fn(u_ext[i + vfv->u_ext_offset], order);
      else 
        prev[i] = NULL;
  else
    for (int i = 0; i < prev_size; i++) 
      prev[i] = NULL;

  Func<double>* v = get_fn(fv, rv, order);
  ExtData<scalar>* ext = init_ext_fns(vfv->ext, rv, order);

  // The actual calculation takes place here.
  scalar res = vfv->value(np, jwt, prev, v, e, ext) * vfv->scaling_factor;

  // Clean up.
  for(int i = 0; i < prev_size; i++)
    if (prev[i] != NULL) { 
      prev[i]->free_fn(); 
      delete prev[i]; 
    }
  delete [] prev;

  if (ext != NULL) {
    ext->free(); 
    delete ext;
  }

  return res;
}

scalar DiscreteProblem::eval_form_adaptive(int order_init, scalar result_init,
                                           WeakForm::VectorFormVol *vfv, Hermes::vector<Solution *> u_ext,
                                           PrecalcShapeset *fv, RefMap *rv)
{
  // Initialize set of all transformable entities.
  std::set<Transformable *> transformable_entities;
  transformable_entities.insert(fv);
  transformable_entities.insert(rv);
  transformable_entities.insert(vfv->ext.begin(), vfv->ext.end());
  transformable_entities.insert(u_ext.begin(), u_ext.end());

  scalar result = 0;
  
  int order_increase = vfv->adapt_order_increase;

  scalar subs_value[4];

#if !defined(H2D_COMPLEX)
  scalar result_current_subelements = 0;
#else
  scalar result_current_subelements = std::complex<double>(0, 0);
#endif
  
  // Clear of the geometry cache for the current active subelement. This needs to be done before AND after
  // as the further division to subelements must be only local.
  this->delete_single_geom_cache(order_init + order_increase);
  for(unsigned int sons_i = 0; sons_i < 4; sons_i++) {
    // Push the transformation to the current son to all functions and reference mappings involved.
    Transformable::push_transforms(transformable_entities, sons_i);

    // The actual calculation.
    subs_value[sons_i] = eval_form_subelement(order_init + order_increase, vfv, 
                                              u_ext, fv, rv);

    // Clear of the geometry cache for the current active subelement.
    this->delete_single_geom_cache(order_init + order_increase);

    result_current_subelements += subs_value[sons_i];

    // Reset the transformation.
    Transformable::pop_transforms(transformable_entities);
  }

  // Tolerance checking.
  // First, if the result is negligible.
  if (std::abs(result_current_subelements) < 1e-6)
    return result_current_subelements;

  // Relative error.
  double rel_error = std::abs(result_current_subelements - result_init) / std::abs(result_current_subelements);

  if (rel_error < vfv->adapt_rel_error_tol)
    // Relative error in tolerance.
    return result_current_subelements;
  else {
#if !defined(H2D_COMPLEX)
    scalar result_recursion = 0;
#else
    scalar result_recursion = std::complex<double>(0, 0);
#endif
    // Relative error exceeds allowed tolerance: Call the function 
    // eval_form_adaptive() in each subelement, with initial values
    // subs_value[sons_i].
    
    // Call the function eval_form_adaptive() recursively in all sons.
    for(unsigned int sons_i = 0; sons_i < 4; sons_i++) {
      // Push the transformation to the current son to all functions and reference mappings involved.
      Transformable::push_transforms(transformable_entities, sons_i);
    
      // Recursion.
      subs_value[sons_i] = eval_form_adaptive(order_init + order_increase, subs_value[sons_i],
                                          vfv, u_ext, fv, rv);

      result_recursion += subs_value[sons_i];

      // Reset the transformation.
      Transformable::pop_transforms(transformable_entities);
    }

    return result_recursion;
  }
}


// Surface matrix forms.

scalar DiscreteProblem::eval_form(WeakForm::MatrixFormSurf *mfs, 
                                  Hermes::vector<Solution *> u_ext, 
                                  PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv, SurfPos* surf_pos)
{
  _F_
  scalar result = 0;

  if (mfs->adapt_eval == false) {
    // Determine the integration order by parsing the form.
    int order = calc_order_matrix_form_surf(mfs, u_ext, fu, fv, ru, rv, surf_pos);
    // Perform non-adaptive numerical quadrature of order "order".
    result = eval_form_subelement(order, mfs, u_ext, fu, fv, ru, rv, surf_pos);
  }
  else {
    // Perform adaptive numerical quadrature starting with order = 2.
    // TODO: The choice of initial order matters a lot for efficiency,
    //       this needs more research.
    int fu_order = fu->get_edge_fn_order(surf_pos->surf_num);
    int fv_order = fv->get_edge_fn_order(surf_pos->surf_num);

    // FIXME - this needs more research.
    int order_init = fu_order + fv_order;

    // Calculate initial value of the form on coarse element.
    scalar result_init = eval_form_subelement(order_init, mfs, u_ext, fu, fv, ru, rv, surf_pos);

    // Calculate the value of the form using adaptive quadrature.    
    result = eval_form_adaptive(order_init, result_init,
                                mfs, u_ext, fu, fv, ru, rv, surf_pos);
  }

  return result;
}
void DiscreteProblem::eval_form(WeakForm::MultiComponentMatrixFormSurf *mfs, 
                                  Hermes::vector<Solution *> u_ext, 
                                  PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv, SurfPos* surf_pos, Hermes::vector<scalar>& result)
{
  _F_
  // Determine the integration order by parsing the form.
  int order = calc_order_matrix_form_surf(mfs, u_ext, fu, fv, ru, rv, surf_pos);
  // Evaluate the form using numerical quadrature of order "order".
  Quad2D* quad = fu->get_quad_2d();

  int eo = quad->get_edge_points(surf_pos->surf_num, order);
  double3* pt = quad->get_points(eo);
  int np = quad->get_num_points(eo);

  // Init geometry and jacobian*weights.
  if (cache_e[eo] == NULL)
  {
    cache_e[eo] = init_geom_surf(ru, surf_pos, eo);
    double3* tan = ru->get_tangent(surf_pos->surf_num, eo);
    cache_jwt[eo] = new double[np];
    for(int i = 0; i < np; i++)
      cache_jwt[eo][i] = pt[i][2] * tan[i][2];
  }
  Geom<double>* e = cache_e[eo];
  double* jwt = cache_jwt[eo];

  // Values of the previous Newton iteration, shape functions and external functions in quadrature points.
  int prev_size = u_ext.size() - mfs->u_ext_offset;
  Func<scalar>** prev = new Func<scalar>*[prev_size];
  if (u_ext != Hermes::vector<Solution *>())
    for (int i = 0; i < prev_size; i++)
      if (u_ext[i + mfs->u_ext_offset] != NULL) 
        prev[i] = init_fn(u_ext[i + mfs->u_ext_offset], eo);
      else 
        prev[i] = NULL;
  else
    for (int i = 0; i < prev_size; i++) 
      prev[i] = NULL;

  Func<double>* u = get_fn(fu, ru, eo);
  Func<double>* v = get_fn(fv, rv, eo);
  ExtData<scalar>* ext = init_ext_fns(mfs->ext, rv, eo);

  // The actual calculation takes place here.
  mfs->value(np, jwt, prev, u, v, e, ext, result);

  for(unsigned int i = 0; i < result.size(); i++)
    result[i] *= mfs->scaling_factor * 0.5;

  // Clean up.
  for(int i = 0; i < prev_size; i++)
    if (prev[i] != NULL) { 
      prev[i]->free_fn(); 
      delete prev[i]; 
    }
  delete [] prev;

  if (ext != NULL) {
    ext->free(); 
    delete ext;
  }
}

int DiscreteProblem::calc_order_matrix_form_surf(WeakForm::MatrixFormSurf *mfs, Hermes::vector<Solution *> u_ext,
                                  PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv, SurfPos* surf_pos)
{
  _F_
  // Order that will be returned.
  int order;

  if(is_fvm)
    order = ru->get_inv_ref_order();
  else {
    int u_ext_length = u_ext.size();      // Number of external solutions.
    int u_ext_offset = mfs->u_ext_offset; // External solutions will start with u_ext[u_ext_offset]
                                          // and there will be only u_ext_length - u_ext_offset of them.
  
    // Increase for multi-valued shape functions.
    int inc = (fu->get_num_components() == 2) ? 1 : 0;

    // Order of solutions from the previous Newton iteration.
    Func<Ord>** oi = new Func<Ord>*[u_ext_length - u_ext_offset];
    if (u_ext != Hermes::vector<Solution *>())
      for(int i = 0; i < u_ext_length - u_ext_offset; i++)
        if (u_ext[i + u_ext_offset] != NULL)
          oi[i] = get_fn_ord(u_ext[i + u_ext_offset]->get_edge_fn_order(surf_pos->surf_num) + inc);
        else
          oi[i] = get_fn_ord(0);
    else
      for(int i = 0; i < u_ext_length - u_ext_offset; i++)
        oi[i] = get_fn_ord(0);

    // Order of shape functions.
    Func<Ord>* ou = get_fn_ord(fu->get_edge_fn_order(surf_pos->surf_num) + inc);
    Func<Ord>* ov = get_fn_ord(fv->get_edge_fn_order(surf_pos->surf_num) + inc);

    // Order of additional external functions.
    ExtData<Ord>* fake_ext = init_ext_fns_ord(mfs->ext, surf_pos->surf_num);

    // Order of geometric attributes (eg. for multiplication of a solution with coordinates, normals, etc.).
    double fake_wt = 1.0;

    // Total order of the matrix form.
    Ord o = mfs->ord(1, &fake_wt, oi, ou, ov, &geom_ord, fake_ext);

    // Increase due to reference map.
    order = ru->get_inv_ref_order();
    order += o.get_order();
    limit_order(order);

    // Cleanup.
    delete [] oi;

    if (fake_ext != NULL) {
      fake_ext->free_ord(); 
      delete fake_ext;
    }
  }
  return order;
}
int DiscreteProblem::calc_order_matrix_form_surf(WeakForm::MultiComponentMatrixFormSurf *mfs, Hermes::vector<Solution *> u_ext,
                                  PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv, SurfPos* surf_pos)
{
  _F_
  // Order that will be returned.
  int order;

  if(is_fvm)
    order = ru->get_inv_ref_order();
  else {
    int u_ext_length = u_ext.size();      // Number of external solutions.
    int u_ext_offset = mfs->u_ext_offset; // External solutions will start with u_ext[u_ext_offset]
                                          // and there will be only u_ext_length - u_ext_offset of them.
  
    // Increase for multi-valued shape functions.
    int inc = (fu->get_num_components() == 2) ? 1 : 0;

    // Order of solutions from the previous Newton iteration.
    Func<Ord>** oi = new Func<Ord>*[u_ext_length - u_ext_offset];
    if (u_ext != Hermes::vector<Solution *>())
      for(int i = 0; i < u_ext_length - u_ext_offset; i++)
        if (u_ext[i + u_ext_offset] != NULL)
          oi[i] = get_fn_ord(u_ext[i + u_ext_offset]->get_edge_fn_order(surf_pos->surf_num) + inc);
        else
          oi[i] = get_fn_ord(0);
    else
      for(int i = 0; i < u_ext_length - u_ext_offset; i++)
        oi[i] = get_fn_ord(0);

    // Order of shape functions.
    Func<Ord>* ou = get_fn_ord(fu->get_edge_fn_order(surf_pos->surf_num) + inc);
    Func<Ord>* ov = get_fn_ord(fv->get_edge_fn_order(surf_pos->surf_num) + inc);

    // Order of additional external functions.
    ExtData<Ord>* fake_ext = init_ext_fns_ord(mfs->ext, surf_pos->surf_num);

    // Order of geometric attributes (eg. for multiplication of a solution with coordinates, normals, etc.).
    double fake_wt = 1.0;

    // Total order of the matrix form.
    Ord o = mfs->ord(1, &fake_wt, oi, ou, ov, &geom_ord, fake_ext);

    // Increase due to reference map.
    order = ru->get_inv_ref_order();
    order += o.get_order();
    limit_order(order);

    // Cleanup.
    delete [] oi;

    if (fake_ext != NULL) {
      fake_ext->free_ord(); 
      delete fake_ext;
    }
  }
  return order;
}

scalar DiscreteProblem::eval_form_subelement(int order, WeakForm::MatrixFormSurf *mfs, Hermes::vector<Solution *> u_ext,
                        PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv, SurfPos* surf_pos)
{
  _F_
  // Evaluate the form using numerical quadrature of order "order".
  Quad2D* quad = fu->get_quad_2d();

  int eo = quad->get_edge_points(surf_pos->surf_num, order);
  double3* pt = quad->get_points(eo);
  int np = quad->get_num_points(eo);

  // Init geometry and jacobian*weights.
  if (cache_e[eo] == NULL)
  {
    cache_e[eo] = init_geom_surf(ru, surf_pos, eo);
    double3* tan = ru->get_tangent(surf_pos->surf_num, eo);
    cache_jwt[eo] = new double[np];
    for(int i = 0; i < np; i++)
      cache_jwt[eo][i] = pt[i][2] * tan[i][2];
  }
  Geom<double>* e = cache_e[eo];
  double* jwt = cache_jwt[eo];

  // Values of the previous Newton iteration, shape functions and external functions in quadrature points.
  int prev_size = u_ext.size() - mfs->u_ext_offset;
  Func<scalar>** prev = new Func<scalar>*[prev_size];
  if (u_ext != Hermes::vector<Solution *>())
    for (int i = 0; i < prev_size; i++)
      if (u_ext[i + mfs->u_ext_offset] != NULL) 
        prev[i] = init_fn(u_ext[i + mfs->u_ext_offset], eo);
      else 
        prev[i] = NULL;
  else
    for (int i = 0; i < prev_size; i++) 
      prev[i] = NULL;

  Func<double>* u = get_fn(fu, ru, eo);
  Func<double>* v = get_fn(fv, rv, eo);
  ExtData<scalar>* ext = init_ext_fns(mfs->ext, rv, eo);

  // The actual calculation takes place here.
  scalar res = mfs->value(np, jwt, prev, u, v, e, ext) * mfs->scaling_factor;

  // Clean up.
  for(int i = 0; i < prev_size; i++)
    if (prev[i] != NULL) { 
      prev[i]->free_fn(); 
      delete prev[i]; 
    }
  delete [] prev;

  if (ext != NULL) {
    ext->free(); 
    delete ext;
  }

  return 0.5 * res; // Edges are parameterized from 0 to 1 while integration weights
                    // are defined in (-1, 1). Thus multiplying with 0.5 to correct
                    // the weights.
}

scalar DiscreteProblem::eval_form_adaptive(int order_init, scalar result_init,
                                           WeakForm::MatrixFormSurf *mfs, Hermes::vector<Solution *> u_ext,
                                           PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv, SurfPos* surf_pos)
{
  // Initialize set of all transformable entities.
  std::set<Transformable *> transformable_entities;
  transformable_entities.insert(fu);
  transformable_entities.insert(fv);
  transformable_entities.insert(ru);
  transformable_entities.insert(rv);
  transformable_entities.insert(mfs->ext.begin(), mfs->ext.end());
  transformable_entities.insert(u_ext.begin(), u_ext.end());

  scalar result = 0;
  
  int order_increase = mfs->adapt_order_increase;

  scalar subs_value[4];

#if !defined(H2D_COMPLEX)
  scalar result_current_subelements = 0;
#else
  scalar result_current_subelements = std::complex<double>(0, 0);
#endif
  
  // Clear of the geometry cache for the current active subelement. This needs to be done before AND after
  // as the further division to subelements must be only local.
  this->delete_single_geom_cache(order_init + order_increase);
  for(unsigned int sons_i = 0; sons_i < 4; sons_i++) {
    // Push the transformation to the current son to all functions and reference mappings involved.
    Transformable::push_transforms(transformable_entities, sons_i);

    // The actual calculation.
    subs_value[sons_i] = eval_form_subelement(order_init + order_increase, mfs, 
                                              u_ext, fu, fv, ru, rv, surf_pos);

    // Clear of the geometry cache for the current active subelement.
    this->delete_single_geom_cache(order_init + order_increase);

    result_current_subelements += subs_value[sons_i];

    // Reset the transformation.
    Transformable::pop_transforms(transformable_entities);
  }

  // Tolerance checking.
  // First, if the result is negligible.
  if (std::abs(result_current_subelements) < 1e-6)
    return result_current_subelements;

  // Relative error.
  double rel_error = std::abs(result_current_subelements - result_init) / std::abs(result_current_subelements);

  if (rel_error < mfs->adapt_rel_error_tol)
    // Relative error in tolerance.
    return result_current_subelements;
  else {
#if !defined(H2D_COMPLEX)
    scalar result_recursion = 0;
#else
    scalar result_recursion = std::complex<double>(0, 0);
#endif
    // Relative error exceeds allowed tolerance: Call the function 
    // eval_form_adaptive() in each subelement, with initial values
    // subs_value[sons_i].
    
    // Call the function eval_form_adaptive() recursively in all sons.
    for(unsigned int sons_i = 0; sons_i < 4; sons_i++) {
      // Push the transformation to the current son to all functions and reference mappings involved.
      Transformable::push_transforms(transformable_entities, sons_i);
    
      // Recursion.
      subs_value[sons_i] = eval_form_adaptive(order_init + order_increase, subs_value[sons_i],
                                          mfs, u_ext, fu, fv, ru, rv, surf_pos);

      result_recursion += subs_value[sons_i];

      // Reset the transformation.
      Transformable::pop_transforms(transformable_entities);
    }

    return result_recursion;
  }
}


// Surface vector forms.

scalar DiscreteProblem::eval_form(WeakForm::VectorFormSurf *vfs, 
                                  Hermes::vector<Solution *> u_ext, 
                                  PrecalcShapeset *fv, RefMap *rv, SurfPos* surf_pos)
{
  _F_
  scalar result = 0;

  if (vfs->adapt_eval == false) {
    // Determine the integration order by parsing the form.
    int order = calc_order_vector_form_surf(vfs, u_ext, fv, rv, surf_pos);
    // Perform non-adaptive numerical quadrature of order "order".
    result = eval_form_subelement(order, vfs, u_ext, fv, rv, surf_pos);
  }
  else {
    // Perform adaptive numerical quadrature starting with order = 2.
    // TODO: The choice of initial order matters a lot for efficiency,
    //       this needs more research.
    int fv_order = fv->get_edge_fn_order(surf_pos->surf_num);

    // FIXME - this needs more research.
    int order_init = fv_order;

    // Calculate initial value of the form on coarse element.
    scalar result_init = eval_form_subelement(order_init, vfs, u_ext, fv, rv, surf_pos);

    // Calculate the value of the form using adaptive quadrature.    
    result = eval_form_adaptive(order_init, result_init,
                                vfs, u_ext, fv, rv, surf_pos);
  }

  return result;
}
void DiscreteProblem::eval_form(WeakForm::MultiComponentVectorFormSurf *vfs, 
                                  Hermes::vector<Solution *> u_ext, 
                                  PrecalcShapeset *fv, RefMap *rv, SurfPos* surf_pos, Hermes::vector<scalar>& result)
{
  _F_
  // Determine the integration order by parsing the form.
  int order = calc_order_vector_form_surf(vfs, u_ext, fv, rv, surf_pos);
  // Evaluate the form using numerical quadrature of order "order".
  Quad2D* quad = fv->get_quad_2d();

  int eo = quad->get_edge_points(surf_pos->surf_num, order);
  double3* pt = quad->get_points(eo);
  int np = quad->get_num_points(eo);

  // Init geometry and jacobian*weights.
  if (cache_e[eo] == NULL) {
    cache_e[eo] = init_geom_surf(rv, surf_pos, eo);
    double3* tan = rv->get_tangent(surf_pos->surf_num, eo);
    cache_jwt[eo] = new double[np];
    for(int i = 0; i < np; i++)
      cache_jwt[eo][i] = pt[i][2] * tan[i][2];
  }
  Geom<double>* e = cache_e[eo];
  double* jwt = cache_jwt[eo];

  // Values of the previous Newton iteration, shape functions and external functions in quadrature points.
  int prev_size = u_ext.size() - vfs->u_ext_offset;
  Func<scalar>** prev = new Func<scalar>*[prev_size];
  if (u_ext != Hermes::vector<Solution *>())
    for (int i = 0; i < prev_size; i++)
      if (u_ext[i + vfs->u_ext_offset] != NULL) 
        prev[i] = init_fn(u_ext[i + vfs->u_ext_offset], eo);
      else 
        prev[i] = NULL;
  else
    for (int i = 0; i < prev_size; i++) 
      prev[i] = NULL;

  Func<double>* v = get_fn(fv, rv, eo);
  ExtData<scalar>* ext = init_ext_fns(vfs->ext, rv, eo);

  // The actual calculation takes place here.
  vfs->value(np, jwt, prev, v, e, ext, result);

  for(unsigned int i = 0; i < result.size(); i++)
    result[i] *= vfs->scaling_factor * 0.5;

  // Clean up.
  for(int i = 0; i < prev_size; i++)
    if (prev[i] != NULL) { 
      prev[i]->free_fn(); 
      delete prev[i]; 
    }
  delete [] prev;

  if (ext != NULL) {
    ext->free(); 
    delete ext;
  }
}

int DiscreteProblem::calc_order_vector_form_surf(WeakForm::VectorFormSurf *vfs, Hermes::vector<Solution *> u_ext,
                                  PrecalcShapeset *fv, RefMap *rv, SurfPos* surf_pos)
{
  _F_
  // Order that will be returned.
  int order;

  if(is_fvm) 
    order = rv->get_inv_ref_order();
  else {
    int u_ext_length = u_ext.size();      // Number of external solutions.
    int u_ext_offset = vfs->u_ext_offset; // External solutions will start with u_ext[u_ext_offset]
                                          // and there will be only u_ext_length - u_ext_offset of them.
  
    // Increase for multi-valued shape functions.
    int inc = (fv->get_num_components() == 2) ? 1 : 0;

    // Order of solutions from the previous Newton iteration.
    Func<Ord>** oi = new Func<Ord>*[u_ext_length - u_ext_offset];
    if (u_ext != Hermes::vector<Solution *>())
      for(int i = 0; i < u_ext_length - u_ext_offset; i++)
        if (u_ext[i + u_ext_offset] != NULL)
          oi[i] = get_fn_ord(u_ext[i]->get_edge_fn_order(surf_pos->surf_num) + inc);
        else
          oi[i] = get_fn_ord(0);
    else
      for(int i = 0; i < u_ext_length - u_ext_offset; i++)
        oi[i] = get_fn_ord(0);

    // Order of the shape function.
    Func<Ord>* ov = get_fn_ord(fv->get_edge_fn_order(surf_pos->surf_num) + inc);

    // Order of additional external functions.
    ExtData<Ord>* fake_ext = init_ext_fns_ord(vfs->ext);

    // Order of geometric attributes (eg. for multiplication of a solution with coordinates, normals, etc.).
    double fake_wt = 1.0;

    // Total order of the vector form.
    Ord o = vfs->ord(1, &fake_wt, oi, ov, &geom_ord, fake_ext);

    // Increase due to reference map.
    order = rv->get_inv_ref_order();
    order += o.get_order();
    limit_order(order);
    
    // Cleanup.
    delete [] oi;

    if (fake_ext != NULL) {
      fake_ext->free_ord(); 
      delete fake_ext;
    }
  }
  return order;
}
int DiscreteProblem::calc_order_vector_form_surf(WeakForm::MultiComponentVectorFormSurf *vfs, Hermes::vector<Solution *> u_ext,
                                  PrecalcShapeset *fv, RefMap *rv, SurfPos* surf_pos)
{
  _F_
  // Order that will be returned.
  int order;

  if(is_fvm) 
    order = rv->get_inv_ref_order();
  else {
    int u_ext_length = u_ext.size();      // Number of external solutions.
    int u_ext_offset = vfs->u_ext_offset; // External solutions will start with u_ext[u_ext_offset]
                                          // and there will be only u_ext_length - u_ext_offset of them.
  
    // Increase for multi-valued shape functions.
    int inc = (fv->get_num_components() == 2) ? 1 : 0;

    // Order of solutions from the previous Newton iteration.
    Func<Ord>** oi = new Func<Ord>*[u_ext_length - u_ext_offset];
    if (u_ext != Hermes::vector<Solution *>())
      for(int i = 0; i < u_ext_length - u_ext_offset; i++)
        if (u_ext[i + u_ext_offset] != NULL)
          oi[i] = get_fn_ord(u_ext[i]->get_edge_fn_order(surf_pos->surf_num) + inc);
        else
          oi[i] = get_fn_ord(0);
    else
      for(int i = 0; i < u_ext_length - u_ext_offset; i++)
        oi[i] = get_fn_ord(0);

    // Order of the shape function.
    Func<Ord>* ov = get_fn_ord(fv->get_edge_fn_order(surf_pos->surf_num) + inc);

    // Order of additional external functions.
    ExtData<Ord>* fake_ext = init_ext_fns_ord(vfs->ext);

    // Order of geometric attributes (eg. for multiplication of a solution with coordinates, normals, etc.).
    double fake_wt = 1.0;

    // Total order of the vector form.
    Ord o = vfs->ord(1, &fake_wt, oi, ov, &geom_ord, fake_ext);

    // Increase due to reference map.
    order = rv->get_inv_ref_order();
    order += o.get_order();
    limit_order(order);
    
    // Cleanup.
    delete [] oi;

    if (fake_ext != NULL) {
      fake_ext->free_ord(); 
      delete fake_ext;
    }
  }
  return order;
}

scalar DiscreteProblem::eval_form_subelement(int order, WeakForm::VectorFormSurf *vfs, Hermes::vector<Solution *> u_ext,
                        PrecalcShapeset *fv, RefMap *rv, SurfPos* surf_pos)
{
  _F_
  // Evaluate the form using numerical quadrature of order "order".
  Quad2D* quad = fv->get_quad_2d();

  int eo = quad->get_edge_points(surf_pos->surf_num, order);
  double3* pt = quad->get_points(eo);
  int np = quad->get_num_points(eo);

  // Init geometry and jacobian*weights.
  if (cache_e[eo] == NULL) {
    cache_e[eo] = init_geom_surf(rv, surf_pos, eo);
    double3* tan = rv->get_tangent(surf_pos->surf_num, eo);
    cache_jwt[eo] = new double[np];
    for(int i = 0; i < np; i++)
      cache_jwt[eo][i] = pt[i][2] * tan[i][2];
  }
  Geom<double>* e = cache_e[eo];
  double* jwt = cache_jwt[eo];

  // Values of the previous Newton iteration, shape functions and external functions in quadrature points.
  int prev_size = u_ext.size() - vfs->u_ext_offset;
  Func<scalar>** prev = new Func<scalar>*[prev_size];
  if (u_ext != Hermes::vector<Solution *>())
    for (int i = 0; i < prev_size; i++)
      if (u_ext[i + vfs->u_ext_offset] != NULL) 
        prev[i] = init_fn(u_ext[i + vfs->u_ext_offset], eo);
      else 
        prev[i] = NULL;
  else
    for (int i = 0; i < prev_size; i++) 
      prev[i] = NULL;

  Func<double>* v = get_fn(fv, rv, eo);
  ExtData<scalar>* ext = init_ext_fns(vfs->ext, rv, eo);

  // The actual calculation takes place here.
  scalar res = vfs->value(np, jwt, prev, v, e, ext) * vfs->scaling_factor;

  // Clean up.
  for(int i = 0; i < prev_size; i++)
    if (prev[i] != NULL) { 
      prev[i]->free_fn(); 
      delete prev[i]; 
    }
  delete [] prev;

  if (ext != NULL) {
    ext->free(); 
    delete ext;
  }

  return 0.5 * res; // Edges are parameterized from 0 to 1 while integration weights
                    // are defined in (-1, 1). Thus multiplying with 0.5 to correct
                    // the weights.
}

scalar DiscreteProblem::eval_form_adaptive(int order_init, scalar result_init,
                                           WeakForm::VectorFormSurf *vfs, Hermes::vector<Solution *> u_ext,
                                           PrecalcShapeset *fv, RefMap *rv, SurfPos* surf_pos)
{
  // Initialize set of all transformable entities.
  std::set<Transformable *> transformable_entities;
  transformable_entities.insert(fv);
  transformable_entities.insert(rv);
  transformable_entities.insert(vfs->ext.begin(), vfs->ext.end());
  transformable_entities.insert(u_ext.begin(), u_ext.end());

  scalar result = 0;
  
  int order_increase = vfs->adapt_order_increase;

  scalar subs_value[4];

#if !defined(H2D_COMPLEX)
  scalar result_current_subelements = 0;
#else
  scalar result_current_subelements = std::complex<double>(0, 0);
#endif
  
  // Clear of the geometry cache for the current active subelement. This needs to be done before AND after
  // as the further division to subelements must be only local.
  this->delete_single_geom_cache(order_init + order_increase);
  for(unsigned int sons_i = 0; sons_i < 4; sons_i++) {
    // Push the transformation to the current son to all functions and reference mappings involved.
    Transformable::push_transforms(transformable_entities, sons_i);

    // The actual calculation.
    subs_value[sons_i] = eval_form_subelement(order_init + order_increase, vfs, 
                                              u_ext, fv, rv, surf_pos);

    // Clear of the geometry cache for the current active subelement.
    this->delete_single_geom_cache(order_init + order_increase);

    result_current_subelements += subs_value[sons_i];

    // Reset the transformation.
    Transformable::pop_transforms(transformable_entities);
  }

  // Tolerance checking.
  // First, if the result is negligible.
  if (std::abs(result_current_subelements) < 1e-6)
    return result_current_subelements;

  // Relative error.
  double rel_error = std::abs(result_current_subelements - result_init) / std::abs(result_current_subelements);

  if (rel_error < vfs->adapt_rel_error_tol)
    // Relative error in tolerance.
    return result_current_subelements;
  else {
#if !defined(H2D_COMPLEX)
    scalar result_recursion = 0;
#else
    scalar result_recursion = std::complex<double>(0, 0);
#endif
    // Relative error exceeds allowed tolerance: Call the function 
    // eval_form_adaptive() in each subelement, with initial values
    // subs_value[sons_i].
    
    // Call the function eval_form_adaptive() recursively in all sons.
    for(unsigned int sons_i = 0; sons_i < 4; sons_i++) {
      // Push the transformation to the current son to all functions and reference mappings involved.
      Transformable::push_transforms(transformable_entities, sons_i);
    
      // Recursion.
      subs_value[sons_i] = eval_form_adaptive(order_init + order_increase, subs_value[sons_i],
                                          vfs, u_ext, fv, rv, surf_pos);

      result_recursion += subs_value[sons_i];

      // Reset the transformation.
      Transformable::pop_transforms(transformable_entities);
    }

    return result_recursion;
  }
}


// DG forms.

int DiscreteProblem::calc_order_dg_matrix_form(WeakForm::MatrixFormSurf *mfs, Hermes::vector<Solution *> u_ext,
                                  PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, SurfPos* surf_pos,
                                  bool neighbor_supp_u, bool neighbor_supp_v, LightArray<NeighborSearch*>& neighbor_searches, int neighbor_index_u)
{
  NeighborSearch* nbs_u = neighbor_searches.get(neighbor_index_u);
  // Order that will be returned.
  int order;
  
  if(this->is_fvm)
    order = ru->get_inv_ref_order();
  else {
    // Order of solutions from the previous Newton iteration.
    int prev_size = u_ext.size() - mfs->u_ext_offset;
    Func<Ord>** oi = new Func<Ord>*[prev_size];
    if (u_ext != Hermes::vector<Solution *>())
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
    DiscontinuousFunc<Ord>* ou = new DiscontinuousFunc<Ord>(get_fn_ord(fu->get_edge_fn_order(surf_pos->surf_num) + inc), neighbor_supp_u);
    DiscontinuousFunc<Ord>* ov = new DiscontinuousFunc<Ord>(get_fn_ord(fv->get_edge_fn_order(surf_pos->surf_num) + inc), neighbor_supp_v);

    // Order of additional external functions.
    ExtData<Ord>* fake_ext = init_ext_fns_ord(mfs->ext, neighbor_searches);

    // Order of geometric attributes (eg. for multiplication of a solution with coordinates, normals, etc.).
    Geom<Ord>* fake_e = new InterfaceGeom<Ord>(&geom_ord, nbs_u->neighb_el->marker, 
      nbs_u->neighb_el->id, nbs_u->neighb_el->get_diameter());
    double fake_wt = 1.0;
    
    // Total order of the matrix form.
    Ord o = mfs->ord(1, &fake_wt, oi, ou, ov, fake_e, fake_ext);

    // Increase due to reference maps.
    order = ru->get_inv_ref_order();

    order += o.get_order();
    limit_order(order);

    // Clean up.
    if (u_ext != Hermes::vector<Solution *>())
      for (int i = 0; i < prev_size; i++)
        if (u_ext[i + mfs->u_ext_offset] != NULL) 
          delete oi[i];
    delete [] oi;
    delete fake_e;
    delete ou;
    delete ov;
    if (fake_ext != NULL) {
      for (int i = 0; i < fake_ext->nf; i++) {
        delete fake_ext->fn[i];
      }
      fake_ext->free_ord(); 
      delete fake_ext;
    }
  }
  return order;
}
int DiscreteProblem::calc_order_dg_matrix_form(WeakForm::MultiComponentMatrixFormSurf *mfs, Hermes::vector<Solution *> u_ext,
                                  PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, SurfPos* surf_pos,
                                  bool neighbor_supp_u, bool neighbor_supp_v, LightArray<NeighborSearch*>& neighbor_searches, int neighbor_index_u)
{
  NeighborSearch* nbs_u = neighbor_searches.get(neighbor_index_u);
  // Order that will be returned.
  int order;
  
  if(this->is_fvm)
    order = ru->get_inv_ref_order();
  else {
    // Order of solutions from the previous Newton iteration.
    int prev_size = u_ext.size() - mfs->u_ext_offset;
    Func<Ord>** oi = new Func<Ord>*[prev_size];
    if (u_ext != Hermes::vector<Solution *>())
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
    DiscontinuousFunc<Ord>* ou = new DiscontinuousFunc<Ord>(get_fn_ord(fu->get_edge_fn_order(surf_pos->surf_num) + inc), neighbor_supp_u);
    DiscontinuousFunc<Ord>* ov = new DiscontinuousFunc<Ord>(get_fn_ord(fv->get_edge_fn_order(surf_pos->surf_num) + inc), neighbor_supp_v);

    // Order of additional external functions.
    ExtData<Ord>* fake_ext = init_ext_fns_ord(mfs->ext, neighbor_searches);

    // Order of geometric attributes (eg. for multiplication of a solution with coordinates, normals, etc.).
    Geom<Ord>* fake_e = new InterfaceGeom<Ord>(&geom_ord, nbs_u->neighb_el->marker, 
      nbs_u->neighb_el->id, nbs_u->neighb_el->get_diameter());
    double fake_wt = 1.0;
    
    // Total order of the matrix form.
    Ord o = mfs->ord(1, &fake_wt, oi, ou, ov, fake_e, fake_ext);

    // Increase due to reference maps.
    order = ru->get_inv_ref_order();

    order += o.get_order();
    limit_order(order);

    // Clean up.
    if (u_ext != Hermes::vector<Solution *>())
      for (int i = 0; i < prev_size; i++)
        if (u_ext[i + mfs->u_ext_offset] != NULL) 
          delete oi[i];
    delete [] oi;
    delete fake_e;
    delete ou;
    delete ov;
    if (fake_ext != NULL) {
      for (int i = 0; i < fake_ext->nf; i++) {
        delete fake_ext->fn[i];
      }
      fake_ext->free_ord(); 
      delete fake_ext;
    }
  }
  return order;
}

scalar DiscreteProblem::eval_dg_form(WeakForm::MatrixFormSurf* mfs, Hermes::vector<Solution *> u_ext,
                                     PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru_central, RefMap *ru_actual, RefMap *rv, 
                                     bool neighbor_supp_u, bool neighbor_supp_v,
                                     SurfPos* surf_pos, LightArray<NeighborSearch*>& neighbor_searches, int neighbor_index_u, int neighbor_index_v)
{
  _F_

  NeighborSearch* nbs_u = neighbor_searches.get(neighbor_index_u);
  NeighborSearch* nbs_v = neighbor_searches.get(neighbor_index_v);

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
  if (cache_e[eo] == NULL) {
    cache_e[eo] = init_geom_surf(ru_central, surf_pos, eo);
    double3* tan = ru_central->get_tangent(surf_pos->surf_num, eo);
    cache_jwt[eo] = new double[np];
    for(int i = 0; i < np; i++)
      cache_jwt[eo][i] = pt[i][2] * tan[i][2];
  }

  Geom<double>* e = new InterfaceGeom<double>(cache_e[eo], nbs_u->neighb_el->marker, 
    nbs_u->neighb_el->id, nbs_u->neighb_el->get_diameter());
  double* jwt = cache_jwt[eo];

  // Values of the previous Newton iteration, shape functions and external functions in quadrature points.
  int prev_size = u_ext.size() - mfs->u_ext_offset;
  Func<scalar>** prev = new Func<scalar>*[prev_size];
  if (u_ext != Hermes::vector<Solution *>())
    for (int i = 0; i < prev_size; i++)
      if (u_ext[i + mfs->u_ext_offset] != NULL) {
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
  
  ExtData<scalar>* ext = init_ext_fns(mfs->ext, neighbor_searches, order);

  scalar res = mfs->value(np, jwt, prev, u, v, e, ext);

  // Clean up.
  for (int i = 0; i < prev_size; i++) {
    if (prev[i] != NULL) {
      prev[i]->free_fn();
      delete prev[i];
    }
  }

  delete [] prev;


  if (ext != NULL) {
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
void DiscreteProblem::eval_dg_form(WeakForm::MultiComponentMatrixFormSurf* mfs, Hermes::vector<Solution *> u_ext,
                                     PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru_central, RefMap *ru_actual, RefMap *rv, 
                                     bool neighbor_supp_u, bool neighbor_supp_v,
                                     SurfPos* surf_pos, LightArray<NeighborSearch*>& neighbor_searches, int neighbor_index_u, int neighbor_index_v, Hermes::vector<scalar>& result)
{
  _F_

  NeighborSearch* nbs_u = neighbor_searches.get(neighbor_index_u);
  NeighborSearch* nbs_v = neighbor_searches.get(neighbor_index_v);

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
  if (cache_e[eo] == NULL) {
    cache_e[eo] = init_geom_surf(ru_central, surf_pos, eo);
    double3* tan = ru_central->get_tangent(surf_pos->surf_num, eo);
    cache_jwt[eo] = new double[np];
    for(int i = 0; i < np; i++)
      cache_jwt[eo][i] = pt[i][2] * tan[i][2];
  }

  Geom<double>* e = new InterfaceGeom<double>(cache_e[eo], nbs_u->neighb_el->marker, 
    nbs_u->neighb_el->id, nbs_u->neighb_el->get_diameter());
  double* jwt = cache_jwt[eo];

  // Values of the previous Newton iteration, shape functions and external functions in quadrature points.
  int prev_size = u_ext.size() - mfs->u_ext_offset;
  Func<scalar>** prev = new Func<scalar>*[prev_size];
  if (u_ext != Hermes::vector<Solution *>())
    for (int i = 0; i < prev_size; i++)
      if (u_ext[i + mfs->u_ext_offset] != NULL) {
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
  
  ExtData<scalar>* ext = init_ext_fns(mfs->ext, neighbor_searches, order);

  mfs->value(np, jwt, prev, u, v, e, ext, result);

  for(unsigned int i = 0; i < result.size(); i++)
    result[i] *= mfs->scaling_factor * 0.5;

  // Clean up.
  for (int i = 0; i < prev_size; i++) {
    if (prev[i] != NULL) {
      prev[i]->free_fn();
      delete prev[i];
    }
  }

  delete [] prev;


  if (ext != NULL) {
    ext->free(); 
    delete ext;
  }

  delete u;
  delete v;
  delete e;
}

int DiscreteProblem::calc_order_dg_vector_form(WeakForm::VectorFormSurf *vfs, Hermes::vector<Solution *> u_ext,
                                  PrecalcShapeset *fv, RefMap *rv, SurfPos* surf_pos,
                                  LightArray<NeighborSearch*>& neighbor_searches, int neighbor_index_v)
{
  NeighborSearch* nbs_v = neighbor_searches.get(neighbor_index_v);
  // Order that will be returned.
  int order;
  int u_ext_length = u_ext.size();      // Number of external solutions.
  int u_ext_offset = vfs->u_ext_offset; // External solutions will start with u_ext[u_ext_offset]
                                        // and there will be only u_ext_length - u_ext_offset of them.
  if(this->is_fvm)
    order = rv->get_inv_ref_order();
  else {
    // Order of solutions from the previous Newton iteration.
    int prev_size = u_ext.size() - vfs->u_ext_offset;
    Func<Ord>** oi = new Func<Ord>*[prev_size];
    if (u_ext != Hermes::vector<Solution *>())
      for (int i = 0; i < prev_size; i++)
        if (u_ext[i + vfs->u_ext_offset] != NULL) 
          oi[i] = init_ext_fn_ord(neighbor_searches.get(u_ext[i]->get_mesh()->get_seq() - min_dg_mesh_seq), u_ext[i]);
        else 
          oi[i] = get_fn_ord(0);
    else
      for (int i = 0; i < prev_size; i++) 
        oi[i] = get_fn_ord(0);

    // Order of the shape function.
    // Determine the integration order.
    int inc = (fv->get_num_components() == 2) ? 1 : 0;
    Func<Ord>* ov = get_fn_ord(fv->get_edge_fn_order(surf_pos->surf_num) + inc);

    // Order of additional external functions.
    ExtData<Ord>* fake_ext = init_ext_fns_ord(vfs->ext, neighbor_searches);

    // Order of geometric attributes (eg. for multiplication of a solution with coordinates, normals, etc.).
    Geom<Ord>* fake_e = new InterfaceGeom<Ord>(&geom_ord,
                        nbs_v->neighb_el->marker, nbs_v->neighb_el->id, nbs_v->neighb_el->get_diameter());
    double fake_wt = 1.0;

    // Total order of the vector form.
    Ord o = vfs->ord(1, &fake_wt, oi, ov, fake_e, fake_ext);

    // Increase due to reference map.
    order = rv->get_inv_ref_order();
    order += o.get_order();
    limit_order(order);

    // Clean up.
    if (u_ext != Hermes::vector<Solution *>())
      for (int i = 0; i < prev_size; i++)
        if (u_ext[i + vfs->u_ext_offset] != NULL) 
          delete oi[i];
    delete [] oi;
    if (fake_ext != NULL) {
      for (int i = 0; i < fake_ext->nf; i++) {
        delete fake_ext->fn[i];
      }
      fake_ext->free_ord(); 
      delete fake_ext;
    }

    delete fake_e;
  }
  return order;
}
int DiscreteProblem::calc_order_dg_vector_form(WeakForm::MultiComponentVectorFormSurf *vfs, Hermes::vector<Solution *> u_ext,
                                  PrecalcShapeset *fv, RefMap *rv, SurfPos* surf_pos,
                                  LightArray<NeighborSearch*>& neighbor_searches, int neighbor_index_v)
{
  NeighborSearch* nbs_v = neighbor_searches.get(neighbor_index_v);
  // Order that will be returned.
  int order;
  int u_ext_length = u_ext.size();      // Number of external solutions.
  int u_ext_offset = vfs->u_ext_offset; // External solutions will start with u_ext[u_ext_offset]
                                        // and there will be only u_ext_length - u_ext_offset of them.
  if(this->is_fvm)
    order = rv->get_inv_ref_order();
  else {
    // Order of solutions from the previous Newton iteration.
    int prev_size = u_ext.size() - vfs->u_ext_offset;
    Func<Ord>** oi = new Func<Ord>*[prev_size];
    if (u_ext != Hermes::vector<Solution *>())
      for (int i = 0; i < prev_size; i++)
        if (u_ext[i + vfs->u_ext_offset] != NULL) 
          oi[i] = init_ext_fn_ord(neighbor_searches.get(u_ext[i]->get_mesh()->get_seq() - min_dg_mesh_seq), u_ext[i]);
        else 
          oi[i] = get_fn_ord(0);
    else
      for (int i = 0; i < prev_size; i++) 
        oi[i] = get_fn_ord(0);

    // Order of the shape function.
    // Determine the integration order.
    int inc = (fv->get_num_components() == 2) ? 1 : 0;
    Func<Ord>* ov = get_fn_ord(fv->get_edge_fn_order(surf_pos->surf_num) + inc);

    // Order of additional external functions.
    ExtData<Ord>* fake_ext = init_ext_fns_ord(vfs->ext, neighbor_searches);

    // Order of geometric attributes (eg. for multiplication of a solution with coordinates, normals, etc.).
    Geom<Ord>* fake_e = new InterfaceGeom<Ord>(&geom_ord,
                        nbs_v->neighb_el->marker, nbs_v->neighb_el->id, nbs_v->neighb_el->get_diameter());
    double fake_wt = 1.0;

    // Total order of the vector form.
    Ord o = vfs->ord(1, &fake_wt, oi, ov, fake_e, fake_ext);

    // Increase due to reference map.
    order = rv->get_inv_ref_order();
    order += o.get_order();
    limit_order(order);

    // Clean up.
    if (u_ext != Hermes::vector<Solution *>())
      for (int i = 0; i < prev_size; i++)
        if (u_ext[i + vfs->u_ext_offset] != NULL) 
          delete oi[i];
    delete [] oi;
    if (fake_ext != NULL) {
      for (int i = 0; i < fake_ext->nf; i++) {
        delete fake_ext->fn[i];
      }
      fake_ext->free_ord(); 
      delete fake_ext;
    }

    delete fake_e;
  }
  return order;
}

scalar DiscreteProblem::eval_dg_form(WeakForm::VectorFormSurf* vfs, Hermes::vector<Solution *> u_ext,
                                     PrecalcShapeset *fv, RefMap *rv, 
                                     SurfPos* surf_pos, LightArray<NeighborSearch*>& neighbor_searches, int neighbor_index_v)
{
  _F_
  NeighborSearch* nbs_v = (neighbor_searches.get(neighbor_index_v));
  int order = calc_order_dg_vector_form(vfs, u_ext, fv, rv, surf_pos, neighbor_searches, neighbor_index_v);
  
  // Evaluate the form using just calculated order.
  Quad2D* quad = fv->get_quad_2d();
  int eo = quad->get_edge_points(surf_pos->surf_num, order);
  int np = quad->get_num_points(eo);
  double3* pt = quad->get_points(eo);

  // A (debug) check.
  assert(surf_pos->surf_num == nbs_v->active_edge);

  // Init geometry and jacobian*weights.
  if (cache_e[eo] == NULL) {
    cache_e[eo] = init_geom_surf(rv, surf_pos, eo);
    double3* tan = rv->get_tangent(surf_pos->surf_num, eo);
    cache_jwt[eo] = new double[np];
    for(int i = 0; i < np; i++)
      cache_jwt[eo][i] = pt[i][2] * tan[i][2];
  }

  Geom<double>* e = new InterfaceGeom<double>(cache_e[eo], nbs_v->neighb_el->marker, 
    nbs_v->neighb_el->id, nbs_v->neighb_el->get_diameter());
  double* jwt = cache_jwt[eo];

  // Values of the previous Newton iteration, shape functions and external functions in quadrature points.
  int prev_size = u_ext.size() - vfs->u_ext_offset;
  Func<scalar>** prev = new Func<scalar>*[prev_size];
  if (u_ext != Hermes::vector<Solution *>())
    for (int i = 0; i < prev_size; i++)
      if (u_ext[i + vfs->u_ext_offset] != NULL) {
        neighbor_searches.get(u_ext[i]->get_mesh()->get_seq() - min_dg_mesh_seq)->set_quad_order(order);
        prev[i]  = neighbor_searches.get(u_ext[i]->get_mesh()->get_seq() - min_dg_mesh_seq)->init_ext_fn(u_ext[i]);
      }
      else prev[i] = NULL;
  else
    for (int i = 0; i < prev_size; i++) 
      prev[i] = NULL;

  Func<double>* v = get_fn(fv, rv, eo);
  ExtData<scalar>* ext = init_ext_fns(vfs->ext, neighbor_searches, order);

  scalar res = vfs->value(np, jwt, prev, v, e, ext);

  // Clean up.
  for (int i = 0; i < prev_size; i++) {
    if (prev[i] != NULL) {
      prev[i]->free_fn(); 
      delete prev[i];
    }
  }

  delete [] prev;

  if (ext != NULL) {
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
void DiscreteProblem::eval_dg_form(WeakForm::MultiComponentVectorFormSurf* vfs, Hermes::vector<Solution *> u_ext,
                                     PrecalcShapeset *fv, RefMap *rv, 
                                     SurfPos* surf_pos, LightArray<NeighborSearch*>& neighbor_searches, int neighbor_index_v, Hermes::vector<scalar>& result)
{
  _F_
  NeighborSearch* nbs_v = (neighbor_searches.get(neighbor_index_v));
  int order = calc_order_dg_vector_form(vfs, u_ext, fv, rv, surf_pos, neighbor_searches, neighbor_index_v);
  
  // Evaluate the form using just calculated order.
  Quad2D* quad = fv->get_quad_2d();
  int eo = quad->get_edge_points(surf_pos->surf_num, order);
  int np = quad->get_num_points(eo);
  double3* pt = quad->get_points(eo);

  // A (debug) check.
  assert(surf_pos->surf_num == nbs_v->active_edge);

  // Init geometry and jacobian*weights.
  if (cache_e[eo] == NULL) {
    cache_e[eo] = init_geom_surf(rv, surf_pos, eo);
    double3* tan = rv->get_tangent(surf_pos->surf_num, eo);
    cache_jwt[eo] = new double[np];
    for(int i = 0; i < np; i++)
      cache_jwt[eo][i] = pt[i][2] * tan[i][2];
  }

  Geom<double>* e = new InterfaceGeom<double>(cache_e[eo], nbs_v->neighb_el->marker, 
    nbs_v->neighb_el->id, nbs_v->neighb_el->get_diameter());
  double* jwt = cache_jwt[eo];

  // Values of the previous Newton iteration, shape functions and external functions in quadrature points.
  int prev_size = u_ext.size() - vfs->u_ext_offset;
  Func<scalar>** prev = new Func<scalar>*[prev_size];
  if (u_ext != Hermes::vector<Solution *>())
    for (int i = 0; i < prev_size; i++)
      if (u_ext[i + vfs->u_ext_offset] != NULL) {
        neighbor_searches.get(u_ext[i]->get_mesh()->get_seq() - min_dg_mesh_seq)->set_quad_order(order);
        prev[i]  = neighbor_searches.get(u_ext[i]->get_mesh()->get_seq() - min_dg_mesh_seq)->init_ext_fn(u_ext[i]);
      }
      else prev[i] = NULL;
  else
    for (int i = 0; i < prev_size; i++) 
      prev[i] = NULL;

  Func<double>* v = get_fn(fv, rv, eo);
  ExtData<scalar>* ext = init_ext_fns(vfs->ext, neighbor_searches, order);

  vfs->value(np, jwt, prev, v, e, ext, result);

  for(unsigned int i = 0; i < result.size(); i++)
    result[i] *= vfs->scaling_factor * 0.5;

  // Clean up.
  for (int i = 0; i < prev_size; i++) {
    if (prev[i] != NULL) {
      prev[i]->free_fn(); 
      delete prev[i];
    }
  }

  delete [] prev;

  if (ext != NULL) {
    ext->free(); 
    delete ext;
  }

  delete e;
}

DiscreteProblem::NeighborNode::NeighborNode(NeighborNode* parent, unsigned int transformation) : parent(parent), transformation(transformation)
{
  left_son = right_son = NULL;
}
DiscreteProblem::NeighborNode::~NeighborNode()
{
  if(left_son != NULL) {
    delete left_son;
    left_son = NULL;
  }
  if(right_son != NULL) {
    delete right_son;
    right_son = NULL;
  }
}
void DiscreteProblem::NeighborNode::set_left_son(DiscreteProblem::NeighborNode* left_son)
{
  this->left_son = left_son;
}
void DiscreteProblem::NeighborNode::set_right_son(DiscreteProblem::NeighborNode* right_son)
{
  this->right_son = right_son;
}
void DiscreteProblem::NeighborNode::set_transformation(unsigned int transformation)
{
  this->transformation = transformation;
}
DiscreteProblem::NeighborNode* DiscreteProblem::NeighborNode::get_left_son()
{
  return left_son;
}
DiscreteProblem::NeighborNode* DiscreteProblem::NeighborNode::get_right_son()
{
  return right_son;
}
unsigned int DiscreteProblem::NeighborNode::get_transformation()
{
  return this->transformation;
}


DiscreteProblem::AssemblingCaches::AssemblingCaches()
{

};

DiscreteProblem::AssemblingCaches::~AssemblingCaches()
{
  _F_
  for (std::map<KeyConst, Func<double>*, CompareConst>::const_iterator it = const_cache_fn_triangles.begin();
       it != const_cache_fn_triangles.end(); it++)
  {
    (it->second)->free_fn(); delete (it->second);
  }
  const_cache_fn_triangles.clear();

  for (std::map<KeyConst, Func<double>*, CompareConst>::const_iterator it = const_cache_fn_quads.begin();
       it != const_cache_fn_quads.end(); it++)
  {
    (it->second)->free_fn(); delete (it->second);
  }
  const_cache_fn_quads.clear();

  for(unsigned int i = 0; i < cache_fn_ord.get_size(); i++)
    if(cache_fn_ord.present(i)) {
      cache_fn_ord.get(i)->free_ord(); 
      delete cache_fn_ord.get(i);
    }
};

double Hermes2D::get_l2_norm(Vector* vec) const 
{
  _F_
  scalar val = 0;
  for (unsigned int i = 0; i < vec->length(); i++) {
    scalar inc = vec->get(i);
    val = val + inc*conj(inc);
  }
  return sqrt(std::abs(val));
}

bool Hermes2D::solve_newton(scalar* coeff_vec, DiscreteProblem* dp, Solver* solver, SparseMatrix* matrix,
                            Vector* rhs, double newton_tol, int newton_max_iter, bool verbose,
                            bool residual_as_function,
                            double damping_coeff, double max_allowed_residual_norm) const
{
  // Prepare solutions for measuring residual norm.
  int num_spaces = dp->get_spaces().size();
  Hermes::vector<Solution*> solutions;
  Hermes::vector<bool> dir_lift_false;
  for (int i=0; i < num_spaces; i++) {
    if (residual_as_function) solutions.push_back(new Solution());
    dir_lift_false.push_back(false);      // No Dirichlet lifts will be considered.
  }

  // The Newton's loop.
  double residual_norm;
  int it = 1;
  while (1)
  {
    // Obtain the number of degrees of freedom.
    int ndof = dp->get_num_dofs();

    // Assemble the residual vector.
    dp->assemble(coeff_vec, NULL, rhs); // NULL = we do not want the Jacobian.

    // Measure the residual norm.
    if (residual_as_function) {
      // Translate the residual vector into a residual function (or multiple functions)
      // in the corresponding finite element space(s) and measure their norm(s) there.
      // This is more meaningful than just measuring the l2-norm of the residual vector,
      // since in the FE space not all components in the residual vector have the same weight.
      // On the other hand, this is slower as it requires global norm calculation, and thus
      // numerical integration over the entire domain. Therefore this option is off by default.
      Solution::vector_to_solutions(rhs, dp->get_spaces(), solutions, dir_lift_false);
      residual_norm = calc_norms(solutions);
    }
    else {
      // Calculate the l2-norm of residual vector, this is the traditional way.
      residual_norm = get_l2_norm(rhs);
    }

    // Info for the user.
    if (it == 1) {
      if (verbose) info("---- Newton initial residual norm: %g", residual_norm);
    }
    else if (verbose) info("---- Newton iter %d, residual norm: %g", it-1, residual_norm);

    // If maximum allowed residual norm is exceeded, fail.
    if (residual_norm > max_allowed_residual_norm) {
      if (verbose) {
        info("Current residual norm: %g", residual_norm);
        info("Maximum allowed residual norm: %g", max_allowed_residual_norm);
        info("Newton solve not successful, returning false.");
      }
      for (unsigned int i = 0; i < solutions.size(); i++)
        delete solutions[i];
      return false;
    }

    // If residual norm is within tolerance, or the maximum number
    // of iteration has been reached, then quit.
    if ((residual_norm < newton_tol || it > newton_max_iter) && it > 1) break;

    // Assemble the Jacobian matrix.
    dp->assemble(coeff_vec, matrix, NULL); // NULL = we do not want the rhs.

    // Multiply the residual vector with -1 since the matrix
    // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
    rhs->change_sign();

    // Solve the linear system.
    if(!solver->solve()) error ("Matrix solver failed.\n");

    // Add \deltaY^{n+1} to Y^n.
    for (int i = 0; i < ndof; i++) coeff_vec[i] += damping_coeff * solver->get_solution()[i];

    it++;
  }

  for (unsigned int i = 0; i < solutions.size(); i++)
    delete solutions[i];

  if (it >= newton_max_iter) {
    if (verbose) info("Maximum allowed number of Newton iterations exceeded, returning false.");
    return false;
  }

  return true;
}

// Perform Picard's iteration.
bool Hermes2D::solve_picard(WeakForm* wf, Space* space, Solution* sln_prev_iter,
                  MatrixSolverType matrix_solver, double picard_tol,
                  int picard_max_iter, bool verbose) const
{
  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Initialize the FE problem.
  bool is_linear = true;
  DiscreteProblem dp(wf, space, is_linear);

  int iter_count = 0;
  while (true) {
    // Assemble the stiffness matrix and right-hand side.
    dp.assemble(matrix, rhs);

    // Solve the linear system and if successful, obtain the solution.
    Solution sln_new;
    if(solver->solve()) Solution::vector_to_solution(solver->get_solution(), space, &sln_new);
    else error ("Matrix solver failed.\n");

    double rel_error = calc_abs_error(sln_prev_iter, &sln_new, HERMES_H1_NORM)
                       / calc_norm(&sln_new, HERMES_H1_NORM) * 100;
    if (verbose) info("---- Picard iter %d, ndof %d, rel. error %g%%",
      iter_count+1, space->get_num_dofs(), rel_error);

    // Stopping criterion.
    if (rel_error < picard_tol) {
      sln_prev_iter->copy(&sln_new);
      delete matrix;
      delete rhs;
      delete solver;
      return true;
    }

    if (iter_count >= picard_max_iter) {
      delete matrix;
      delete rhs;
      delete solver;
      if (verbose) info("Maximum allowed number of Picard iterations exceeded, returning false.");
      return false;
    }

    // Saving solution for the next iteration;
    sln_prev_iter->copy(&sln_new);

    iter_count++;
  }
}
