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

#define HERMES_REPORT_INFO
#define HERMES_REPORT_WARN

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

std::map<DiscreteProblem::SurfVectorFormsKey, double*, DiscreteProblem::SurfVectorFormsKeyCompare>
DiscreteProblem::surf_forms_cache =
    *new std::map<DiscreteProblem::SurfVectorFormsKey, double*, DiscreteProblem::SurfVectorFormsKeyCompare>();

std::map<DiscreteProblem::VolVectorFormsKey, double*, DiscreteProblem::VolVectorFormsKeyCompare>
DiscreteProblem::vol_forms_cache =
    *new std::map<DiscreteProblem::VolVectorFormsKey, double*, DiscreteProblem::VolVectorFormsKeyCompare>();

DiscreteProblem::SurfVectorFormsKey DiscreteProblem::surf_forms_key =
  DiscreteProblem::SurfVectorFormsKey(NULL, 0, 0, 0, 0);

DiscreteProblem::VolVectorFormsKey DiscreteProblem::vol_forms_key =
  DiscreteProblem::VolVectorFormsKey(NULL, 0, 0);

void DiscreteProblem::empty_form_caches()
{
  std::map<DiscreteProblem::SurfVectorFormsKey, double*, DiscreteProblem::SurfVectorFormsKeyCompare>::iterator its;
  for(its = DiscreteProblem::surf_forms_cache.begin(); its != DiscreteProblem::surf_forms_cache.end(); its++)
    delete [] (*its).second;

  // This maybe is not needed.
  DiscreteProblem::surf_forms_cache.clear();

  std::map<DiscreteProblem::VolVectorFormsKey, double*, DiscreteProblem::VolVectorFormsKeyCompare>::iterator itv;
  for(itv = DiscreteProblem::vol_forms_cache.begin(); itv != DiscreteProblem::vol_forms_cache.end(); itv++)
    delete [] (*itv).second;

  // This maybe is not needed.
  DiscreteProblem::vol_forms_cache.clear();
};

DiscreteProblem::DiscreteProblem(WeakForm* wf, Hermes::vector<Space *> spaces, bool is_linear) :
  wf(wf), is_linear(is_linear), wf_seq(-1), spaces(spaces)
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
  this->pss = new PrecalcShapeset*[wf->get_neq()];
  for (unsigned int i = 0; i < wf->get_neq(); i++) this->pss[i] = NULL;
  this->num_user_pss = 0;
  for (unsigned int i = 0; i < wf->get_neq(); i++){
    Shapeset *shapeset = spaces[i]->get_shapeset();
    if (shapeset == NULL) error("Internal in DiscreteProblem::init_spaces().");
    PrecalcShapeset *p = new PrecalcShapeset(shapeset);
    if (p == NULL) error("New PrecalcShapeset could not be allocated in DiscreteProblem::init_spaces().");
    this->pss[i] = p;
    this->num_user_pss++;
  }

  // Create global enumeration of dof and fill the ndof variable.
  this->ndof = Space::assign_dofs(this->spaces);

  // Update the weak formulation with the user-supplied string markers
  // according to the conversion table contained in the mesh.
  this->wf->update_markers_acc_to_conversion(spaces[0]->get_mesh()->markers_conversion);

  // There is a special function that sets a DiscreteProblem to be FVM.
  // Purpose is that this constructor looks cleaner and is simpler.
  this->is_fvm = false;

  vector_valued_forms = false;
}

DiscreteProblem::~DiscreteProblem()
{
  _F_
  free();
  if (sp_seq != NULL) delete [] sp_seq;
  for(int i = 0; i < num_user_pss; i++)
    delete pss[i];
}

void DiscreteProblem::free()
{
  _F_
  struct_changed = values_changed = true;
  memset(sp_seq, -1, sizeof(int) * wf->get_neq());
  wf_seq = -1;
}

int DiscreteProblem::get_num_dofs()
{
  _F_
  ndof = 0;
  for (unsigned int i = 0; i < wf->get_neq(); i++)
    ndof += Space::get_num_dofs(spaces[i]);
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

//// matrix structure precalculation ///////////////////////////////////////////////////////////////

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

//// matrix creation ///////////////////////////////////////////////////////////////////////////////

// This functions is identical in H2D and H3D.
void DiscreteProblem::create_sparse_structure(SparseMatrix* mat, Vector* rhs, bool rhsonly,
                             bool force_diagonal_blocks, Table* block_weights)
{
  _F_

  if (is_up_to_date())
  {
    if (!rhsonly && mat != NULL)
    {
      verbose("Reusing matrix sparse structure.");
      mat->zero();
    }
    if (rhs != NULL) rhs->zero();
    return;
  }

  // For DG, the sparse structure is different as we have to account for over-edge calculations.
  bool is_DG = false;
  for(unsigned int i = 0; i < this->wf->mfsurf.size(); i++) {
    if(this->wf->mfsurf[i].area == H2D_DG_INNER_EDGE) {
      is_DG = true;
      break;
    }
  }
  for(unsigned int i = 0; i < this->wf->vfsurf.size(); i++) {
    if(this->wf->vfsurf[i].area == H2D_DG_INNER_EDGE) {
      is_DG = true;
      break;
    }
  }

  int ndof = get_num_dofs();

  if (mat != NULL)  // mat may be NULL when assembling the rhs for NOX
  {
    // Spaces have changed: create the matrix from scratch.
    mat->free();
    mat->prealloc(ndof);

    AUTOLA_CL(AsmList, al, wf->get_neq());
    AUTOLA_OR(Mesh*, meshes, wf->get_neq());
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
            ns.set_active_edge(ed, false);
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

        // Deallocation an array of arrays of neighboring elements for every mesh x edge.
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
  for (unsigned int i = 0; i < wf->get_neq(); i++)
    sp_seq[i] = spaces[i]->get_seq();

  wf_seq = wf->get_seq();

  struct_changed = true;
  have_matrix = true;
}

//// assembly //////////////////////////////////////////////////////////////////////////////////////

// Light version for linear problems.
// The Table is here for optional weighting of matrix blocks in systems.
void DiscreteProblem::assemble(SparseMatrix* mat, Vector* rhs, bool rhsonly,
                               bool force_diagonal_blocks, Table* block_weights)
{
  _F_
  scalar* coeff_vec = NULL;
  bool add_dir_lift = false;
  assemble(coeff_vec, mat, rhs, rhsonly, force_diagonal_blocks, add_dir_lift, block_weights);
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

void DiscreteProblem::assemble(scalar* coeff_vec, SparseMatrix* mat, Vector* rhs, bool rhsonly,
                               bool force_diagonal_blocks, bool add_dir_lift, Table* block_weights)
{
  _F_
  // Sanity checks.
  assemble_sanity_checks(block_weights);

  // Creating matrix sparse structure.
  create_sparse_structure(mat, rhs, rhsonly, force_diagonal_blocks, block_weights);
 
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
  if (mat != NULL) 
    get_matrix_buffer(9);
  
  // Create assembling stages.
  std::vector<WeakForm::Stage> stages = std::vector<WeakForm::Stage>();
  wf->get_stages(spaces, u_ext, stages, rhsonly);


  // Loop through all assembling stages -- the purpose of this is increased performance
  // in multi-mesh calculations, where, e.g., only the right hand side uses two meshes.
  // In such a case, the matrix forms are assembled over one mesh, and only the rhs
  // traverses through the union mesh. On the other hand, if you don't use multi-mesh
  // at all, there will always be only one stage in which all forms are assembled as usual.
  for (unsigned ss = 0; ss < stages.size(); ss++) 
    // Assemble one stage. One stage is a collection of functions, and meshes that can not be further minimized.
    // E.g. if a linear form uses two external solution, each of which is defined on a different mesh, and different to the
    // mesh of the current test function, then the stage would have three meshes.
    // By stage functions, all functions are meant: shape function (their precalculated values), and mesh functions.
    assemble_one_stage(stages[ss], mat, rhs, rhsonly, force_diagonal_blocks, block_weights, spss, refmap, u_ext);
  
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
      SparseMatrix* mat, Vector* rhs, bool rhsonly, bool force_diagonal_blocks, Table* block_weights,
        Hermes::vector<PrecalcShapeset *>& spss, Hermes::vector<RefMap *>& refmap, Hermes::vector<Solution *>& u_ext)
{
  _F_
  // Boundary flags. bnd[i] == true := i-th edge of the current Element is a boundary edge.
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


  // Loop through all assembling states.
  // Assemble each one.
  Element** e;
  while ((e = trav.get_next_state(bnd, surf_pos)) != NULL)
    // One state is a collection of (virtual) elements sharing the same physical location on (possibly) different meshes.
    // This is then the same element of the virtual union mesh. The proper sub-element mappings to all the functions of
    // this stage is supplied by the function Traverse::get_next_state() called in the while loop.
    assemble_one_state(stage, mat, rhs, rhsonly, force_diagonal_blocks, block_weights, spss, refmap, u_ext, e, bnd, surf_pos, trav.get_base());

  if (mat != NULL) mat->finish();
  if (rhs != NULL) rhs->finish();
  trav.finish();
  
}

Element* DiscreteProblem::init_state(WeakForm::Stage& stage, Hermes::vector<PrecalcShapeset *>& spss, 
  Hermes::vector<RefMap *>& refmap, Element** e, Hermes::vector<bool>& isempty, Hermes::vector<AsmList *>& al)
{
  _F_
  // Find a non-NULL e[i].
  Element* e0;
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
    e[i]->visited = true;
  }
  return e0;
}

void DiscreteProblem::assemble_one_state(WeakForm::Stage& stage, 
      SparseMatrix* mat, Vector* rhs, bool rhsonly, bool force_diagonal_blocks, Table* block_weights,
        Hermes::vector<PrecalcShapeset *>& spss, Hermes::vector<RefMap *>& refmap, Hermes::vector<Solution *>& u_ext, Element** e, 
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

  // Element usage flag: iempty[i] == true := The current state does not posses an active element in the i-th space.
  Hermes::vector<bool> isempty;
  for(unsigned int i = 0; i < wf->get_neq(); i++) 
    isempty.push_back(false);
    
  // Initialize the state, return a non-NULL element, if no such Element found, return.
  Element* rep_element = init_state(stage, spss, refmap, e, isempty, al);
  if(rep_element == NULL)
    return;

  init_cache();

  /// Assemble volume matrix forms.
  assemble_volume_matrix_forms(stage, mat, rhs, rhsonly, force_diagonal_blocks, block_weights, spss, refmap, u_ext, isempty, rep_element->marker, al);

  /// Assemble volume vector forms.
  if (rhs != NULL)
    assemble_volume_vector_forms(stage, mat, rhs, rhsonly, force_diagonal_blocks, block_weights, spss, refmap, u_ext, isempty, rep_element->marker, al);

  // Assemble surface integrals now: loop through surfaces of the element.
  for (int isurf = 0; isurf < e[0]->get_num_surf(); isurf++)
   assemble_surface_integrals(stage, mat, rhs, rhsonly, force_diagonal_blocks, block_weights, spss, refmap, u_ext, isempty, 
   surf_pos[isurf].marker, al, bnd[isurf], surf_pos[isurf], nat, isurf, e, trav_base, rep_element);
  
  // Delete assembly lists.
  for(unsigned int i = 0; i < wf->get_neq(); i++) 
    delete al[i];

  delete_cache();
}

void DiscreteProblem::assemble_volume_matrix_forms(WeakForm::Stage& stage, 
      SparseMatrix* mat, Vector* rhs, bool rhsonly, bool force_diagonal_blocks, Table* block_weights,
       Hermes::vector<PrecalcShapeset *>& spss, Hermes::vector<RefMap *>& refmap, Hermes::vector<Solution *>& u_ext, 
       Hermes::vector<bool>& isempty, int marker, Hermes::vector<AsmList *>& al)
{
  _F_
  if (mat != NULL) {
    for (unsigned ww = 0; ww < stage.mfvol.size(); ww++) {
      WeakForm::MatrixFormVol* mfv = stage.mfvol[ww];
      int m = mfv->i;
      int n = mfv->j;
      if (isempty[m] || isempty[n])
        continue;
      if (fabs(mfv->scaling_factor) < 1e-12)
        continue;
      if (mfv->area != HERMES_ANY && !this->wf->is_in_area(marker, mfv->area))
        continue;

      // If a block scaling table is provided, and if the scaling coefficient
      // A_mn for this block is zero, then the form does not need to be assembled.
      scalar block_scaling_coeff = 1.;
      if (block_weights != NULL) {
        if (fabs(block_weights->get_A(m, n)) < 1e-12)
          continue;
        block_scaling_coeff = block_weights->get_A(m, n);
      }
      bool tra = (m != n) && (mfv->sym != 0);
      bool sym = (m == n) && (mfv->sym == 1);

      // Assemble the local stiffness matrix for the form mfv.
      scalar **local_stiffness_matrix;
      if(rhsonly == false)
        local_stiffness_matrix = get_matrix_buffer(std::max(al[m]->cnt, al[n]->cnt));

      for (unsigned int i = 0; i < al[m]->cnt; i++) {
        if (!tra && al[m]->dof[i] < 0) 
          continue;
        spss[m]->set_active_shape(al[m]->idx[i]);
        
        // Unsymmetric block .
        if (!sym) {
          for (unsigned int j = 0; j < al[n]->cnt; j++) {
            
            pss[n]->set_active_shape(al[n]->idx[j]);
            
            if (al[n]->dof[j] < 0) {
              // Linear problems only: Subtracting Dirichlet lift contribution from the RHS:
              if (rhs != NULL && this->is_linear) {
                scalar val = eval_form(mfv, u_ext, pss[n], spss[m], refmap[n],
                        refmap[m]) * al[n]->coef[j] * al[m]->coef[i];
                if(al[m]->dof[i] >= 0)
                  rhs->add(al[m]->dof[i], -val);
              }
            }
            else if (rhsonly == false) {
              scalar val = block_scaling_coeff * eval_form(mfv, u_ext, pss[n], spss[m], refmap[n],
                      refmap[m]) * al[n]->coef[j] * al[m]->coef[i];
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
                scalar val = eval_form(mfv, u_ext, pss[n], spss[m], refmap[n],
                        refmap[m]) * al[n]->coef[j] * al[m]->coef[i];
                if(al[m]->dof[i] >= 0)
                  rhs->add(al[m]->dof[i], -val);
              }
            }
            else if (rhsonly == false) {
              scalar val = block_scaling_coeff * eval_form(mfv, u_ext, pss[n], spss[m], refmap[n],
                      refmap[m]) * al[n]->coef[j] * al[m]->coef[i];
              local_stiffness_matrix[i][j] = local_stiffness_matrix[j][i] = val;
            }
          }
        }
      }

      // Insert the local stiffness matrix into the global one.
      if (rhsonly == false)
        mat->add(al[m]->cnt, al[n]->cnt, local_stiffness_matrix, al[m]->dof, al[n]->dof);

      // Insert also the off-diagonal (anti-)symmetric block, if required.
      if (tra) {
        if (mfv->sym < 0)
          chsgn(local_stiffness_matrix, al[m]->cnt, al[n]->cnt);

        transpose(local_stiffness_matrix, al[m]->cnt, al[n]->cnt);

        if (rhsonly == false)
          mat->add(al[n]->cnt, al[m]->cnt, local_stiffness_matrix, al[n]->dof, al[m]->dof);

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
}

void DiscreteProblem::assemble_volume_vector_forms(WeakForm::Stage& stage, 
      SparseMatrix* mat, Vector* rhs, bool rhsonly, bool force_diagonal_blocks, Table* block_weights,
       Hermes::vector<PrecalcShapeset *>& spss, Hermes::vector<RefMap *>& refmap, Hermes::vector<Solution *>& u_ext, 
       Hermes::vector<bool>& isempty, int marker, Hermes::vector<AsmList *>& al)
{
  _F_
  for (unsigned int ww = 0; ww < stage.vfvol.size(); ww++) {
    WeakForm::VectorFormVol* vfv = stage.vfvol[ww];
    int m = vfv->i;
    if (isempty[vfv->i]) 
      continue;
    if (fabs(vfv->scaling_factor) < 1e-12) 
      continue;
    if (vfv->area != HERMES_ANY && !this->wf->is_in_area(marker, vfv->area)) 
      continue;

    for (unsigned int i = 0; i < al[m]->cnt; i++) {
      if (al[m]->dof[i] < 0) 
        continue;
      
      spss[m]->set_active_shape(al[m]->idx[i]);
      
      if(vector_valued_forms) {
        vol_forms_key = VolVectorFormsKey(vfv->fn, spss[m]->get_active_element()->id, al[m]->idx[i]);
        if(vol_forms_cache[vol_forms_key] == NULL) 
          rhs->add(al[m]->dof[i], eval_form(vfv, u_ext, spss[m], refmap[m]) * al[m]->coef[i]);
        else
          rhs->add(al[m]->dof[i], vol_forms_cache[vol_forms_key][m]);
      }
      else {
        scalar val = eval_form(vfv, u_ext, spss[m], refmap[m]) * al[m]->coef[i];
        rhs->add(al[m]->dof[i], val);
      }
    }
  }
}

void DiscreteProblem::assemble_surface_integrals(WeakForm::Stage& stage, 
      SparseMatrix* mat, Vector* rhs, bool rhsonly, bool force_diagonal_blocks, Table* block_weights,
       Hermes::vector<PrecalcShapeset *>& spss, Hermes::vector<RefMap *>& refmap, Hermes::vector<Solution *>& u_ext, 
       Hermes::vector<bool>& isempty, int marker, Hermes::vector<AsmList *>& al, bool bnd, SurfPos& surf_pos, 
       Hermes::vector<bool>& nat, int isurf, Element** e, Element* trav_base, Element* rep_element)
{
  _F_
  // Obtain the list of shape functions which are nonzero on this surface.
  for (unsigned int i = 0; i < stage.idx.size(); i++) {
    int j = stage.idx[i];
    if (isempty[j])
      continue;
    // For inner edges (with marker == 0), bc_types should not be called,
    // for them it is not important what value (true/false) is set, as it
    // is not read anywhere.
    if(marker > 0)
      nat[j] = (spaces[j]->bc_types->get_type(marker) == BC_NATURAL);
    spaces[j]->get_boundary_assembly_list(e[i], isurf, al[j]);
  }

  // Assemble boundary edges: 
  if(bnd == 1) {
    if (mat != NULL)
      assemble_surface_matrix_forms(stage, mat, rhs, rhsonly, force_diagonal_blocks, block_weights, spss, refmap, u_ext, isempty, 
        marker, al, bnd, surf_pos, nat, isurf, e, trav_base);
    if (rhs != NULL)
      assemble_surface_vector_forms(stage, mat, rhs, rhsonly, force_diagonal_blocks, block_weights, spss, refmap, u_ext, isempty, 
        marker, al, bnd, surf_pos, nat, isurf, e, trav_base);
  }
  // Assemble inner edges (in discontinuous Galerkin discretization): 
  else {
    // Initialize the NeighborSearches.
    std::map<unsigned int, NeighborSearch> neighbor_searches = init_neighbors(stage, isurf);
    
    if(neighbor_searches.empty())
      error("Empty neighbor_searches in DiscreteProblem::assemble_surface_integrals().");

    // Create a multimesh tree;
    DiscreteProblem::NeighborNode* root = new DiscreteProblem::NeighborNode(NULL, 0);
    build_multimesh_tree(root, neighbor_searches);
    
    // Update all NeighborSearches according to the multimesh tree.
    // After this, all NeighborSearches in neighbor_searches should have the same count of neighbors and proper set of transformations
    // for the central and the neighbor element(s) alike.
    for(std::map<unsigned int, NeighborSearch>::iterator it = neighbor_searches.begin(); it != neighbor_searches.end(); it++)
      update_neighbor_search(it->second, root);

    // Check that every NeighborSearch has the same number of neighbor elements.
    unsigned int num_neighbors = neighbor_searches.begin()->second.get_num_neighbors();
    for(std::map<unsigned int, NeighborSearch>::iterator it = neighbor_searches.begin(); it != neighbor_searches.end(); it++)
      if(it->second.get_num_neighbors() != num_neighbors)
        error("Num_neighbors of different NeighborSearches not matching in DiscreteProblem::assemble_surface_integrals().");

    // Create neighbor psss, refmaps.
    Hermes::vector<PrecalcShapeset *> npss;
    Hermes::vector<PrecalcShapeset *> nspss;
    Hermes::vector<RefMap *> nrefmap;

    // Initialize neighbor precalc shapesets and refmaps.      
    for (unsigned int i = 0; i < stage.idx.size(); i++) {
      npss.push_back(new PrecalcShapeset(pss[i]->get_shapeset()));
      npss[i]->set_quad_2d(&g_quad_2d_std);
      nspss.push_back(new PrecalcShapeset(npss[i]));
      nspss[i]->set_quad_2d(&g_quad_2d_std);
      nrefmap.push_back(new RefMap());
      nrefmap[i]->set_quad_2d(&g_quad_2d_std);
    }

    for(unsigned int neighbor_i = 0; neighbor_i < num_neighbors; neighbor_i++) {
      // Set the active segment in all NeighborSearches.

      // Push all the necessary transformations to all functions of this stage.
      // The important thing is that the transformations to the current subelement are already there.
      // Also store the current neighbor element and neighbor edge in neighb_el, neighbor_edge.
      for(unsigned int fns_i = 0; fns_i < stage.fns.size(); fns_i++)
        for(unsigned int trf_i = 0; trf_i < neighbor_searches.at(stage.meshes[fns_i]->get_seq()).central_n_trans[neighbor_i]; trf_i++) {
          neighbor_searches.at(stage.meshes[fns_i]->get_seq()).original_central_el_transform = stage.fns[fns_i]->get_transform();
          neighbor_searches.at(stage.meshes[fns_i]->get_seq()).active_segment = neighbor_i;
          stage.fns[fns_i]->push_transform(neighbor_searches.at(stage.meshes[fns_i]->get_seq()).central_transformations[neighbor_i][trf_i]);
        }
      // For neighbor psss.
      for(unsigned int idx_i = 0; idx_i < stage.idx.size(); idx_i++) {
        npss[idx_i]->set_active_element((*neighbor_searches.at(stage.meshes[idx_i]->get_seq()).get_neighbors())[neighbor_i]);
        for(unsigned int trf_i = 0; trf_i < neighbor_searches.at(stage.meshes[idx_i]->get_seq()).neighbor_n_trans[neighbor_i]; trf_i++)
          npss[idx_i]->push_transform(neighbor_searches.at(stage.meshes[idx_i]->get_seq()).neighbor_transformations[neighbor_i][trf_i]);
      }

      // Also push the transformations to the slave psss and refmaps.
      for (unsigned int i = 0; i < stage.idx.size(); i++) {
        if(isempty[stage.idx[i]])
          continue;
        spss[stage.idx[i]]->set_master_transform();
        refmap[stage.idx[i]]->force_transform(pss[stage.idx[i]]->get_transform(), pss[stage.idx[i]]->get_ctm());

        // Neighbor.
        nspss[stage.idx[i]]->set_active_element(npss[stage.idx[i]]->get_active_element());
        nspss[stage.idx[i]]->set_master_transform();
        nrefmap[stage.idx[i]]->set_active_element(npss[stage.idx[i]]->get_active_element());
        nrefmap[stage.idx[i]]->force_transform(npss[stage.idx[i]]->get_transform(), npss[stage.idx[i]]->get_ctm());
      }


      // The computation takes place here.
    if (mat != NULL)
        assemble_DG_matrix_forms(stage, mat, rhs, rhsonly, force_diagonal_blocks, block_weights, spss, refmap, npss, nspss, nrefmap, neighbor_searches, u_ext, isempty, 
        marker, al, bnd, surf_pos, nat, isurf, e, trav_base, rep_element);
    if (rhs != NULL)
        assemble_DG_vector_forms(stage, mat, rhs, rhsonly, force_diagonal_blocks, block_weights, spss, refmap, neighbor_searches, u_ext, isempty, 
        marker, al, bnd, surf_pos, nat, isurf, e, trav_base, rep_element);


      // This is just cleaning after ourselves.
      // Clear the transformations from the RefMaps and all functions.
      for(unsigned int fns_i = 0; fns_i < stage.fns.size(); fns_i++)
        for(unsigned int trf_i = 0; trf_i < neighbor_searches.at(stage.meshes[fns_i]->get_seq()).central_n_trans[neighbor_i]; trf_i++) {
          stage.fns[fns_i]->set_transform(neighbor_searches.at(stage.meshes[fns_i]->get_seq()).original_central_el_transform);
        }

      // Also clear the transformations from the slave psss and refmaps.
      for (unsigned int i = 0; i < stage.idx.size(); i++) {
        if(isempty[stage.idx[i]])
          continue;
        spss[stage.idx[i]]->set_master_transform();
        refmap[stage.idx[i]]->force_transform(pss[stage.idx[i]]->get_transform(), pss[stage.idx[i]]->get_ctm());
      }

       // Deinitialize neighbor pss's, refmaps.
      for(std::vector<PrecalcShapeset *>::iterator it = nspss.begin(); it != nspss.end(); it++)
        delete *it;
      for(std::vector<PrecalcShapeset *>::iterator it = npss.begin(); it != npss.end(); it++)
        delete *it;
      for(std::vector<RefMap *>::iterator it = nrefmap.begin(); it != nrefmap.end(); it++)
        delete *it;
    }
  }
}

std::map<unsigned int, NeighborSearch> DiscreteProblem::init_neighbors(const WeakForm::Stage& stage, const int& isurf)
{
  // Setup a mapping of mesh sequence numbers to current subelement transformations.
  std::map<unsigned int, uint64_t> initial_sub_idxs;

  // Initialize the NeighborSearches.
  std::map<unsigned int, NeighborSearch> neighbor_searches;
  for(unsigned int i = 0; i < stage.meshes.size(); i++) {
    neighbor_searches.insert(std::pair<unsigned int, NeighborSearch>(stage.meshes[i]->get_seq(), NeighborSearch(stage.fns[i]->get_active_element(), stage.meshes[i])));
    initial_sub_idxs.insert(std::pair<unsigned int, uint64_t>(stage.meshes[i]->get_seq(), stage.fns[i]->get_transform()));
  }

  // Calculate respective neighbors.
  for(std::map<unsigned int, NeighborSearch>::iterator it = neighbor_searches.begin(); it != neighbor_searches.end(); it++)
    it->second.set_active_edge_multimesh(isurf, initial_sub_idxs[it->first]);

  // Clear the initial_sub_idxs from the central element transformations of NeighborSearches with multiple neighbors.
  for(std::map<unsigned int, NeighborSearch>::iterator it = neighbor_searches.begin(); it != neighbor_searches.end(); it++)
    it->second.clear_initial_sub_idx(initial_sub_idxs[it->first]);

  return neighbor_searches;
}


void DiscreteProblem::build_multimesh_tree(DiscreteProblem::NeighborNode* root, const std::map<unsigned int, NeighborSearch>& neighbor_searches)
{
  for(std::map<unsigned int, NeighborSearch>::const_iterator it = neighbor_searches.begin(); it != neighbor_searches.end(); it++) {
    if(it->second.get_num_neighbors() == 1)
      break;
    for(unsigned int i = 0; i < it->second.get_num_neighbors(); i++)
      insert_into_multimesh_tree(root, it->second.central_transformations[i], it->second.central_n_trans[i]);
  }
}

void DiscreteProblem::insert_into_multimesh_tree(NeighborNode* node, int* transformations, unsigned int transformation_count)
{
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

/* This function is not used. It may be used when the implementation changes. Will be deleted when found not necessary.*/
Hermes::vector<Hermes::vector<int>*> DiscreteProblem::get_multimesh_neighbors_transformations(DiscreteProblem::NeighborNode* multimesh_tree)
{
  // Initialize the vector.
  Hermes::vector<Hermes::vector<int>*> running_transformations;
  // Prepare the first neighbor's vector.
  running_transformations.push_back(new Hermes::vector<int>);
  // Fill the vector.
  traverse_multimesh_tree(multimesh_tree, running_transformations);
  return running_transformations;
}

/* This function is not used. It may be used when the implementation changes. Will be deleted when found not necessary.*/
void DiscreteProblem::traverse_multimesh_tree(DiscreteProblem::NeighborNode* node, Hermes::vector<Hermes::vector<int>*>& running_transformations)
{
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
    Hermes::vector<int>* new_neighbor_transformations = new Hermes::vector<int>;
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

void DiscreteProblem::update_neighbor_search(NeighborSearch& ns, NeighborNode* multimesh_tree)
{
  // This has to be done, because we pass ns by reference and the number of neighbors is changing.
  unsigned int num_neighbors = ns.get_num_neighbors();

  for(unsigned int i = 0; i < num_neighbors; i++) {
    // Find the node corresponding to this neighbor in the tree.
    NeighborNode* node = find_node(ns.central_transformations[i], ns.central_n_trans[i], multimesh_tree);
    // Update the NeighborSearch.
    update_ns_subtree(ns, node, i);
  }
}

DiscreteProblem::NeighborNode* DiscreteProblem::find_node(int* transformations, int transformation_count, DiscreteProblem::NeighborNode* node)
{
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
    
void DiscreteProblem::update_ns_subtree(NeighborSearch& ns, DiscreteProblem::NeighborNode* node, unsigned int ith_neighbor)
{
  // No subtree => no work.
  // Also check the assertion that if one son is null, then the other too.
  if(node->get_left_son() == NULL) {
    if(node->get_right_son() != NULL)
      error("Only one son (right) not null in DiscreteProblem::update_ns_subtree.");
    return;
  }
  // Check the opposite and proceed.
  else
    if(node->get_right_son() == NULL)
      error("Only one son (left) not null in DiscreteProblem::update_ns_subtree.");

  // Key part.
  // Begin with storing the info about the current neighbor.
  Element* neighbor = ns.neighbors[ith_neighbor];
  NeighborSearch::NeighborEdgeInfo edge_info = ns.neighbor_edges[ith_neighbor];

  // Initialize the vector for central transformations.
  Hermes::vector<Hermes::vector<int>*> running_central_transformations;
  // Prepare the first new neighbor's vector. Push back the current transformations (in case of GO_DOWN neighborhood).
  running_central_transformations.push_back(new Hermes::vector<int>);
  for(unsigned int i = 0; i < ns.central_n_trans[ith_neighbor]; i++)
    running_central_transformations.back()->push_back(ns.central_transformations[ith_neighbor][i]);

  // Initialize the vector for neighbor transformations.
  Hermes::vector<Hermes::vector<int>*> running_neighbor_transformations;
  // Prepare the first new neighbor's vector. Push back the current transformations (in case of GO_UP/NO_TRF neighborhood).
  running_neighbor_transformations.push_back(new Hermes::vector<int>);
  for(unsigned int i = 0; i < ns.neighbor_n_trans[ith_neighbor]; i++)
    running_neighbor_transformations.back()->push_back(ns.neighbor_transformations[ith_neighbor][i]);

  // Delete the current neighbor.
  ns.delete_neighbor(ith_neighbor);

  // Move down the subtree.
  if(node->get_left_son() != NULL)
    traverse_multimesh_subtree(node->get_left_son(), running_central_transformations, running_neighbor_transformations, edge_info, ns.active_edge, ns.central_el->get_mode());
  if(node->get_right_son() != NULL)
    traverse_multimesh_subtree(node->get_right_son(), running_central_transformations, running_neighbor_transformations, edge_info, ns.active_edge, ns.central_el->get_mode());

  // Delete the last neighbors' info (this is a dead end, caused by the function traverse_multimesh_subtree.
  running_central_transformations.pop_back();
  running_neighbor_transformations.pop_back();

  // Insert new neighbors.
  for(unsigned int i = 0; i < running_central_transformations.size(); i++) {
    ns.neighbors.insert(ns.neighbors.begin() + ith_neighbor + i, neighbor);
    ns.neighbor_edges.insert(ns.neighbor_edges.begin() + ith_neighbor + i, edge_info);
    ns.central_n_trans.insert(ns.central_n_trans.begin() + ith_neighbor + i, running_central_transformations[i]->size());
    ns.neighbor_n_trans.insert(ns.neighbor_n_trans.begin() + ith_neighbor + i, running_neighbor_transformations[i]->size());
    ns.central_transformations.insert(ns.central_transformations.begin() + ith_neighbor + i, new int[NeighborSearch::max_n_trans]);
    for(unsigned int ii = 0; ii < ns.central_n_trans[ith_neighbor + i]; ii++)
      ns.central_transformations[ith_neighbor + i][ii] = (*running_central_transformations[i])[ii];
    ns.neighbor_transformations.insert(ns.neighbor_transformations.begin() + ith_neighbor + i, new int[NeighborSearch::max_n_trans]);
    for(unsigned int ii = 0; ii < ns.neighbor_n_trans[ith_neighbor + i]; ii++)
      ns.neighbor_transformations[ith_neighbor + i][ii] = (*running_neighbor_transformations[i])[ii];
    ns.n_neighbors++;
  }
}

void DiscreteProblem::traverse_multimesh_subtree(DiscreteProblem::NeighborNode* node, Hermes::vector<Hermes::vector<int>*>& running_central_transformations,
      Hermes::vector<Hermes::vector<int>*>& running_neighbor_transformations, const NeighborSearch::NeighborEdgeInfo& edge_info, const int& active_edge, const int& mode)
{
  // If we are in a leaf.
  if(node->get_left_son() == NULL && node->get_right_son() == NULL) {
    // Create vectors for the new neighbor.
    Hermes::vector<int>* new_neighbor_central_transformations = new Hermes::vector<int>;
    Hermes::vector<int>* new_neighbor_neighbor_transformations = new Hermes::vector<int>;
    
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
      if ((active_edge == 0 && node->get_transformation() == 1) ||
          (active_edge == 1 && node->get_transformation() == 2) ||
          (active_edge == 2 && node->get_transformation() == 3))
        running_neighbor_transformations.back()->push_back((!edge_info.orientation ? edge_info.local_num_of_edge : (edge_info.local_num_of_edge + 1) % 3));
			else
        running_neighbor_transformations.back()->push_back((edge_info.orientation ? edge_info.local_num_of_edge : (edge_info.local_num_of_edge + 1) % 3));
	    // Quads.
	  else
      if ((active_edge == 0 && (node->get_transformation() == 1 || node->get_transformation() == 7)) ||
          (active_edge == 1 && (node->get_transformation() == 2 || node->get_transformation() == 5)) ||
          (active_edge == 2 && (node->get_transformation() == 3 || node->get_transformation() == 0)) ||
          (active_edge == 3 && (node->get_transformation() == 4 || node->get_transformation() == 6)))
        running_neighbor_transformations.back()->push_back((!edge_info.orientation ? edge_info.local_num_of_edge : (edge_info.local_num_of_edge + 1) % 4));
      else
        running_neighbor_transformations.back()->push_back((edge_info.orientation ? edge_info.local_num_of_edge : (edge_info.local_num_of_edge + 1) % 4));

    // Make the new_neighbor_neighbor_transformations the current running neighbor transformation.
    running_central_transformations.push_back(new_neighbor_neighbor_transformations);
    
    return;
  }
  else {
    // Insert this leaf into the current running central transformation, thus complete it.
    running_central_transformations.back()->push_back(node->get_transformation());

    // Insert appropriate info from this leaf into the current running neighbor transformation, thus complete it.
    // Triangles.
    if(mode == HERMES_MODE_TRIANGLE)
      if ((active_edge == 0 && node->get_transformation() == 1) ||
          (active_edge == 1 && node->get_transformation() == 2) ||
          (active_edge == 2 && node->get_transformation() == 3))
        running_neighbor_transformations.back()->push_back((!edge_info.orientation ? edge_info.local_num_of_edge : (edge_info.local_num_of_edge + 1) % 3));
			else
        running_neighbor_transformations.back()->push_back((edge_info.orientation ? edge_info.local_num_of_edge : (edge_info.local_num_of_edge + 1) % 3));
	  // Quads.
	  else
      if ((active_edge == 0 && (node->get_transformation() == 1 || node->get_transformation() == 7)) ||
          (active_edge == 1 && (node->get_transformation() == 2 || node->get_transformation() == 5)) ||
          (active_edge == 2 && (node->get_transformation() == 3 || node->get_transformation() == 0)) ||
          (active_edge == 3 && (node->get_transformation() == 4 || node->get_transformation() == 6)))
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
      SparseMatrix* mat, Vector* rhs, bool rhsonly, bool force_diagonal_blocks, Table* block_weights,
       Hermes::vector<PrecalcShapeset *>& spss, Hermes::vector<RefMap *>& refmap, Hermes::vector<Solution *>& u_ext, 
       Hermes::vector<bool>& isempty, int marker, Hermes::vector<AsmList *>& al, bool bnd, SurfPos& surf_pos, Hermes::vector<bool>& nat, 
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
        && !wf->is_in_area(marker, mfs->area)) continue;

    // If a block scaling table is provided, and if the scaling coefficient
    // A_mn for this block is zero, then the form does not need to be assembled.
    scalar block_scaling_coeff = 1.;
    if (block_weights != NULL) {
      if (fabs(block_weights->get_A(m, n)) < 1e-12) continue;
      block_scaling_coeff = block_weights->get_A(m, n);
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
            scalar val = eval_form(mfs, u_ext, pss[n], spss[m], refmap[n],
                    refmap[m], &surf_pos) * al[n]->coef[j] * al[m]->coef[i];
            rhs->add(al[m]->dof[i], -val);
          }
        }
        else if (rhsonly == false) {
          scalar val = block_scaling_coeff * eval_form(mfs, u_ext, pss[n], spss[m], refmap[n],
                  refmap[m], &surf_pos) * al[n]->coef[j] * al[m]->coef[i];
          local_stiffness_matrix[i][j] = val;
        }
      }
    }
    if (rhsonly == false)
      mat->add(al[m]->cnt, al[n]->cnt, local_stiffness_matrix, al[m]->dof, al[n]->dof);
  }
}

void DiscreteProblem::assemble_surface_vector_forms(WeakForm::Stage& stage, 
      SparseMatrix* mat, Vector* rhs, bool rhsonly, bool force_diagonal_blocks, Table* block_weights,
       Hermes::vector<PrecalcShapeset *>& spss, Hermes::vector<RefMap *>& refmap, Hermes::vector<Solution *>& u_ext, 
       Hermes::vector<bool>& isempty, int marker, Hermes::vector<AsmList *>& al, bool bnd, SurfPos& surf_pos, Hermes::vector<bool>& nat, 
       int isurf, Element** e, Element* trav_base)
{
  _F_
  for (unsigned int ww = 0; ww < stage.vfsurf.size(); ww++) {
      WeakForm::VectorFormSurf* vfs = stage.vfsurf[ww];
      int m = vfs->i;
      if (isempty[m]) continue;
      if (fabs(vfs->scaling_factor) < 1e-12) continue;
      if (vfs->area == H2D_DG_INNER_EDGE) continue;
      if (vfs->area != HERMES_ANY && vfs->area != H2D_DG_BOUNDARY_EDGE 
          && !wf->is_in_area(marker, vfs->area)) continue;

      if (vfs->area == HERMES_ANY && !nat[m]) continue;

      surf_pos.base = trav_base;
      surf_pos.space_v = spaces[m];

      for (unsigned int i = 0; i < al[m]->cnt; i++) {
        if (al[m]->dof[i] < 0) continue;
        spss[m]->set_active_shape(al[m]->idx[i]);

        if (vector_valued_forms) {
          surf_forms_key = SurfVectorFormsKey(vfs->fn, spss[m]->get_active_element()->id, isurf, al[m]->idx[i],
              spss[m]->get_transform());
          if(surf_forms_cache[surf_forms_key] == NULL)
            rhs->add(al[m]->dof[i], eval_form(vfs, u_ext, spss[m], refmap[m], 
                      &surf_pos) * al[m]->coef[i]);
          else
            rhs->add(al[m]->dof[i], 0.5 * surf_forms_cache[surf_forms_key][m]);
        }
        else {
          scalar val = eval_form(vfs, u_ext, spss[m], refmap[m], &surf_pos) * al[m]->coef[i];
          rhs->add(al[m]->dof[i], val);
        }
      }
      }
}

void DiscreteProblem::assemble_DG_matrix_forms(WeakForm::Stage& stage, 
      SparseMatrix* mat, Vector* rhs, bool rhsonly, bool force_diagonal_blocks, Table* block_weights,
       Hermes::vector<PrecalcShapeset *>& spss, Hermes::vector<RefMap *>& refmap, Hermes::vector<PrecalcShapeset *>& npss, 
       Hermes::vector<PrecalcShapeset *>& nspss, Hermes::vector<RefMap *>& nrefmap, std::map<unsigned int, NeighborSearch>& neighbor_searches, Hermes::vector<Solution *>& u_ext, 
       Hermes::vector<bool>& isempty, int marker, Hermes::vector<AsmList *>& al, bool bnd, SurfPos& surf_pos, Hermes::vector<bool>& nat, 
       int isurf, Element** e, Element* trav_base, Element* rep_element)
{
  _F_

  // If the active segment has already been processed (when the neighbor element was assembled), it is skipped.
  // We test all neighbor searches, because in the case of intra-element edge, the neighboring (the same as central) element
  // will be marked as visited, even though the edge was not calculated.
  bool processed = true;
  for(std::map<unsigned int, NeighborSearch>::iterator it = neighbor_searches.begin(); it != neighbor_searches.end(); it++)
    if(!it->second.neighbors[it->second.active_segment]->visited)
      processed = false;
  if(processed)
    return;

  for (unsigned int ww = 0; ww < stage.mfsurf.size(); ww++) {
            WeakForm::MatrixFormSurf* mfs = stage.mfsurf[ww];
            int m = mfs->i;
            int n = mfs->j;

            if (isempty[m] || isempty[n]) continue;
            if (fabs(mfs->scaling_factor) < 1e-12) continue;
            if (mfs->area != H2D_DG_INNER_EDGE) continue;

            surf_pos.base = trav_base;
            surf_pos.space_v = spaces[m];
            surf_pos.space_u = spaces[n];

              // Create the extended shapeset on the union of the central element and its current neighbor.
    NeighborSearch::ExtendedShapeset* ext_asmlist_u = neighbor_searches.at(stage.meshes[n]->get_seq()).create_extended_asmlist(spaces[n], al[n]);
    NeighborSearch::ExtendedShapeset* ext_asmlist_v = neighbor_searches.at(stage.meshes[m]->get_seq()).create_extended_asmlist(spaces[m], al[m]);

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
        if (!ext_asmlist_v->has_support_on_neighbor(i)) {
          pss[n]->set_active_shape(ext_asmlist_u->central_al->idx[i]);
          fu = pss[n];
          ru = refmap[n];
          support_neigh_u = false;
        }
        else {
          npss[n]->set_active_shape(ext_asmlist_u->neighbor_al->idx[i - ext_asmlist_u->central_al->cnt]);
          fu = npss[n];
          ru = nrefmap[n];
          support_neigh_u = true;
        }

        if (ext_asmlist_u->dof[j] < 0) {
                    if (rhs != NULL && this->is_linear) {
                      // Evaluate the form with the activated discontinuous shape functions.
            scalar val = eval_dg_form(mfs, u_ext, fu, fv, refmap[n], ru, rv, support_neigh_u, support_neigh_v, &surf_pos, neighbor_searches, stage.meshes[n]->get_seq(), stage.meshes[m]->get_seq())
              * support_neigh_u ? ext_asmlist_u->neighbor_al->coef[i - ext_asmlist_u->central_al->cnt]: ext_asmlist_u->central_al->coef[i]
              * support_neigh_v ? ext_asmlist_v->neighbor_al->coef[i - ext_asmlist_v->central_al->cnt]: ext_asmlist_v->central_al->coef[i];
                      // Add the contribution to the global dof index.
            rhs->add(ext_asmlist_v->dof[i], -val);
                    }
                  }
                  else if (rhsonly == false) {
          scalar val = eval_dg_form(mfs, u_ext, fu, fv, refmap[n], ru, rv, support_neigh_u, support_neigh_v, &surf_pos, neighbor_searches, stage.meshes[n]->get_seq(), stage.meshes[m]->get_seq())
            * support_neigh_u ? ext_asmlist_u->neighbor_al->coef[i - ext_asmlist_u->central_al->cnt]: ext_asmlist_u->central_al->coef[i]
            * support_neigh_v ? ext_asmlist_v->neighbor_al->coef[i - ext_asmlist_v->central_al->cnt]: ext_asmlist_v->central_al->coef[i];
                    local_stiffness_matrix[i][j] = val;
                  }
                }
              if (rhsonly == false)
        mat->add(ext_asmlist_u->cnt, ext_asmlist_v->cnt, local_stiffness_matrix,
                ext_asmlist_v->dof, ext_asmlist_u->dof);
            }
          }
}

void DiscreteProblem::assemble_DG_vector_forms(WeakForm::Stage& stage, 
      SparseMatrix* mat, Vector* rhs, bool rhsonly, bool force_diagonal_blocks, Table* block_weights,
       Hermes::vector<PrecalcShapeset *>& spss, Hermes::vector<RefMap *>& refmap, std::map<unsigned int, NeighborSearch>& neighbor_searches, Hermes::vector<Solution *>& u_ext, 
       Hermes::vector<bool>& isempty, int marker, Hermes::vector<AsmList *>& al, bool bnd, SurfPos& surf_pos, Hermes::vector<bool>& nat, 
       int isurf, Element** e, Element* trav_base, Element* rep_element)
{
  _F_
  for (unsigned int ww = 0; ww < stage.vfsurf.size(); ww++) {
            WeakForm::VectorFormSurf* vfs = stage.vfsurf[ww];
            int m = vfs->i;
            if (isempty[m]) continue;
            if (fabs(vfs->scaling_factor) < 1e-12) continue;
            if (vfs->area != H2D_DG_INNER_EDGE) continue;
            
    // Here we use the standard pss, possibly just transformed by NeighborSearch.
              for (unsigned int i = 0; i < al[m]->cnt; i++) {
                if (al[m]->dof[i] < 0) continue;
      spss[m]->set_active_shape(al[m]->idx[i]);
                if(vector_valued_forms) {
        surf_forms_key = SurfVectorFormsKey(vfs->fn, spss[m]->get_active_element()->id, isurf,
            al[m]->idx[i], spss[m]->get_transform());
                  if(surf_forms_cache[surf_forms_key] == NULL)
          rhs->add(al[m]->dof[i], eval_dg_form(vfs, u_ext, spss[m], refmap[m], &surf_pos, neighbor_searches, stage.meshes[m]->get_seq()) * al[m]->coef[i]);
                  else
                    rhs->add(al[m]->dof[i], 0.5 * surf_forms_cache[surf_forms_key][m]);
                }
      else
        rhs->add(al[m]->dof[i], eval_dg_form(vfs, u_ext, spss[m], refmap[m], &surf_pos, neighbor_searches, stage.meshes[m]->get_seq()) * al[m]->coef[i]);
            }
          }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Initialize integration order for external functions
ExtData<Ord>* DiscreteProblem::init_ext_fns_ord(Hermes::vector<MeshFunction *> &ext)
{
  _F_
  ExtData<Ord>* fake_ext = new ExtData<Ord>;
  fake_ext->nf = ext.size();
  Func<Ord>** fake_ext_fn = new Func<Ord>*[fake_ext->nf];
  for (int i = 0; i < fake_ext->nf; i++)
    fake_ext_fn[i] = get_fn_ord(ext[i]->get_fn_order());
  fake_ext->fn = fake_ext_fn;

  return fake_ext;
}

// Initialize external functions (obtain values, derivatives,...)
ExtData<scalar>* DiscreteProblem::init_ext_fns(Hermes::vector<MeshFunction *> &ext, RefMap *rm, const int order)
{
  _F_
  ExtData<scalar>* ext_data = new ExtData<scalar>;
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
ExtData<Ord>* DiscreteProblem::init_ext_fns_ord(Hermes::vector<MeshFunction *> &ext, int edge)
{
  _F_
  ExtData<Ord>* fake_ext = new ExtData<Ord>;
  fake_ext->nf = ext.size();
  Func<Ord>** fake_ext_fn = new Func<Ord>*[fake_ext->nf];
  for (int i = 0; i < fake_ext->nf; i++)
    fake_ext_fn[i] = get_fn_ord(ext[i]->get_edge_fn_order(edge));
  fake_ext->fn = fake_ext_fn;

  return fake_ext;
}

// Initialize discontinuous external functions (obtain values, derivatives,... on both sides of the
// supplied NeighborSearch's active edge).
ExtData<scalar>* DiscreteProblem::init_ext_fns(Hermes::vector<MeshFunction *> &ext, std::map<unsigned int, NeighborSearch>& neighbor_searches)
{
  Func<scalar>** ext_fns = new Func<scalar>*[ext.size()];
  for(unsigned int j = 0; j < ext.size(); j++)
    ext_fns[j] = neighbor_searches.at(ext[j]->get_mesh()->get_seq()).init_ext_fn(ext[j]);

  ExtData<scalar>* ext_data = new ExtData<scalar>;
  ext_data->fn = ext_fns;
  ext_data->nf = ext.size();

  return ext_data;
}

// Initialize integration order for discontinuous external functions.
ExtData<Ord>* DiscreteProblem::init_ext_fns_ord(Hermes::vector<MeshFunction *> &ext, std::map<unsigned int, NeighborSearch>& neighbor_searches)
{
  Func<Ord>** fake_ext_fns = new Func<Ord>*[ext.size()];
  for (unsigned int j = 0; j < ext.size(); j++)
    fake_ext_fns[j] = neighbor_searches.at(ext[j]->get_mesh()->get_seq()).init_ext_fn_ord(ext[j]);

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
    AssemblingCaches::KeyNonConst key(256 - fu->get_active_shape(), order, fu->get_transform(), fu->get_shapeset()->get_id());
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

//  Evaluation of forms  ///////////////////////////////////////////////////////////////////////

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
    Geom<Ord>* fake_e = init_geom_ord();

    // Total order of the matrix form.
    Ord o = mfv->ord(1, &fake_wt, oi, ou, ov, fake_e, fake_ext);

    // Increase due to reference map.
    order = ru->get_inv_ref_order();
    order += o.get_order();
    limit_order_nowarn(order);
    
    // Cleanup.
    delete [] oi;
    delete fake_e;
  }
  return order;
}

// Actual evaluation of volume matrix form (calculates integral)
scalar DiscreteProblem::eval_form(WeakForm::MatrixFormVol *mfv, Hermes::vector<Solution *> u_ext,
                                  PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv)
{
  _F_
  // Determine the integration order.
  int order = calc_order_matrix_form_vol(mfv, u_ext, fu, fv, ru, rv);

  // Evaluate the form using the quadrature of the just calculated order.
  Quad2D* quad = fu->get_quad_2d();
  double3* pt = quad->get_points(order);
  int np = quad->get_num_points(order);

  // Init geometry and jacobian*weights.
  if (cache_e[order] == NULL)
  {
    cache_e[order] = init_geom_vol(ru, order);
    double* jac;
    if(ru->is_jacobian_const()) {
      jac = new double[np];
      double const_jacobian = ru->get_const_jacobian();
      for(int i = 0; i < np; i++) jac[i] = const_jacobian;
    }
    else jac = ru->get_jacobian(order);
    cache_jwt[order] = new double[np];
    for(int i = 0; i < np; i++) cache_jwt[order][i] = pt[i][2] * jac[i];
    if(ru->is_jacobian_const()) delete [] jac;
  }
  Geom<double>* e = cache_e[order];
  double* jwt = cache_jwt[order];

  // Values of the previous Newton iteration, shape functions and external functions in quadrature points.
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
  scalar res = mfv->fn(np, jwt, prev, u, v, e, ext) * mfv->scaling_factor;

  // Clean up.
  for(int i = 0; i < prev_size; i++)
    if (prev[i] != NULL) { 
      prev[i]->free_fn(); 
      delete prev[i]; 
    }
  delete [] prev;

  if (ext != NULL) 
  {
    ext->free(); 
    delete ext;
  }

  return res;
}

int DiscreteProblem::calc_order_vector_form_vol(WeakForm::VectorFormVol *vfv, Hermes::vector<Solution *> u_ext,
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

    // Order of geometric attributes (eg. for multiplication of a solution with coordinates, normals, etc.).
    double fake_wt = 1.0;
    Geom<Ord>* fake_e = init_geom_ord();

    // Total order of the vector form.
    Ord o = vfv->ord(1, &fake_wt, oi, ov, fake_e, fake_ext);

    // Increase due to reference map.
    order = rv->get_inv_ref_order();
    order += o.get_order();
    limit_order_nowarn(order);
    
    // Cleanup.
    delete [] oi;
  }
  return order;
}

// Actual evaluation of volume vector form (calculates integral)
scalar DiscreteProblem::eval_form(WeakForm::VectorFormVol *vfv, 
                                  Hermes::vector<Solution *> u_ext, 
                                  PrecalcShapeset *fv, RefMap *rv)
{
  _F_
  // Determine the integration order.
  int order = calc_order_vector_form_vol(vfv, u_ext, fv, rv);;

  // Evaluate the form using the quadrature of the just calculated order.
  Quad2D* quad = fv->get_quad_2d();
  double3* pt = quad->get_points(order);
  int np = quad->get_num_points(order);

  // Init geometry and jacobian*weights.
  if (cache_e[order] == NULL)
  {
    cache_e[order] = init_geom_vol(rv, order);
    double* jac;
    if(rv->is_jacobian_const()) {
      jac = new double[np];
      double const_jacobian = rv->get_const_jacobian();
      for(int i = 0; i < np; i++)
        jac[i] = const_jacobian;
    }
    else
      jac = rv->get_jacobian(order);
    cache_jwt[order] = new double[np];
    for(int i = 0; i < np; i++)
      cache_jwt[order][i] = pt[i][2] * jac[i];
    if(rv->is_jacobian_const())
      delete [] jac;
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
  scalar res = vfv->fn(np, jwt, prev, v, e, ext) * vfv->scaling_factor;

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
    Geom<Ord>* fake_e = init_geom_ord();

    // Total order of the matrix form.
    Ord o = mfs->ord(1, &fake_wt, oi, ou, ov, fake_e, fake_ext);

    // Increase due to reference map.
    order = ru->get_inv_ref_order();
    order += o.get_order();
    limit_order_nowarn(order);

    // Cleanup.
    delete [] oi;
  }
  return order;
}

// Actual evaluation of surface matrix forms (calculates integral)
scalar DiscreteProblem::eval_form(WeakForm::MatrixFormSurf *mfs, Hermes::vector<Solution *> u_ext,
                        PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru, RefMap *rv, SurfPos* surf_pos)
{
  _F_
  // Determine the integration order.
  int order = calc_order_matrix_form_surf(mfs, u_ext, fu, fv, ru, rv, surf_pos);
  
  // Evaluate the form using the quadrature of the just calculated order.
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
  scalar res = mfs->fn(np, jwt, prev, u, v, e, ext) * mfs->scaling_factor;

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
    Geom<Ord>* fake_e = init_geom_ord();

    // Total order of the vector form.
    Ord o = vfs->ord(1, &fake_wt, oi, ov, fake_e, fake_ext);

    // Increase due to reference map.
    order = rv->get_inv_ref_order();
    order += o.get_order();
    limit_order_nowarn(order);
    
    // Cleanup.
    delete [] oi;
  }
  return order;
}
// Actual evaluation of surface vector form (calculates integral)
scalar DiscreteProblem::eval_form(WeakForm::VectorFormSurf *vfs, Hermes::vector<Solution *> u_ext,
                        PrecalcShapeset *fv, RefMap *rv, SurfPos* surf_pos)
{
  _F_
  // Determine the integration order.
  int order = calc_order_vector_form_surf(vfs, u_ext, fv, rv, surf_pos);

  // Evaluate the form using the quadrature of the just calculated order.
  Quad2D* quad = fv->get_quad_2d();

  int eo = quad->get_edge_points(surf_pos->surf_num, order);
  double3* pt = quad->get_points(eo);
  int np = quad->get_num_points(eo);

  // Init geometry and jacobian*weights.
  if (cache_e[eo] == NULL)
  {
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
  scalar res = vfs->fn(np, jwt, prev, v, e, ext) * vfs->scaling_factor;

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

scalar DiscreteProblem::eval_dg_form(WeakForm::MatrixFormSurf* mfs, Hermes::vector<Solution *> u_ext,
                                     PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *ru_central, RefMap *ru_actual, RefMap *rv, 
                                     bool neighbor_supp_u, bool neighbor_supp_v,
                                     SurfPos* surf_pos, std::map<unsigned int, NeighborSearch>& neighbor_searches, int neighbor_index_u, int neighbor_index_v)
{
  _F_

  NeighborSearch* nbs_u = &(neighbor_searches.at(neighbor_index_u));
  NeighborSearch* nbs_v = &(neighbor_searches.at(neighbor_index_v));

  int order;
  int u_ext_length = u_ext.size();      // Number of external solutions.
  int u_ext_offset = mfs->u_ext_offset; // External solutions will start with u_ext[u_ext_offset]
                                        // and there will be only u_ext_length - u_ext_offset of them.
  if(this->is_fvm)
    order = ru_actual->get_inv_ref_order();
  else {
    // Order of solutions from the previous Newton iteration.
    AUTOLA_OR(Func<Ord>*, oi, u_ext_length - u_ext_offset);
    if (u_ext != Hermes::vector<Solution *>())
      for (int i = u_ext_offset; i < u_ext_length; i++)
        if (u_ext[i] != NULL) 
          oi[i - u_ext_offset] = neighbor_searches.at(u_ext[i]->get_mesh()->get_seq()).init_ext_fn_ord(u_ext[i]);
        else 
          oi[i - u_ext_offset] = init_fn_ord(0);
    else
      for (int i = u_ext_offset; i < u_ext_length; i++) 
        oi[i - u_ext_offset] = init_fn_ord(0);

    // Order of shape functions.
    int inc = (fv->get_num_components() == 2) ? 1 : 0;
    DiscontinuousFunc<Ord>* ou = new DiscontinuousFunc<Ord>(get_fn_ord(fu->get_edge_fn_order(surf_pos->surf_num) + inc), neighbor_supp_u);
    DiscontinuousFunc<Ord>* ov = new DiscontinuousFunc<Ord>(get_fn_ord(fv->get_edge_fn_order(surf_pos->surf_num) + inc), neighbor_supp_v);

    // Order of additional external functions.
    ExtData<Ord>* fake_ext = init_ext_fns_ord(mfs->ext, neighbor_searches);

    // Order of geometric attributes (eg. for multiplication of a solution with coordinates, normals, etc.).
    Geom<Ord>* fake_e = new InterfaceGeom<Ord>(init_geom_ord(), nbs_u->neighb_el->marker, 
        nbs_u->neighb_el->id, nbs_u->neighb_el->get_diameter());
    double fake_wt = 1.0;

    // Total order of the matrix form.
    Ord o = mfs->ord(1, &fake_wt, oi, ou, ov, fake_e, fake_ext);

    // Increase due to reference maps.
    order = ru_actual->get_inv_ref_order();

    order += o.get_order();
    limit_order(order);

    // Clean up.
    for (int i = u_ext_offset; i < u_ext_length; i++) {
      if (oi[i - u_ext_offset] != NULL) { oi[i - u_ext_offset]->free_ord(); delete oi[i - u_ext_offset]; }
    }
    if (ou != NULL) {
      ou->free_ord(); delete ou;
    }
    if (ov != NULL) {
      ov->free_ord(); delete ov;
    }
  }

  // Evaluate the form using just calculated order.
  for(std::map<unsigned int, NeighborSearch>::iterator it = neighbor_searches.begin(); it != neighbor_searches.end(); it++)
    it->second.set_quad_order(order);

  Quad2D* quad = fu->get_quad_2d();
  int eo = quad->get_edge_points(surf_pos->surf_num, order);
  int np = quad->get_num_points(eo);
  double3* pt = quad->get_points(eo);

  // Init geometry and jacobian*weights.
  // TODO: A cache.
  assert(surf_pos->surf_num = neighbor_searches.at(neighbor_index_u).active_edge);

  Geom<double>* e = new InterfaceGeom<double>(init_geom_surf(ru_central, surf_pos, eo), nbs_u->neighb_el->marker, 
    nbs_u->neighb_el->id, nbs_u->neighb_el->get_diameter());

  double3* tan = ru_central->get_tangent(surf_pos->surf_num, eo);

  double *jwt = new double[np];
  for(int i = 0; i < np; i++)
    jwt[i] = pt[i][2] * tan[i][2];

  // Values of the previous Newton iteration, shape functions and external functions in quadrature points.
  AUTOLA_OR(Func<scalar>*, prev, u_ext_length - u_ext_offset);
  if (u_ext != Hermes::vector<Solution *>()) {
    for (int i = u_ext_offset; i < u_ext_length; i++) {
      if (u_ext[i] != NULL) prev[i - u_ext_offset]  = neighbor_searches.at(u_ext[i]->get_mesh()->get_seq()).init_ext_fn(u_ext[i]);
      else prev[i - u_ext_offset] = NULL;
    }
  }
  else {
    for (int i = u_ext_offset; i < u_ext_length; i++) prev[i - u_ext_offset] = NULL;
  }

  // Values of the previous Newton iteration, shape functions and external functions in quadrature points.
  DiscontinuousFunc<double>* u = new DiscontinuousFunc<double>(get_fn(fu, ru_actual, nbs_u->get_quad_eo(neighbor_supp_u)),
    neighbor_supp_u, nbs_u->neighbor_edges[nbs_u->active_segment].orientation);
  DiscontinuousFunc<double>* v = new DiscontinuousFunc<double>(get_fn(fv, rv, nbs_v->get_quad_eo(neighbor_supp_v)),
    neighbor_supp_v, nbs_v->neighbor_edges[nbs_v->active_segment].orientation);
  
  ExtData<scalar>* ext = init_ext_fns(mfs->ext, neighbor_searches);

  scalar res = mfs->fn(nbs_v->get_quad_np(), jwt, prev, u, v, e, ext);

  // Clean up.
  for (int i = u_ext_offset; i < u_ext_length; i++) {
    if (prev[i - u_ext_offset] != NULL) {
      prev[i - u_ext_offset]->free_fn(); delete prev[i - u_ext_offset];
    }
  }
  if (ext != NULL) {ext->free(); delete ext;}

  // Delete the DiscontinuousFunctions. This does not clear their component functions (DiscontinuousFunc::free_fn()
  // must be called in order to do that) as they are contained in cache_fn and may be used by another form - they
  // will be cleared in DiscreteProblem::delete_cache.
  delete u;
  delete v;

  // Scaling.
  res *= mfs->scaling_factor;

  return 0.5 * res; // Edges are parametrized from 0 to 1 while integration weights
                    // are defined in (-1, 1). Thus multiplying with 0.5 to correct
                    // the weights.
}


// Actual evaluation of surface linear form, just in case using information from neighbors.
// Used only for inner edges.
scalar DiscreteProblem::eval_dg_form(WeakForm::VectorFormSurf* vfs, Hermes::vector<Solution *> u_ext,
                                     PrecalcShapeset *fv, RefMap *rv, 
                                     SurfPos* surf_pos, std::map<unsigned int, NeighborSearch>& neighbor_searches, int neighbor_index_v)
{
  _F_
  NeighborSearch* nbs_v = &(neighbor_searches.at(neighbor_index_v));
  int order;
  int u_ext_length = u_ext.size();      // Number of external solutions.
  int u_ext_offset = vfs->u_ext_offset; // External solutions will start with u_ext[u_ext_offset]
                                        // and there will be only u_ext_length - u_ext_offset of them.
  if(this->is_fvm)
    order = rv->get_inv_ref_order();
  else {
    // Order of solutions from the previous Newton iteration.
    AUTOLA_OR(Func<Ord>*, oi, u_ext_length - u_ext_offset);
    if (u_ext != Hermes::vector<Solution *>())
      for (int i = u_ext_offset; i < u_ext_length; i++)
        if (u_ext[i] != NULL) 
          oi[i - u_ext_offset] = neighbor_searches.at(u_ext[i]->get_mesh()->get_seq()).init_ext_fn_ord(u_ext[i]);
        else 
          oi[i - u_ext_offset] = init_fn_ord(0);
    else
      for (int i = u_ext_offset; i < u_ext_length; i++) 
        oi[i - u_ext_offset] = init_fn_ord(0);

    // Order of the shape function.
    // Determine the integration order.
    int inc = (fv->get_num_components() == 2) ? 1 : 0;
    Func<Ord>* ov = get_fn_ord(fv->get_edge_fn_order(surf_pos->surf_num) + inc);

    // Order of additional external functions.
    ExtData<Ord>* fake_ext = init_ext_fns_ord(vfs->ext, neighbor_searches);

    // Order of geometric attributes (eg. for multiplication of a solution with coordinates, normals, etc.).
    Geom<Ord>* fake_e = new InterfaceGeom<Ord>(init_geom_ord(),
                        nbs_v->neighb_el->marker, nbs_v->neighb_el->id, nbs_v->neighb_el->get_diameter());
    double fake_wt = 1.0;

    // Total order of the vector form.
    Ord o = vfs->ord(1, &fake_wt, oi, ov, fake_e, fake_ext);

    // Increase due to reference map.
    order = rv->get_inv_ref_order();
    order += o.get_order();
    limit_order(order);

    // Clean up.
    for (int i = u_ext_offset; i < u_ext_length; i++) {
      if (oi[i - u_ext_offset] != NULL) { oi[i - u_ext_offset]->free_ord(); delete oi[i - u_ext_offset]; }
    }
    if (ov != NULL) {
      ov->free_ord(); delete ov;
    }
  }

  // Evaluate the form.
  for(std::map<unsigned int, NeighborSearch>::iterator it = neighbor_searches.begin(); it != neighbor_searches.end(); it++)
    it->second.set_quad_order(order);

  Quad2D* quad = fv->get_quad_2d();
  int eo = quad->get_edge_points(surf_pos->surf_num, order);
  int np = quad->get_num_points(eo);
  double3* pt = quad->get_points(eo);

  // Init geometry and jacobian*weights.
  // TODO: A cache.
  assert(surf_pos->surf_num = nbs_v->active_edge);

  Geom<double>* e = new InterfaceGeom<double>(init_geom_surf(rv, surf_pos, eo), nbs_v->neighb_el->marker, 
    nbs_v->neighb_el->id, nbs_v->neighb_el->get_diameter());

  double3* tan = rv->get_tangent(surf_pos->surf_num, eo);

  double *jwt = new double[np];
  for(int i = 0; i < np; i++)
    jwt[i] = pt[i][2] * tan[i][2];

  // Values of the previous Newton iteration, shape functions and external functions in quadrature points.
  AUTOLA_OR(Func<scalar>*, prev, u_ext_length - u_ext_offset);
  if (u_ext != Hermes::vector<Solution *>()) {
    for (int i = u_ext_offset; i < u_ext_length; i++) {
      if (u_ext[i] != NULL) prev[i - u_ext_offset]  = neighbor_searches.at(u_ext[i]->get_mesh()->get_seq()).init_ext_fn(u_ext[i]);
      else prev[i - u_ext_offset] = NULL;
    }
  }
  else {
    for (int i = u_ext_offset; i < u_ext_length; i++) prev[i - u_ext_offset] = NULL;
  }

  Func<double>* v = get_fn(fv, rv, nbs_v->get_quad_eo());
  ExtData<scalar>* ext = init_ext_fns(vfs->ext, neighbor_searches);

  scalar res = vfs->fn(np, jwt, prev, v, e, ext);

  // Clean up.
  for (int i = u_ext_offset; i < u_ext_length; i++) {
    if (prev[i - u_ext_offset] != NULL) {
      prev[i - u_ext_offset]->free_fn(); delete prev[i - u_ext_offset];
    }
  }
  if (ext != NULL) {ext->free(); delete ext;}

  // Scaling.
  res *= vfs->scaling_factor;

  return 0.5 * res; // Edges are parametrized from 0 to 1 while integration weights
                    // are defined in (-1, 1). Thus multiplying with 0.5 to correct
                    // the weights.
}


DiscreteProblem::NeighborNode::NeighborNode(NeighborNode* parent, unsigned int transformation) : parent(parent), transformation(transformation)
{
  left_son = right_son = NULL;
}
DiscreteProblem::NeighborNode::~NeighborNode()
{
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

double get_l2_norm(Vector* vec)
{
  _F_
  scalar val = 0;
  for (unsigned int i = 0; i < vec->length(); i++) {
    scalar inc = vec->get(i);
    val = val + inc*conj(inc);
  }
  return sqrt(std::abs(val));
}

// Performs uniform global refinement of a FE space.
Hermes::vector<Space *> * construct_refined_spaces(Hermes::vector<Space *> coarse, int order_increase)
{
  _F_
  Hermes::vector<Space *> * ref_spaces = new Hermes::vector<Space *>;
  for (unsigned int i = 0; i < coarse.size(); i++)
  {
    Mesh* ref_mesh = new Mesh;
    ref_mesh->copy(coarse[i]->get_mesh());
    ref_mesh->refine_all_elements();
    ref_spaces->push_back(coarse[i]->dup(ref_mesh, order_increase));
  }

  return ref_spaces;
}

// Light version for a single space.
Space* construct_refined_space(Space* coarse, int order_increase)
{
  _F_
  Mesh* ref_mesh = new Mesh;
  ref_mesh->copy(coarse->get_mesh());
  ref_mesh->refine_all_elements();
  Space* ref_space = coarse->dup(ref_mesh, order_increase);

  return ref_space;
}

bool solve_newton(scalar* coeff_vec, DiscreteProblem* dp, Solver* solver, SparseMatrix* matrix,
                  Vector* rhs, double newton_tol, int newton_max_iter, bool verbose,
                  bool residual_as_function,
                  double damping_coeff, double max_allowed_residual_norm)
{
  // Prepare solutions for measuring residual norm.
  int num_spaces = dp->get_num_spaces();
  Hermes::vector<Solution*> solutions;
  Hermes::vector<bool> dir_lift_false;
  for (int i=0; i < num_spaces; i++) {
    solutions.push_back(new Solution());
    dir_lift_false.push_back(false);      // No Dirichlet lifts will be considered.
  }

  // The Newton's loop.
  double residual_norm;
  int it = 1;
  while (1)
  {
    // Obtain the number of degrees of freedom.
    int ndof = dp->get_num_dofs();

    // Assemble the Jacobian matrix and residual vector.
    dp->assemble(coeff_vec, matrix, rhs, false);

    // Multiply the residual vector with -1 since the matrix
    // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
    rhs->change_sign();

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
    if (verbose) info("---- Newton iter %d, ndof %d, residual norm %g", it, ndof, residual_norm);

    // If maximum allowed residual norm is exceeded, fail.
    if (residual_norm > max_allowed_residual_norm) {
      if (verbose) {
        info("Current residual norm: %g", residual_norm);
        info("Maximum allowed residual norm: %g", max_allowed_residual_norm);
        info("Newton solve not successful, returning false.");
      }
      return false;
    }

    // If residual norm is within tolerance, or the maximum number
    // of iteration has been reached, then quit.
    if ((residual_norm < newton_tol || it > newton_max_iter) && it > 1) break;

    // Solve the linear system.
    if(!solver->solve()) error ("Matrix solver failed.\n");

    // Add \deltaY^{n+1} to Y^n.
    for (int i = 0; i < ndof; i++) coeff_vec[i] += damping_coeff * solver->get_solution()[i];

    it++;
  }

  if (it >= newton_max_iter) {
    if (verbose) info("Maximum allowed number of Newton iterations exceeded, returning false.");
    return false;
  }

  return true;
}

// Perform Picard's iteration.
bool solve_picard(WeakForm* wf, Space* space, Solution* sln_prev_iter,
                  MatrixSolverType matrix_solver, double picard_tol,
                  int picard_max_iter, bool verbose)
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
                 iter_count+1, Space::get_num_dofs(space), rel_error);

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