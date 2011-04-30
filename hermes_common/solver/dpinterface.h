#ifndef DPINTERFACE_H
#define DPINTERFACE_H

//#include "../matrix.h"
//#include "../vector.h"
//#include "../tables.h"

class SparseMatrix;
class Vector;
class Table;

/// Minimalistic DiscreteProblem interface required by NoxProblemInterface.
class DiscreteProblemInterface
{
public:
  /// Get the number of unknowns.
  virtual int get_num_dofs() = 0;
  
  /// Get info about presence of a matrix.
  virtual bool is_matrix_free() = 0;
  
  /// Preassembling.
  /// Precalculate matrix sparse structure.
  /// If force_diagonal_block == true, then (zero) matrix
  /// antries are created in diagonal blocks even if corresponding matrix weak
  /// forms do not exist. This is useful if the matrix is later to be merged with
  /// a matrix that has nonzeros in these blocks. The Table serves for optional
  /// weighting of matrix blocks in systems.
  virtual void create_sparse_structure(SparseMatrix* mat, Vector* rhs = NULL,
                                       bool force_diagonal_blocks = false, Table* block_weights = NULL) = 0;

  /// Assembling.
  /// General assembling procedure for nonlinear problems. coeff_vec is the
  /// previous Newton vector. If force_diagonal_block == true, then (zero) matrix
  /// antries are created in diagonal blocks even if corresponding matrix weak
  /// forms do not exist. This is useful if the matrix is later to be merged with
  /// a matrix that has nonzeros in these blocks. The Table serves for optional
  /// weighting of matrix blocks in systems. The parameter add_dir_lift decides 
  /// whether Dirichlet lift will be added while coeff_vec is converted into 
  /// Solutions.
  virtual void assemble(scalar* coeff_vec, SparseMatrix* mat, Vector* rhs = NULL,
             bool force_diagonal_blocks = false, bool add_dir_lift = true, Table* block_weights = NULL) = 0; 
                
  virtual void invalidate_matrix() = 0;
  
};

#endif // DPINTERFACE_H
