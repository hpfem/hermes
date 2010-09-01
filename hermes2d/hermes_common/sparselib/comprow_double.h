/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*             ********   ***                                 SparseLib++    */
/*          *******  **  ***       ***      ***                              */
/*           *****      ***     ******** ********                            */
/*            *****    ***     ******** ********              R. Pozo        */
/*       **  *******  ***   **   ***      ***                 K. Remington   */
/*        ********   ********                                 A. Lumsdaine   */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                                                                           */
/*                                                                           */
/*                     SparseLib++ : Sparse Matrix Library                   */
/*                                                                           */
/*               National Institute of Standards and Technology              */
/*                        University of Notre Dame                           */
/*              Authors: R. Pozo, K. Remington, A. Lumsdaine                 */
/*                                                                           */
/*                                 NOTICE                                    */
/*                                                                           */
/* Permission to use, copy, modify, and distribute this software and         */
/* its documentation for any purpose and without fee is hereby granted       */
/* provided that the above notice appear in all copies and supporting        */
/* documentation.                                                            */
/*                                                                           */
/* Neither the Institutions (National Institute of Standards and Technology, */
/* University of Notre Dame) nor the Authors make any representations about  */
/* the suitability of this software for any purpose.  This software is       */
/* provided ``as is'' without expressed or implied warranty.                 */
/*                                                                           */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*             Compressed row sparse matrix (0-based, Fortran)           */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


#ifndef CompRow_Mat_double_H
#define CompRow_Mat_double_H


#include "vecdefs.h"
#include VECTOR_H

class CompCol_Mat_double;
class Coord_Mat_double;

class CompRow_Mat_double {

private:
       VECTOR_double val_;       // data values (nz_ elements)
       VECTOR_int    rowptr_;    // row_ptr (dim_[0]+1 elements)
       VECTOR_int    colind_;    // col_ind  (nz_ elements)

       int base_;                 // index base: offset of first element
       int nz_;                   // number of nonzeros
       int dim_[2];               // number of rows, cols
  
public:
       CompRow_Mat_double(void);
       CompRow_Mat_double(const CompRow_Mat_double &S);
       CompRow_Mat_double(const CompCol_Mat_double &C);
       CompRow_Mat_double(const Coord_Mat_double &CO);
       CompRow_Mat_double(int M, int N, int nz, double *val, int *r, 
                               int *c, int base=0);
       CompRow_Mat_double(int M, int N, int nz, const VECTOR_double &val, 
                               const VECTOR_int &r, const VECTOR_int &c,
                               int base=0);
      ~CompRow_Mat_double() {};
    
/*******************************/
/*  Access and info functions  */
/*******************************/

       double&      val(int i) { return val_(i); }
       int&         row_ptr(int i) { return rowptr_(i); }
       int&         col_ind(int i) { return colind_(i);}

       const double&      val(int i) const { return val_(i); }
       const int&         row_ptr(int i) const { return rowptr_(i); }
       const int&         col_ind(int i) const { return colind_(i);}

       int          dim(int i) const {return dim_[i];};
       int          size(int i) const {return dim_[i];};
       int          NumNonzeros() const {return nz_;};
       int          base() const {return base_;}

/*******************************/
/*  Assignment  operator       */
/*******************************/

       CompRow_Mat_double& operator=(const CompRow_Mat_double &R);
       CompRow_Mat_double& newsize(int M, int N, int nz);

/***********************************/
/*  General access function (slow) */
/***********************************/

       double       operator() (int i, int j) const;        
       double&      set(int i, int j);

/***********************************/
/*  Matrix/Vector multiply         */
/***********************************/

       VECTOR_double operator*(const VECTOR_double &x) const;
       VECTOR_double trans_mult(const VECTOR_double &x) const;

};

std::ostream& operator << (std::ostream & os, const CompRow_Mat_double & mat);

#endif  
/* CompRow_Mat_double_H */

