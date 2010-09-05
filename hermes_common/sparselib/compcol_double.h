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
/*        Compressed column sparse matrix  (O-based, Fortran)            */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


#ifndef CompCol_Mat_double_H
#define CompCol_Mat_double_H

#include "vecdefs.h"
#include VECTOR_H

class CompRow_Mat_double;
class Coord_Mat_double;

class CompCol_Mat_double {

private:
       VECTOR_double     val_;       // data values (nz_ elements)
       VECTOR_int    rowind_;    // row_ind (nz_ elements)
       VECTOR_int    colptr_;    // col_ptr (dim_[1]+1 elements)

       int base_;                 // index base: offset of first element
       int nz_;                   // number of nonzeros
       int dim_[2];               // number of rows, cols
  
public:
       CompCol_Mat_double(void);
       CompCol_Mat_double(const CompCol_Mat_double &S);
       CompCol_Mat_double(const CompRow_Mat_double &R);
       CompCol_Mat_double(const Coord_Mat_double &CO);
       CompCol_Mat_double(int M, int N, int nz, double *val, int *r, 
                               int *c, int base=0);
       CompCol_Mat_double(int M, int N, int nz, const VECTOR_double &val, 
                               const VECTOR_int &r, const VECTOR_int &c, 
                               int base=0);
      ~CompCol_Mat_double() {};
        
/*******************************/
/*  Access and info functions  */
/*******************************/

       double&      val(int i) { return val_(i); }
       int&         row_ind(int i) { return rowind_(i); }
       int&         col_ptr(int i) { return colptr_(i);}

       const double&      val(int i) const { return val_(i); }
       const int&         row_ind(int i) const { return rowind_(i); }
       const int&         col_ptr(int i) const { return colptr_(i);}

       int          dim(int i) const {return dim_[i];};
       int          size(int i) const {return dim_[i];};
       int          NumNonzeros() const {return nz_;};
       int          base() const {return base_;}

/*******************************/
/*  Assignment  operator       */
/*******************************/

       CompCol_Mat_double& operator=(const CompCol_Mat_double &C);
       CompCol_Mat_double& newsize(int M, int N, int nz);

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

std::ostream& operator << (std::ostream & os, const CompCol_Mat_double & mat);
void readHB(const char *c, CompCol_Mat_double &M);

#endif  
/* CompCol_Mat_double_H */

