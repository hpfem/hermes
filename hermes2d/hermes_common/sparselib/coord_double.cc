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
/*              Coordinate sparse matrix (0-based)                       */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include <iostream>
#include <stdlib.h>
#include "compcol_double.h"
#include "comprow_double.h"
#include "coord_double.h"

#include "spblas.h"

/*****************************/
/*  Constructor(s)           */
/*****************************/

Coord_Mat_double::Coord_Mat_double(void) : 
        val_(0), rowind_(0), colind_(0), base_(0), nz_(0)
{
        dim_[0] = 0;
        dim_[1] = 0;
}

/*****************************/
/*  Copy constructor         */
/*****************************/

Coord_Mat_double::Coord_Mat_double(const Coord_Mat_double &S) : 
        val_(S.val_), rowind_(S.rowind_), colind_(S.colind_), base_(S.base_), 
        nz_(S.nz_)
{
        dim_[0] = S.dim_[0];
        dim_[1] = S.dim_[1];
}

/***********************************/
/*  Construct from storage vectors */
/***********************************/

Coord_Mat_double::Coord_Mat_double(int M, int N, int nz, double *val, int *r, 
                                 int *c, int base) : 
        val_(val, nz), rowind_(r, nz), colind_(c, nz), base_(base), nz_(nz)
{
        dim_[0] = M;
        dim_[1] = N;
}


/**********************************************************************/
/*  Construct a Coord_Mat_double from a CompRow_Mat_double              */
/*  (i.e. convert compressed row storage to coordinate storage)       */
/**********************************************************************/

Coord_Mat_double::Coord_Mat_double(const CompRow_Mat_double &R) :
        val_(R.NumNonzeros()), rowind_(R.NumNonzeros()), 
        colind_(R.NumNonzeros()), base_(R.base()), nz_(R.NumNonzeros())
{


        dim_[0] = R.dim(0);
        dim_[1] = R.dim(1);


        int count = 0;
        int i, j;
//      Loop through rows...
        for (i=1;i<=dim_[0];i++)
        {
           for (j=count;j<R.row_ptr(i);j++)
           {
              val_[count] = R.val(count);
              colind_[count] = R.col_ind(count);
              rowind_[count] = i-1;
              count++;
           }
       }
}
          

/**********************************************************************/
/*  Construct a Coord_Mat_double from a CompCol_Mat_double              */
/*  (i.e. convert compressed column storage to coordinate storage)    */
/**********************************************************************/
              
Coord_Mat_double::Coord_Mat_double(const CompCol_Mat_double &C) :
        val_(C.NumNonzeros()), rowind_(C.NumNonzeros()), 
        colind_(C.NumNonzeros()), base_(C.base()),  nz_(C.NumNonzeros())
{

        dim_[0] = C.dim(0);
        dim_[1] = C.dim(1);

        int count = 0;
        int i, j;
//      Loop through columns...
        for (j=1;j<=dim_[1];j++)
        {
          for (i=count;i<C.col_ptr(j);i++)
          {
              val_[count] = C.val(count);
              rowind_[count] = C.row_ind(count);
              colind_[count] = j-1;
              count++;
           }
         }
}


/***************************/
/* Assignment operator...  */
/***************************/

Coord_Mat_double& Coord_Mat_double::operator=(const Coord_Mat_double &C)
{
         dim_[0] = C.dim_[0];
         dim_[1] = C.dim_[1];
         base_ = C.base_;
         nz_ = C.nz_;
         val_ = C.val_;
         rowind_ = C.rowind_;
         colind_ = C.colind_;
         return *this;
}

/***************************/
/* newsize()               */
/***************************/

Coord_Mat_double& Coord_Mat_double::newsize(int M, int N, int nz)
{
         dim_[0] = M;
         dim_[1] = N;
         nz_ = nz;
         val_.newsize(nz);
         rowind_.newsize(nz);
         colind_.newsize(nz);
         return *this;
}


/*********************/
/*   Array access    */
/*********************/

double Coord_Mat_double::operator()(int i, int j)  const
{
         for (int t=0; t<nz_; t++)
            if (rowind_(t) == i && colind_(t) == j) return val_(t);
         if (i < dim_[0] && j < dim_[1]) return 0.0;
         std::cerr << "Array accessing exception -- out of bounds." << "\n";
         exit(1);
         return 0;          // return to suppress compiler warning message
}


double& Coord_Mat_double::set(int i, int j)
{
         for (int t=0; t<nz_; t++)
            if (rowind_(t) == i && colind_(t) == j) return val_(t);
         std::cerr << "Array element (" << i << "," << j ;
         std::cerr << ") not in sparse structure -- cannot assign." << "\n";
         exit(1);
         return val_(0);    // return to suppress compiler warning message
}

/*************/
/*   I/O     */
/*************/
using namespace std;
ostream& operator << (ostream & os, const Coord_Mat_double & mat)
{
         int nnz = mat.NumNonzeros();
         int rowp1, colp1;
         int M = mat.dim(0);
         int N = mat.dim(1);
         int flag = 0;


        std::ios::fmtflags olda = os.setf(ios::right,ios::adjustfield);
        std::ios::fmtflags oldf = os.setf(ios::scientific,ios::floatfield);

         int oldp = os.precision(12);

//       Loop through Nonzeros
         for (int j = 0; j < nnz ; j++) 
         {
            rowp1 = mat.row_ind(j) +1;
            colp1 = mat.col_ind(j) +1;
            if ( rowp1 == M && colp1 == N ) flag = 1;
            os.width(14);
            os <<  rowp1 ; os << "    " ;
            os.width(14);
            os <<  colp1 ; os << "    " ;
            os.width(20);
            os <<  mat.val(j) << "\n";
         }
         if (flag == 0) 
         {
            os.width(14);
            os <<  M ; os << "    " ;
            os.width(14);
            os <<  N ; os << "    " ;
            os.width(20);
            os <<  mat(M-1,N-1) << "\n";
         }
         os.setf(olda,ios::adjustfield);
         os.setf(oldf,ios::floatfield);
         os.precision(oldp);
         return os;
}

/***************************************/
/* Matrix-MV_Vector multiplication...  */
/***************************************/

VECTOR_double Coord_Mat_double::operator*(const VECTOR_double &x) const
{
         int M = dim(0);
         int N = dim(1);

//       Check for compatible dimensions:
         if (x.size() != N) 
         {
            std::cerr << "Error in Coord Matvec -- incompatible dimensions." 
                 << "\n";
            exit(1);
            return x;        // return to suppress compiler warning message.
         }
  
         VECTOR_double result(M, 0.0);
         VECTOR_double work(M);
  
         int descra[9];
         descra[0] = 0;
         descra[1] = 0;
         descra[2] = 0;
         descra[4] = 0;

         F77NAME(dcoomm) (0, M, 1, N, 1.0,
                  descra, &val_(0), &rowind_(0), &colind_(0), nz_,
                  &x(1), N, 1.0, &result(1), M,
                  &work(0), M);

         return result;
}

/*************************************************/
/* Matrix-Transpose-MV_Vector multiplication...  */
/*************************************************/

VECTOR_double Coord_Mat_double::trans_mult(const VECTOR_double &x) const
{
         int M = dim(0);
         int N = dim(1);

//       Check for compatible dimensions:
         if (x.size() != M) 
         {
            std::cerr << "Error in Coord trans_mult -- incompatible dimensions."  
                 << "\n";
            exit(1);
            return x;        // return to suppress compiler warning message.
         }

         VECTOR_double result(N, 0.0);
         VECTOR_double work(N);
 
         int descra[9];
         descra[0] = 0;
         descra[1] = 0;
         descra[2] = 0;
         descra[4] = 0;
  
         F77NAME(dcoomm) (1, N, 1, M, 1.0,
                          descra, &val_(0), &rowind_(0), &colind_(0), nz_,
                          &x(1), M, 1.0, &result(1), N,
                          &work(0), N);
         return result;
}

