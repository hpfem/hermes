
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                                                                           */
/*                                                                           */
/*                   MV++ Numerical Matrix/Vector C++ Library                */
/*                             MV++ Version 1.5                              */
/*                                                                           */
/*                                  R. Pozo                                  */
/*               National Institute of Standards and Technology              */
/*                                                                           */
/*                                  NOTICE                                   */
/*                                                                           */
/* Permission to use, copy, modify, and distribute this software and         */
/* its documentation for any purpose and without fee is hereby granted       */
/* provided that this permission notice appear in all copies and             */
/* supporting documentation.                                                 */
/*                                                                           */
/* Neither the Institution (National Institute of Standards and Technology)  */
/* nor the author makes any representations about the suitability of this    */
/* software for any purpose.  This software is provided ``as is''without     */
/* expressed or implied warranty.                                            */
/*                                                                           */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

//
//          Basic matrix class (int precision)
//

#include <iostream>                                 
#include "mvmi.h"

 int MV_ColMat_int::dim(int i) const 
{
    if (i==0) return dim0_;
    if (i==1) return dim1_;
    else
    {
     std::cerr << "Called MV_ColMat::dim(" << i << ")  must be 0 or 1 " << "\n";
     exit(1);
    }

    // never should be here, but many compilers warn about not
    // returning a value
    return 0;
}

// NOTE: null construct have ref_ flag turned OFF, otherwise, we can
//          never reset the dim of matrix....
MV_ColMat_int::MV_ColMat_int()  : v_(), dim0_(0), dim1_(0) , lda_(0), ref_(0){}
                                                                

MV_ColMat_int::MV_ColMat_int( int m,  int n) : v_(m*n),
        dim0_(m), dim1_(n), lda_(m), ref_(0) {}

MV_ColMat_int::MV_ColMat_int( int m,  int n, const int &s) : v_(m*n),
        dim0_(m), dim1_(n), lda_(m), ref_(0) 
{
    operator=(s);
}

// operators and member functions






MV_ColMat_int& MV_ColMat_int::operator=(const int & s) 
{
    int M = dim(0);
    int N = dim(1);

    if (lda_ == M)      // if continuous, then just assign as a ?
        v_ =  s;        // single long vector.

    else
    {
    // this should run much faster than the just accessing each (i,j)
    // element individually 
    //

        MV_VecIndex I(0,M-1);
        for (int j=0; j<N; j++)
        {
            v_(I) = s;
            I += lda_;
        }
    }

    return *this;
}

MV_ColMat_int& MV_ColMat_int::newsize( int M,  int N)
{
    v_.newsize(M*N);
    dim0_ = M;
    dim1_ = N;
    lda_ = M;

    return *this;
}

MV_ColMat_int& MV_ColMat_int::operator=(const MV_ColMat_int & m) 
{

    int lM = dim0_;     // left hand arg  (this)
    int lN = dim1_;

    int rM = m.dim0_;   // right hand arg (m)
    int rN = m.dim1_;


    // if the left-hand side is a matrix reference, the we copy the
    // elements of m *into* the region specfied by the reference.
    // i.e. inject().

    if (ref_)
    {
        // check conformance,       
        if (lM != rM  || lN != rN)      
        {
            std::cerr << "MV_ColMatRef::operator=  non-conformant assignment.\n";
            exit(1);
        }
    }
    else
    {
        newsize(rM,rN);
    }

    // at this point the left hand and right hand sides are conformant

    // this should run much faster than the just accessing each (i,j)
    // element individually 

    // if both sides are contigous, then just copy as one vector
    if ( lM == lda_ && rM == m.lda_)
    {
        MV_VecIndex I(0,rM*rN-1);
        v_(I) = m.v_(I);
    }
    else
    {
        // slower way...

        MV_VecIndex I(0,rM-1);
        MV_VecIndex K(0,rM-1);
        for (int j=0; j<rN; j++)
        {
            v_(I) = m.v_(K);
            I += lda_;
            K += m.lda_;
        }
    }

    return *this;   
}

MV_ColMat_int::MV_ColMat_int(const MV_ColMat_int & m) : 
        v_(m.dim0_*m.dim1_), dim0_(m.dim0_),
        dim1_(m.dim1_), lda_(m.dim0_), ref_(0)
{

    int M = m.dim0_;
    int N = m.dim1_;

    // this should run much faster than the just accessing each (i,j)
    // element individually 

    MV_VecIndex I(0,M-1);
    MV_VecIndex K(0,M-1);
    for (int j=0; j<N; j++)
    {
        v_(I) = m.v_(K);
        I += lda_;
        K += m.lda_;
    }
}



MV_ColMat_int::MV_ColMat_int(int* d,  int m,  int n) :
    v_(m*n), dim0_(m), dim1_(n), lda_(m), ref_(0)
{
    int mn = m*n;

    // d is contiguous, so just copy 1-d vector
    for (int i=0; i< mn; i++)
            v_[i] = d[i];
}


MV_ColMat_int::MV_ColMat_int(int* d,  int m,  int n, 
         int lda) :
    v_(m*n), dim0_(m), dim1_(n), lda_(lda), ref_(0)
{
    for (int j=0; j< n; j++)
        for (int i=0; i<m; i++)
            operator()(i,j) = d[j*lda + i];   // could be made faster!!
}


MV_ColMat_int MV_ColMat_int::operator()(const MV_VecIndex &I, const MV_VecIndex &J)
{
    // check that index is not out of bounds
    //
    if (I.end() >= dim0_  || J.end() >= dim1_)
    {
        std::cerr << "Matrix index: (" << I.start() << ":" << I.end()  
             << "," << J.start() << ":" << J.end()   
             << ") not a subset of (0:" << dim0_ - 1 << ", 0:" 
             << dim1_-1 << ") " << "\n";
        exit(1);
    }

    // this automatically returns a reference
    // 
    return MV_ColMat_int(&v_[J.start()*lda_ + I.start()], 
            I.end() - I.start() + 1, 
            J.end() - J.start() + 1, lda_, MV_Matrix_::ref);
}

const MV_ColMat_int MV_ColMat_int::operator()(const MV_VecIndex &I, 
    const MV_VecIndex &J) const
{

    std::cerr << "Const operator()(MV_VecIndex, MV_VecIndex) called " << "\n";

    // check that index is not out of bounds
    //
    if (I.end() >= dim0_  || J.end() >= dim1_)
    {
        std::cerr << "Matrix index: (" << I.start() << ":" << I.end()  
             << "," << J.start() << ":" << J.end()   
             << ") not a subset of (0:" << dim0_ - 1 << ", 0:" 
             << dim1_-1 << ") " << "\n";
        exit(1);
    }

    // this automatically returns a reference.  we need to 
    // "cast away" constness here, so the &v_[] arg will
    // not cause a compiler error.
    //
    MV_ColMat_int *t =  (MV_ColMat_int*) this;
    return MV_ColMat_int(&(t->v_[J.start()*lda_ + I.start()]), 
            I.end() - I.start() + 1, 
            J.end() - J.start() + 1, lda_, MV_Matrix_::ref);
}

MV_ColMat_int::~MV_ColMat_int() {}

std::ostream&   operator<<(std::ostream& s, const MV_ColMat_int& V)
{
    int M = V.dim(0);
    int N = V.dim(1);

    for (int i=0; i<M; i++)
    {
        for (int j=0; j<N; j++)
            s << V(i,j) << " " ;
        s << "\n";
    }
    
    return s;
}



