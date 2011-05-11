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

#include "nurbs.h"

/// Constructor for the NURBS class.
NurbsCurve::NurbsCurve(int deg, int CPtMaxIndx, double Px[], double Py[], double w[], int KntMaxIndx, double Uvals[]) {

  p = deg;
  n = CPtMaxIndx;
  m = KntMaxIndx;

  U = new double[m+1];

  for (int j = 0; j < m+1; j++)
    U[j] = Uvals[j];

  P = new HPoint2D[n+1];

  for (int i = 0; i < n+1; i++) {
    P[i].wx = Px[i]*w[i];
    P[i].wy = Py[i]*w[i];
    P[i].w = w[i];

  }
}

/// Destructor for the NURBS class.
  NurbsCurve::~NurbsCurve() {

  delete [] U;
  delete [] P;

 }

///  Determines the knot span idex.
///  Algorithm A2.1 in The NURBS book p68.
///  Parameter:  
///    u - Control parameter. Determines where to evaluate knot span index.
///  Returns: 
///    The knot span index while \f$ u < U[mid] \lor u >= U[mid+1] \f$. 
///    That is \f$i\f$ if \f$ u_i <= u < u_{i+1} \f$.
///                                                                     
  int NurbsCurve::FindSpan(double u) const {
 
  if ( u >= U[n+1]) 
    return n;
  if ( u<= U[p]) 
    return p;

  int low = p;
  int high = n+1;
  int mid = (int) (low + high)/2;

    while (u < U[mid] || u >= U[mid+1]) {
      if ( u < U[mid] ) 
        high = mid;
      else
        low = mid;
      mid = (int) (low + high)/2;
    } 

  return mid;

  }

/// Computes the nonvanishing basis functions.
/// Algorithm A2.2 of The NURBS book, p70.
/// Based on Eq. 2.5 and Cox-deBoor-Mansfield reccurence algorithm.
/// Parameters: 
///   span - Specifies the span at which we compute basis functions.
///   u - Control parameter value. 
/// Output: 
///   N - Non-zero basis functions in the array \f$ N[0],...N[p] \f$ of length p+1.
///   This array has to be allocated in the calling function.
  void NurbsCurve::BasisFuns(int span, double u, double N[]) const {

  int i = span;
  double *left = new double[2*(p+1)];
  double *right = &left[p+1];

  double temp, saved;

  N[0] = 1.0;

  for(int j=1; j <= p; j++) {
    left[j] = u-U[i+1-j];
    right[j] = U[i+j]-u;
    saved = 0.0;
    for(int r = 0; r < j; r++) {
      temp = N[r]/(right[r+1]+left[j-r]);
      N[r] = saved+right[r+1] * temp;
      saved = left[j-r] * temp;
    }
    N[j] = saved;
  }

  delete [] left;

}

/// Computes the nonvanishing basis functions.
/// Based on Eq. 2.5 and Cox-deBoor-Mansfield reccurence algorithm.
/// Algorithm A2.2 of The NURBS book, p70. 
/// Parameters: 
///   span - Specifies the span at which we compute basis functions. 
///   u - Control parameter value. 
///   deg - Degree of the basis function we compute.  
///   U - The knot vector.
/// Output: 
///   N - Non-zero basis functions in the array \f$ N[0],...N[deg] \f$ of length deg+1.
///       This array has to be allocated in the calling function. 
///       To relate N to the basis functions: 
///       \f$ N_{all}[span -deg +i] = N[i] \f$, where \f$ i=0... deg \f$.
  void NurbsCurve::nurbsBasisFuns(int span, double u, int deg, double U[], double N[]) const {

  int i = span;
  double *left = new double[2*(deg+1)];
  double *right = &left[deg+1];

  double temp, saved;

  N[0] = 1.0;

  for(int j=1; j <= deg; j++) {
    left[j] = u-U[i+1-j];
    right[j] = U[i+j]-u;
    saved = 0.0;
    for(int r = 0; r < j; r++) {
      temp = N[r]/(right[r+1]+left[j-r]);
      N[r] = saved+right[r+1] * temp;
      saved = left[j-r] * temp;
    }
    N[j] = saved;
  }

  delete [] left;

}

/// Computes the basis functions of the NURBS curve and their derivatives at parametric value "u". 
/// Based on Eq. 2.10 and Algorithm A2.3 in The NURBS book, p72. 
/// The result is stored in the ders matrix, where ders is of 
/// size [d+1][p+1] and the derivative \f$ N'_i(u) = ders(1,i=span-p+j) \f$ where \f$ j = 0...p+1 \f$.
/// Parameters:
///   d - The degree of the derivation. 
///   u - Control parameter value.
///   span - Specifies the span at which we compute basis functions and their derivatives. 
/// Output:
///   ders - A [d+1][p+1] matrix containing the derivatives of the curve.
///   This array has to be allocated in the calling function.
  void NurbsCurve::DersBasisFuns(int d, double u, int span, Matrix<double> &ders) const {

  double left[p+1];
  double right[p+1];
  

  double ndu[p+1][p+1];

  double saved,temp;
  int j, r;

  ndu[0][0] = 1.0;
  for(j=1; j<= p;j++) {
    left[j] = u-U[span+1-j];
    right[j] = U[span+j]-u;
    saved = 0.0;
    
    for(r=0; r<j; r++) {
      // Lower triangle.
      ndu[j][r] = right[r+1]+left[j-r];
      temp = ndu[r][j-1]/ndu[j][r];
      // Upper triangle.
      ndu[r][j] = saved+right[r+1] * temp;
      saved = left[j-r] * temp;
    }

    ndu[j][j] = saved;
  }

  for(j=0; j<=p; j++) // Load the basis functions.
    ders[0][j] = ndu[j][p];

  // Compute the derivatives Eq. 2.10 in The NURBS book.

  double a[2][p+1];

  for(r=0; r<=p; r++) {
    int s1,s2;
    s1 = 0; 
    s2 = 1; // Alternate rows in array a.
    a[0][0] = 1.0;
    // Compute the kth derivative.
    for(int k=1; k<=d; k++) {
      double d;
      int rk,pk,j1,j2;
      d = 0.0;
      rk = r-k; 
      pk = p-k;

      if(r>=k) {
	a[s2][0] = a[s1][0]/ndu[pk+1][rk];
	d = a[s2][0]*ndu[rk][pk];
      }

      if(rk>=-1) {
	j1 = 1;
      }
      else {
	j1 = -rk;
      }

      if(r-1 <= pk) {
	j2 = k-1;
      }
      else {
	j2 = p-r;
      }

      for(j=j1; j<=j2; j++) {
	a[s2][j] = (a[s1][j]-a[s1][j-1])/ndu[pk+1][rk+j];
	d += a[s2][j]*ndu[rk+j][pk];
      }
      
      if(r<=pk) {
	a[s2][k] = -a[s1][k-1]/ndu[pk+1][r];
	d += a[s2][k]*ndu[r][pk];
      }
      ders[k][r] = d;
      j = s1; 
      s1 = s2; 
      s2 = j; // Switch rows.
    }
  }

  // Multiply through by the correct factors.
  r = p;
  for(int k=1; k<=d; k++) {
    for(j=0; j<=p; j++)
      ders[k][j] *= r;
    r *= p-k;
  }

}

/// Computes the i-th basis function of degree p of the NURBS curve at parametric value u.
/// Algorithm A2.4 of The NURBS book, p74. 
/// Parameters:
///   u - Control parameter value.
///   i - Specifies the span at which we compute particular basis function.
///   p - the degree to which the basis function is computed
//  Returns:
///  \f$ N_{i,p}(u) \f$ - The value of i-th basis function of degree p at u.  
  double NurbsCurve::SingleBasisFun(double u, int i, int p) const {

  double Nip;
  double saved,Uleft,Uright,temp;

  if((i==0 && u == U[0]) || (i == m-p-1 && u==U[m])){
    Nip = 1.0;
    return Nip;
  }
  if(u < U[i] || u >= U[i+p+1]) {
    Nip = 0.0;
    return Nip;
  }

  double N[p+1];
  
  int j;
  // Inialize zeroth degree functions.
  for(j = 0; j <= p; j++) {
    if(u >= U[i+j] && u < U[i+j+1]) 
      N[j] = 1.0;
    else
      N[j] = 0.0;
  }

  // Compute triangular table.
  for(int k=1; k<=p ; k++){
    if(N[0] == 0.0)
      saved = 0.0;
    else
      saved = ((u-U[i])*N[0])/(U[i+k]-U[i]);
    for(j = 0;j < p-k+1; j++) {
      Uleft = U[i+j+1];
      Uright = U[i+j+k+1];

      if(N[j+1]==0.0) {
	N[j] = saved;
	saved = 0.0;
      }
      else {
	temp = N[j+1]/(Uright-Uleft);
	N[j] = saved+(Uright-u)*temp;
	saved = (u-Uleft)*temp;
      }
    }
  }
  Nip = N[0];

  return Nip;  
}

/// Simple modification BasisFuns function enable calculation of 
/// all nonzero basis functions of all degrees from 0 up to p.
/// Parameters: 
///   span - Specifies the span at which we compute basis functions. 
///   u - Control parameter value.
/// Output:
///   A matrix N, where N[j][i] is the j-th basis function of i-th degree. 
///   \f[ N_{span-i+j, i}(u), \; 0 <= i <= p, \; 0 <= j <= i \f].
  void NurbsCurve::AllBasisFuns(int span, double u, Matrix<double> &N) const {

 for (int i = 0; i<=p; i++) {
   double *Nb = new double[p+1];
   nurbsBasisFuns(span, u, p, U, Nb);
   for (int j = 0; j<=p; j++) {
     N[j][i] = Nb[j];
   }
   delete [] Nb;
 }

}

/// Computes point on NURBS curve. 
/// Algorithm A4.1 The NURBS book p124. 
/// Parameter: 
///   u - Parametric value we evaluate the curve at.
/// Returns: 
///   C(u) - a point of type Point2D It is found by projection 
///   from Cw a point in three-dimensional space with coordinates \f$ {wx, wy, w} \f$.
  Point2D NurbsCurve::CurvePoint(double u) const {

  int span = FindSpan(u);
  double N[p+1];
  BasisFuns(span, u, N);

  HPoint2D Cw;
  Cw.wx = 0.; Cw.wy = 0.; Cw.w = 0.;

  for (int i=0; i<=p; i++) {
    Cw.wx += N[i]*P[span-p+i].wx; 
    Cw.wy += N[i]*P[span-p+i].wy;
    Cw.w  += N[i]*P[span-p+i].w; 
  }

  // Project it to normal 2D space.
  Point2D C;
  C.x = Cw.wx/Cw.w;
  C.y = Cw.wy/Cw.w;

  return C;

}

/// Computes point on NURBS curve.
/// Algorithm A4.1 The NURBS book p124.  
/// Parameter: 
///   u - Parametric value we evaluate the curve at.
///   span - Specifies the span at which we compute basis functions.
/// Returns: 
///   C(u) - a point of type Point2D. It is found by projection 
///   from Cw, a point in three-dimensional space with coordinates \f$ {wx, wy, w} \f$.
  Point2D NurbsCurve::CurvePoint(double u, int span, double N[]) const {

  HPoint2D Cw;
  Cw.wx = 0.; Cw.wy = 0.; Cw.w = 0.;

  for (int i=0; i<=p; i++) {
    Cw.wx += N[i]*P[span-p+i].wx; 
    Cw.wy += N[i]*P[span-p+i].wy;
    Cw.w  += N[i]*P[span-p+i].w; 
  }

  // Project it to normal 2D space.
  Point2D C;
  C.x = Cw.wx/Cw.w;
  C.y = Cw.wy/Cw.w;

  return C;

}

/// Computes point on rational B-spline curve.
/// Algorithm A4.1 The NURBS book p124.  
/// Parameter: 
///   u - Parametric value we evaluate the curve at.
///   span - Specifies the span at which we compute basis functions.
/// Returns: 
///   Cw(u) - a point of type HPoint2D. 
///   It is a point in three-dimensional space with homogenous coordinates \f$ {wx, wy, w} \f$.
  HPoint2D NurbsCurve::CurvePointH(double u, int span) const {

  double N[p+1];
  BasisFuns(span, u, N);

  HPoint2D Cw;
  Cw.wx = 0.; Cw.wy = 0.; Cw.w = 0.;

  for (int i=0; i<=p; i++) {
    Cw.wx += N[i]*P[span-p+i].wx; 
    Cw.wy += N[i]*P[span-p+i].wy;
    Cw.w  += N[i]*P[span-p+i].w; 
  }

  return Cw;

}

/// Calculate first derivatives of the NURBS curve at end points.
/// Eq 4.9 and Eq 4.10 on p126 in The NURBS book.
/// Parameters: 
///   Cd - An array with two elements of type Point2D.
///   It should be initialized in the calling routine with "Point2D Cd[2];",
///   before the function is called.
  void NurbsCurve::FirstD_at_ends(struct Point2D Cd[]) const {
 
 Cd[0].x = double(p)/U[p+1]*(P[1].w/P[0].w)*(P[1].wx/P[1].w-P[0].wx/P[0].w);
 Cd[0].y = double(p)/U[p+1]*(P[1].w/P[0].w)*(P[1].wy/P[1].w-P[0].wy/P[0].w);

 Cd[1].x = double(p)/(1-U[m-p-1])*(P[n-1].w/P[n].w)*(P[n].wx/P[n].w-P[n-1].wx/P[n-1].w);
 Cd[1].y = double(p)/(1-U[m-p-1])*(P[n-1].w/P[n].w)*(P[n].wy/P[n].w-P[n-1].wy/P[n-1].w);
 
 }

/// Computes the derivative of degree d in the homonegeous domain.
/// of the curve at parametric value u.
/// Algorithm A3.2 on p93 in The NURBS book.
/// Parameters:
///   u - Control parameter value.
///   d - Degree of the derivative.
/// Output:
///   ders - An array of size [d+1], containing the derivatives of the curve at u. 
///   This array has to be allocated in the calling function.
  void NurbsCurve::deriveAtH(double u, int d, struct HPoint2D ders[]) const {

  int du = (d < p) ? d : p;

  int span = FindSpan(u);

  Matrix<double> nders(du+1,p+1);
  DersBasisFuns(du, u, span, nders);

  for(int k=0; k<=du; k++) {
    ders[k].wx = 0.0; ders[k].wy = 0.0; ders[k].w = 0.0;
    for(int j=0; j<=p; j++) {
      ders[k].wx += nders[k][j]*P[span-p+j].wx; // Essential eq.
      ders[k].wy += nders[k][j]*P[span-p+j].wy;
      ders[k].w += nders[k][j]*P[span-p+j].w;
    }
  }

}


/// Computes the derivative of degree d in the homonegeous domain. 
/// of the curve at parametric value u.
/// Algorithm A3.2 on p93 in The NURBS book.
/// Parameters:
///   u - Control parameter value.
///   d - Degree of the derivative.
///   span - Specifies the span at which we compute basis functions.
/// Output:
///   ders - An array of size [d+1], with elements of type HPoint2D, containing the derivatives of the curve at u.
///   This array has to be allocated in the calling function.
  void NurbsCurve::deriveAtH(double u, int d, int span, struct HPoint2D ders[]) const {

  int du = (d < p) ? d : p;

  Matrix<double> nders(du+1,p+1);
  DersBasisFuns(du, u, span, nders);

  for(int k=0; k<=du; k++) {
    ders[k].wx = 0.0; ders[k].wy = 0.0; ders[k].w = 0.0;
    for(int j=0; j<=p; j++) {
      ders[k].wx += nders[k][j]*P[span-p+j].wx; // Essential eq.
      ders[k].wy += nders[k][j]*P[span-p+j].wy;
      ders[k].w += nders[k][j]*P[span-p+j].w;
    }
  }

}

/// Computes control points of all derivative curves up to and icluding k-th derivative \f$(d <= p)\f$.
/// Alghorithm A3.3 in The NURBS book.
/// Parameters: 
/// d - Degree of the derivative. 
/// r1, r2 - Auxilliary integers - they take values span-p and span respectively. 
///  With \f$ 0 <= k <= d \f$ and \f$ r1 <= i <= r2-k \f$. If \f$r1 = 0\f$ and \f$r2 = n\f$, all control points are computed.
/// Output: 
///   A matrix PK of size [d+1][r+1] where PK[k][i] is the i-th control point of k-th derivative curve.
/// Mathematical background: 
/// k-th derivative curve's control point $P_i$  is:
/// \f[
///  P_i = \left\{
///  \begin{array}{l l}
///    P_i & \quad \text{if $k = 0$}\\
///    \frac{p-k+1}{u[i+p+1]-u[i+k]}(P_{i+1}^{(k-1)} - P_{i}^{(k-1)}) & \quad \text{if $k > 0$}\\
///  \end{array} \right.
/// \f]
///
  void NurbsCurve::DerivCpts(int d, int r1, int r2, Matrix<HPoint2D> &PK) const {

    int r = r2 - r1;

    for (int i=0; i<=r; i++) { 
       PK[0][i].wx = P[r1 + i].wx; 
       PK[0][i].wy = P[r1 + i].wy; 
       PK[0][i].w  = P[r1 + i].w;
    }   

    for (int k=1; k<=d; k++) {
      int tmp = p-k+1;
      for (int i=0; i<=r-k; i++) {
        PK[k][i].wx = double(tmp)*(PK[k-1][i+1].wx - PK[k-1][i].wx)/(U[r1+i+p+1]-U[r1+i+k] ); 
        PK[k][i].wy = double(tmp)*(PK[k-1][i+1].wy - PK[k-1][i].wy)/(U[r1+i+p+1]-U[r1+i+k] );
        PK[k][i].w = double(tmp)*(PK[k-1][i+1].w - PK[k-1][i].w)/(U[r1+i+p+1]-U[r1+i+k] );
      }
    }

  }

 
/// Computes curve derivatives using Eq. 3.8. p97.
/// Alghorithm A3.4 of the NURBS book.
/// Parameters: 
///   u - Control parameter value where we evaluate derivatives. 
///   d - Degree of the derivative. 
/// Output: 
///   ders - An array of size [d+1], with elements of type HPoint2D, containing the derivatives of the curve at u. 
///   This array has to be allocated in the calling function.
  void NurbsCurve::deriveAtH2(double u, int d, struct HPoint2D ders[]) const {

   int du = ( d < p ) ? d : p;
   for (int k = p+1; k <= d; k++) 
     { ders[k].wx = 0.0; ders[k].wy = 0.0; ders[k].w = 0.0; }
    
   int span = FindSpan(u);

   Matrix<double> N(span+1, p+1);
   AllBasisFuns(span, u, N); 

   Matrix<HPoint2D> PK(d+1,p+1);
 
   DerivCpts(du, span-p, span, PK); 

     for (int k=0; k<=du; k++) {
       ders[k].wx = 0.0; ders[k].wy = 0.0; ders[k].w = 0.0;
       for (int j=0; j<=p-k; j++) {
         ders[k].wx += N[j][p-k]*PK[k][j].wx; // Essential eq.
         ders[k].wy += N[j][p-k]*PK[k][j].wy;
         ders[k].w  += N[j][p-k]*PK[k][j].w;
       }
     }
   
 } 

/// Sets up a matrix containing binomial coefficients Bin.
/// \f[ \tt{Bin}[i][j] = \left( \begin{array}{c}i \\ j\end{array} \right)\f]
/// The binomial coefficients are defined as follows:
/// \f[ \left(\begin{array}{c}   n \\ k \end{array} \right)= \frac{ n!}{k!(n-k)!} \mbox{for $0\leq k \leq n$} \f]
/// and the following relationship applies: 
/// \f[ \left(\begin{array}{c} n+1 \\ k \end{array} \right) = 
/// \left(\begin{array}{c} n \\ k \end{array} \right) +
/// \left(\begin{array}{c} n \\ k-1 \end{array} \right) \f]
/// Parameters:
///   d - Number of rows and cols of "Bin" matrix.
///   Bin - The binomial matrix of size [d][d] to fill.
  void binomialCoef(int d, Matrix<double> &Bin) {
  int n,k;
  // Setup the first line.
  Bin[0][0] = 1.0 ;
  for(k=d-1;k>0;--k)
    Bin[0][k] = 0.0 ;
  // Setup the other lines.
  for(n=0; n<d-1; n++){
    Bin[n+1][0] = 1.0;
    for(k=1;k<d;k++)
      if(n+1<k)
	Bin[n][k] = 0.0;
      else
	Bin[n+1][k] = Bin[n][k] + Bin[n][k-1];
  }
}

/// Computes the derivative at the parameter u. 
/// It uses derivatives from \f$ C^w(u) \f$ to compute derivatives in \f$ C(u) \f$.
/// Functions deriveAtH, and deriveAtH2 work for B-splines. For NURBS we call this function. 
/// Algorithm A4.2 The NURBS book p127.  
/// Parameters:
///   u - Parametric value at which the derivative is computed.
///   d - Degree of derivation.
/// Output:
///   ders - An array of size [d+1] containing the derivatives of the point at u.
///   This array has to be allocated in the calling function.
void NurbsCurve::deriveAt(double u, int d, struct Point2D ders[]) const {

  HPoint2D dersW[d+1];
  deriveAtH(u,d,dersW);
  Point2D v;
  int i,k;


  Matrix<double> Bin(d+1,d+1);  
  binomialCoef(d+1, Bin);


  for(k=0; k<=d; k++){
    v.x = dersW[k].wx;
    v.y = dersW[k].wy;

    for(i=k; i>0;--i) {
      v.x -= (Bin[k][i]*dersW[i].w)*ders[k-i].x;
      v.y -= (Bin[k][i]*dersW[i].w)*ders[k-i].y;
    }
    ders[k].x = v.x; 
    ders[k].y = v.y;
    // Projection.
    ders[k].x /= dersW[0].w; 
    ders[k].y /= dersW[0].w;
  }

}

/// Sometimes we need just first derivative...
/// This function computes the first derivative in the homogenous space.
/// See Eq. 3.4 pp 94 in The NURBS book.
/// Parameter:
///  u - Control parameter value at which we evaluate first derivative.
/// Returns:
/// Cd - The first derivative - a point in homogenous space HPoint2D.
///
/// Note we look for Basis Functions of (p-1)-th degree, where p is the degree of our NURBS curve,
/// and then 
/// There is also an alternative way to solve this problem. We use different knot vector
/// \f[
/// U=\{\underbrace{a,\ldots,a}_{p},u_{p+1},\ldots,u_{m-p-1},\underbrace{b,\ldots,b}_{p} \}, and look for p-th degree basis functions.
/// \f]
/// Check out Eq. 3.6 in The NURBS book for that.
 HPoint2D NurbsCurve::firstD(double u) const {

  int span = FindSpan(u); 

  double N[p];
  nurbsBasisFuns(span, u, p-1,U,N);

  HPoint2D Cd;
     Cd.wx = 0.0; Cd.wx = 0.0; Cd.w = 0.0;
  HPoint2D Qi;

  for(int i=0; i<=p-1; i++) {
    int j = span-p+i;
 
    Qi.wx = (P[j+1].wx-P[j].wx); 
    Qi.wy = (P[j+1].wy-P[j].wy); 
    Qi.w  = (P[j+1].w-P[j].w);

    Qi.wx *= double(p)/(U[j+p+1]-U[j+1]); 
    Qi.wy *= double(p)/(U[j+p+1]-U[j+1]);
    Qi.w  *= double(p)/(U[j+p+1]-U[j+1]);

    Cd.wx += N[i]*Qi.wx;
    Cd.wy += N[i]*Qi.wy;
    Cd.w  += N[i]*Qi.w; 
    
  }
  
  return Cd;
}

/// Sometimes we need just first derivative...
/// This function computes the first derivative in the homogenous space.
/// See Eq. 3.4 pp 94 in The NURBS book.
/// Parameter:
///  u - Control parameter value at which we evaluate first derivative.
///  span - Specifies the span determined by value of u.
/// Returns:
/// Cd - The first derivative - a point in homogenous space HPoint2D.
///
/// Note we look for Basis Functions of (p-1)-th degree, where p is the degree of our NURBS curve,
/// and then 
/// There is also an alternative way to solve this problem. We use different knot vector
/// \f[
/// U=\{\underbrace{a,\ldots,a}_{p},u_{p+1},\ldots,u_{m-p-1},\underbrace{b,\ldots,b}_{p} \}, and look for p-th degree basis functions.
/// \f]
/// Check out Eq. 3.6 in The NURBS book for that.
 HPoint2D NurbsCurve::firstD(double u, int span) const {

  double N[p];
  nurbsBasisFuns(span, u, p-1,U,N);

  HPoint2D Cd;
     Cd.wx = 0.0; Cd.wx = 0.0; Cd.w = 0.0;
  HPoint2D Qi;

  for(int i=0; i<=p-1; i++) {
    int j = span-p+i;
 
    Qi.wx = (P[j+1].wx-P[j].wx); 
    Qi.wy = (P[j+1].wy-P[j].wy); 
    Qi.w  = (P[j+1].w-P[j].w);

    Qi.wx *= double(p)/(U[j+p+1]-U[j+1]); 
    Qi.wy *= double(p)/(U[j+p+1]-U[j+1]);
    Qi.w  *= double(p)/(U[j+p+1]-U[j+1]);

    Cd.wx += N[i]*Qi.wx;
    Cd.wy += N[i]*Qi.wy;
    Cd.w  += N[i]*Qi.w; 
  }
 
  return Cd;
}

/// Computes the first derivative in the normal space.  
/// Parameter: 
///   u - parameter which determines where to compute the derivative
/// Returns: 
///  Cp - The first derivative in normal space.
Point2D NurbsCurve::firstDn(double u) const {

  int span = FindSpan(u); 

  HPoint2D Cd; 
  Cd = firstD(u);

  Point2D Cp;
  Cp.x = Cd.wx; Cp.y = Cd.wy;
  double w = Cd.w;
  Cd = CurvePointH(u,span); 

  Cp.x -= w*project(Cd).x;
  Cp.y -= w*project(Cd).y;

  Cp.x /= Cd.w;
  Cp.y /= Cd.w;

  return Cp;
}

/// Computes normal of the curve at position determined by the value of parameter u.
/// Parameter:
///  u - Control parameter value.
/// Returns:
///  Cn - normal at parametric value u.
Point2D NurbsCurve::normal_at(double u) const {

  Point2D Cn;
  Point2D Cd = firstDn(u);

  Cn.x = -Cd.y;
  Cn.y = Cd.x;

  return Cn;

}

/// Computes normal of the curve at end. 
/// Parameter:
///    which_end: It can be 0 or 1 : parametric beggining or the end of the curve.
Point2D NurbsCurve::normal_at_end(int which_end) const {

  Point2D Cn; 

  Point2D Cd[2];
  FirstD_at_ends(Cd);
  
  // This checks is it at beggining or at the end of the curve.
  if (which_end == 0) {
    Cn.x = -Cd[0].y;
    Cn.y = Cd[0].x;

  } else {
    Cn.x = -Cd[1].y;
    Cn.y = Cd[1].x;   
  }

  return Cn;

}

/// Projects a HPoint2D point, that is a point in homogenous space to normal space.
Point2D project(HPoint2D &Cw) {

  Point2D C;
  
  C.x = Cw.wx/Cw.w;
  C.y = Cw.wy/Cw.w;

  return C;

}

