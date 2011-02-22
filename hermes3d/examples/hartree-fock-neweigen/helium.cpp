const int P=2;
const int Nnuc=1;// number of nuclei
const double E0=  -4.0000000000000000;// number of nuclei
double C[Nnuc]= {   1.0000000000000000}; // coefficients in cusp factor
double Z[Nnuc]= {   2.0000000000000000}; // nuclear charges
double Rnuc[Nnuc][3]={{   0.0000000000000000,   0.0000000000000000,   0.0000000000000000}};// coordinates of the two nuclei

// the following functions are needed for the cusp factor  
double PI=3.1415926535897931;
double const sigma=0.5;
double const sigma2=sigma*sigma;
double gc=1.0/pow(sigma,3)/pow(2*PI,1.5);
double ri(int i,double x,double y,double z){
  double dx,dy,dz;
  dx=x-Rnuc[i][0];
  dy=y-Rnuc[i][1];
  dz=z-Rnuc[i][2];
  return sqrt(dx*dx+dy*dy+dz*dz);
}// this function gives the distances from the point (x,y,z)  to the i-th nucleus
scalar gaussian_pot(double x,double y,double z,scalar &dx, scalar &dy, scalar &dz)
{
   double tmp=0.0;
   double ri1;
   for (int i = 0; i < Nnuc; i++) {
   ri1=ri(i,x,y,z);
   tmp=tmp+2*Z[i]/ri1*erf(ri1/sigma/sqrt(2.0));
}
   return tmp;
}// this function provides the potential for the  gaussian charge distribution around all nuclei compensating the gaussian counter charges used
// for the solution of the poisson equation 

scalar gaussian_cdist(double x,double y,double z,scalar &dx, scalar &dy, scalar &dz)
{
   double tmp=0.0;
   double ri1;
   for (int i = 0; i < Nnuc; i++) {
   ri1=ri(i,x,y,z);
   tmp=tmp-2*Z[i]*gc*exp(-ri1*ri1/2.0/sigma2);
}
   return tmp;
}// this function the provides the gaussian charge distribution around all nuclei used to simplify the solution of the poisson equation

double f(double x,double y,double z){
   double tmp=1.0;
   for (int i = 0; i < Nnuc; i++) tmp=tmp+C[i]*exp(-2.0*Z[i]*ri(i,x,y,z));
   return tmp;
}// this is the factor that satisfies the cusp conditions at the nuclei
double laplacef(double x,double y,double z){
   double tmp=0.0 ;
   for (int i = 0; i < Nnuc; i++)   tmp=tmp+4.0*(C[i]*Z[i]*Z[i]*exp(-2.0*Z[i]*ri(i,x,y,z))-Z[i]/ri(i,x,y,z)*C[i]*exp(-2.0*Z[i]*ri(i,x,y,z)));
   return tmp;
}// this is laplace of f
// The last two functions have dummy argument to make them fit with the ExactSolution Class!
scalar pot(double x, double y, double z,scalar &dx, scalar &dy, scalar &dz){
  scalar  tmp=-laplacef(x,y,z)/f(x,y,z);
  for (int i = 0; i < Nnuc; i++)   tmp=tmp-2.0*Z[i]/ri(i,x,y,z);
  return tmp;
}//this is the potential with the singularities at the nuclei removed and the 1/r replaced by a kink 
scalar  wfun (double x,double y, double z, scalar &dx, scalar &dy, scalar &dz){
  scalar tmp=f(x,y,z)*f(x,y,z);
  return tmp;
}//this is the weight function f**2

