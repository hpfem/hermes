#!/usr/bin/env python
from numpy import *
import os,sys
def remove_duplicate(list):
    tmp={}
    for i in range(len(list)):
        tmp[list[i]]=i
    outlist=[]
    for x in tmp.keys():
        outlist.append(x)
    outlist.sort()
    return outlist

xmax=input()
Nnuc=input()
x,y,z,Z=[],[],[],[]
for i in range(Nnuc):
    ZZ,xx,yy,zz=input()
    Z.append(ZZ)
    x.append(xx)
    y.append(yy)
    z.append(zz)
meshname=input()
r=input()
P=input()
E0=input()
def refine(g,n):
    gcopy=g[:]
    s=arange(1,n)/float(n)
    for i in range(1,len(gcopy)):
        i1=g.index(gcopy[i])
        gleft=g[i1-1]
        gright=g[i1]
        if gleft == -xmax:
            t=1-(1-s)**2
        elif gright == xmax:
            t=s**2
        else:
            t=0.5-0.5*cos(s*pi)
        g=g[:i1]+(g[i1-1]+t*(g[i1]-g[i1-1])).tolist()+g[i1:]
    return g    
x.sort()
y.sort()
z.sort()
x=array(x)
y=array(y)
z=array(z)
# create initial hex grid with nuclei as edge nodes 
xi=(-xmax,)+tuple(x)+(xmax,)
yi=(-xmax,)+tuple(y)+(xmax,)
zi=(-xmax,)+tuple(z)+(xmax,)
xi=list(xi)
yi=list(yi)
zi=list(zi)
xi=remove_duplicate(xi)
yi=remove_duplicate(yi)
zi=remove_duplicate(zi)
xi=refine(xi,r)
yi=refine(yi,r)
zi=refine(zi,r)
xs=("["+(len(xi)-1)*"%21.16f,"+"%21.16f ]") % tuple(xi)
ys=("["+(len(yi)-1)*"%21.16f,"+"%21.16f ]") % tuple(yi)
zs=("["+(len(zi)-1)*"%21.16f,"+"%21.16f ]") % tuple(zi)
command='./mkmesh3d.py ' +'"'+xs+'"    "'+ys+'"   "' +zs+'" '+'> %s' % meshname
os.system(command)
R=zeros((Nnuc,3),"d")
R[:,0]=x
R[:,1]=y
R[:,2]=z

def norm(v):
    return sum(v*v)**0.5

b=zeros((Nnuc,Nnuc),"d")
for i in range(Nnuc):
    for j in range(Nnuc):
        if i != j:
            rij=norm(R[i]-R[j])
            b[i,j]=-exp(-2*Z[j]*rij)
        else:
            b[i,j]=1.0
d=ones(Nnuc,"d")
C=linalg.solve(b,d)
R.shape=3*Nnuc
print "const int P=%d;" % P
print "const int Nnuc=%d;// number of nuclei" % Nnuc 
print "const double E0=%21.16f;// number of nuclei" % E0
S=((Nnuc-1)*"%21.16f, "+"%21.16f") % tuple(C)
print "double C[Nnuc]= {"+S+"}; // coefficients in cusp factor" 
S=((Nnuc-1)*"%21.16f, "+"%21.16f") % tuple(Z)
print "double Z[Nnuc]= {"+S+"}; // nuclear charges"
S=("{"+(Nnuc-1)*"{%21.16f,%21.16f,%21.16f},"+"{%21.16f,%21.16f,%21.16f}};" ) % tuple(R)
print "double Rnuc[Nnuc][3]=" +S+"// coordinates of the two nuclei"
print """
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
"""
