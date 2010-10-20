from numpy import *
import sys
from pysparse import jdsym,spmatrix,itsolvers,precon
matfiles=sys.argv[1:3]
tau=eval(sys.argv[3])
kmax=eval(sys.argv[4])
H=spmatrix.ll_mat_from_mtx(matfiles[0])
U=spmatrix.ll_mat_from_mtx(matfiles[1])
shape=H.shape
T=H.copy()
T.shift(-tau,U)
K=precon.ssor(T.to_sss(),1.0,1)
A=H.to_sss()
M=U.to_sss()
k_conv,lmbd,Q,it,itall=jdsym.jdsym(A,M,K,kmax,tau,1e-12,1000,itsolvers.minres)
NEIG=len(lmbd)
for lam in lmbd:
    print lam,lam-int(lam)
eivecfile=open("eivecs.dat","w")
N=len(Q[:,0])
print >> eivecfile,N
print >> eivecfile,NEIG
for ieig in range(len(lmbd)):
    eivec=Q[:,ieig]
    for val in eivec:
        print >> eivecfile, val
eivecfile.close()
    
