from numpy import *
import sys
from pysparse import jdsym, spmatrix, itsolvers, precon
matfiles = sys.argv[1:5]
target_value = eval(sys.argv[3])
eigenval_num = eval(sys.argv[4])
jdtol = eval(sys.argv[5])
max_iter = eval(sys.argv[6])
mat_left = spmatrix.ll_mat_from_mtx(matfiles[0])
mat_right = spmatrix.ll_mat_from_mtx(matfiles[1])
shape = mat_left.shape
T = mat_left.copy()
T.shift(-target_value, mat_right)
K = precon.ssor(T.to_sss(), 1.0, 1) # K is preconditioner.
A = mat_left.to_sss()
M = mat_right.to_sss()
k_conv, lmbd, Q, it, itall = jdsym.jdsym(A, M, K, eigenval_num, target_value, jdtol, max_iter, itsolvers.minres)
NEIG = len(lmbd)
for lam in lmbd:
    print lam, lam - int(lam)
eivecfile = open("eivecs.dat", "w")
N = len(Q[:,0])
print >> eivecfile, N
print >> eivecfile, NEIG
for ieig in range(len(lmbd)):
    eivec = Q[:,ieig]
    for val in eivec:
        print >> eivecfile, val
eivecfile.close()
    
