# Compile with:
#
# cython --cplus eigen.pyx

import warnings

def solve_eig_numpy(A, B):
    """
    A, B .... scipy sparse matrices

    Uses numpy to solve the A*x = lambda*B*x eigenproblem.
    """
    from numpy import array, dot
    from numpy.linalg import inv, eig
    A = A.todense()
    B = B.todense()

    print "inverting"
    M = dot(inv(B), A)
    print "solving"
    w, v = eig(M)
    print "sorting the eigenvalues"

    r = []
    for i in range(len(w)):
        vec = v[:, i]
        r.append((w[i], vec))
    r.sort(key=lambda x: x[0])
    print "eigenvalues:"
    eigs = []
    for w, vec in r:
        if w > 0:
            break
        print w
        eigs.append(vec)
    return r

def solve_eig_scipy(A, B):
    """
    A, B .... scipy sparse matrices

    Uses scipy to solve the A*x = lambda*B*x eigenproblem.
    """
    from numpy import array
    from scipy.linalg import eigh
    A = A.todense()
    B = B.todense()

    print "solving"
    w, v = eigh(A, B)
    print "sorting the eigenvalues"

    r = []
    for i in range(len(w)):
        vec = v[:, i]
        E = w[i]
        r.append((E, vec))
    r.sort(key=lambda x: x[0])
    print "eigenvalues:"
    eigs = []
    for w, vec in r:
        #if w > 0:
        #    break
        #print w
        eigs.append(vec)
    return r

def convert_mat(mtx):
    """
    Converts a scipy matrix "mtx" to a pysparse matrix.
    """
    from pysparse import spmatrix
    mtx = mtx.tocsr()
    A = spmatrix.ll_mat(*mtx.shape)
    for i in xrange( mtx.indptr.shape[0] - 1 ):
        ii = slice( mtx.indptr[i], mtx.indptr[i+1] )
        n_in_row = ii.stop - ii.start
        A.update_add_at( mtx.data[ii], [i] * n_in_row, mtx.indices[ii] )
    return A

def solve_eig_pysparse(A, B, n_eigs=4, target_value=-1, verbose=False,
        jdtol=1e-6, max_iter=150):
    """
    Solves the generalized eigenvalue problem.

    A, B ..... scipy matrices
    n_eigs ... number of eigenvalues to solve for

    returns a list of (lmbd, vec), where lmbd is the eigenvalue and vec is the
        eigenvector
    """
    # Suppress the pysparse numpy warnings:
    warnings.simplefilter("ignore", DeprecationWarning)
    from pysparse import jdsym, precon, itsolvers
    if verbose:
        print "converting to pysparse"
    n = A.shape[0]
    A = convert_mat(A)
    B = convert_mat(B)
    if verbose:
        print "solving (%d x %d)" % (n, n)
    Atau = A.copy()
    tau = target_value
    Atau.shift(-tau, B)
    K = precon.jacobi(Atau)
    A = A.to_sss()
    B = B.to_sss()
    kconv, lmbd, Q, it, it_in = jdsym.jdsym(A, B, K, n_eigs, tau, jdtol,
            max_iter, itsolvers.qmrs)
    if verbose:
        print "number of converged eigenvalues:", kconv

    r = []
    for i in range(len(lmbd)):
        vec = Q[:, i]
        r.append((lmbd[i], vec))
    r.sort(key=lambda x: x[0])
    #print "eigenvalues:"
    #eigs = []
    #for w, vec in r:
    #    if w > 0:
    #        break
    #    print w
    #    eigs.append(vec)
    return r
