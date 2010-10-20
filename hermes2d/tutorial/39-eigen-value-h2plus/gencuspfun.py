from numpy import *
def cuspfun(Z,R):
    """
    factory for cusp function satisfying cusp condition at N Nuclei
    arguments:
    Z[i],i=0,N-1: charges of N nuclei
    R[i],i=0,N=1: coordinates of N nuclei
    """
    def norm(v):
        return sum(v*v)**0.5
    N=len(Z)
    b=zeros((N,N),"d")
    for i in range(N):
        for j in range(N):
            if i != j:
                rij=norm(R[i]-R[j])
                b[i,j]=-exp(-2*Z[j]*rij)
            else:
                b[i,j]=1.0
    #print b
    d=ones(N,"d")
    ci=linalg.solve(b,d)
    print "ci=%24.15f %24.15f" % tuple(ci)
    def ri(i,x,y,z):
        return sqrt( (x-R[i][0])**2 + (y-R[i][1])**2 + (z-R[i][2])**2 )
    def f(x,y,z):
        tmp=1.0
        for i in range(N):
           tmp=tmp+ci[i]*exp(-2*Z[i]*ri(i,x,y,z))
        return tmp
    def laplacef(x,y,z):
        tmp=0*x
        for i in range(N):
           tmp=tmp+4*Z[i]**2*ci[i]*exp(-2*Z[i]*ri(i,x,y,z))
        for i in range(N):
           tmp=tmp-4*Z[i]/ri(i,x,y,z)*ci[i]*exp(-2*Z[i]*ri(i,x,y,z))           
        return tmp

    def pot(x,y,z):
        tmp=0*x
        for i in range(N):
            tmp=tmp-2*Z[i]/ri(i,x,y,z)
        return tmp-laplacef(x,y,z)/f(x,y,z)
    return f,pot


