# import libraries
import numpy, pylab
from pylab import *

# plot DOF convergence graph
axis('equal')
pylab.title("Error convergence")
pylab.xlabel("Degrees of freedom")
pylab.ylabel("Error [%]")
data = numpy.loadtxt("conv_dof_h1_m.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="h-FEM (p=1, multi)")
data = numpy.loadtxt("conv_dof_h1_s.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="h-FEM (p=1, single)")
legend()

# finalize
show()
