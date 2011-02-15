# import libraries
import numpy, pylab
from pylab import *

# plot DOF convergence graph
pylab.title("Error convergence")
pylab.xlabel("Degrees of freedom")
pylab.ylabel("Error [%]")
axis('equal')
data = numpy.loadtxt("conv_dof_est.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, '-s', label="error (est)")
legend()

# initialize new window
pylab.figure()

# plot CPU convergence graph
pylab.title("Error convergence")
pylab.xlabel("CPU time (s)")
pylab.ylabel("Error [%]")
axis('equal')
data = numpy.loadtxt("conv_cpu_est.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, '-s', label="error (est)")
legend()

# finalize
show()
