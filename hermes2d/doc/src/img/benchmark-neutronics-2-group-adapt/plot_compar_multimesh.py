# import libraries
import numpy, pylab
from pylab import *

# name of the configuration
name = "hp_iso"

# plot DOF convergence graph
axis('equal')
pylab.title("Error convergence")
pylab.xlabel("Degrees of freedom")
pylab.ylabel("Error [%]")
data = numpy.loadtxt("conv_dof_"+name+"_multi.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="error (multi)")
data = numpy.loadtxt("conv_dof_"+name+"_single.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="error (single)")

legend()

# initialize new window
pylab.figure()

axis('equal')
pylab.title("Error convergence")
pylab.xlabel("CPU time [s]")
pylab.ylabel("Error [%]")
data = numpy.loadtxt("conv_cpu_"+name+"_multi.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="error (multi)")
data = numpy.loadtxt("conv_cpu_"+name+"_single.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="error (single)")
legend()


# finalize
show()
