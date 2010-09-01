# import libraries
import numpy, pylab
from pylab import *

# plot DOF convergence graph
axis('equal')
pylab.title("Error convergence")
pylab.xlabel("CPU time")
pylab.ylabel("Error [%]")
data = numpy.loadtxt("conv_cpu_hp.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="hp-FEM (multi)")
data = numpy.loadtxt("conv_cpu_hp_single.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="hp-FEM (single)")
legend()

# finalize
show()
