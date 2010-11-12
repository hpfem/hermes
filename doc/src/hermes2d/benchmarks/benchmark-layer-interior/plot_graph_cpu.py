# import libraries
import numpy, pylab
from pylab import *

# plot DOF convergence graph
axis('equal')
pylab.title("Error convergence")
pylab.xlabel("CPU time")
pylab.ylabel("Error [%]")
data = numpy.loadtxt("conv_cpu_exact_h1.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="h-FEM (p=1)")
data = numpy.loadtxt("conv_cpu_exact_h2.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="h-FEM (p=2)")
data = numpy.loadtxt("conv_cpu_exact_hp.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="hp-FEM")
legend()

# finalize
show()
