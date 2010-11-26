# import libraries
import numpy, pylab
from pylab import *

# plot DOF convergence graph
axis('equal')
pylab.title("Error convergence")
pylab.xlabel("CPU time")
pylab.ylabel("Error [%]")
data = numpy.loadtxt("conv_cpu_hp_iso.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="hp-FEM (iso)")
data = numpy.loadtxt("conv_cpu_hp_aniso.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="hp-FEM (aniso)")
legend()

# finalize
show()
