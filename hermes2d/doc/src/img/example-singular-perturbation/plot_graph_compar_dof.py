# import libraries
import numpy, pylab
from pylab import *

# plot DOF convergence graph
axis('equal')
pylab.title("Error convergence")
pylab.xlabel("Degrees of freedom")
pylab.ylabel("Error [%]")
data = numpy.loadtxt("conv_dof_h1_aniso.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="h1-FEM (aniso)")
data = numpy.loadtxt("conv_dof_h2_aniso.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="h2-FEM (aniso)")
data = numpy.loadtxt("conv_dof_hp_aniso.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="hp-FEM (aniso)")
legend()

# finalize
show()
