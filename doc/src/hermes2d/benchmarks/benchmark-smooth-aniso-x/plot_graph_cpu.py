# import libraries
import numpy, pylab
from pylab import *

# plot DOF convergence graph
axis('equal')
pylab.title("Error convergence")
pylab.xlabel("CPU time")
pylab.ylabel("Error [%]")
data = numpy.loadtxt("conv_cpu_h1_aniso.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="h1-FEM (aniso)")
data = numpy.loadtxt("conv_cpu_h1_iso.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="h1-FEM (iso)")
data = numpy.loadtxt("conv_cpu_h2_aniso.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="h2-FEM (aniso)")
data = numpy.loadtxt("conv_cpu_h2_iso.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="h2-FEM (iso)")
data = numpy.loadtxt("conv_cpu_hp_aniso.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="hp-FEM (aniso)")
data = numpy.loadtxt("conv_cpu_hp_iso.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="hp-FEM (iso)")
legend()

# finalize
show()
