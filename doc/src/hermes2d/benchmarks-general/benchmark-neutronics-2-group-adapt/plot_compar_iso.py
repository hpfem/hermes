# import libraries
import numpy, pylab
from pylab import *

# plot DOF convergence graph
axis('equal')
pylab.title("Error convergence")
pylab.xlabel("Degrees of freedom")
pylab.ylabel("Error [%]")
data = numpy.loadtxt("conv_dof_h1_1_iso_multi.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="h-FEM (p=1)")
data = numpy.loadtxt("conv_dof_h2_2_iso_multi.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="h-FEM (p=2)")
data = numpy.loadtxt("conv_dof_hp_iso_multi.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="hp-FEM")
legend()

# initialize new window
pylab.figure()

# plot CPU convergence graph
axis('equal')
pylab.title("Error convergence")
pylab.xlabel("CPU time [s]")
pylab.ylabel("Error [%]")
data = numpy.loadtxt("conv_cpu_h1_1_iso_multi.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="h-FEM (p=1)")
data = numpy.loadtxt("conv_cpu_h2_2_iso_multi.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="h-FEM (p=2)")
data = numpy.loadtxt("conv_cpu_hp_iso_multi.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="hp-FEM")
legend()

# finalize
show()
