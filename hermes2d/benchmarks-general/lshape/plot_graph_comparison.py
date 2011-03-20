# This is a utility that produces comparison plots 
# of DOF and CPU convergence. To use it, your directory 
# must contain Hermes convergence files 
# conv_dof_exact_h1.dat, conv_dof_exact_h2.dat
# and conv_dof_exact_hp.dat

# import libraries
import numpy, pylab
from pylab import *

# plot DOF convergence graph
pylab.yscale("log")
pylab.title("Error convergence")
pylab.xlabel("Degrees of freedom")
pylab.ylabel("Error [%]")
axis('equal')
data = numpy.loadtxt("conv_dof_exact_h1.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, '-s', label="h-FEM (p=1)")
data = numpy.loadtxt("conv_dof_exact_h2.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, '-s', label="h-FEM (p=2)")
data = numpy.loadtxt("conv_dof_exact_hp.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, '-s', label="hp-FEM")
legend()

# initialize new window
pylab.figure()

# plot CPU convergence graph
pylab.yscale("log")
pylab.title("Error convergence")
pylab.xlabel("CPU time (s)")
pylab.ylabel("Error [%]")
axis('equal')
data = numpy.loadtxt("conv_cpu_exact_h1.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, '-s', label="h-FEM (p=1)")
data = numpy.loadtxt("conv_cpu_exact_h2.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, '-s', label="h-FEM (p=2)")
data = numpy.loadtxt("conv_cpu_exact_hp.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, '-s', label="hp-FEM")
legend()

# finalize
show()
