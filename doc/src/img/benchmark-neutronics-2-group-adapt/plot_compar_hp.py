# import libraries
import numpy, pylab
from pylab import *

# plot DOF convergence graph
axis('equal')
pylab.title("Error convergence")
pylab.xlabel("Degrees of freedom")
pylab.ylabel("Error [%]")
data = numpy.loadtxt("conv_dof_hp_iso_multi.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="hp-FEM (iso)")
data = numpy.loadtxt("conv_dof_hp_anisoh_multi.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="hp-FEM (aniso-h)")
data = numpy.loadtxt("conv_dof_hp_anisop_multi.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="hp-FEM (aniso-p)")
data = numpy.loadtxt("conv_dof_hp_aniso_multi.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="hp-FEM (aniso-hp)")
legend()

# initialize new window
pylab.figure()

# plot CPU convergence graph
axis('equal')
pylab.title("Error convergence")
pylab.xlabel("CPU time [s]")
pylab.ylabel("Error [%]")
data = numpy.loadtxt("conv_cpu_hp_iso_multi.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="hp-FEM (iso)")
data = numpy.loadtxt("conv_cpu_hp_anisoh_multi.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="hp-FEM (aniso-h)")
data = numpy.loadtxt("conv_cpu_hp_anisop_multi.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="hp-FEM (aniso-p)")
data = numpy.loadtxt("conv_cpu_hp_aniso_multi.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="hp-FEM (aniso-hp)")
legend()

# finalize
show()
