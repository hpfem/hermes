# import libraries
import numpy, pylab
from pylab import *

# plot DOF convergence graph
axis('equal')
pylab.title("Error convergence")
pylab.xlabel("Degrees of freedom")
pylab.ylabel("Error [%]")
data = numpy.loadtxt("conv_dof_exact_hp_iso.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="hp-FEM (iso)")
data = numpy.loadtxt("conv_dof_exact_hp_anisoh.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="hp-FEM (aniso h)")
data = numpy.loadtxt("conv_dof_exact_hp_aniso.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="hp-FEM (aniso hp)")

legend()

# initialize new window
pylab.figure()

axis('equal')
pylab.title("Error convergence")
pylab.xlabel("CPU time (s)")
pylab.ylabel("Error [%]")
data = numpy.loadtxt("conv_cpu_exact_hp_iso.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="hp-FEM (iso)")
data = numpy.loadtxt("conv_cpu_exact_hp_anisoh.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="hp-FEM (aniso h)")
data = numpy.loadtxt("conv_cpu_exact_hp_aniso.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="hp-FEM (aniso hp)")

legend()


# finalize
show()
