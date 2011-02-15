# import libraries
import numpy, pylab
from pylab import *

# plot DOF convergence graph
pylab.title("Error convergence")
pylab.xlabel("Degrees of freedom")
pylab.ylabel("Error [%]")
axis('equal')
data = numpy.loadtxt("conv_dof.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, '-s', label="error (est)")
legend()

# initialize new window
pylab.figure()

# plot CPU convergence graph
pylab.title("Error convergence")
pylab.xlabel("CPU time (s)")
pylab.ylabel("Error [%]")
axis('equal')
data = numpy.loadtxt("conv_cpu.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, '-s', label="error (est)")
legend()

# initialize new window
pylab.figure()

# plot DOF convergence graph w.r.t. keff
pylab.title("K_eff error convergence")
pylab.xlabel("Degrees of freedom")
pylab.ylabel("Error [pcm]")
axis('equal')
data = numpy.loadtxt("conv_dof_keff.dat")
x = data[:, 0]
y = data[:, 1]
semilogx(x, y, '-s', label="k_eff error")
legend()

# initialize new window
pylab.figure()

# plot CPU convergence graph w.r.t. keff 
pylab.title("K_eff error convergence")
pylab.xlabel("CPU time (s)")
pylab.ylabel("Error [pcm]")
axis('equal')
data = numpy.loadtxt("conv_cpu_keff.dat")
x = data[:, 0]
y = data[:, 1]
semilogx(x, y, '-s', label="k_eff error")
legend()

# finalize
show()
