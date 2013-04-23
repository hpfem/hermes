# import libraries
import numpy, pylab
from pylab import *

# plot DOF convergence graph
pylab.title("Error convergence")
pylab.xlabel("Degrees of freedom")
pylab.ylabel("Error [%]")
axis('equal')
pylab.grid(True)
ax = pylab.gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(20)
data = numpy.loadtxt("conv_dof_est.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, '-s', label="error (est)")
legend()

# finalize
show()
