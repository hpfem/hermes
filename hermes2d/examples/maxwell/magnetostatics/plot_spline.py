# import libraries
import numpy, pylab
from pylab import *

# plot DOF convergence graph
#pylab.yscale("log")
#pylab.title("Error convergence")
#pylab.xlabel("Degrees of freedom")
#pylab.ylabel("Error [%]")
#axis('equal')
data = numpy.loadtxt("spline.dat")
x = data[:, 0]
y = data[:, 1]
plot(x, y, '-o', label="cubic spline")
data = numpy.loadtxt("spline_der.dat")
x = data[:, 0]
y = data[:, 1]
plot(x, y, '-*', label="derivative")

legend()

# finalize
show()
