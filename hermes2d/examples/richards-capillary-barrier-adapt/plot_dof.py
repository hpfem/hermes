# import libraries
import numpy, pylab
from pylab import *

# plot DOF convergence graph
pylab.title("Problem size history")
pylab.xlabel("Physical time (days)")
pylab.ylabel("Number of degrees of freedom")
data = numpy.loadtxt("time_dof.dat")
x = data[:, 0]
y = data[:, 1]
plot(x, y, "-s", label="problem size")

legend()

# initialize new window
#pylab.figure()



# finalize
show()
