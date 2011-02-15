# import libraries
import numpy, pylab
from pylab import *

# plot DOF convergence graph
pylab.title("Time step history")
pylab.xlabel("Physical time (days)")
pylab.ylabel("Time step size")
data = numpy.loadtxt("time_step_history.dat")
x = data[:, 0]
y = data[:, 1]
plot(x, y, "-s", label="dof")

legend()

# initialize new window
#pylab.figure()



# finalize
show()
