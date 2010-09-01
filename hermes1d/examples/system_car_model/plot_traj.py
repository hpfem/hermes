from pylab import plot, show, legend
import numpy
data = numpy.loadtxt("trajectory.gp")
x = data[:, 0]
y = data[:, 1]
plot(x, y, ".", label="trajectory")
legend()
show()
