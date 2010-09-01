from pylab import plot, show, legend
import numpy
data = numpy.loadtxt("solution.gp_2")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="velocity")
data = numpy.loadtxt("solution.gp_3")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="steering angle")
legend()
show()
