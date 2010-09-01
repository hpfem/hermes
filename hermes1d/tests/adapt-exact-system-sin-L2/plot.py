from pylab import plot, show, legend
import numpy
data = numpy.loadtxt("solution.gp_0")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="0")
data = numpy.loadtxt("solution.gp_1")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="1")
data = numpy.loadtxt("solution_ref.gp_0")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="0 ref")
data = numpy.loadtxt("solution_ref.gp_1")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="1 ref")
legend()
show()
