from pylab import plot, show, legend
import numpy
data = numpy.loadtxt("solution_ref.gp")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="reference")
data = numpy.loadtxt("solution.gp")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="coarse")
legend()
show()
