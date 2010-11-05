from pylab import plot, show, legend
import numpy
data = numpy.loadtxt("solution.dat")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="solution")
legend()
show()
