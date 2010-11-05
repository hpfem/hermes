from pylab import plot, show, legend
import numpy
data = numpy.loadtxt("space_ref.dat")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="space_ref")
data = numpy.loadtxt("space.dat")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="space")
legend()
show()
