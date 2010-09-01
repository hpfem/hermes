from pylab import plot, show, legend
import numpy
data = numpy.loadtxt("reach.gp")
x = data[:, 0]
y = data[:, 1]
plot(x, y, ".")#, label="reach")
legend()
show()
