from pylab import plot, show
import numpy
data = numpy.loadtxt("solution_320W.gp")
x = data[:, 0]
y = data[:, 1]
plot(x, y)
show()
