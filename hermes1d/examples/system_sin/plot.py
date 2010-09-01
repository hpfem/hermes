try:
    from jsplot import plot, show
except ImportError:
    from pylab import plot, show
import numpy
data = numpy.loadtxt("solution.gp_0")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="0")
data = numpy.loadtxt("solution.gp_1")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="1")
show()
