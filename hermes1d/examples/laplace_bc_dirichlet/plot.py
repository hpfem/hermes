try:
    from jsplot import plot, show
except ImportError:
    from pylab import plot, show
import numpy
data = numpy.loadtxt("solution.gp")
x = data[:, 0]
y = data[:, 1]
plot(x, y)
show()
