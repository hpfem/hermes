try:
    from jsplot import plot, show, legend
except ImportError:
    from pylab import plot, show, legend
import numpy
data = numpy.loadtxt("solution_ref.gp")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="FTR array")
data = numpy.loadtxt("solution.gp")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="solution")
legend()
show()
