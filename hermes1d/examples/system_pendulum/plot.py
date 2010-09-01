try:
    from jsplot import plot, show, grid, legend
except ImportError:
    from pylab import plot, show, grid, legend
import numpy
data = numpy.loadtxt("solution.gp_0")
x = data[:, 0]
y = data[:, 1]
label = "Pendulum angle"
plot(x, y, label=label)
legend()
grid(True)
data = numpy.loadtxt("solution.gp_1")
x = data[:, 0]
y = data[:, 1]
label = "Angular velocity"
plot(x, y, label=label)
legend()
grid(True)
show()
