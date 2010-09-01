from pylab import plot, show, legend
import numpy
data = numpy.loadtxt("solution.gp_0")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="bas-0")
data = numpy.loadtxt("solution.gp_1")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="bas-1")
data = numpy.loadtxt("solution.gp_2")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="bas-2")
data = numpy.loadtxt("solution.gp_3")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="bas-3")
"""
data = numpy.loadtxt("solution_ref.gp_0")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="ref-0")
data = numpy.loadtxt("solution_ref.gp_1")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="ref-1")
data = numpy.loadtxt("solution_ref.gp_2")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="ref-2")
data = numpy.loadtxt("solution_ref.gp_3")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="ref-3")
"""
legend()
show()
