from pylab import plot, show, legend
import numpy
data = numpy.loadtxt("solution.gp_0")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="0")
data = numpy.loadtxt("solution.gp_1")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="1")
data = numpy.loadtxt("solution.gp_2")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="2")
data = numpy.loadtxt("solution.gp_3")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="3")
"""
data = numpy.loadtxt("solution_ref.gp_0")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="0 ref")
data = numpy.loadtxt("solution_ref.gp_1")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="1 ref")
data = numpy.loadtxt("solution_ref.gp_2")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="2 ref")
data = numpy.loadtxt("solution_ref.gp_3")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="3 ref")
"""
legend()
show()
