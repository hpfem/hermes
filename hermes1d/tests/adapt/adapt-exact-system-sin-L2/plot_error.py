from pylab import plot, show, legend
import numpy
data = numpy.loadtxt("error_est.gp")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="error (est)")
data = numpy.loadtxt("error_exact.gp")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="error (exact)")
legend()
show()
