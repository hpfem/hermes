from pylab import plot, show, legend
import numpy
data = numpy.loadtxt("error_est.gp")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="max FTR error")
data = numpy.loadtxt("error_exact.gp")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="exact error")
legend()
show()
