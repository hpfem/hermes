# import libraries
import numpy, pylab
from pylab import *

data = numpy.loadtxt("spline.dat")
x = data[:, 0]
y = data[:, 1]
plot(x, y, '-o', label="cubic spline")
data = numpy.loadtxt("spline_der.dat")
x = data[:, 0]
y = data[:, 1]
plot(x, y, '-*', label="derivative")

legend()

# finalize
show()
