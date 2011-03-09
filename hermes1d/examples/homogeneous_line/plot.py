from pylab import plot, show, legend, subplot, ylim
import numpy as np
import scipy as sc
data = np.loadtxt("solution.gp_0")
x1 = data[:, 0]
y1 = data[:, 1]
subplot(4,1,1)
plot(x1, y1, label="0")

data = np.loadtxt("solution.gp_1")
x2 = data[:, 0]
y2 = data[:, 1]
subplot(4,1,2)
plot(x1, y2, label="1")

data = np.loadtxt("solution.gp_2")
x3 = data[:, 0]
y3 = data[:, 1]
subplot(4,1,3)
plot(x3, y3, label="2")

data = np.loadtxt("solution.gp_3")
x4 = data[:, 0]
y4 = data[:, 1]
subplot(4,1,4)
plot(x4, y4, label="3")
show()
