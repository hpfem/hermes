# import libraries
import numpy, pylab
from pylab import *

# initialize new window
pylab.figure(figsize=(9, 7), dpi=80)

# plot DOF convergence graph
axis('equal')
pylab.title("Error convergence")
pylab.xlabel("Degrees of freedom")
pylab.ylabel("Relative error")
data = numpy.loadtxt("dgh0/conv_dof_ex_dg-upwind_h0_iso.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="h-DGM (p=0)")
data = numpy.loadtxt("dgh1/conv_dof_ex_dg-upwind_h1_iso.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="h-DGM (p=1)")
data = numpy.loadtxt("dghp/conv_dof_ex_dg-upwind_hp_iso.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="hp-DGM")
data = numpy.loadtxt("supgh1/conv_dof_ex_supg_h1_iso.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="h-SUPGM (p=1)")
data = numpy.loadtxt("supgh2/conv_dof_ex_supg_h2_iso.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="h-SUPGM (p=2)")
data = numpy.loadtxt("supghp/conv_dof_ex_supg_hp_iso.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="hp-SUPGM")
legend()

# initialize new window
pylab.figure(figsize=(9, 7), dpi=80)

# plot CPU convergence graph
axis('equal')
pylab.title("Error convergence")
pylab.xlabel("CPU time [s]")
pylab.ylabel("Relative error")
data = numpy.loadtxt("dgh0/conv_cpu_ex_dg-upwind_h0_iso.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="h-DGM (p=0)")
data = numpy.loadtxt("dgh1/conv_cpu_ex_dg-upwind_h1_iso.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="h-DGM (p=1)")
data = numpy.loadtxt("dghp/conv_cpu_ex_dg-upwind_hp_iso.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="hp-DGM")
data = numpy.loadtxt("supgh1/conv_cpu_ex_supg_h1_iso.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="h-SUPGM (p=1)")
data = numpy.loadtxt("supgh2/conv_cpu_ex_supg_h2_iso.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="h-SUPGM (p=2)")
data = numpy.loadtxt("supghp/conv_cpu_ex_supg_hp_iso.dat")
x = data[:, 0]
y = data[:, 1]
loglog(x, y, "-s", label="hp-SUPGM")
legend()

# finalize
show()
