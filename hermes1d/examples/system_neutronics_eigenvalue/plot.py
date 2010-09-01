import matplotlib.pyplot as plt
import numpy as np
import sys

fig = plt.figure()

# one axes for each group

ax1 = fig.add_subplot(211)
ax1.grid(True)
ax1.axhline(0, color='black', lw=2)
ax2 = fig.add_subplot(212, sharex=ax1)
ax2.grid(True)
ax2.axhline(0, color='black', lw=2)

# computed solution

# plot solutions corresponding to power method iterations specified as arguments
# (enter 'final' as an argument to plot the converged solution) 
for arg in sys.argv[1:]:
	sarg = str(arg)
	
	if sarg == 'final':
		fname = "solution"
	else:
		fname = "solution_"+sarg
	
	# group 1
	data = np.loadtxt(fname+".gp_0")
	x = data[:, 0]
	y = data[:, 1]
	y0 = y[0]	#for normalization of the analytic solution
	ax1.plot(x,y,label='iteration '+sarg)
	
	# group 2
	data = np.loadtxt(fname+".gp_1")
	x = data[:, 0]
	y = data[:, 1]
	ax2.plot(x,y,label='iteration '+sarg)

# analytic solution

# group 1

# region 1	
x1 = np.arange(0, 40, 0.05)
y1 = 0.65259*np.cos(0.023596*x1)-0.0012912*np.cosh(0.11331*x1)
# region 2
x2 = np.arange(40, 70, 0.05)
y2 = 0.12628*np.sinh(0.055596*(70-x2))
# normalization with the same condition as the computed solution
A = y0/y1[0]
# plot normalized solution in both regions
ax1.plot(np.concatenate((x1,x2)), A*np.concatenate((y1,y2)), 'k--', label='reference')

# group 2

# region 1
x1 = np.arange(0, 40, 0.05)
y1 = 0.25577*np.cos(0.023596*x1)+0.0013523*np.cosh(0.11331*x1)
# region 2
x2 = np.arange(40, 70, 0.05)
y2 = np.sinh(0.0212*(70-x2))-0.18263*np.sinh(0.055596*(70-x2))
ax2.plot(np.concatenate((x1,x2)), A*np.concatenate((y1,y2)), 'k--', label='reference')

plt.axes(ax1)
plt.legend()
plt.axes(ax2)
plt.legend()

plt.show()

