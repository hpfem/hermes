import matplotlib.pyplot as plt
import numpy as np

# material data
Q  = [0.0, 1.5, 1.8, 1.5, 1.8, 1.8, 1.5]
D1 = 7*[1.2]
D2 = 7*[0.4]
S1 = 7*[0.03]
S2 = [0.1, 0.2, 0.25, 0.2, 0.25, 0.25, 0.2]
S12= [0.02] + 6*[0.015]
nSf1 = [0.005] + 6*[0.0075]
nSf2 = 7*[0.1]

fig = plt.figure()

# one axes for each group

ax1 = fig.add_subplot(211)
ax1.grid(True)
ax1.axhline(0, color='black', lw=2)
ax2 = fig.add_subplot(212, sharex=ax1)
ax2.grid(True)
ax2.axhline(0, color='black', lw=2)

# computed solution

# group 1
data = np.loadtxt("solution.gp_0")
x = data[:, 0]
y = data[:, 1]
ax1.plot(x,y,label='approximate')
# group 2
data = np.loadtxt("solution.gp_1")
x = data[:, 0]
y = data[:, 1]
ax2.plot(x,y,label='approximate')

# analytic solution (valid only with certain distance from interfaces)

phi1 = np.array([]) # group 1
phi2 = np.array([]) # group 2
for i in range(0,7,1):
	phi1loc = np.array(100*[ Q[i]/(S1[i] - nSf1[i] - nSf2[i]*S12[i]/S2[i]) ])
	phi1 = np.append(phi1, phi1loc)
	phi2 = np.append(phi2, S12[i]/S2[i]*phi1loc)

x = np.arange(0, 700, 1)	
ax1.plot(x, phi1, 'k--', label='reference')
ax2.plot(x, phi2, 'k--', label='reference')


plt.axes(ax1)
plt.legend()
plt.axes(ax2)
plt.legend()

plt.show()

