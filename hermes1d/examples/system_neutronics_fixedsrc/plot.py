import matplotlib.pyplot as plt
import numpy as np
from numpy import cosh 

# material data
a = 80.0
Q = 1.5
D1 = 1.2
D2 = 0.4
S1 = 0.03
S2 = 0.1
S12 = 0.02
L1 = np.sqrt(D1/S1)
L2 = np.sqrt(D2/S2)

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

# analytic solution

x = np.arange(0, a, 0.05)

# group 1
phi1 = Q/S1 * (1 - cosh(x/L1)/cosh(a/L1))
ax1.plot(x, phi1, 'k--', label='reference')

# group 2
phi2 = Q*S12/(S1*S2) * (1 - L1**2/(L1**2 - L2**2)*cosh(x/L1)/cosh(a/L1) + L2**2/(L1**2-L2**2)*cosh(x/L2)/cosh(a/L2))
ax2.plot(x, phi2, 'k--', label='reference')

plt.axes(ax1)
plt.legend()
plt.axes(ax2)
plt.legend()

plt.show()

