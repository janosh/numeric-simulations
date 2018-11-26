import matplotlib as mpl
font = {'family' : 'serif',
'size'   : 12}
mpl.rc('font', **font)

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

rho = np.loadtxt("rho.txt")
rho_xconst = rho[len(rho)/2]
rho_yconst = rho[:,len(rho[0])/2]

Lx = 3.0
Ly = 1.5
x = np.linspace(0, Lx, len(rho))
y = np.linspace(0, Ly, len(rho[0]))

plt.figure()
plt.plot(x, rho_yconst, '-^', label=r'$\rho(x)$', linewidth=2)
plt.ylabel(r"$\rho$", fontsize=20)
plt.xlabel("$x$",fontsize=20)
plt.xlim(min(x), max(x))
plt.ylim(0, 1.1 * max(rho_yconst))
plt.legend(loc='upper left')
plt.savefig('rho-yconst.pdf', bbox_inches='tight', transparent=True)

plt.figure()
plt.plot(y, rho_xconst, '-^', label=r'$\rho(y)$', linewidth=2)
plt.ylabel(r"$\rho$", fontsize=20)
plt.xlabel("$y$",fontsize=20)
plt.xlim(min(y), max(y))
plt.ylim(0, 1.1 * max(rho_xconst))
plt.legend(loc='upper left')
plt.savefig('rho-xconst.pdf', bbox_inches='tight', transparent=True)
