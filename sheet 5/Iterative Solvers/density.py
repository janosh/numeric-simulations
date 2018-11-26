import matplotlib as mpl
font = {'family' : 'serif',
        'size'   : 12}
mpl.rc('font', **font)

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def fct(x, y):
  return 10 * np.exp(-50 * ((x - 1.0/2)**2 + (y - 1.0/2)**2))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
x = y = np.arange(0, 1, 0.005)
X, Y = np.meshgrid(x, y)
zs = np.array([fct(x,y) for x,y in zip(np.ravel(X), np.ravel(Y))])
Z = zs.reshape(X.shape)

ax.plot_surface(X, Y, Z)

ax.set_xlabel("$x$",fontsize=20)
ax.set_ylabel("$y$",fontsize=20)
ax.set_zlabel(r"$\rho(\mathbf{r})$",fontsize=20)

plt.savefig('density.pdf', bbox_inches='tight', transparent=True)
