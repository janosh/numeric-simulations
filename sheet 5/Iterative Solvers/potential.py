import matplotlib as mpl
font = {'family' : 'serif',
        'size'   : 12}
mpl.rc('font', **font)

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

solvers = ["Jacobi/phi_jac", "Gauss-Seidel/phi_gs", "Multigrid/phi_mg"]

for i in solvers:
    i = np.fromfile(i + '.dat', dtype=float, count=-1, sep='')
    i = i.reshape((np.sqrt(len(i)), np.sqrt(len(i))))
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = '3d')
    x = y = np.arange(0, 1, 1.0/np.sqrt(len(i)))
    x, y = np.meshgrid(x, y)
    ax.plot_surface(x, y, i)
    ax.set_xlabel("$x$",fontsize=20)
    ax.set_ylabel("$y$",fontsize=20)
    ax.set_zlabel(r"$\phi(\mathbf{r})$",fontsize=20)
    plt.savefig(i + '.pdf', bbox_inches='tight', transparent=True)
