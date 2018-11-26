import matplotlib as mpl
font = {'family' : 'serif',
        'size'   : 16}
mpl.rc('font', **font)

import numpy as np
jacobian = np.loadtxt('Jacobi/res_jacobi.txt')
gaussseidel = np.loadtxt('Gauss-Seidel/res_gausseidel.txt')
#cb_gaussseidel = np.loadtxt('Chessboard Gauss-Seidel/res_gausseidel_cb.txt')
multigrid = np.loadtxt('Multigrid/res_multigrid.txt')
steps = np.arange(0, len(jacobian), 1)

import matplotlib.pyplot as plt

plt.plot(steps, jacobian, label='$S_\mathrm{Jac}$')
plt.plot(steps, gaussseidel, label='$S_\mathrm{GS}$')
#plt.plot(steps, cb_gaussseidel, label='$S_\mathrm{cbGS}$')
plt.plot(steps, multigrid, label='$S_\mathrm{mult}$')

plt.yscale('log')
plt.xlabel("Step $N$",fontsize=20)
plt.ylabel("Residual $S$",fontsize=20)
plt.ylim(0.95 * min(min(jacobian),min(gaussseidel),min(multigrid)), 1.05 * max(max(jacobian),max(gaussseidel),max(multigrid)))
plt.legend(loc = 'upper right')

from matplotlib.ticker import FormatStrFormatter
plt.gca().yaxis.set_minor_formatter(FormatStrFormatter("%.2f"))
plt.gca().yaxis.set_major_formatter(FormatStrFormatter("%.2f"))

plt.savefig('res-decay.pdf', bbox_inches='tight', transparent=True)
