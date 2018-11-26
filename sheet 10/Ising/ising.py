import matplotlib as mpl
font = {'family' : 'serif','size' : 12}
mpl.rc('font', **font)

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

magnetization_data = np.loadtxt('magnetization.txt')
betas = magnetization_data[:,0]
'''
for beta in betas:
    lattice = np.loadtxt('lattice%02d.txt' % int(10 * beta))
    plt.matshow(lattice)

    plt.xlabel('$x$',fontsize=20)
    plt.ylabel('$y$',fontsize=20)

    plt.savefig('lattice%02d.pdf' % int(10 * beta), bbox_inches='tight', transparent=True)
'''
temp = 1/betas
magnetization = magnetization_data[:,1]

plt.figure()
plt.plot(temp, magnetization, label=r'$\langle|M(T)|\rangle$', linewidth = 2)
plt.axvline(1.0/np.log(1 + np.sqrt(2)), color='red', label=r'$T = 1/\beta_c$', linewidth = 2)
plt.ylabel(r'$\langle|M(T)|\rangle$', fontsize=20)
plt.xlabel('$T$',fontsize=20)
plt.legend(loc='upper right')
plt.savefig('magnetization.pdf', bbox_inches='tight', transparent=True)
