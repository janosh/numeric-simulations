import matplotlib as mpl
font = {'family' : 'serif',
        'size'   : 16}
mpl.rc('font', **font)

import numpy as np
radii = np.fromfile('radii.dat', dtype=float, count=-1, sep='')
forces = np.fromfile('forces.dat', dtype=float, count=-1, sep='')

import matplotlib.pyplot as plt
plt.figure()

def f(r):
    return 2/r

r = np.arange(min(radii), 10, 1)

plt.scatter(radii, forces, label='$a_\mathrm{data}(r)$')
plt.plot(r, f(r), label='$a(r) = 2/r$')
plt.axvline(1.0/256, label='$L/N_\mathrm{grid}$', color = 'r')

plt.xlabel("$r$",fontsize=20)
plt.ylabel("$a(r)$",fontsize=20)
plt.yscale('log')
plt.xscale('log')
plt.xlim(min(radii), 1.1 * max(radii))
plt.legend(loc='upper right')

plt.savefig('force-law.pdf', bbox_inches='tight', transparent=True)
