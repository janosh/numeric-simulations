import numpy as np
data_a = np.fromfile('error.dat', dtype=float, count=-1, sep='')

data_a1 = data_a[::2]
data_a2 = data_a[1::2]

import matplotlib as mpl
font = {'family' : 'serif',
        'size'   : 16}
mpl.rc('font', **font)

import matplotlib.pyplot as plt
fig = plt.figure()

plt.plot(data_a1, data_a2, label='$\eta_E(t)$', linewidth=2)

plt.xlabel("$t \, [a.u.]$",fontsize=20)
plt.ylabel("$\eta_E(t)$",fontsize=20)
plt.legend(loc='upper right')

plt.savefig('error.pdf', bbox_inches='tight', transparent=True)
