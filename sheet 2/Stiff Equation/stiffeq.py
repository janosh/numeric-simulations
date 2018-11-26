import matplotlib as mpl
font = {'family' : 'serif',
        'size'   : 16}
mpl.rc('font', **font)

import numpy as np
data_a = np.fromfile('temp.dat', dtype=float, count=-1, sep='')
data_c = np.fromfile('adaptive_temp.dat', dtype=float, count=-1, sep='')

data_a1 = data_a[::2]
data_a2 = data_a[1::2]
data_c1 = data_c[::2]
data_c2 = data_c[1::2]

magn = int(np.log10(max(data_c1)))

import matplotlib.pyplot as plt
plt.figure()

plt.plot(data_a1 / 10**magn, data_a2, label='$T_{(a)}(t)$', linewidth=2)
plt.plot(data_c1 / 10**magn, data_c2, label='$T_{(c)}(t)$', linewidth=2, color='r', linestyle = 'dashed')

plt.yscale('log')
plt.xlabel("$t[10^{%d}s]$" % magn,fontsize=20)
plt.ylabel("$T\,[K]$",fontsize=20)
plt.xlim(right=max(data_c1) / 10**magn)
plt.legend(loc='upper right')

plt.savefig('stiffeq.pdf', bbox_inches='tight', transparent=True)
