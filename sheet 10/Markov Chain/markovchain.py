import matplotlib as mpl
font = {'family' : 'serif',
        'size'   : 16}
mpl.rc('font', **font)

import numpy as np
histogram = np.loadtxt('markovchain.txt')

def p(x):
	return np.exp(-(x + 2 * np.cos(x)**2)**2)/1.104
x = np.linspace(min(histogram), max(histogram), 1000)

import matplotlib.pyplot as plt
plt.figure()
plt.hist(histogram, bins=np.arange(min(histogram), max(histogram), 0.02), label='$p_\mathrm{hist}(x)$', normed = 1, linewidth = 0.1)
plt.plot(x, p(x), label='$p(x)$', color = 'red', linewidth = 2)
plt.ylabel('$p(x)$', fontsize=20)
plt.xlabel('$x$',fontsize=20)
plt.xlim(min(histogram), max(histogram))
plt.legend(loc='upper right')
plt.savefig('markovchain_histogram.pdf', bbox_inches='tight', transparent=True)
