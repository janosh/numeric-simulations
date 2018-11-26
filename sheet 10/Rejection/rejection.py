import matplotlib as mpl
font = {'family' : 'serif',
        'size'   : 16}
mpl.rc('font', **font)

import numpy as np
histogram = np.loadtxt("rejection_sample.txt")

p_0 = 1/12.8136
def p(x):
	return p_0/((x - 2)**4 + np.sin(x - 3)**8)
x = np.linspace(0.0, 5.0, 500)

import matplotlib.pyplot as plt
plt.figure()
plt.hist(histogram, 100, label='$p_\mathrm{uni}(x)$', normed = 1, log=True)
plt.plot(x, p(x), label='$p(x)$', color = 'red', linewidth = 3)
plt.ylabel("$p(x)$", fontsize=20)
plt.xlabel("$x$",fontsize=20)
plt.legend(loc='upper right')
plt.savefig('rejection_histogram.pdf', bbox_inches='tight', transparent=True)

aux_histogram = np.loadtxt("rejection_aux_sample.txt")

plt.figure()
plt.hist(aux_histogram, 100, label='$p_\mathrm{aux}(x)$', normed = 1, log=True)
plt.plot(x, p(x), label='$p(x)$', color = 'red', linewidth = 3)
plt.ylabel("$p(x)$", fontsize=20)
plt.xlabel("$x$",fontsize=20)
plt.legend(loc='upper right')
plt.savefig('rejection_aux_histogram.pdf', bbox_inches='tight', transparent=True)


#plots the envelope function f(x)
def f(p1, p2):
	x = np.linspace(p1[0], p2[0], 2)
	return (p2[1] - p1[1])/(p2[0] - p1[0]) * (x - p1[0]) + p1[1]

p0 = [0, 0.01]; p1 = [1.8, 0.15]; p2 = [2.35, 2.5]; p3 = [3.0, 0.1]; p4 = [5.0, 0.002]

plt.figure()
plt.plot(x, p(x), label='$p(x)$', color = 'red', linewidth = 2)
plt.plot([p0[0],p1[0]], f(p0, p1), label='$f(x)$', color = 'blue', linewidth = 2)
plt.plot([p1[0],p2[0]], f(p1, p2), color = 'blue', linewidth = 2)
plt.plot([p2[0],p3[0]], f(p2, p3), color = 'blue', linewidth = 2)
plt.plot([p3[0],p4[0]], f(p3, p4), color = 'blue', linewidth = 2)
plt.ylabel("$p(x)$", fontsize=20)
plt.xlabel("$x$",fontsize=20)
plt.legend(loc='upper right')
plt.savefig('rejection_envelope.pdf', bbox_inches='tight', transparent=True)
