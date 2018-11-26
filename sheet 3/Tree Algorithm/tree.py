import matplotlib as mpl
font = {'family' : 'serif',
        'size'   : 16}
mpl.rc('font', **font)

import numpy as np
data_N = np.fromfile('N.dat', dtype=np.int32, count=-1, sep='')
data_tree = np.fromfile('tree.dat', dtype=float, count=-1, sep='')
data_dir = np.fromfile('direct.dat', dtype=float, count=-1, sep='')

import matplotlib.pyplot as plt
plt.figure()

log_N = np.log(data_N)
log_tree = np.log(data_tree)
log_dir = np.log(data_dir)

slope_tree, intercept_tree = np.polyfit(log_N, log_tree, 1)
slope_dir, intercept_dir = np.polyfit(log_N, log_dir, 1)
print "m_tree = %g\nb_tree = %g\nm_dir = %g\nb_dir = %g" % (slope_tree, intercept_tree, slope_dir, intercept_dir)

def tree(x):
    return np.exp(slope_tree * np.log(x) + intercept_tree)
def direct(x):
    return np.exp(slope_dir * np.log(x) + intercept_dir)

points = np.arange(1, 1e10, 1e3)
plt.plot(data_N, data_dir, label='$t_\mathrm{dir}(N)$', linewidth=4)
plt.plot(points, direct(points), '--', color = 'b', label='$t_\mathrm{dir,fit}(N)$', linewidth=2)
plt.plot(data_N, data_tree, label='$t_\mathrm{tree}(N)$', linewidth=4)
plt.plot(points, tree(points), '--', color = 'g', label='$t_\mathrm{tree,fit}(N)$', linewidth=2)


plt.xlabel("$N$",fontsize=20)
plt.ylabel("$t(N) \, [s]$",fontsize=20)
plt.yscale('log')
plt.xscale('log')
plt.xlim(min(points),max(points))
plt.legend(loc='lower right')

plt.savefig('projection.pdf', bbox_inches='tight', transparent=True)
