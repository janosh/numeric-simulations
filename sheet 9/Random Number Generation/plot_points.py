import matplotlib as mpl
font = {'family' : 'serif',
        'size'   : 16}
mpl.rc('font', **font)

import numpy as np
randu = np.loadtxt("randu.txt")
randu_x = randu[:,0]
randu_y = randu[:,1]

import matplotlib.pyplot as plt
plt.figure()
plt.scatter(randu_x, randu_y)
plt.ylabel("$y$", fontsize=20)
plt.xlabel("$x$",fontsize=20)
plt.xlim(0,1)
plt.ylim(0,1)
plt.savefig('randu.pdf', bbox_inches='tight', transparent=True)


drand = np.loadtxt("drand.txt")
drand_x = drand[:,0]
drand_y = drand[:,1]

plt.figure()
plt.scatter(drand_x, drand_y)
plt.ylabel("$y$", fontsize=20)
plt.xlabel("$x$",fontsize=20)
plt.xlim(0,1)
plt.ylim(0,1)
plt.savefig('drand.pdf', bbox_inches='tight', transparent=True)


zoom_randu = np.loadtxt("randu_zoom.txt")
zoom_randu_x = zoom_randu[:,0]
zoom_randu_y = zoom_randu[:,1]

plt.figure()
plt.scatter(zoom_randu_x, zoom_randu_y)
plt.ylabel("$y$", fontsize=20)
plt.xlabel("$x$",fontsize=20)
plt.xlim(0.2,0.2005)
plt.ylim(0.3,0.3005)
plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
plt.savefig('randu-zoom.pdf', bbox_inches='tight', transparent=True)


drand_zoom = np.loadtxt("drand_zoom.txt")
drand_zoom_x = drand_zoom[:,0]
drand_zoom_y = drand_zoom[:,1]

plt.figure()
plt.scatter(drand_zoom_x, drand_zoom_y)
plt.ylabel("$y$", fontsize=20)
plt.xlabel("$x$",fontsize=20)
plt.xlim(0.2,0.2005)
plt.ylim(0.3,0.3005)
plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
plt.savefig('drand-zoom.pdf', bbox_inches='tight', transparent=True)
