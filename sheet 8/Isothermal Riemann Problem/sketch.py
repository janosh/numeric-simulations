import matplotlib as mpl
font = {'family' : 'serif',
        'size'   : 16}
mpl.rc('font', **font)

import numpy as np
x = np.linspace(-10,10,100)

def rho(x):
    if x < -4:
        return 1.0  # add .0 to first return value to make whole function return floats instead of integers
    elif -4 <= x < -3.5:
        return 4 * (x + 4) + 1
    elif -3.5 <= x < 3.5:
        return 3
    elif 3.5 <= x < 4:
        return -4 * (x - 4) + 1
    elif 4 <= x:
        return 1

def u(x):
    if x < -4:
        return 1.0
    elif -4 <= x < -3.5:
        return -2 * (x + 4) + 1
    elif -3.5 <= x < 3.5:
        return 0
    elif 3.5 <= x < 4:
        return -2 * (x - 4) - 1
    elif 4 <= x:
        return -1

vrho = np.vectorize(rho)
vu = np.vectorize(u)


import matplotlib.pyplot as plt
plt.figure()
plt.plot(x, vrho(x), label=r'$\rho(x)$', linewidth=2)
plt.axvline(-3.75, label='shock', color = 'r', ls='dashed')
plt.axvline(3.75, color = 'r', ls='dashed')
plt.ylabel(r"$\rho/\rho_0$", fontsize=20)
plt.ylim(0,4)
plt.xlabel("$x$",fontsize=20)
plt.legend(loc='upper left')
plt.savefig('rho-sketch.pdf', bbox_inches='tight', transparent=True)

plt.figure()
plt.plot(x, vu(x), label='$u(x)$', linewidth=2)
plt.axvline(-3.75, label='shock', color = 'r', ls='dashed')
plt.axvline(3.75, color = 'r', ls='dashed')
plt.ylabel("$u/u_0$", fontsize=20)
plt.ylim(-2,2)
plt.xlabel("$x$",fontsize=20)
plt.legend(loc='upper left')
plt.savefig('u-sketch.pdf', bbox_inches='tight', transparent=True)
