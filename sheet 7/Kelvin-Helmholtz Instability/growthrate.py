import matplotlib as mpl
font = {'family' : 'serif',
        'size'   : 16}
mpl.rc('font', **font)

import numpy as np
time = np.loadtxt('Output 512/kh.hst', dtype=float, usecols=[0])
Ekiny = np.loadtxt('Output 512/kh.hst', dtype=float, usecols=[8])
timefit = time[7:29]
Ekinyfit = Ekiny[7:29]

from scipy.optimize import curve_fit
k = 4 * np.pi
u1 = 0.3
u2 = -0.3
rho1 = 1.0
rho2 = 2.0
omega = k * abs(u1 - u2) * np.sqrt(rho1 * rho2)/(rho1 + rho2)
def growthline(t, a):
    return a * np.exp(omega * t)

params, cov = curve_fit(growthline, timefit, Ekinyfit)

import matplotlib.pyplot as plt
plt.figure()
plt.plot(time, Ekiny, label='$E_\mathrm{kin,y}$')
plt.plot(time, growthline(time, *params), label='$E_\mathrm{fit} \propto \, e^{\omega t}$')
plt.yscale('log')
plt.ylabel("$E_\mathrm{kin}$",fontsize=20)
plt.xlabel("$t$",fontsize=20)
plt.legend(loc='upper left')
plt.savefig('growth-rate.pdf', bbox_inches='tight', transparent=True)
