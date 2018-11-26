import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


font = {'family' : 'serif',
        'size'   : 16}
mpl.rc('font', **font)


# calculation of the heat capacities at T = 75 K and 400 K
relaxed_Etot_at70K = np.loadtxt('70K.txt', skiprows=50000, usecols=[3])
avr_relaxed_Etot_at70K = np.average(relaxed_Etot_at70K)
relaxed_Etot_at80K = np.loadtxt('80K.txt', skiprows=50000, usecols=[3])
avr_relaxed_Etot_at80K = np.average(relaxed_Etot_at80K)

relaxed_Etot_at395K = np.loadtxt('395K.txt', skiprows=50000, usecols=[3])
avr_relaxed_Etot_at395K = np.average(relaxed_Etot_at395K)
relaxed_Etot_at405K = np.loadtxt('405K.txt', skiprows=50000, usecols=[3])
avr_relaxed_Etot_at405K = np.average(relaxed_Etot_at405K)

epsilon = 120.0
Ntot = 512.0
Delta_T = 10.0
heat_capa = epsilon * (avr_relaxed_Etot_at80K - avr_relaxed_Etot_at70K) / (Ntot * Delta_T)
print "C_V(T=75K) = ", heat_capa
heat_capa = epsilon * (avr_relaxed_Etot_at405K - avr_relaxed_Etot_at395K) / (Ntot * Delta_T)
print "C_V(T=400K) = ", heat_capa


# plotting the total energy at T = 80 K
Etot80K = np.loadtxt('80K_novscale.txt', usecols=[3])/512.0 # we divide each entry by the number of particles
step = np.arange(0, len(Etot80K), 1)
dt = 0.01

plt.figure()

plt.plot(dt * step, Etot80K, label='$E_\mathrm{tot}(t)$')

plt.xlabel("$t$",fontsize=20)
plt.ylabel("$E(t)$",fontsize=20)
plt.ylim(0.995 * min(Etot80K), 1.005 * max(Etot80K))
plt.legend(loc='lower right')

plt.savefig('Etot_at80K.pdf', bbox_inches='tight', transparent=True)


# plotting all energies at T = 30 K
Ekin = np.loadtxt('30K.txt', usecols=[1])/512.0 # we divide each entry by the number of particles
Epot = np.loadtxt('30K.txt', usecols=[2])/512.0
Etot30K = np.loadtxt('30K.txt', usecols=[3])/512.0

plt.figure()

plt.plot(dt * step, Ekin, label='$E_\mathrm{kin}(t)$')
plt.plot(dt * step, Epot, label='$E_\mathrm{pot}(t)$')
plt.plot(dt * step, Etot30K, label='$E_\mathrm{tot}(t)$')

plt.xlabel("$t$",fontsize=20)
plt.ylabel("$E(t)$",fontsize=20)
plt.legend(loc='lower left')

plt.savefig('energies_at30K.pdf', bbox_inches='tight', transparent=True)
