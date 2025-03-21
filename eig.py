import numpy as np
import matplotlib.pyplot as plt

def parabolic_distance(x, y):
    alpha = (2 * np.pi**2 - 3) / (12 * np.pi)
    return -x + alpha * (y - 1)**2

w = 4
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 18

fig, ax = plt.subplots(1,1, figsize=(w, w))

X = np.linspace(-0.1, 0, 200)
Y = np.linspace(0, 2, 200)
X, Y = np.meshgrid(X, Y)
Z = parabolic_distance(X, Y)
cs = ax.contour(X, Y, Z, [0.02, 0.04, 0.06], colors='k')

for i, w in enumerate([10, 20, 30]):
    omega = w * np.pi
    
    x, y = np.loadtxt(f"solution/eig_{w}.txt").T
    x = x / omega
    y = y / omega

    d = parabolic_distance(x, y)
    j = np.argmin(d)

    ax.plot(x, y, '.', c=f"C{i}", ms=3, label=r'$\omega=%d\pi$'%w)
    ax.plot(x[j], y[j], 'd', c=f"C{i}", ms=10, markerfacecolor="None", markeredgewidth=2)

ax.axvline(0, 0, 100, ls='--', c='k')#, label='imaginary axis')
# ax.set_xlabel(r"$\Re\{\lambda/\omega\}$")
# ax.set_ylabel(r"$\Im\{\lambda/\omega\}$")
ax.set_ylim([0, 5])
ax.set_xlim([None, 0.005])
ax.legend()
ax.clabel(cs, cs.levels, inline=True, fontsize=10)
fig.savefig('plots/fd_eigs.pdf')