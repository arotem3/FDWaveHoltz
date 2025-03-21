#%% setup python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker
import sklearn.linear_model as sklm
from numpy import pi

plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 14

#%% load simulation data
X = pd.read_csv("solution/convergence_rate1d.txt")
X['omega'] *= pi

logw = np.log(X['omega'].to_numpy()).reshape(-1, 1)

#%% omega vs rate (first iteration)
fig, ax = plt.subplots(1,1,figsize=(4, 4))

logmu = np.log(1-X['mu'])
lm = sklm.LinearRegression()
lm.fit(logw, logmu)
bmu = lm.coef_[0]

print("mu ~ %.2f omega^(%.2f)"%(np.exp(lm.intercept_), bmu))

p1, = ax.loglog(X['omega'],1 - X['e'], 'o-', ms=4, label=r'$1 - \frac{\|e_h^{(1)}\|}{\|e_h^{(0)}\|}$')
p2, = ax.loglog(X['omega'],1 - X['mu'], 'd-', ms=4, label=r'$1 - \frac{\|\mu_h^{(1)}\|}{\|\mu_h^{(0)}\|}$')
p3, = ax.loglog(X['omega'],1 - X['max beta'], '-s', ms=4, label=r'$1 - \rho(\mathcal{S}_h)$')

l1, = ax.loglog(X['omega'], 15 / X['omega']**2, '--', c='C0', label=r"$O(\omega^{-2})$")
l2, = ax.loglog(X['omega'], 1.5 * X['omega']**bmu, '--', c='C1', label=r"$O(\omega^{%0.2f})$"%bmu)
l3, = ax.loglog(X['omega'], 0.7 * X['omega']**-0.72, '--', c='C2', label=r"$O(\omega^{-0.72})$")

ax.set_ylim([1e-3, 2e-1])

ax.xaxis.set_minor_locator(ticker.NullLocator())
ax.xaxis.set_major_locator(ticker.FixedLocator([10*pi, 20*pi, 30*pi]))
ax.xaxis.set_major_formatter(ticker.FixedFormatter([r"$10\pi$", r"$20\pi$", r"$30\pi$"]))

leg1 = ax.legend(handles=[p1, l1], loc='lower left')
ax.add_artist(leg1)
leg2 = ax.legend(handles=[p2, l2], loc='upper right')
ax.add_artist(leg2)
ax.legend(handles=[p3, l3], loc='center right')
fig.tight_layout()
fig.savefig("plots/omega_v_first_rate.pdf")
# fig.show()

#%% omega vs rate (avg)
fig, ax = plt.subplots(1,1,figsize=(4,4))

logmu = np.log(1-X['avg rate (e)'])
lm = sklm.LinearRegression()
lm.fit(logw, logmu)
bmu = lm.coef_[0]

print("avg rate (e) ~ %.2f omega^(%.2f)"%(np.exp(lm.intercept_), bmu))

ax.loglog(X['omega'],1 - X['avg rate (e)'], 'o-', ms=4, label=r'$1 - \bar{r}_e$')
ax.loglog(X['omega'],1 - X['avg rate (mu)'], 'd-', ms=4, label=r'$1 - \bar{r}_\mu$')
ax.loglog(X['omega'], 0.82 * X['omega']**(bmu), '--k', label=r"$O(\omega^{%0.2f})$"%bmu)

ax.set_ylim([1e-2, 1e-1])

ax.xaxis.set_minor_locator(ticker.NullLocator())
ax.xaxis.set_major_locator(ticker.FixedLocator([10*pi, 20*pi, 30*pi]))
ax.xaxis.set_major_formatter(ticker.FixedFormatter([r"$10\pi$", r"$20\pi$", r"$30\pi$"]))

ax.yaxis.set_minor_locator(ticker.NullLocator())
ax.yaxis.set_major_locator(ticker.FixedLocator([1e-2,1e-1]))

ax.legend()
fig.tight_layout()
fig.savefig("plots/omega_v_avg_rate.pdf")
# fig.show()

#%% omega vs parabolic distance
plt.rcParams['font.size'] = 12
fig, ax = plt.subplots(1,1,figsize=(4, 4))

logd = np.log(X['min parabolic distance'])
lm = sklm.LinearRegression()
lm.fit(logw, logd)
b = lm.coef_[0]

print("min distance ~ %.2f omega^(%.2f)"%(np.exp(lm.intercept_), b))

p1, = ax.loglog(X['omega'], X['min parabolic distance'], 'o-', label=r'$\varepsilon^\star$')
p2, = ax.loglog(X['omega'], 1-X['max beta'], 's-', label=r'$1 - \rho(\mathcal{S}_h)$')

l1, = ax.loglog(X['omega'], 0.3 * X['omega']**b, '--', c="C0", label=r"$O(\omega^{%.2f})$"%b)
l2, = ax.loglog(X['omega'], 0.7 * X['omega']**-0.72, '--', c="C1", label=r"$O(\omega^{-0.72})$")

ax.xaxis.set_minor_locator(ticker.NullLocator())
ax.xaxis.set_major_locator(ticker.FixedLocator([10*pi, 20*pi, 30*pi]))
ax.xaxis.set_major_formatter(ticker.FixedFormatter([r"$10\pi$", r"$20\pi$", r"$30\pi$"]))

leg1 = ax.legend(handles=[p1, l1], loc='lower left')
ax.add_artist(leg1)
ax.legend(handles=[p2, l2], loc='upper right')
fig.tight_layout()
fig.savefig('plots/omega_v_parabolic_distance.pdf')
fig.show()

#%% omega vs kappa
fig, ax = plt.subplots(1,1, figsize=(4,4))

logk = np.log(X['kappa'])
lm = sklm.LinearRegression()
lm.fit(logw, logk)
bk = lm.coef_[0]

print("kappa ~ %.2f omega^(%.2f)"%(np.exp(lm.intercept_), bk))

ax.loglog(X['omega'], X['kappa'], 'o-', label=r"$\kappa$")
ax.loglog(X['omega'], 0.3 * X['omega']**bk, '--k', label=r"$O(\omega^{%.2f})$"%bk)

ax.xaxis.set_minor_locator(ticker.NullLocator())
ax.xaxis.set_major_locator(ticker.FixedLocator([10*pi, 20*pi, 30*pi]))
ax.xaxis.set_major_formatter(ticker.FixedFormatter([r"$10\pi$", r"$20\pi$", r"$30\pi$"]))

ax.legend()
fig.tight_layout()
fig.savefig('plots/fd1d_kappa.pdf')
# plt.show()

# %%
