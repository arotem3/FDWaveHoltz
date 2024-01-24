#%% setup python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker
import sklearn.linear_model as sklm
from numpy import pi

plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 20

#%% load simulation data
X = pd.read_csv("solution/convergence_rate1d.txt")
X['omega'] *= pi

logw = np.log(X['omega'].to_numpy()).reshape(-1, 1)

#%% omega vs rate (first iteration)
fig, ax = plt.subplots(1,1,figsize=(8,8))

logmu = np.log(1-X['mu'])
lm = sklm.LinearRegression()
lm.fit(logw, logmu)
bmu = lm.coef_[0]

print("mu ~ %.2f omega^(%.2f)"%(np.exp(lm.intercept_), bmu))

ax.loglog(X['omega'],1 - X['e'], 'o-', label=r'$1 - \frac{\|e_h^{(1)}\|}{\|e_h^{(0)}\|}$')
ax.loglog(X['omega'],1 - X['mu'], 'd-', label=r'$1 - \frac{\|\mu_h^{(1)}\|}{\|\mu_h^{(0)}\|}$')
ax.loglog(X['omega'],1 - X['max beta'], '-s', label=r'$1 - \max_\lambda |\beta(\lambda)|$')

ax.loglog(X['omega'], 15 / X['omega']**2, '--', c='C0', label=r"$O(\omega^{-2})$")
ax.loglog(X['omega'], 1.5 * X['omega']**bmu, '--', c='C1', label=r"$O(\omega^{%0.2f})$"%bmu)
ax.loglog(X['omega'], 0.7 * X['omega']**-0.72, '--', c='C2', label=r"$O(\omega^{-0.72})$")

ax.set_ylim([1e-3, 1e-1])

ax.xaxis.set_minor_locator(ticker.NullLocator())
ax.xaxis.set_major_locator(ticker.FixedLocator([10*pi, 20*pi, 30*pi]))
ax.xaxis.set_major_formatter(ticker.FixedFormatter([r"$10\pi$", r"$20\pi$", r"$30\pi$"]))

ax.legend(ncol=2)
fig.tight_layout()
fig.savefig("plots/omega_v_first_rate.pdf")
# fig.show()

#%% omega vs rate (avg)
fig, ax = plt.subplots(1,1,figsize=(8,8))

logmu = np.log(1-X['avg rate (e)'])
lm = sklm.LinearRegression()
lm.fit(logw, logmu)
bmu = lm.coef_[0]

print("avg rate (e) ~ %.2f omega^(%.2f)"%(np.exp(lm.intercept_), bmu))

ax.loglog(X['omega'],1 - X['avg rate (e)'], 'o-', label=r'$1 - \mathrm{avg}\frac{\|e_h^{(n+1)}\|}{\|e_h^{(n)}\|}$')
ax.loglog(X['omega'],1 - X['avg rate (mu)'], 'd-', label=r'$1 - \mathrm{avg}\frac{\|\mu_h^{(n+1)}\|}{\|\mu_h^{(n)}\|}$')
ax.loglog(X['omega'], 0.82 * X['omega']**(bmu), '--k', label=r"$O(\omega^{%0.2f})$"%bmu)

ax.set_ylim([1e-2, 1e-1])

ax.xaxis.set_minor_locator(ticker.NullLocator())
ax.xaxis.set_major_locator(ticker.FixedLocator([10*pi, 20*pi, 30*pi]))
ax.xaxis.set_major_formatter(ticker.FixedFormatter([r"$10\pi$", r"$20\pi$", r"$30\pi$"]))

ax.yaxis.set_minor_locator(ticker.NullLocator())
ax.yaxis.set_major_locator(ticker.FixedLocator([1e-2,1e-1]))

ax.legend(ncol=2)
fig.tight_layout()
# fig.savefig("plots/omega_v_avg_rate.pdf")
fig.show()

#%% omega vs parabolic distance
fig, ax = plt.subplots(1,1,figsize=(8,8))

logd = np.log(X['min parabolic distance'])
lm = sklm.LinearRegression()
lm.fit(logw, logd)
b = lm.coef_[0]

print("min distance ~ %.2f omega^(%.2f)"%(np.exp(lm.intercept_), b))

ax.loglog(X['omega'], X['min parabolic distance'], 'o-', label=r'$\min_\lambda \{ -x + \alpha (y - 1)^2 \}$')
ax.loglog(X['omega'], 1-X['max beta'], 's-', label=r'$1 - \max_\lambda |\beta(\lambda)|$')

ax.loglog(X['omega'], 0.3 * X['omega']**b, '--', c="C0", label=r"$O(\omega^{%.2f})$"%b)
ax.loglog(X['omega'], 0.7 * X['omega']**-0.72, '--', c="C1", label=r"$O(\omega^{-0.72})$")

ax.xaxis.set_minor_locator(ticker.NullLocator())
ax.xaxis.set_major_locator(ticker.FixedLocator([10*pi, 20*pi, 30*pi]))
ax.xaxis.set_major_formatter(ticker.FixedFormatter([r"$10\pi$", r"$20\pi$", r"$30\pi$"]))

ax.legend(ncols=2, loc='lower left')
fig.tight_layout()
fig.savefig('plots/omega_v_parabolic_distance.pdf')
# fig.show()

print(bmu**2 / abs(b))

#%% omega vs kappa
fig, ax = plt.subplots(1,1, figsize=(8,8))

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
# fig.savefig('plots/fd1d_kappa.pdf')
plt.show()

# %%
