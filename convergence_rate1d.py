import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker
import sklearn.linear_model as sklm
from numpy import pi

plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 20

fig, ax = plt.subplots(1,1,figsize=(8,8))

X = pd.read_csv("solution/convergence_rate1d.txt")
X['omega'] *= pi

logw = np.log(X['omega'].to_numpy()).reshape(-1, 1)
logmu = np.log(1-X['mu'])
lm = sklm.LinearRegression()
lm.fit(logw, logmu)
bmu = lm.coef_[0]

print("mu ~ %.2f omega^(%.2f)"%(np.exp(lm.intercept_), bmu))

ax.loglog(X['omega'],1 - X['e'], 'o-', label=r'$1 - \frac{\|e_h^{(1)}\|}{\|e_h^{(0)}\|}$')
ax.loglog(X['omega'],1 - X['mu'], 'd-', label=r'$1 - \frac{\|\mu_h^{(1)}\|}{\|\mu_h^{(0)}\|}$')
ax.loglog(X['omega'], 1.5 * X['omega']**bmu, '--k', label=r"$O(\omega^{%0.2f})$"%bmu)
ax.loglog(X['omega'], 15 / X['omega']**2, '--r', label=r"$O(\omega^{-2})$")

ax.xaxis.set_minor_locator(ticker.NullLocator())
ax.xaxis.set_major_locator(ticker.FixedLocator([10*pi, 20*pi, 30*pi]))
ax.xaxis.set_major_formatter(ticker.FixedFormatter([r"$10\pi$", r"$20\pi$", r"$30\pi$"]))

ax.legend(ncol=2)
fig.tight_layout()
fig.savefig("plots/omega_v_first_rate.pdf")

fig, ax = plt.subplots(1,1,figsize=(8,8))

logd = np.log(X['min parabolic distance'])
lm = sklm.LinearRegression()
lm.fit(logw, logd)
b = lm.coef_[0]

print("min distance ~ %.2f omega^(%.2f)"%(np.exp(lm.intercept_), b))

ax.loglog(X['omega'], X['min parabolic distance'], 'o-', label=r'$\min_\lambda \{ -x + \alpha (y - 1)^2 \}$')
ax.loglog(X['omega'], 0.3 * X['omega']**b, '--k', label=r"$O(\omega^{%.2f})$"%b)

ax.xaxis.set_minor_locator(ticker.NullLocator())
ax.xaxis.set_major_locator(ticker.FixedLocator([10*pi, 20*pi, 30*pi]))
ax.xaxis.set_major_formatter(ticker.FixedFormatter([r"$10\pi$", r"$20\pi$", r"$30\pi$"]))

ax.legend()
fig.tight_layout()
fig.savefig('plots/omega_v_parabolic_distance.pdf')

print(bmu**2 / abs(b))

fig, ax = plt.subplots(1,1, figsize=(8,8))

logk = np.log(X['kappa'])
lm = sklm.LinearRegression()
lm.fit(logw, logk)
bk = lm.coef_[0]

print("kappa ~ %.2f omega^(%.2f)"%(np.exp(lm.intercept_), bk))
print(bk)

ax.loglog(X['omega'], X['kappa'], 'o-', label=r"$\kappa$")
ax.loglog(X['omega'], 0.3 * X['omega']**bk, '--k', label=r"$O(\omega^{%.2f})$"%bk)

ax.xaxis.set_minor_locator(ticker.NullLocator())
ax.xaxis.set_major_locator(ticker.FixedLocator([10*pi, 20*pi, 30*pi]))
ax.xaxis.set_major_formatter(ticker.FixedFormatter([r"$10\pi$", r"$20\pi$", r"$30\pi$"]))

ax.legend()
fig.tight_layout()
fig.savefig('plots/fd1d_kappa.pdf')
