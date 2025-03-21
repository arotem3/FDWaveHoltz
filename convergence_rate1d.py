#%% setup python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker
import sklearn.linear_model as sklm
from numpy import pi

plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 14

data_dir = "build/solution"

#%% load simulation data
X = pd.read_csv(data_dir + "/convergence_rate1d.txt")
X['omega'] *= pi

logw = np.log(X['omega'].to_numpy()).reshape(-1, 1)

#%% omega vs parabolic distance
logrho = np.log(1-X['max beta'])
lm = sklm.LinearRegression()
lm.fit(logw, logrho)
a0 = np.exp(lm.intercept_)
b0 = lm.coef_[0]

print("ρ ~ %.2f × ω^(%.2f)"%(a0, b0))

logd = np.log(X['min parabolic distance'])
lm = sklm.LinearRegression()
lm.fit(logw, logd)
a1 = np.exp(lm.intercept_)
b1 = lm.coef_[0]

print("ε ~ %.2f × ω^(%.2f)"%(a1, b1))

plt.rcParams['font.size'] = 12
fig, ax = plt.subplots(1,1,figsize=(4, 4))

p1, = ax.loglog(X['omega'], X['min parabolic distance'], 'o-', label=r'$\varepsilon^\star$')
p2, = ax.loglog(X['omega'], 1-X['max beta'], 's-', label=r'$1 - \rho(\mathcal{S}_h)$')
p3, = ax.loglog(X['omega'], 1-X['rho estimate'], 'd-', ms=4, label=r'$1 - \hat{\rho}$')

l1, = ax.loglog(X['omega'], 1.1 * a1 * X['omega']**b1, '--', c="C0", label=r"$O(\omega^{%.2f})$"%b1)
l2, = ax.loglog(X['omega'], 0.9 * a0 * X['omega']**b0, '--', c="C1", label=r"$O(\omega^{%.2f})$"%b0)

ax.xaxis.set_minor_locator(ticker.NullLocator())
ax.xaxis.set_major_locator(ticker.FixedLocator([10*pi, 20*pi, 30*pi]))
ax.xaxis.set_major_formatter(ticker.FixedFormatter([r"$10\pi$", r"$20\pi$", r"$30\pi$"]))

leg1 = ax.legend(handles=[p1, l1], loc='lower left')
ax.add_artist(leg1)
ax.legend(handles=[p2, p3, l2], loc='upper right')
fig.tight_layout()
fig.savefig('plots/fd1d_parabolic_distance.pdf')
plt.show()

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

# %% omega vs. number of iterations

logit = np.log(X['FP#'])
lm = sklm.LinearRegression()
lm.fit(logw, logit)
a0 = np.exp(lm.intercept_)
b0 = lm.coef_[0]

print("FP# ~ %.2f × ω^(%.2f)"%(a0, b0))

logit = np.log(X['GMRES#'])
lm = sklm.LinearRegression()
lm.fit(logw, logit)
a1 = np.exp(lm.intercept_)
b1 = lm.coef_[0]

print("GMRES# ~ %.2f × ω^(%.2f)"%(a1, b1))

fig, ax = plt.subplots(1,1, figsize=(4,4))

p1, = ax.loglog(X['omega'], X['FP#'], 'o-', label=r"fixed-point iterations")
p2, = ax.loglog(X['omega'], X['GMRES#'], 's-', label=r"GMRES iterations")

l1, = ax.loglog(X['omega'], 1.1 * a0 * X['omega']**b0, '--', c="C0", label=r"$O(\omega^{%.2f})$"%b0)
l2, = ax.loglog(X['omega'], 1.1 * a1 * X['omega']**b1, '--', c="C1", label=r"$O(\omega^{%.2f})$"%b1)

ax.xaxis.set_minor_locator(ticker.NullLocator())
ax.xaxis.set_major_locator(ticker.FixedLocator([10*pi, 20*pi, 30*pi]))
ax.xaxis.set_major_formatter(ticker.FixedFormatter([r"$10\pi$", r"$20\pi$", r"$30\pi$"]))

ax.yaxis.set_minor_locator(ticker.NullLocator())
ax.yaxis.set_major_locator(ticker.FixedLocator([100, 200, 300, 400, 500, 600]))
ax.yaxis.set_major_formatter(ticker.FixedFormatter([r"$100$", r"$200$", r"$300$", r"$400$", r"$500$", r"$600$"]))

leg1 = ax.legend(handles=[p1, l1], loc='upper left')
ax.add_artist(leg1)
ax.legend(handles=[p2, l2], loc='lower right')
fig.tight_layout()
fig.savefig('plots/fd1d_omega_v_iters.pdf')
plt.show()
