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

#%%

def power_fit(X, y):
    lm = sklm.LinearRegression()
    lm.fit(np.log(np.array(X)).reshape(-1,1), np.log(np.array(y)))
    a = np.exp(lm.intercept_)
    b = lm.coef_[0]
    print(f"y ~ %.2f * x^%.2f"%(a,b))
    return a,b

#%% load simulation data
X = pd.read_csv(data_dir + "/convergence_rate2d.txt")
X.sort_values('omega', inplace=True)

#%% omega vs iter

fig, ax = plt.subplots(1,1, figsize=(4, 4))

a, b = power_fit(X['omega'], X['FP#'])

p1, = ax.loglog(X['omega'], X['FP#'], '-o', c='C0', label=r'fixed-point iterations')
l1, = ax.loglog(X['omega'], 1.1 * a * X['omega']**b, '--', c='C0', label=r"$O(\omega^{%.2f})$"%b)

a, b = power_fit(X['omega'], X['GMRES#'])

p2, = ax.loglog(X['omega'], X['GMRES#'], '-s', c='C1', label=r'GMRES iterations')
l2, = ax.loglog(X['omega'], 1.1 * a * X['omega']**b, '--', c='C1', label=r"$O(\omega^{%.2f})$"%b)

ax.xaxis.set_minor_locator(ticker.NullLocator())
ax.xaxis.set_major_locator(ticker.FixedLocator([10*pi, 20*pi, 30*pi]))
ax.xaxis.set_major_formatter(ticker.FixedFormatter([r"$10\pi$", r"$20\pi$", r"$30\pi$"]))

ax.yaxis.set_minor_locator(ticker.NullLocator())
ax.yaxis.set_major_locator(ticker.FixedLocator([200, 300, 400, 500, 600]))
ax.yaxis.set_major_formatter(ticker.FixedFormatter([r"$200$", r"$300$", r"$400$", r"$500$", r"$600$"]))

leg1 = ax.legend(handles=[p1, l1], loc='upper left')
ax.add_artist(leg1)
ax.legend(handles=[p2, l2], loc='lower right')

fig.tight_layout()
fig.savefig("plots/fd2d_iters.pdf")

# %% omega vs rho

a, b = power_fit(X['omega'], 1-X['rho'])

fig, ax = plt.subplots(1,1, figsize=(4, 4))

p1, = ax.loglog(X['omega'], 1-X['rho'], '-o', c='C0', label=r'$1 - \hat{\rho}$')
l1, = ax.loglog(X['omega'], 1.1 * a * X['omega']**b, '--', c='C0', label=r"$O(\omega^{%.2f})$"%b)

ax.xaxis.set_minor_locator(ticker.NullLocator())
ax.xaxis.set_major_locator(ticker.FixedLocator([10*pi, 20*pi, 30*pi]))
ax.xaxis.set_major_formatter(ticker.FixedFormatter([r"$10\pi$", r"$20\pi$", r"$30\pi$"]))

ax.legend()
fig.tight_layout()
fig.savefig("plots/fd2d_rho.pdf")
plt.show()