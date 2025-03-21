#%% setup python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker
import sklearn.linear_model as sklm
from numpy import pi

plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 16

#%%

def power_fit(X, y):
    lm = sklm.LinearRegression()
    lm.fit(np.log(np.array(X)).reshape(-1,1), np.log(np.array(y)))
    a = np.exp(lm.intercept_)
    b = lm.coef_[0]
    print(f"y ~ %.2f * x^%.2f"%(a,b))
    return a,b

#%% load simulation data
X = pd.read_csv("solution/convergence_rate2d.txt")
X.sort_values('omega', inplace=True)

#%% omega vs iter

fig, ax = plt.subplots(1,1, figsize=(4, 4))

a, b = power_fit(X['omega'], X['iter'])

ax.loglog(X['omega'], X['iter'], '-o', label=r'$\mathrm{Iterations}$')
fit=ax.loglog(X['omega'], 1.1 * a * X['omega']**b, '--', c='C0', label=r"$O(\omega^{%.2f})$"%b)

ax.xaxis.set_minor_locator(ticker.NullLocator())
ax.xaxis.set_major_locator(ticker.FixedLocator([10*pi, 20*pi, 30*pi]))
ax.xaxis.set_major_formatter(ticker.FixedFormatter([r"$10\pi$", r"$20\pi$", r"$30\pi$"]))

ax.yaxis.set_minor_locator(ticker.NullLocator())
ax.yaxis.set_major_locator(ticker.FixedLocator([300, 400, 500, 600]))
ax.yaxis.set_major_formatter(ticker.FixedFormatter([r"$300$", r"$400$", r"$500$", r"$600$"]))

ax.legend(handles=fit)

fig.tight_layout()
fig.savefig("plots/fd2d_omega_v_iters.pdf")

# %%

