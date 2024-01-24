#%% setup python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker
import sklearn.linear_model as sklm
from numpy import pi

plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 20

#%% iteration vs rel error
fig, ax = plt.subplots(1,1, figsize=(12,8))

n = np.arange(1, 1000+1)
for k, w in enumerate([10,20,30]):
    X = pd.read_csv(f"solution/iters{w}.txt")

    ax.semilogy(n, X['e'], '-', c=f"C{k}", label=r'$\frac{\| e_h^{(n)} \| }{\| e_h^{(0)} \| } (\omega = %d \pi)$'%w)
    ax.semilogy(n, X['mu'], '--', c=f"C{k}", label=r'$\frac{\| \mu_h^{(n)} \| }{\| \mu_h^{(0)} \| } (\omega = %d \pi)$'%w)

# ax.set_xlim([1,3])
# ax.set_ylim([7e-1,1])
ax.legend(ncols=3)
fig.tight_layout()
# fig.savefig("plots/iters_v_rel_error.pdf")
# %%
