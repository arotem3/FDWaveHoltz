#%% setup python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker
import sklearn.linear_model as sklm
from numpy import pi

plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 16

data_dir = "build/solution"

#%% iteration vs rel error
fig, ax = plt.subplots(1,1, figsize=(5, 2*5/3))

locs = ['lower left', 'lower right', 'upper right']
for k, w in enumerate([10, 20, 30]):
    X = pd.read_csv(f"{data_dir}/iters{w}.txt")
    Y = np.loadtxt(f"{data_dir}/iters{w}_gmres.txt")

    n = np.arange(1, len(X)+1)

    p1, = ax.semilogy(n, X['e'], '-', c=f"C{k}", label=r'$\frac{\| e_h^{(n)} \| }{\| e_h^{(0)} \| } (\omega = %d \pi)$'%w)
    l1, = ax.semilogy(n, X['mu'], '--', c=f"C{k}", label=r'$\frac{\| \mu_h^{(n)} \| }{\| \mu_h^{(0)} \| } (\omega = %d \pi)$'%w)

    n = np.arange(1, len(Y)+1)
    p2, = ax.semilogy(n, Y/Y[0], ':', c=f"C{k}", label=r'GMRES rel. res. $(\omega = %d \pi)$'%w)

    # leg = ax.legend(handles=[p1, l1], loc=locs[k])
    # ax.add_artist(leg)

# ax.set_xlim([1,3])
# ax.set_ylim([7e-1,1])
# ax.legend(ncols=3)
fig.tight_layout()
fig.savefig("plots/fd1d_rel_error.pdf")
plt.show()
# %%
