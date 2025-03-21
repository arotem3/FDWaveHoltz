#%% setup python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker
import sklearn.linear_model as sklm
from numpy import pi

plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 20

data_dir = "build/solution"

#%% iteration vs rel error
fig, ax = plt.subplots(1,1, figsize=(5, 4))

for k, w in enumerate([10, 20, 30]):
    X = np.loadtxt(f"{data_dir}/iter2d_{w}.txt")
    Y = np.loadtxt(f"{data_dir}/gmres2d_{w}.txt")

    n = np.arange(1, len(X)+1)
    ax.semilogy(n, X, '-', c=f"C{k}", label=r'$r^{(n)} (\omega = %d \pi)$'%w)

    n = np.arange(1, len(Y)+1)
    ax.semilogy(n, Y/Y[0], ':', c=f"C{k}", label=r'GMRES rel. res. $(\omega = %d \pi)$'%w)
    

# ax.set_xlim([1,3])
# ax.set_ylim([7e-1,1])
# ax.legend(ncols=3)
fig.tight_layout()
fig.savefig("plots/fd2d_rel_error.pdf")
plt.show()
# %%
