import numpy as np
import matplotlib.pyplot as plt

x = np.fromfile("solution/x")
n = len(x)

fig, axes = plt.subplots(4, 4, figsize=(8,8))

nt = 8000
skip = 10

K = (nt * skip // 16) // skip
for i in range(4):
    for j in range(4):
        k = (j + 4*i) * K
        u = np.fromfile("solution/u%06d"%k).reshape(n, n)
        v = np.abs(u).max()

        # axes[i, j].contour(x, x, u, 20, cmap='seismic', vmin=-v, vmax=v)
        axes[i, j].imshow(u, cmap='seismic', vmin=-v, vmax=v, origin='lower')
        axes[i, j].set_xticks([])
        axes[i, j].set_yticks([])
plt.tight_layout()
plt.show()