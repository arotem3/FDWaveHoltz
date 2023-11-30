import numpy as np
import matplotlib.pyplot as plt

x = np.fromfile("solution/x")

fig, axes = plt.subplots(4, 4, figsize=(8,8))

nt = 8000
skip = 10

K = (nt * skip // 16) // skip
for i in range(4):
    for j in range(4):
        k = (j + 4*i) * K
        u = np.fromfile("solution/u%06d"%k)

        axes[i, j].plot(x, u)
        axes[i, j].set_ylim([-0.15, 0.15])
        axes[i, j].set_xticks([])
        axes[i, j].set_yticks([])
plt.tight_layout()
plt.show()