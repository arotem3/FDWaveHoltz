import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import LightSource
import matplotlib.tri as mtri

def load_u(it):
    p = np.fromfile("solution/u%06d"%it)
    p = np.reshape(p, (500,500), "F")
    return p

# --- ANIMATION ---
fig, ax = plt.subplots(1, 1, figsize=(20,2))

nframes = 40000
skip = 20
offset = 10000
nf = 1 + (nframes-offset) // skip

u = load_u(offset)[50:100]

im = ax.imshow(u, extent=[-1, 1, -1, -0.8], cmap='seismic', vmin=-0.2, vmax=0.2, origin='lower')

ax.set_xticks([])
ax.set_yticks([])
ax.axis('equal')
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.spines["bottom"].set_visible(False)
fig.tight_layout()

def update_fig(it):
    it = it * skip + offset
    print(f"{it:>5} / {nframes}", end='\r')
    p = load_u(it)[50:100]

    im.set_data(p)

    return im,

FFwriter = animation.FFMpegWriter(fps=60)
ani = animation.FuncAnimation(fig, update_fig, frames=nf, blit=True)
ani.save('plots/wave.mov', writer=FFwriter)
print('\ndone')
