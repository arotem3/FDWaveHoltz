import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import LightSource
import matplotlib.tri as mtri

def load_u(it):
    p = np.fromfile("build/solution/u%06d"%it)
    return p

# --- ANIMATION ---
fig, ax = plt.subplots(1, 1, figsize=(20,2))

nframes = 40000
skip = 20
offset = 20000
nf = 1 + (nframes-offset) // skip

u = load_u(offset)
u_max = np.abs(u).max()

im, = ax.plot(u)
ax.axvline(0, c='k')

ax.set_xticks([])
ax.set_yticks([])
ax.set_ylim([-1.2*u_max, 1.2*u_max])
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.spines["bottom"].set_visible(False)
fig.tight_layout()

def update_fig(it):
    it = it * skip + offset
    print(f"{it:>5} / {nframes}", end='\r')
    p = load_u(it)

    im.set_ydata(p)

    return im,

FFwriter = animation.FFMpegWriter(fps=60)
ani = animation.FuncAnimation(fig, update_fig, frames=nf, blit=True)
ani.save('plots/wave1d_30.mov', writer=FFwriter)
print('\ndone')
