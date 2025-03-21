import numpy as np
import matplotlib.pyplot as plt

for w in [10, 20, 30]:
    u = np.fromfile(f'solution/u{w}')
    n = int(np.sqrt(u.shape[0] // 2))
    u = np.reshape(u, (n, n, 2), 'F')
    E = np.hypot(u[...,0], u[...,1]).T

    x = np.linspace(-1, 1, n)

    # fig, axes = plt.subplots(1, 2, figsize=(16,8))
    # axes[0].contourf(x, x, u[...,0], vmin=-0.15, vmax=0.15, cmap='seismic')
    # axes[1].contourf(x, x, u[...,1], vmin=-0.15, vmax=0.15, cmap='seismic')

    # axes[1].set_yticks([])

    # fig.tight_layout()

    fig, ax = plt.subplots(1, 1, figsize=(8,8))
    # ax.contourf(x, x, E, 50, cmap='Blues')
    ax.imshow(E, cmap='magma', origin='lower', interpolation='bicubic', extent=[-1,1,-1,1])
    ax.set_xticks([])
    ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(f"plots/fd2d_sol{w}.pdf")
    # plt.show()