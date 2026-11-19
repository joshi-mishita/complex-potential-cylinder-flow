"""
Step 5 - Spinning cylinder (Magnus effect).
Adds circulation Gamma to the potential; stagnation points move and
the flow becomes asymmetric, which produces lift.
Saves: figures/fig4_magnus_effect.png
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from flow_model import complex_potential, stagnation_points, make_grid, U_INF, RADIUS

os.makedirs("figures", exist_ok=True)

X, Y, Z = make_grid(ylim=(-3, 2), ny=500)
inside = np.abs(Z) < RADIUS
cases = [0.0, 2 * np.pi, 4 * np.pi, 6 * np.pi]   # Gamma / (4 pi U a) = 0, 0.5, 1, 1.5

fig, axes = plt.subplots(2, 2, figsize=(11, 9))
for ax, g in zip(axes.ravel(), cases):
    psi = complex_potential(Z, gamma=g).imag
    psi[inside] = np.nan
    ax.contour(X, Y, psi, levels=40, colors="tab:blue", linewidths=0.8)
    ax.add_patch(plt.Circle((0, 0), RADIUS, color="0.6", zorder=3))
    sp = stagnation_points(gamma=g)
    sp = sp[np.abs(sp) >= RADIUS - 1e-9]   # keep only points in the flow (outside the body)
    ax.plot(sp.real, sp.imag, "o", color="crimson", ms=7, zorder=4)
    ratio = g / (4 * np.pi * U_INF * RADIUS)
    ax.set_title(f"Γ = {g:.2f} m²/s   (Γ / 4πUa = {ratio:.1f})")
    ax.set_aspect("equal"); ax.set_xlim(-3, 3); ax.set_ylim(-3, 2)
    print(f"Gamma = {g:7.3f}: stagnation points = " +
          ", ".join(f"{p.real:+.3f}{p.imag:+.3f}i" for p in sp))
fig.suptitle("Magnus effect: circulation moves the stagnation points (red) and creates lift")
fig.tight_layout()
fig.savefig("figures/fig4_magnus_effect.png", dpi=200)
print("saved figures/fig4_magnus_effect.png")
