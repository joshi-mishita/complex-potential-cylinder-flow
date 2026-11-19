"""
Step 2 - Streamlines (psi = const) and equipotentials (phi = const)
around the cylinder, plus the speed field |f'(z)|.
Saves: figures/fig1_streamlines.png, figures/fig2_speed_field.png
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from flow_model import complex_potential, complex_velocity, make_grid, RADIUS

os.makedirs("figures", exist_ok=True)

X, Y, Z = make_grid()
inside = np.abs(Z) < RADIUS
F = complex_potential(Z)
phi, psi = F.real.copy(), F.imag.copy()
phi[inside] = np.nan
psi[inside] = np.nan

# Figure 1: streamlines and equipotentials
fig, ax = plt.subplots(figsize=(9, 6))
ax.contour(X, Y, psi, levels=np.linspace(-3, 3, 31), colors="tab:blue", linewidths=1)
ax.contour(X, Y, phi, levels=np.linspace(-6, 6, 31), colors="tab:orange",
           linestyles="dashed", linewidths=0.8)
ax.add_patch(plt.Circle((0, 0), RADIUS, color="0.6", zorder=3))
ax.plot([-RADIUS, RADIUS], [0, 0], "o", color="crimson", ms=8, zorder=4,
        label="Stagnation points z = ±a")
ax.plot([], [], color="tab:blue", label="Streamlines  ψ = const")
ax.plot([], [], color="tab:orange", ls="--", label="Equipotentials  φ = const")
ax.set_aspect("equal")
ax.set_xlabel("x (m)"); ax.set_ylabel("y (m)")
ax.set_title("Ideal flow past a cylinder:  f(z) = U(z + a²/z),  U = 1 m/s, a = 1 m")
ax.legend(loc="upper right", fontsize=9, framealpha=0.95)
fig.tight_layout()
fig.savefig("figures/fig1_streamlines.png", dpi=200)
print("saved figures/fig1_streamlines.png")

# Figure 2: speed field
speed = np.abs(complex_velocity(Z))
speed[inside] = np.nan
fig, ax = plt.subplots(figsize=(9, 6))
pc = ax.pcolormesh(X, Y, speed, cmap="viridis", shading="auto", vmin=0, vmax=2)
ax.contour(X, Y, psi, levels=np.linspace(-3, 3, 21), colors="white", linewidths=0.6)
ax.add_patch(plt.Circle((0, 0), RADIUS, color="0.6", zorder=3))
fig.colorbar(pc, ax=ax, label="Speed |f′(z)|  (m/s)")
ax.set_aspect("equal")
ax.set_xlabel("x (m)"); ax.set_ylabel("y (m)")
ax.set_title("Speed field: zero at the stagnation points, 2U at the top and bottom")
fig.tight_layout()
fig.savefig("figures/fig2_speed_field.png", dpi=200)
print("saved figures/fig2_speed_field.png")
print(f"max speed in field = {np.nanmax(speed):.4f} m/s (theory: 2U = 2.0000)")
