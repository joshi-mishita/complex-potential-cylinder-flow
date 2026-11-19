"""
Step 6 - Joukowski transformation: cylinder flow -> airfoil flow.

    zeta = z + c^2 / z     (conformal wherever z != ±c)

A circle centred at mu = (-0.08 + 0.08i) passing through z = c = 1 is
mapped to a cambered airfoil. The circulation is fixed by the Kutta
condition (smooth flow off the sharp trailing edge).
Saves: figures/fig5_joukowski_airfoil.png
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import flow_model  # noqa: F401  (applies shared plot settings)

os.makedirs("figures", exist_ok=True)

U, c = 1.0, 1.0
alpha = np.radians(8)                 # angle of attack
mu = -0.08 + 0.08j                    # centre of the circle
R = abs(c - mu)                       # circle passes through z = c
beta = np.angle(c - mu)
gamma = 4 * np.pi * U * R * np.sin(alpha - beta)   # Kutta condition


def potential(z):
    """Flow past circle centred at mu, free stream at angle alpha, circulation gamma."""
    s = (z - mu) * np.exp(-1j * alpha)
    return U * (s + R**2 / s) + 1j * gamma / (2 * np.pi) * np.log(s)


def joukowski(z):
    return z + c**2 / z


# polar grid outside the circle
r = np.linspace(R * 1.0001, 4.5, 300)
t = np.linspace(0, 2 * np.pi, 600)
RR, TT = np.meshgrid(r, t)
Zc = mu + RR * np.exp(1j * TT)
psi = potential(Zc).imag
Zeta = joukowski(Zc)

circle = mu + R * np.exp(1j * t)
airfoil = joukowski(circle)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.8))
for ax, W, body, title in [
    (ax1, Zc, circle, "z-plane: flow past a circle"),
    (ax2, Zeta, airfoil, "ζ-plane: same flow past a Joukowski airfoil"),
]:
    ax.contour(W.real, W.imag, psi, levels=50, colors="tab:blue", linewidths=0.7)
    ax.fill(body.real, body.imag, color="0.6", zorder=3)
    ax.set_aspect("equal"); ax.set_title(title)
ax1.set_xlim(-3, 3); ax1.set_ylim(-2.2, 2.2)
ax2.set_xlim(-3.5, 3.5); ax2.set_ylim(-2.2, 2.2)
fig.suptitle(f"Conformal map ζ = z + c²/z   (angle of attack {np.degrees(alpha):.0f}°, "
             f"Kutta circulation Γ = {gamma:.3f} m²/s)")
fig.tight_layout()
fig.savefig("figures/fig5_joukowski_airfoil.png", dpi=200)

chord = airfoil.real.max() - airfoil.real.min()
cl = 2 * gamma / (U * chord)          # lift coefficient per unit span, L = rho U Gamma
print(f"circle radius R = {R:.4f}, centre mu = {mu}")
print(f"Kutta circulation Gamma = {gamma:.4f} m^2/s")
print(f"airfoil chord = {chord:.4f} m, lift coefficient C_L = {cl:.4f}")
print("saved figures/fig5_joukowski_airfoil.png")
