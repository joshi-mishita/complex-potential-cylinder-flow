"""
Step 3 - Surface speed and pressure coefficient on |z| = a.
Compares the numerical value |f'(a e^{i theta})| with the exact result
q = 2U|sin theta| and Cp = 1 - 4 sin^2 theta.
Saves: figures/fig3_surface_pressure.png
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from flow_model import complex_velocity, pressure_coefficient, U_INF, RADIUS

os.makedirs("figures", exist_ok=True)

theta = np.linspace(0, 2 * np.pi, 721)
z_surface = RADIUS * np.exp(1j * theta) * (1 + 1e-12)
q_num = np.abs(complex_velocity(z_surface))
q_exact = 2 * U_INF * np.abs(np.sin(theta))
cp_num = pressure_coefficient(q_num)
cp_exact = 1 - 4 * np.sin(theta) ** 2

print(f"max |q_num - q_exact|   = {np.max(np.abs(q_num - q_exact)):.2e}")
print(f"max |Cp_num - Cp_exact| = {np.max(np.abs(cp_num - cp_exact)):.2e}")
print("\n theta(deg)   speed q/U     Cp")
for deg in (0, 30, 45, 60, 90, 120, 150, 180):
    t = np.radians(deg)
    q = 2 * abs(np.sin(t))
    print(f"   {deg:5d}     {q:7.4f}   {1 - q**2:8.4f}")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
deg = np.degrees(theta)
ax1.plot(deg, q_exact, color="tab:blue", lw=2, label="Exact 2U|sin θ|")
ax1.plot(deg[::30], q_num[::30], "o", color="black", ms=4, label="Numerical |f′(z)|")
ax1.set_xlabel("θ (degrees)"); ax1.set_ylabel("Surface speed q / U")
ax1.set_title("Speed on the cylinder surface")
ax1.set_xticks(range(0, 361, 45)); ax1.grid(alpha=0.3); ax1.legend()

ax2.plot(deg, cp_exact, color="tab:red", lw=2, label="Exact 1 − 4 sin²θ")
ax2.plot(deg[::30], cp_num[::30], "o", color="black", ms=4, label="Numerical")
ax2.axhline(0, color="0.5", lw=0.8)
ax2.set_xlabel("θ (degrees)"); ax2.set_ylabel("Pressure coefficient Cp")
ax2.set_title("Pressure on the cylinder surface")
ax2.set_xticks(range(0, 361, 45)); ax2.grid(alpha=0.3); ax2.legend()
fig.tight_layout()
fig.savefig("figures/fig3_surface_pressure.png", dpi=200)
print("\nsaved figures/fig3_surface_pressure.png")
