"""
Step 4 - Forces on the cylinder by numerical integration of surface pressure.

    Drag D = -∫_0^{2π} p(θ) cos θ · a dθ
    Lift L = -∫_0^{2π} p(θ) sin θ · a dθ

The integrals are evaluated with Riemann sums (left, midpoint) and the
trapezoidal rule to show convergence. Expected:
    Gamma = 0  ->  D = 0, L = 0         (d'Alembert's paradox)
    Gamma != 0 ->  D = 0, L = rho U Gamma  (Kutta-Joukowski theorem)
"""

import numpy as np
from flow_model import surface_speed_exact, U_INF, RADIUS, RHO


def surface_pressure(theta, gamma):
    """Gauge pressure p - p_inf = 0.5 rho (U^2 - q^2)  (Bernoulli)."""
    q = surface_speed_exact(theta, gamma=gamma)
    return 0.5 * RHO * (U_INF**2 - q**2)


def integrate(fun, n, rule):
    h = 2 * np.pi / n
    k = np.arange(n)
    if rule == "left":
        return h * np.sum(fun(k * h))
    if rule == "midpoint":
        return h * np.sum(fun((k + 0.5) * h))
    if rule == "trapezoid":
        t = np.linspace(0, 2 * np.pi, n + 1)
        y = fun(t)
        return h * (np.sum(y) - 0.5 * (y[0] + y[-1]))
    raise ValueError(rule)


def forces(gamma, n=2000, rule="trapezoid"):
    drag = -integrate(lambda t: surface_pressure(t, gamma) * np.cos(t) * RADIUS, n, rule)
    lift = -integrate(lambda t: surface_pressure(t, gamma) * np.sin(t) * RADIUS, n, rule)
    return drag, lift


print("Convergence of the lift integral (Gamma = 4 m^2/s)")
print(f"exact rho*U*Gamma = {RHO * U_INF * 4:.6f} N/m")
print("     n        left        midpoint     trapezoid")
for n in (4, 8, 16, 64, 256):
    vals = [forces(4.0, n, r)[1] for r in ("left", "midpoint", "trapezoid")]
    print(f"{n:6d}   " + "   ".join(f"{v:10.6f}" for v in vals))

print("\nForces for several circulations (n = 2000, trapezoid rule)")
print(" Gamma (m^2/s)   Drag (N/m)     Lift (N/m)    rho*U*Gamma (N/m)")
for g in (0.0, 2.0, 4.0, 8.0, 4 * np.pi):
    d, l = forces(g)
    print(f"   {g:8.4f}    {d:11.2e}    {l:11.6f}     {RHO * U_INF * g:11.6f}")
