"""
flow_model.py
-------------
Core model for ideal (potential) flow around a circular cylinder.

Complex potential (with optional circulation Gamma):
    f(z) = U (z + a^2 / z) + (i Gamma / 2 pi) log z
    f = phi + i psi
        phi : velocity potential  (velocity = grad phi)
        psi : stream function     (streamlines are psi = constant)

Complex velocity:
    w(z) = f'(z) = U (1 - a^2 / z^2) + i Gamma / (2 pi z) = u - i v

Real and Complex Analysis (26B11MA318) - PBL, JIIT Noida.
"""

import numpy as np
import matplotlib

# draw negative-level contours solid (not dashed) in all figures
matplotlib.rcParams["contour.negative_linestyle"] = "solid"

# Default parameters (SI-style units)
U_INF = 1.0   # free-stream speed  (m/s)
RADIUS = 1.0  # cylinder radius a   (m)
RHO = 1.225   # air density         (kg/m^3)


def complex_potential(z, U=U_INF, a=RADIUS, gamma=0.0):
    """f(z) = U(z + a^2/z) + (i*gamma/2pi) log z."""
    z = np.asarray(z, dtype=complex)
    f = U * (z + a**2 / z)
    if gamma != 0.0:
        f = f + 1j * gamma / (2 * np.pi) * np.log(z)
    return f


def complex_velocity(z, U=U_INF, a=RADIUS, gamma=0.0):
    """w(z) = f'(z) = u - i v."""
    z = np.asarray(z, dtype=complex)
    return U * (1 - a**2 / z**2) + 1j * gamma / (2 * np.pi * z)


def velocity_components(z, U=U_INF, a=RADIUS, gamma=0.0):
    """Return (u, v) from w = u - i v."""
    w = complex_velocity(z, U, a, gamma)
    return w.real, -w.imag


def surface_speed_exact(theta, U=U_INF, a=RADIUS, gamma=0.0):
    """Tangential speed on |z| = a:  q_theta = -2U sin(theta) - gamma/(2 pi a)."""
    return -2 * U * np.sin(theta) - gamma / (2 * np.pi * a)


def pressure_coefficient(speed, U=U_INF):
    """Bernoulli: Cp = 1 - (q / U)^2."""
    return 1 - (speed / U) ** 2


def stagnation_points(U=U_INF, a=RADIUS, gamma=0.0):
    """Roots of f'(z) = 0, i.e. U z^2 + (i gamma / 2pi) z - U a^2 = 0."""
    return np.roots([U, 1j * gamma / (2 * np.pi), -U * a**2])


def make_grid(xlim=(-3, 3), ylim=(-2, 2), nx=600, ny=400):
    """Cartesian grid as a complex array z = x + i y."""
    x = np.linspace(*xlim, nx)
    y = np.linspace(*ylim, ny)
    X, Y = np.meshgrid(x, y)
    return X, Y, X + 1j * Y
