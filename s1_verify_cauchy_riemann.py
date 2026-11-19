"""
Step 1 - Symbolic verification (SymPy)
  * phi and psi satisfy the Cauchy-Riemann equations
  * phi and psi are harmonic (Laplace's equation)
  * the circle |z| = a is a streamline (psi = 0 there)
  * f'(z) and the stagnation points
"""

import sympy as sp

x, y = sp.symbols("x y", real=True)
U, a = sp.symbols("U a", positive=True)
z = x + sp.I * y

f = U * (z + a**2 / z)
phi = sp.simplify(sp.re(sp.expand_complex(f)))
psi = sp.simplify(sp.im(sp.expand_complex(f)))

print("=" * 64)
print("Complex potential  f(z) = U (z + a^2/z)")
print("=" * 64)
print("phi(x, y) =", phi)
print("psi(x, y) =", psi)

# Cauchy-Riemann equations
cr1 = sp.simplify(sp.diff(phi, x) - sp.diff(psi, y))
cr2 = sp.simplify(sp.diff(phi, y) + sp.diff(psi, x))
print("\nCauchy-Riemann equations")
print("  phi_x - psi_y =", cr1)
print("  phi_y + psi_x =", cr2)

# Laplace's equation
lap_phi = sp.simplify(sp.diff(phi, x, 2) + sp.diff(phi, y, 2))
lap_psi = sp.simplify(sp.diff(psi, x, 2) + sp.diff(psi, y, 2))
print("\nHarmonicity")
print("  Laplacian(phi) =", lap_phi)
print("  Laplacian(psi) =", lap_psi)

# Orthogonality of the two families: grad phi . grad psi = 0
dot = sp.simplify(sp.diff(phi, x) * sp.diff(psi, x) + sp.diff(phi, y) * sp.diff(psi, y))
print("\nOrthogonality  grad(phi) . grad(psi) =", dot)

# Circle is a streamline: polar form
r, t = sp.symbols("r theta", positive=True)
psi_polar = sp.simplify(psi.subs({x: r * sp.cos(t), y: r * sp.sin(t)}))
print("\nPolar form  psi(r, theta) =", sp.factor(psi_polar))
print("  psi on r = a           =", sp.simplify(psi_polar.subs(r, a)))

# Complex velocity and stagnation points
Z = sp.symbols("z")
fz = U * (Z + a**2 / Z)
w = sp.diff(fz, Z)
print("\nComplex velocity  f'(z) =", sp.simplify(w))
print("Stagnation points (f'(z) = 0):", sp.solve(sp.Eq(w, 0), Z))
print("Far field  lim_{z->oo} f'(z) =", sp.limit(w, Z, sp.oo))

ok = all(e == 0 for e in (cr1, cr2, lap_phi, lap_psi, dot))
print("\nRESULT:", "all checks passed" if ok else "a check FAILED")
