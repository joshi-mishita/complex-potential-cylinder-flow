"""Run every step of the project in order and save all figures to ./figures."""

import os
import runpy

os.chdir(os.path.dirname(os.path.abspath(__file__)))

STEPS = [
    "s1_verify_cauchy_riemann.py",
    "s2_streamlines.py",
    "s3_surface_pressure.py",
    "s4_forces.py",
    "s5_magnus_effect.py",
    "s6_joukowski_airfoil.py",
]

for step in STEPS:
    print("\n" + "#" * 70 + f"\n# {step}\n" + "#" * 70)
    runpy.run_path(step, run_name="__main__")

print("\nDone. All figures are in the 'figures' folder.")
