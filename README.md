# Ideal Airflow Around a Circular Cylinder Using Complex Potential

PBL project for **Real and Complex Analysis (26B11MA318)**, B.Tech Mathematics and Computing, JIIT Noida.

The flow is modelled with the complex potential

    f(z) = U (z + a²/z) + (iΓ / 2π) log z

and the code checks the theory (Cauchy–Riemann equations, harmonic functions), plots the flow, computes the surface pressure and forces, and maps the cylinder to an airfoil with the Joukowski transformation.

## Files

| File | What it does |
| --- | --- |
| `flow_model.py` | Shared model: complex potential, velocity, pressure, stagnation points |
| `s1_verify_cauchy_riemann.py` | SymPy proof: Cauchy–Riemann, Laplace's equation, circle is a streamline |
| `s2_streamlines.py` | Streamlines, equipotentials and speed field (Figures 1–2) |
| `s3_surface_pressure.py` | Surface speed and pressure coefficient Cp (Figure 3) |
| `s4_forces.py` | Drag and lift by numerical integration (d'Alembert, Kutta–Joukowski) |
| `s5_magnus_effect.py` | Spinning cylinder, moving stagnation points (Figure 4) |
| `s6_joukowski_airfoil.py` | Conformal map from cylinder to airfoil (Figure 5) |
| `run_all.py` | Runs all six steps in order |

## How to run in VS Code

1. Install **Python 3.10 or newer** from https://www.python.org/downloads/
   (on Windows, tick "Add Python to PATH" during installation).
2. Install **VS Code** and, inside it, the **Python** extension by Microsoft
   (Extensions panel: Ctrl+Shift+X, search "Python").
3. Unzip `cylinder_flow_pbl.zip`. In VS Code choose **File → Open Folder…** and open the `cylinder_flow_pbl` folder.
4. Open a terminal: **Terminal → New Terminal** (Ctrl+`).
5. (Recommended) create a virtual environment:
   - Windows: `python -m venv .venv` then `.venv\Scripts\activate`
   - macOS / Linux: `python3 -m venv .venv` then `source .venv/bin/activate`

   If VS Code asks "Select this environment for the workspace?", click **Yes**.
6. Install the libraries:

       pip install -r requirements.txt

7. Run everything:

       python run_all.py

   (use `python3` on macOS/Linux). Results print in the terminal and the
   figures are saved in the `figures` folder.
8. To run a single step, open its file (e.g. `s2_streamlines.py`) and press the
   **▶ Run Python File** button at the top right, or type `python s2_streamlines.py`.

### Common problems

- `python is not recognized` → reinstall Python with "Add to PATH" ticked, or use `py run_all.py` on Windows.
- `ModuleNotFoundError: No module named 'numpy'` → step 6 was not run in the same environment; run `pip install -r requirements.txt` again in the VS Code terminal.
- PowerShell blocks `activate` → run `Set-ExecutionPolicy -Scope CurrentUser RemoteSigned` once, then activate again.
- Figure windows do not open → that is expected; the scripts save PNG files to `figures/` instead.

## Changing the experiment

Edit the constants at the top of `flow_model.py` (`U_INF`, `RADIUS`, `RHO`), the list `cases` in
`s5_magnus_effect.py`, or `alpha` and `mu` in `s6_joukowski_airfoil.py`, then run again.
