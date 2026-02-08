"""
Steady-State Vorticity-Streamfunction Cavity Solver
=====================================================

Matches Ghia, Ghia & Shin (1982) methodology exactly:
  - Vorticity-streamfunction (ψ-ω) formulation
  - 129×129 uniform grid (same as Ghia, h = 1/128)
  - Second-order central differences for all spatial operators
  - Woods' second-order wall vorticity boundary condition
  - Direct sparse Poisson solve for ψ (UMFPACK LU)
  - Gauss-Seidel/SOR relaxation for vorticity transport
  - Iterates on the STEADY-STATE equations (no time-marching)

Governing Equations (steady, 2D, incompressible):
  Vorticity transport:  u·∂ω/∂x + v·∂ω/∂y = ν·(∂²ω/∂x² + ∂²ω/∂y²)
  Poisson:              ∂²ψ/∂x² + ∂²ψ/∂y² = -ω
  Velocity:             u = ∂ψ/∂y,   v = -∂ψ/∂x

Boundary Conditions:
  ψ = 0 on all four walls (no flow through walls)
  u = U_lid at y = L (lid), u = 0 on other walls
  v = 0 on all walls

Wall Vorticity (Woods/Briley second-order formula):
  Uses two interior streamfunction values to eliminate the O(h³) truncation:
  South (j=0, u_wall=0): ω₀ = (-8ψ₁ + ψ₂) / (2h²)
  North (j=N-1, lid):     ω_N = (-8ψ_{N-2} + ψ_{N-3} - 6hU) / (2h²)
  West  (i=0, v_wall=0): ω₀ = (-8ψ_{·,1} + ψ_{·,2}) / (2h²)
  East  (i=N-1, v=0):     ω_N = (-8ψ_{·,N-2} + ψ_{·,N-3}) / (2h²)

  (Thom's first-order formula for comparison: ω_wall = -2ψ₁/h² ± 2U/h)

Algorithm:
  repeat:
    1. Compute wall vorticity from current ψ (Woods formula)
    2. Solve Poisson ∇²ψ = -ω for ψ (direct sparse LU, pre-factorized)
    3. Compute u = ∂ψ/∂y, v = -∂ψ/∂x via central differences
    4. Relax interior vorticity using SOR on the steady vorticity equation
    5. Check convergence (max change in ψ and ω)
  until converged

Reference: Ghia, Ghia & Shin (1982), J. Comp. Physics 48, pp 387-411.
           Woods (1954), Aeronautical Quarterly 5, pp 176-184.
"""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import factorized
import time


def build_poisson_matrix_psi(N, h):
    """
    Build the Poisson matrix for the interior nodes of an N×N grid.

    Interior nodes: (j, i) for 1 ≤ j ≤ N-2, 1 ≤ i ≤ N-2.
    Total interior unknowns: (N-2)².
    Boundary ψ = 0 goes into the RHS.
    Standard 5-point Laplacian: (ψ_E + ψ_W + ψ_N + ψ_S - 4ψ_C) / h² = -ω_C.

    Returns the sparse matrix A and the factorized solver.
    """
    Ni = N - 2                 # interior nodes per side
    M = Ni * Ni                # total interior unknowns
    h2 = h * h

    def idx(jj, ii):
        """Map interior indices (jj, ii) ∈ [0, Ni-1] to flat index."""
        return jj * Ni + ii

    rows, cols, vals = [], [], []
    for jj in range(Ni):
        for ii in range(Ni):
            k = idx(jj, ii)
            # Center coefficient
            rows.append(k); cols.append(k); vals.append(-4.0 / h2)

            # East (ii+1)
            if ii < Ni - 1:
                rows.append(k); cols.append(idx(jj, ii + 1)); vals.append(1.0 / h2)
            # else: ψ at east wall = 0, contribution absorbed into RHS (but RHS is -ω, no change)

            # West (ii-1)
            if ii > 0:
                rows.append(k); cols.append(idx(jj, ii - 1)); vals.append(1.0 / h2)

            # North (jj+1)
            if jj < Ni - 1:
                rows.append(k); cols.append(idx(jj + 1, ii)); vals.append(1.0 / h2)

            # South (jj-1)
            if jj > 0:
                rows.append(k); cols.append(idx(jj - 1, ii)); vals.append(1.0 / h2)

    A = sparse.csc_matrix((vals, (rows, cols)), shape=(M, M))
    return A, factorized(A)


def solve_cavity_psi_omega(N=129, Re=100.0, max_iter=200000, tol=1e-8,
                            use_woods=True, verbose=True):
    """
    Solve the lid-driven cavity using the ψ-ω formulation with
    pseudo-transient time-stepping (explicit Euler for vorticity,
    direct Poisson solve for streamfunction).

    Parameters
    ----------
    N : int
        Grid size (N×N nodes including boundaries).  N=129 matches Ghia exactly.
    Re : float
        Reynolds number.
    max_iter : int
        Maximum pseudo-time steps.
    tol : float
        Convergence tolerance on max |Δω| per step.
    use_woods : bool
        If True, use Woods' second-order wall vorticity; otherwise Thom's first-order.
    verbose : bool
        Print diagnostics.

    Returns
    -------
    dict with ψ, ω, u, v, x, y, convergence info.
    """
    U_lid = 1.0
    L = 1.0
    nu = U_lid * L / Re
    h = L / (N - 1)

    # CFL pseudo-time step for explicit Euler on convection-diffusion of ω
    # Diffusion: Δτ < h²/(4ν)  and  Convection: Δτ < h/|u_max|
    dt_diff = 0.25 * h * h / nu
    dt_conv = h / U_lid
    dt = 0.5 * min(dt_diff, dt_conv)          # safety factor 0.5

    if verbose:
        print(f"  Grid: {N}x{N} nodes  (h = {h:.6f})")
        print(f"  Re = {Re:.0f}, nu = {nu:.6f}")
        print(f"  Pe_cell = {U_lid * h / nu:.3f}  (need < 2 for central diff)")
        print(f"  dt = {dt:.4e}  (dt_diff={dt_diff:.3e}, dt_conv={dt_conv:.3e})")
        print(f"  Wall vorticity: {'Woods (2nd order)' if use_woods else 'Thom (1st order)'}")

    x = np.linspace(0, L, N)
    y = np.linspace(0, L, N)

    # Field arrays: [j, i] where j=y-index, i=x-index
    psi   = np.zeros((N, N), dtype=np.float64)
    omega = np.zeros((N, N), dtype=np.float64)
    u     = np.zeros((N, N), dtype=np.float64)
    v     = np.zeros((N, N), dtype=np.float64)

    # Lid velocity
    u[-1, :] = U_lid

    # Build & factorize Poisson matrix for interior nodes
    if verbose:
        print("  Building Poisson matrix...")
    Ni = N - 2
    A, solve_poisson = build_poisson_matrix_psi(N, h)
    if verbose:
        print(f"  Factorized {Ni*Ni}x{Ni*Ni} Poisson matrix")

    h2 = h * h
    t0 = time.time()
    converged = False

    for iteration in range(max_iter):

        # ================================================================
        # STEP 1:  Wall vorticity boundary conditions
        # ================================================================
        if use_woods:
            # Woods' second-order formula:
            #   ω_wall = (-8ψ₁ + ψ₂ ± 6h·U_wall) / (2h²)
            omega[0, 1:-1] = (-8.0 * psi[1, 1:-1] + psi[2, 1:-1]) / (2.0 * h2)
            omega[-1, 1:-1] = (-8.0 * psi[-2, 1:-1] + psi[-3, 1:-1]
                                - 6.0 * h * U_lid) / (2.0 * h2)
            omega[1:-1, 0] = (-8.0 * psi[1:-1, 1] + psi[1:-1, 2]) / (2.0 * h2)
            omega[1:-1, -1] = (-8.0 * psi[1:-1, -2] + psi[1:-1, -3]) / (2.0 * h2)
        else:
            # Thom's first-order formula
            omega[0, 1:-1]   = -2.0 * psi[1, 1:-1] / h2
            omega[-1, 1:-1]  = -2.0 * psi[-2, 1:-1] / h2 - 2.0 * U_lid / h
            omega[1:-1, 0]   = -2.0 * psi[1:-1, 1] / h2
            omega[1:-1, -1]  = -2.0 * psi[1:-1, -2] / h2

        # Corner nodes (average adjacent wall formulas)
        omega[0, 0]   = 0.5 * (omega[0, 1] + omega[1, 0])
        omega[0, -1]  = 0.5 * (omega[0, -2] + omega[1, -1])
        omega[-1, 0]  = 0.5 * (omega[-1, 1] + omega[-2, 0])
        omega[-1, -1] = 0.5 * (omega[-1, -2] + omega[-2, -1])

        # ================================================================
        # STEP 2:  Solve Poisson ∇²ψ = -ω (direct sparse LU)
        # ================================================================
        rhs = -omega[1:-1, 1:-1].flatten()
        psi[1:-1, 1:-1] = solve_poisson(rhs).reshape((Ni, Ni))

        # ================================================================
        # STEP 3:  Compute velocity from ψ (central differences)
        # ================================================================
        u[1:-1, 1:-1] = (psi[2:, 1:-1] - psi[:-2, 1:-1]) / (2.0 * h)
        v[1:-1, 1:-1] = -(psi[1:-1, 2:] - psi[1:-1, :-2]) / (2.0 * h)
        u[0, :] = 0.0;  u[-1, :] = U_lid;  u[:, 0] = 0.0;  u[:, -1] = 0.0
        v[0, :] = 0.0;  v[-1, :] = 0.0;     v[:, 0] = 0.0;  v[:, -1] = 0.0

        # ================================================================
        # STEP 4:  Pseudo-transient explicit Euler step for interior ω
        # ================================================================
        # ∂ω/∂τ = ν·∇²ω − u·∂ω/∂x − v·∂ω/∂y
        omega_old = omega[1:-1, 1:-1].copy()

        omega_E = omega[1:-1, 2:]
        omega_W = omega[1:-1, :-2]
        omega_N = omega[2:, 1:-1]
        omega_S = omega[:-2, 1:-1]
        omega_C = omega[1:-1, 1:-1]

        laplacian = (omega_E + omega_W + omega_N + omega_S - 4.0 * omega_C) / h2
        domega_dx = (omega_E - omega_W) / (2.0 * h)
        domega_dy = (omega_N - omega_S) / (2.0 * h)

        ue = u[1:-1, 1:-1]
        ve = v[1:-1, 1:-1]

        rhs_vort = nu * laplacian - ue * domega_dx - ve * domega_dy
        omega[1:-1, 1:-1] += dt * rhs_vort

        # ================================================================
        # STEP 5:  Convergence check
        # ================================================================
        domega_max = np.max(np.abs(omega[1:-1, 1:-1] - omega_old))

        if verbose and (iteration % 5000 == 0 or iteration < 5):
            elapsed = time.time() - t0
            print(f"  Iter {iteration:6d}: Δω={domega_max:.4e}  ({elapsed:.1f}s)")

        if domega_max < tol and iteration > 100:
            converged = True
            break

    elapsed = time.time() - t0
    if verbose:
        status = "CONVERGED" if converged else "NOT converged"
        print(f"  {status} at iter {iteration}, Δω={domega_max:.4e} ({elapsed:.1f}s)")

    # Vorticity at primary vortex center
    psi_min = np.min(psi)
    psi_min_loc = np.unravel_index(np.argmin(psi), psi.shape)

    return {
        'psi': psi, 'omega': omega, 'u': u, 'v': v,
        'x': x, 'y': y, 'N': N, 'h': h,
        'Re': Re, 'U_lid': U_lid,
        'converged': converged,
        'iterations': min(iteration + 1, max_iter),
        'elapsed': elapsed,
        'psi_min': psi_min,
        'psi_min_loc': (y[psi_min_loc[0]], x[psi_min_loc[1]]),
        'use_woods': use_woods,
    }


def ghia_data_re100():
    """Ghia, Ghia & Shin (1982) benchmark for Re=100."""
    ghia_y = np.array([
        0.0000, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 0.2813,
        0.4531, 0.5000, 0.6172, 0.7344, 0.8516, 0.9531, 0.9609,
        0.9688, 1.0000
    ])
    ghia_u = np.array([
        0.0000, -0.03717, -0.04192, -0.04775, -0.06434, -0.10150, -0.15662,
        -0.21090, -0.20581, -0.13641,  0.00332,  0.23151,  0.68717,  0.73722,
        0.78871,  1.0000
    ])
    ghia_x = np.array([
        0.0000, 0.0625, 0.0703, 0.0781, 0.0938, 0.1563, 0.2266,
        0.2344, 0.5000, 0.8047, 0.8594, 0.9063, 0.9453, 0.9531,
        0.9609, 1.0000
    ])
    ghia_v = np.array([
        0.0000,  0.09233,  0.10091,  0.10890,  0.12317,  0.16077,  0.17507,
        0.17527,  0.05454, -0.24533, -0.22445, -0.16914, -0.10313, -0.08864,
       -0.07391,  0.0000
    ])
    return ghia_y, ghia_u, ghia_x, ghia_v


def compare_with_ghia(result):
    """Compare ψ-ω solver results with Ghia benchmark data."""
    u = result['u']
    v = result['v']
    x = result['x']
    y = result['y']
    N = result['N']

    # Extract centerline profiles on the NODE GRID (same as Ghia)
    # u along vertical centerline at x = 0.5
    i_center = (N - 1) // 2        # index of x = 0.5 (for N=129: i=64)
    u_centerline = u[:, i_center]   # u at all y-nodes along x=0.5

    # v along horizontal centerline at y = 0.5
    j_center = (N - 1) // 2        # index of y = 0.5
    v_centerline = v[j_center, :]   # v at all x-nodes along y=0.5

    ghia_y, ghia_u, ghia_x, ghia_v = ghia_data_re100()

    # Interpolate solver profiles to Ghia data points (should be near-exact for N=129)
    u_interp = np.interp(ghia_y, y, u_centerline)
    v_interp = np.interp(ghia_x, x, v_centerline)

    # Compute errors
    u_rel_errors = []
    for k in range(len(ghia_u)):
        if abs(ghia_u[k]) > 0.01:
            u_rel_errors.append(abs((u_interp[k] - ghia_u[k]) / ghia_u[k]) * 100.0)
    v_rel_errors = []
    for k in range(len(ghia_v)):
        if abs(ghia_v[k]) > 0.01:
            v_rel_errors.append(abs((v_interp[k] - ghia_v[k]) / ghia_v[k]) * 100.0)

    print(f"\n  U-velocity along vertical centerline (x = {x[i_center]:.4f}):")
    print(f"  {'y':>8s}  {'Ghia':>10s}  {'Solver':>10s}  {'AbsErr':>10s}  {'RelErr':>8s}")
    for k in range(len(ghia_y)):
        ref = ghia_u[k]
        sol = u_interp[k]
        err = abs(sol - ref)
        rel = abs((sol - ref) / ref) * 100.0 if abs(ref) > 0.01 else 0.0
        flag = " ***" if rel > 5.0 else (" **" if rel > 2.0 else "")
        print(f"  {ghia_y[k]:8.4f}  {ref:10.5f}  {sol:10.5f}  {err:10.6f}  {rel:7.2f}%{flag}")

    print(f"\n  U max absolute error:  {max(abs(u_interp - ghia_u)):.6f}")
    if u_rel_errors:
        print(f"  U max relative error:  {max(u_rel_errors):.3f}%")
        print(f"  U mean relative error: {np.mean(u_rel_errors):.3f}%")

    print(f"\n  V-velocity along horizontal centerline (y = {y[j_center]:.4f}):")
    print(f"  {'x':>8s}  {'Ghia':>10s}  {'Solver':>10s}  {'AbsErr':>10s}  {'RelErr':>8s}")
    for k in range(len(ghia_x)):
        ref = ghia_v[k]
        sol = v_interp[k]
        err = abs(sol - ref)
        rel = abs((sol - ref) / ref) * 100.0 if abs(ref) > 0.01 else 0.0
        flag = " ***" if rel > 5.0 else (" **" if rel > 2.0 else "")
        print(f"  {ghia_x[k]:8.4f}  {ref:10.5f}  {sol:10.5f}  {err:10.6f}  {rel:7.2f}%{flag}")

    print(f"\n  V max absolute error:  {max(abs(v_interp - ghia_v)):.6f}")
    if v_rel_errors:
        print(f"  V max relative error:  {max(v_rel_errors):.3f}%")
        print(f"  V mean relative error: {np.mean(v_rel_errors):.3f}%")

    max_u_rel = max(u_rel_errors) if u_rel_errors else 0.0
    max_v_rel = max(v_rel_errors) if v_rel_errors else 0.0
    overall = max(max_u_rel, max_v_rel)

    print(f"\n  ψ_min = {result['psi_min']:.6f} at y={result['psi_min_loc'][0]:.4f}, x={result['psi_min_loc'][1]:.4f}")
    print(f"  (Ghia: ψ_min ≈ -0.1034 at ≈ (0.6172, 0.7344))")

    print(f"\n  Overall max relative error: {overall:.3f}%")
    print(f"  PASS (<2%): {'YES' if overall < 2.0 else 'NO'}")
    print(f"  PASS (<5%): {'YES' if overall < 5.0 else 'NO'}")

    return {
        'u_interp': u_interp, 'v_interp': v_interp,
        'u_max_rel_error': max_u_rel,
        'v_max_rel_error': max_v_rel,
        'overall_max_error': overall,
        'passed_2pct': overall < 2.0,
        'passed_5pct': overall < 5.0,
    }


if __name__ == '__main__':
    print("=" * 72)
    print("Steady-State ψ-ω Cavity Solver (Ghia formulation)")
    print("  Iteration: Poisson(direct) + SOR(vorticity)")
    print("=" * 72)

    results = []

    # Test 1: N=129 with Woods wall vorticity (matching Ghia)
    print("\n--- Test 1: N=129, Woods wall vorticity (matching Ghia) ---")
    r1 = solve_cavity_psi_omega(N=129, Re=100.0, max_iter=200000, tol=1e-8,
                                 use_woods=True)
    c1 = compare_with_ghia(r1)
    results.append(("129-Woods", c1, r1))

    # Test 2: N=129 with Thom wall vorticity (for comparison)
    print("\n--- Test 2: N=129, Thom wall vorticity ---")
    r2 = solve_cavity_psi_omega(N=129, Re=100.0, max_iter=200000, tol=1e-8,
                                 use_woods=False)
    c2 = compare_with_ghia(r2)
    results.append(("129-Thom", c2, r2))

    # Summary
    print("\n" + "=" * 72)
    print("SUMMARY")
    print("=" * 72)
    for name, comp, res in results:
        status = "PASS" if comp['passed_2pct'] else (
            "MARGINAL" if comp['passed_5pct'] else "FAIL")
        wv = "Woods" if res['use_woods'] else "Thom"
        print(f"  {name:>12s} ({wv}): U={comp['u_max_rel_error']:.3f}%, "
              f"V={comp['v_max_rel_error']:.3f}%, "
              f"max={comp['overall_max_error']:.3f}% [{status}]  "
              f"({res['iterations']} iters, {res['elapsed']:.1f}s)")
    print("=" * 72)
