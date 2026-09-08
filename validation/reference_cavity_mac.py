"""
MAC Grid Chorin Projection Solver for Lid-Driven Cavity
=========================================================

Uses a staggered (MAC) grid to eliminate checkerboard pressure oscillations
that plague collocated-grid projection methods. On the MAC grid, the discrete
divergence and gradient operators are the EXACT adjoint (transpose) of each
other, which guarantees that the pressure projection produces a EXACTLY
divergence-free velocity field (to machine precision with a direct solve).

Grid Layout (ny pressure cells in y, nx in x):
  p[j, i]   at cell centers:      x=(i+0.5)*dx, y=(j+0.5)*dy   shape (ny, nx)
  u[j, i]   at vertical faces:    x=i*dx,       y=(j+0.5)*dy   shape (ny, nx+1)
  v[j, i]   at horizontal faces:  x=(i+0.5)*dx, y=j*dy          shape (ny+1, nx)

  Domain: [0, L] x [0, L],  dx = L/nx,  dy = L/ny

  j=0 → south wall (y=0)         j=ny → north boundary (y=L, lid)
  i=0 → west face/wall (x=0)     i=nx → east face/wall (x=L)

Boundary Conditions (lid-driven cavity):
  u[:, 0] = 0       (west wall)
  u[:, nx] = 0      (east wall)
  v[0, :] = 0       (south wall)
  v[ny, :] = 0      (north wall / lid, vertical component)
  u at y=0:  no-slip → u_wall = 0
  u at y=L:  lid     → u_wall = U_lid
  v at x=0:  no-slip → v_wall = 0
  v at x=L:  no-slip → v_wall = 0

Ghost Cells (QUADRATIC, for second-order accuracy at walls):
  Linear ghosts (u_g = 2*u_wall - u_interior) give first-order Laplacians near
  walls because the wall is at 0.5*h from the first interior point.  Quadratic
  extrapolation through u_wall at 0.5*h and two interior points restores the
  correct Laplacian coefficient (4/3 of the naive value → exact u'').

  South ghost for u:  u_g = (8/3)*0       - 2*u[0,:]  + (1/3)*u[1,:]  = -2*u[0,:] + u[1,:]/3
  North ghost for u:  u_g = (8/3)*U_lid   - 2*u[-1,:] + (1/3)*u[-2,:]
  West ghost for v:   v_g = (8/3)*0       - 2*v[:,0]  + (1/3)*v[:,1]  = -2*v[:,0] + v[:,1]/3
  East ghost for v:   v_g = (8/3)*0       - 2*v[:,-1] + (1/3)*v[:,-2] = -2*v[:,-1]+ v[:,-2]/3

  CFL Note: Quadratic ghost makes the near-wall center coefficient -4/h² (vs -2/h²
  for interior), so the diffusive CFL tightens from h²/(4ν) to h²/(6ν).

Algorithm (Chorin's Projection / Fractional Step):
  1. Predictor:  u* = u^n + dt * (nu * Lap(u^n) - (u^n . grad) u^n)
  2. Apply wall BCs to u*
  3. Pressure:   Lap(p) = (1/dt) * div(u*)   [sparse LU direct solve]
  4. Corrector:  u^{n+1} = u* - dt * grad(p) [staggered gradient]
  5. Wall face velocities remain unchanged

Key Property:
  div(u^{n+1}) = div(u*) - dt * div(grad(p))
               = div(u*) - dt * Lap(p)           [on MAC grid: div(grad) = Lap !]
               = div(u*) - div(u*)                [from the Poisson equation]
               = 0                                [EXACT, to machine precision]

  This consistency is IMPOSSIBLE on collocated grids with central differences,
  where div(grad(p)) != Lap_5pt(p), producing O(1) divergence errors.

Spatial:  2nd-order central differences (stable for Pe_cell < 2);
          for Re=100: Pe_cell = U*dx/nu = 100*dx, needs dx < 0.02 → nx >= 50
Temporal: Explicit Euler (steady-state target; final answer is dt-independent)
Pressure: scipy.sparse.linalg.factorized (UMFPACK-based direct LU)

Reference: Ghia, Ghia & Shin (1982), J. Comp. Physics 48, pp 387-411.
"""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import factorized
import time
import sys


# ─────────────────────────── Pressure Poisson Matrix ───────────────────────────

def build_poisson_matrix(nx, ny, dx, dy):
    """
    Build the pressure Poisson matrix for nx*ny cell-centered unknowns.

    Interior cells: standard 5-point Laplacian.
    Boundary cells: Neumann dP/dn=0 via ghost-cell elimination
        (ghost = interior neighbor → effective one-sided stencil).
    Reference pressure: p[0,0] = 0 (removes null-space of all-Neumann system).

    Returns a factorized (LU-decomposed) solver function.
    """
    N = nx * ny
    dx2 = dx * dx
    dy2 = dy * dy

    def idx(j, i):
        return j * nx + i

    rows, cols, vals = [], [], []

    for j in range(ny):
        for i in range(nx):
            k = idx(j, i)

            # Reference pressure: fix p[0,0] = 0
            if j == 0 and i == 0:
                rows.append(k); cols.append(k); vals.append(1.0)
                continue

            center = 0.0

            # --- x-direction ---
            if i > 0:                                       # west neighbor
                rows.append(k); cols.append(idx(j, i - 1)); vals.append(1.0 / dx2)
                center -= 1.0 / dx2
            else:
                # Neumann: ghost p = p[j,0] → (p_E - p_C)/dx2
                # No off-diagonal entry; center gets -1/dx2 from east side only
                pass  # center NOT modified; the missing neighbor is handled below

            if i < nx - 1:                                  # east neighbor
                rows.append(k); cols.append(idx(j, i + 1)); vals.append(1.0 / dx2)
                center -= 1.0 / dx2
            else:
                pass

            # Neumann boundaries contribute an additional -1/dx2 that is absorbed
            # because the ghost = center value. With ghost substituted:
            #   (p_neighbor + p_ghost - 2*p_center)/dx2 = (p_neighbor - p_center)/dx2
            # So each EXISTING off-diagonal gets +1/dx2 and center gets -1/dx2.
            # For a MISSING side (Neumann), the pair cancels so neither the
            # off-diagonal nor the center includes that direction — BUT the
            # standard Laplacian *always* has -2/dx2 from each direction.
            # Neumann removes one of the two neighbors, reducing center from
            # -2/dx2 to -1/dx2 for that direction.  This is exactly what the
            # if/else above achieves: each existing neighbor adds -1/dx2 to center.

            # --- y-direction ---
            if j > 0:                                       # south neighbor
                rows.append(k); cols.append(idx(j - 1, i)); vals.append(1.0 / dy2)
                center -= 1.0 / dy2

            if j < ny - 1:                                  # north neighbor
                rows.append(k); cols.append(idx(j + 1, i)); vals.append(1.0 / dy2)
                center -= 1.0 / dy2

            rows.append(k); cols.append(k); vals.append(center)

    A = sparse.csc_matrix((vals, (rows, cols)), shape=(N, N))
    return A


# ─────────────────────────── MAC Cavity Solver ─────────────────────────────────

def solve_cavity_mac(nx=64, ny=64, Re=100.0, max_iter=200000,
                     tol=1e-7, verbose=True):
    """
    Solve the lid-driven cavity flow on a MAC (staggered) grid using
    Chorin's projection method with explicit Euler time integration.

    Parameters
    ----------
    nx, ny : int
        Number of pressure cells in x and y.  Total grid points for
        comparison with Ghia: nx+1 vertical-face nodes, ny+1 horizontal-face nodes.
    Re : float
        Reynolds number based on lid velocity and cavity size.
    max_iter : int
        Maximum number of time steps.
    tol : float
        Convergence tolerance on max velocity change per step.
    verbose : bool
        Print iteration diagnostics.

    Returns
    -------
    dict with u, v, p arrays, grid coordinates, convergence info.
    """
    U_lid = 1.0
    L = 1.0
    rho = 1.0
    nu = U_lid * L / Re       # kinematic viscosity

    dx = L / nx
    dy = L / ny

    # Cell-Peclet stability check for central differences
    Pe_cell = U_lid * max(dx, dy) / nu
    if Pe_cell > 2.0 and verbose:
        print(f"  WARNING: Pe_cell = {Pe_cell:.2f} > 2; central differences may oscillate!")

    # CFL time step (explicit Euler)
    # With quadratic ghost cells, the near-wall Laplacian center coefficient
    # is -4/h² (vs -2/h² interior), tightening the diffusive CFL to h²/(6ν).
    h_min = min(dx, dy)
    dt_diff = h_min**2 / (6.0 * nu)            # 2D diffusive limit with quad ghost
    dt_conv = h_min / U_lid                     # convective limit
    safety = 0.9
    dt = safety * min(dt_diff, dt_conv)

    if verbose:
        print(f"  Grid: {nx}x{ny} cells  (dx={dx:.6f}, dy={dy:.6f})")
        print(f"  Re={Re:.0f}, nu={nu:.6f}, Pe_cell={Pe_cell:.3f}")
        print(f"  dt={dt:.6e}  (dt_diff={dt_diff:.3e}, dt_conv={dt_conv:.3e})")

    # ── Field arrays ──
    u = np.zeros((ny, nx + 1), dtype=np.float64)    # u at vertical faces
    v = np.zeros((ny + 1, nx), dtype=np.float64)     # v at horizontal faces
    p = np.zeros((ny, nx), dtype=np.float64)          # pressure at cell centers

    # Initial BC: lid velocity
    # u at y=L is not directly stored; it's imposed via ghost cells.
    # Wall faces are already zero.

    # ── Build & factorize pressure Poisson matrix ──
    if verbose:
        print("  Building pressure Poisson matrix...")
    A = build_poisson_matrix(nx, ny, dx, dy)
    solve_poisson = factorized(A)
    if verbose:
        print(f"  Factorized {nx*ny}x{nx*ny} Poisson matrix")

    t0 = time.time()
    converged = False

    for iteration in range(max_iter):
        u_old = u.copy()
        v_old = v.copy()

        # ================================================================
        # STEP 1:  Predictor — tentative velocity u*, v* (no pressure)
        # ================================================================

        # ── u-momentum at interior vertical faces: i = 1..nx-1, j = 0..ny-1 ──
        # Build ghost rows for u (shape ny+2 × nx+1), QUADRATIC extrapolation
        u_ext = np.empty((ny + 2, nx + 1), dtype=np.float64)
        u_ext[1:-1, :] = u
        # Quadratic ghost: u_g = (8/3)*u_wall - 2*u_1 + (1/3)*u_2
        u_ext[0, :] = -2.0 * u[0, :] + (1.0/3.0) * u[1, :]       # south (u_wall=0)
        u_ext[-1, :] = ((8.0/3.0) * U_lid
                        - 2.0 * u[-1, :]
                        + (1.0/3.0) * u[-2, :])                    # north (u_wall=U_lid)

        # Slices for interior u faces: rows 1..ny in u_ext → j=0..ny-1 of u
        # columns 1:-1 of u_ext → columns 1..nx-1 of u (skip wall faces)
        u_c = u_ext[1:-1, 1:-1]          # (ny, nx-1)  center
        u_e = u_ext[1:-1, 2:]            # (ny, nx-1)  east
        u_w = u_ext[1:-1, :-2]           # (ny, nx-1)  west
        u_n = u_ext[2:, 1:-1]            # (ny, nx-1)  north (uses ghost at j=ny)
        u_s = u_ext[:-2, 1:-1]           # (ny, nx-1)  south (uses ghost at j=-1)

        # Diffusion of u
        diff_u = nu * ((u_e - 2.0*u_c + u_w) / (dx*dx) +
                        (u_n - 2.0*u_c + u_s) / (dy*dy))

        # Convection of u:  u * du/dx + v_interp * du/dy  (non-conservative)
        du_dx = (u_e - u_w) / (2.0 * dx)
        du_dy = (u_n - u_s) / (2.0 * dy)

        # Interpolate v to u-face locations:
        # v_at_u[j, i'] = 0.25*(v[j,i-1] + v[j,i] + v[j+1,i-1] + v[j+1,i])
        # where i' = i-1 maps u interior index i (1..nx-1) to v index i-1 (0..nx-2)
        # and v index i (1..nx-1).
        # v[j, i-1] for j=0..ny-1, i=1..nx-1 :  v[:-1, :-1]   (ny, nx-1)
        # v[j, i]   for j=0..ny-1, i=1..nx-1 :  v[:-1, 1:]    (ny, nx-1)
        # v[j+1,i-1]                           :  v[1:, :-1]    (ny, nx-1)
        # v[j+1,i]                             :  v[1:, 1:]     (ny, nx-1)
        v_at_u = 0.25 * (v[:-1, :-1] + v[:-1, 1:] + v[1:, :-1] + v[1:, 1:])

        conv_u = u_c * du_dx + v_at_u * du_dy

        # Update tentative u at interior faces
        u[:, 1:-1] = u_c + dt * (diff_u - conv_u)

        # ── v-momentum at interior horizontal faces: j = 1..ny-1, i = 0..nx-1 ──
        # Build ghost columns for v (shape ny+1 × nx+2), QUADRATIC extrapolation
        v_ext = np.empty((ny + 1, nx + 2), dtype=np.float64)
        v_ext[:, 1:-1] = v
        # Quadratic ghost: v_g = (8/3)*v_wall - 2*v_1 + (1/3)*v_2
        v_ext[:, 0] = -2.0 * v[:, 0] + (1.0/3.0) * v[:, 1]       # west (v_wall=0)
        v_ext[:, -1] = -2.0 * v[:, -1] + (1.0/3.0) * v[:, -2]    # east (v_wall=0)

        # Slices for interior v faces: rows 1:-1 → j=1..ny-1, cols 1:-1 → i=0..nx-1
        v_c = v_ext[1:-1, 1:-1]          # (ny-1, nx)
        v_e = v_ext[1:-1, 2:]            # (ny-1, nx)
        v_w = v_ext[1:-1, :-2]           # (ny-1, nx)
        v_n = v_ext[2:, 1:-1]            # (ny-1, nx)
        v_s = v_ext[:-2, 1:-1]           # (ny-1, nx)

        diff_v = nu * ((v_e - 2.0*v_c + v_w) / (dx*dx) +
                        (v_n - 2.0*v_c + v_s) / (dy*dy))

        dv_dx = (v_e - v_w) / (2.0 * dx)
        dv_dy = (v_n - v_s) / (2.0 * dy)

        # Interpolate u to v-face locations:
        # u_at_v[j', i] = 0.25*(u[j-1,i] + u[j-1,i+1] + u[j,i] + u[j,i+1])
        # where j' = j-1 maps v interior index j (1..ny-1) to u index j-1 (0..ny-2)
        # and u index j (1..ny-1).
        # u[j-1, i]   for j=1..ny-1, i=0..nx-1 :  u[:-1, :-1]   (ny-1, nx)
        # u[j-1, i+1] for j=1..ny-1, i=0..nx-1 :  u[:-1, 1:]    (ny-1, nx)
        # u[j, i]                                :  u[1:, :-1]    (ny-1, nx)
        # u[j, i+1]                              :  u[1:, 1:]     (ny-1, nx)
        u_at_v = 0.25 * (u[:-1, :-1] + u[:-1, 1:] + u[1:, :-1] + u[1:, 1:])

        conv_v = u_at_v * dv_dx + v_c * dv_dy

        v[1:-1, :] = v_c + dt * (diff_v - conv_v)

        # ================================================================
        # STEP 2:  Apply wall BCs to tentative velocity
        # ================================================================
        u[:, 0] = 0.0           # west
        u[:, nx] = 0.0          # east
        v[0, :] = 0.0           # south
        v[ny, :] = 0.0          # north (lid: no vertical penetration)

        # ================================================================
        # STEP 3:  Pressure Poisson equation
        # ================================================================
        # div(u*) at cell centers using staggered differences:
        #   div = (u*[j,i+1] - u*[j,i]) / dx + (v*[j+1,i] - v*[j,i]) / dy
        div_ustar = ((u[:, 1:] - u[:, :-1]) / dx +
                     (v[1:, :] - v[:-1, :]) / dy)        # shape (ny, nx)

        rhs = div_ustar / dt                              # shape (ny, nx)
        rhs_flat = rhs.flatten()
        rhs_flat[0] = 0.0                                 # reference p[0,0] = 0

        p_flat = solve_poisson(rhs_flat)
        p = p_flat.reshape((ny, nx))

        # ================================================================
        # STEP 4:  Velocity correction (projection)
        # ================================================================
        # grad(p) at u-faces: dp/dx = (p[j,i] - p[j,i-1]) / dx  for i=1..nx-1
        u[:, 1:-1] -= dt * (p[:, 1:] - p[:, :-1]) / dx

        # grad(p) at v-faces: dp/dy = (p[j,i] - p[j-1,i]) / dy  for j=1..ny-1
        v[1:-1, :] -= dt * (p[1:, :] - p[:-1, :]) / dy

        # Wall BCs unchanged (u[:,0], u[:,nx], v[0,:], v[ny,:])

        # ================================================================
        # STEP 5:  Convergence check
        # ================================================================
        du_max = np.max(np.abs(u - u_old))
        dv_max = np.max(np.abs(v - v_old))
        vel_change = max(du_max, dv_max)

        if verbose and (iteration % 5000 == 0 or iteration < 5):
            # Divergence of corrected velocity (should be ~ eps_machine)
            div_corrected = ((u[:, 1:] - u[:, :-1]) / dx +
                             (v[1:, :] - v[:-1, :]) / dy)
            div_max = np.max(np.abs(div_corrected))
            elapsed = time.time() - t0
            print(f"  Iter {iteration:6d}: Δvel={vel_change:.4e}  "
                  f"|div|={div_max:.2e}  ({elapsed:.1f}s)")

        if vel_change < tol and iteration > 100:
            converged = True
            break

    elapsed = time.time() - t0
    if verbose:
        div_final = np.max(np.abs(
            (u[:, 1:] - u[:, :-1]) / dx + (v[1:, :] - v[:-1, :]) / dy))
        status = "CONVERGED" if converged else "NOT converged"
        print(f"  {status} at iter {iteration}, Δvel={vel_change:.4e}, "
              f"|div|={div_final:.2e} ({elapsed:.1f}s)")

    # ── Build coordinate arrays for output ──
    x_u = np.linspace(0, L, nx + 1)                    # u-face x-coords
    y_u = (np.arange(ny) + 0.5) * dy                    # u-face y-coords (cell centers)

    x_v = (np.arange(nx) + 0.5) * dx                    # v-face x-coords (cell centers)
    y_v = np.linspace(0, L, ny + 1)                      # v-face y-coords

    return {
        'u': u, 'v': v, 'p': p,
        'x_u': x_u, 'y_u': y_u,
        'x_v': x_v, 'y_v': y_v,
        'nx': nx, 'ny': ny, 'dx': dx, 'dy': dy,
        'Re': Re, 'U_lid': U_lid,
        'converged': converged,
        'iterations': min(iteration + 1, max_iter),
        'elapsed': elapsed,
        'dt': dt,
    }


# ─────────────────────────── Ghia Benchmark Data ──────────────────────────────

def ghia_data_re100():
    """
    Ghia, Ghia & Shin (1982) benchmark for Re=100 lid-driven cavity.
    Returns (y_positions, u_values) along vertical centerline x=0.5,
    and (x_positions, v_values) along horizontal centerline y=0.5.
    """
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


# ─────────────────────────── Comparison with Ghia ─────────────────────────────

def compare_with_ghia(result):
    """
    Extract velocity profiles along centerlines and compare with Ghia et al.

    On the MAC grid:
      u along x=0.5:  u[:, nx//2] at y-locations y_u = (j+0.5)*dy
                       with u=0 at y=0 and u=U_lid at y=1 added as endpoints
      v along y=0.5:  v[ny//2, :] at x-locations x_v = (i+0.5)*dx
                       with v=0 at x=0 and v=0 at x=1 added as endpoints
    """
    u = result['u']
    v = result['v']
    nx = result['nx']
    ny = result['ny']
    dy = result['dy']
    dx = result['dx']
    U_lid = result['U_lid']

    # ── u at vertical centerline x = 0.5 ──
    # u-face at x = (nx//2)*dx.  For even nx: exactly 0.5.
    i_center = nx // 2
    x_actual_u = i_center * dx
    u_centerline = u[:, i_center]                        # shape (ny,)
    y_profile_u = (np.arange(ny) + 0.5) * dy            # cell-center y

    # Add wall boundary values for complete profile
    y_full = np.concatenate([[0.0], y_profile_u, [1.0]])
    u_full = np.concatenate([[0.0], u_centerline, [U_lid]])

    # ── v at horizontal centerline y = 0.5 ──
    j_center = ny // 2
    y_actual_v = j_center * dy
    v_centerline = v[j_center, :]                        # shape (nx,)
    x_profile_v = (np.arange(nx) + 0.5) * dx            # cell-center x

    x_full = np.concatenate([[0.0], x_profile_v, [1.0]])
    v_full = np.concatenate([[0.0], v_centerline, [0.0]])

    # ── Get Ghia reference data ──
    ghia_y, ghia_u, ghia_x, ghia_v = ghia_data_re100()

    # ── Interpolate solver profiles to Ghia data points ──
    u_interp = np.interp(ghia_y, y_full, u_full)
    v_interp = np.interp(ghia_x, x_full, v_full)

    # ── Compute errors ──
    u_abs_err = np.abs(u_interp - ghia_u)
    v_abs_err = np.abs(v_interp - ghia_v)

    u_rel_errors = []
    for k in range(len(ghia_u)):
        if abs(ghia_u[k]) > 0.01:
            u_rel_errors.append(abs((u_interp[k] - ghia_u[k]) / ghia_u[k]) * 100.0)
    v_rel_errors = []
    for k in range(len(ghia_v)):
        if abs(ghia_v[k]) > 0.01:
            v_rel_errors.append(abs((v_interp[k] - ghia_v[k]) / ghia_v[k]) * 100.0)

    # ── Print detailed comparison ──
    print(f"\n  U-velocity along vertical centerline (x = {x_actual_u:.4f}):")
    print(f"  {'y':>8s}  {'Ghia':>10s}  {'Solver':>10s}  {'AbsErr':>10s}  {'RelErr':>8s}")
    for k in range(len(ghia_y)):
        ref = ghia_u[k]
        sol = u_interp[k]
        err = abs(sol - ref)
        rel = abs((sol - ref) / ref) * 100.0 if abs(ref) > 0.01 else 0.0
        flag = " ***" if rel > 5.0 else (" **" if rel > 2.0 else "")
        print(f"  {ghia_y[k]:8.4f}  {ref:10.5f}  {sol:10.5f}  {err:10.6f}  {rel:7.2f}%{flag}")

    print(f"\n  U max absolute error:  {np.max(u_abs_err):.6f}")
    if u_rel_errors:
        print(f"  U max relative error:  {max(u_rel_errors):.3f}%")
        print(f"  U mean relative error: {np.mean(u_rel_errors):.3f}%")

    print(f"\n  V-velocity along horizontal centerline (y = {y_actual_v:.4f}):")
    print(f"  {'x':>8s}  {'Ghia':>10s}  {'Solver':>10s}  {'AbsErr':>10s}  {'RelErr':>8s}")
    for k in range(len(ghia_x)):
        ref = ghia_v[k]
        sol = v_interp[k]
        err = abs(sol - ref)
        rel = abs((sol - ref) / ref) * 100.0 if abs(ref) > 0.01 else 0.0
        flag = " ***" if rel > 5.0 else (" **" if rel > 2.0 else "")
        print(f"  {ghia_x[k]:8.4f}  {ref:10.5f}  {sol:10.5f}  {err:10.6f}  {rel:7.2f}%{flag}")

    print(f"\n  V max absolute error:  {np.max(v_abs_err):.6f}")
    if v_rel_errors:
        print(f"  V max relative error:  {max(v_rel_errors):.3f}%")
        print(f"  V mean relative error: {np.mean(v_rel_errors):.3f}%")

    max_u_rel = max(u_rel_errors) if u_rel_errors else 0.0
    max_v_rel = max(v_rel_errors) if v_rel_errors else 0.0
    overall = max(max_u_rel, max_v_rel)

    print(f"\n  Overall max relative error: {overall:.3f}%")
    print(f"  PASS (<2%): {'YES' if overall < 2.0 else 'NO'}")
    print(f"  PASS (<5%): {'YES' if overall < 5.0 else 'NO'}")

    return {
        'u_interp': u_interp,
        'v_interp': v_interp,
        'u_max_rel_error': max_u_rel,
        'v_max_rel_error': max_v_rel,
        'overall_max_error': overall,
        'passed_2pct': overall < 2.0,
        'passed_5pct': overall < 5.0,
    }


# ─────────────────────────── Main ─────────────────────────────────────────────

if __name__ == '__main__':
    print("=" * 72)
    print("MAC Grid Chorin Projection — Lid-Driven Cavity Solver")
    print("  Staggered grid | Direct sparse LU | Central differences")
    print("=" * 72)

    results = []

    # ── Test 1: 64x64 cells (quick accuracy check) ──
    print("\n--- Test 1: 64x64 cells (Re = 100) ---")
    r64 = solve_cavity_mac(nx=64, ny=64, Re=100.0, max_iter=200000, tol=1e-7)
    c64 = compare_with_ghia(r64)
    results.append(("64x64", c64, r64))

    # ── Test 2: 128x128 cells (match Ghia grid resolution) ──
    print("\n--- Test 2: 128x128 cells (Re = 100) ---")
    r128 = solve_cavity_mac(nx=128, ny=128, Re=100.0, max_iter=400000, tol=1e-7)
    c128 = compare_with_ghia(r128)
    results.append(("128x128", c128, r128))

    # ── Summary ──
    print("\n" + "=" * 72)
    print("SUMMARY")
    print("=" * 72)
    for name, comp, res in results:
        status = "PASS" if comp['passed_2pct'] else (
            "MARGINAL" if comp['passed_5pct'] else "FAIL")
        print(f"  {name:>10s}: U={comp['u_max_rel_error']:.3f}%, "
              f"V={comp['v_max_rel_error']:.3f}%, "
              f"max={comp['overall_max_error']:.3f}% [{status}]  "
              f"({res['iterations']} iters, {res['elapsed']:.1f}s)")

    # Convergence order (if both grids available)
    if len(results) >= 2:
        err_coarse = results[0][1]['overall_max_error']
        err_fine   = results[1][1]['overall_max_error']
        if err_fine > 0 and err_coarse > 0:
            order = np.log2(err_coarse / err_fine)
            print(f"\n  Convergence order (64→128): {order:.2f}  (expected ~2.0)")

    print("=" * 72)
