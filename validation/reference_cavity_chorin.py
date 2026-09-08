"""
Lid-driven cavity solver using Chorin's Projection (fractional step) method.
This directly computes (u, v, p) without a streamfunction, avoiding wall 
vorticity BC issues that limit the ω-ψ approach to first-order at walls.

Algorithm:
  1. Predictor: u* = u^n + dt*(-u·∇u + ν∇²u)  (no pressure)
  2. Pressure Poisson: ∇²p = (ρ/dt)*∇·u*
  3. Corrector: u^{n+1} = u* - (dt/ρ)*∇p

Spatial: 2nd-order central differences (convection + diffusion)
Temporal: Explicit Euler
Pressure: scipy sparse direct solve (factorized LU) — solved EXACTLY each step
Grid: Collocated uniform nx × ny including boundaries

Reference: Ghia, Ghia & Shin (1982), Re=100

Array convention: arr[j, i] where j=y-direction (rows), i=x-direction (cols)
  j=0 is south wall, j=ny-1 is north wall (lid)
  i=0 is west wall, i=nx-1 is east wall
"""
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import factorized
import time
import sys


def build_pressure_poisson_matrix(nx, ny, dx, dy):
    """
    Build the pressure Poisson matrix for ALL grid nodes.
    
    Interior nodes: standard 5-point Laplacian
    Top wall (j=ny-1): Dirichlet p=0 (identity equation)
    Bottom wall (j=0): Neumann dp/dy=0 → p[0,i] = p[1,i]
    Left wall (i=0): Neumann dp/dx=0 → p[j,0] = p[j,1]
    Right wall (i=nx-1): Neumann dp/dx=0 → p[j,nx-1] = p[j,nx-2]
    """
    N = nx * ny
    dx2 = dx * dx
    dy2 = dy * dy
    
    def idx(j, i):
        return j * nx + i
    
    rows = []
    cols = []
    vals = []
    
    for j in range(ny):
        for i in range(nx):
            k = idx(j, i)
            
            if j == ny - 1:
                # Top wall: Dirichlet p = 0
                rows.append(k); cols.append(k); vals.append(1.0)
                
            elif j == 0:
                # Bottom wall: Neumann dp/dy = 0 → p[0,i] = p[1,i]
                rows.append(k); cols.append(k); vals.append(1.0)
                rows.append(k); cols.append(idx(1, i)); vals.append(-1.0)
                
            elif i == 0:
                # Left wall: Neumann dp/dx = 0 → p[j,0] = p[j,1]
                rows.append(k); cols.append(k); vals.append(1.0)
                rows.append(k); cols.append(idx(j, 1)); vals.append(-1.0)
                
            elif i == nx - 1:
                # Right wall: Neumann dp/dx = 0 → p[j,nx-1] = p[j,nx-2]
                rows.append(k); cols.append(k); vals.append(1.0)
                rows.append(k); cols.append(idx(j, nx - 2)); vals.append(-1.0)
                
            else:
                # Interior node: 5-point Laplacian
                center = -2.0 / dx2 - 2.0 / dy2
                rows.append(k); cols.append(k); vals.append(center)
                
                # East: (j, i+1)
                rows.append(k); cols.append(idx(j, i + 1)); vals.append(1.0 / dx2)
                # West: (j, i-1)
                rows.append(k); cols.append(idx(j, i - 1)); vals.append(1.0 / dx2)
                # North: (j+1, i)
                rows.append(k); cols.append(idx(j + 1, i)); vals.append(1.0 / dy2)
                # South: (j-1, i)
                rows.append(k); cols.append(idx(j - 1, i)); vals.append(1.0 / dy2)
    
    A = sparse.csc_matrix((vals, (rows, cols)), shape=(N, N))
    return A


def solve_cavity_chorin(nx=41, ny=41, Re=100.0, max_iter=100000,
                         tol=1e-7, verbose=True):
    """
    Solve lid-driven cavity using Chorin's projection method.
    
    Parameters
    ----------
    nx, ny : int
        Total grid nodes in x, y (including boundaries)
    Re : float
        Reynolds number
    max_iter : int
        Maximum time steps
    tol : float
        Convergence tolerance on L-infinity velocity change
    verbose : bool
        Print progress
    
    Returns
    -------
    dict with u, v, p, x, y, converged, etc.
    """
    U_lid = 1.0
    L = 1.0
    rho = 1.0
    nu = U_lid * L / Re
    
    dx = L / (nx - 1)
    dy = L / (ny - 1)
    
    x = np.linspace(0, L, nx)
    y = np.linspace(0, L, ny)
    
    # CFL-based time step for explicit scheme with central differences
    # Diffusion: dt < dx^2 / (4*nu) for 2D
    # Convection: dt < dx / |u_max|
    dt_diff = 0.25 * min(dx, dy)**2 / nu
    dt_conv = 0.5 * min(dx, dy) / U_lid
    dt = min(dt_diff, dt_conv) * 0.5  # safety factor 0.5
    
    if verbose:
        print(f"  Grid: {nx}x{ny}, dx={dx:.6f}, dy={dy:.6f}")
        print(f"  Re={Re}, nu={nu:.6f}, dt={dt:.6e}")
        print(f"  dt_diff={dt_diff:.6e}, dt_conv={dt_conv:.6e}")
    
    # Fields [j, i] where j=y, i=x
    u = np.zeros((ny, nx))
    v = np.zeros((ny, nx))
    p = np.zeros((ny, nx))
    
    # Lid velocity BC
    u[-1, :] = U_lid
    
    # Build and factorize pressure Poisson matrix (once)
    if verbose:
        print("  Building pressure Poisson matrix...")
    A_poisson = build_pressure_poisson_matrix(nx, ny, dx, dy)
    solve_poisson = factorized(A_poisson)
    if verbose:
        print(f"  Poisson matrix factorized: {nx*ny}x{nx*ny}")
    
    t0 = time.time()
    converged = False
    
    for iteration in range(max_iter):
        un = u.copy()
        vn = v.copy()
        
        # =============================================================
        # Step 1: Predictor — compute tentative velocity (no pressure)
        # =============================================================
        # u* = u^n + dt * (-u·∂u/∂x - v·∂u/∂y + ν∇²u)
        # v* = v^n + dt * (-u·∂v/∂x - v·∂v/∂y + ν∇²v)
        # Central differences for both convection and diffusion
        
        # Diffusion terms (central, 2nd order)
        d2u_dx2 = (un[1:-1, 2:] - 2.0*un[1:-1, 1:-1] + un[1:-1, 0:-2]) / (dx*dx)
        d2u_dy2 = (un[2:, 1:-1] - 2.0*un[1:-1, 1:-1] + un[0:-2, 1:-1]) / (dy*dy)
        
        d2v_dx2 = (vn[1:-1, 2:] - 2.0*vn[1:-1, 1:-1] + vn[1:-1, 0:-2]) / (dx*dx)
        d2v_dy2 = (vn[2:, 1:-1] - 2.0*vn[1:-1, 1:-1] + vn[0:-2, 1:-1]) / (dy*dy)
        
        # Convection terms (central, 2nd order)
        du_dx = (un[1:-1, 2:] - un[1:-1, 0:-2]) / (2.0 * dx)
        du_dy = (un[2:, 1:-1] - un[0:-2, 1:-1]) / (2.0 * dy)
        
        dv_dx = (vn[1:-1, 2:] - vn[1:-1, 0:-2]) / (2.0 * dx)
        dv_dy = (vn[2:, 1:-1] - vn[0:-2, 1:-1]) / (2.0 * dy)
        
        u_int = un[1:-1, 1:-1]
        v_int = vn[1:-1, 1:-1]
        
        # Tentative velocity (no pressure gradient)
        u[1:-1, 1:-1] = (un[1:-1, 1:-1] + 
                          dt * (-u_int * du_dx - v_int * du_dy + 
                                nu * (d2u_dx2 + d2u_dy2)))
        
        v[1:-1, 1:-1] = (vn[1:-1, 1:-1] + 
                          dt * (-u_int * dv_dx - v_int * dv_dy + 
                                nu * (d2v_dx2 + d2v_dy2)))
        
        # Apply BCs to tentative velocity
        u[0, :] = 0.0     # south
        u[-1, :] = U_lid   # north (lid)
        u[:, 0] = 0.0      # west
        u[:, -1] = 0.0     # east
        v[0, :] = 0.0      # south
        v[-1, :] = 0.0     # north
        v[:, 0] = 0.0      # west
        v[:, -1] = 0.0     # east
        
        # =============================================================
        # Step 2: Pressure Poisson equation
        # =============================================================
        # ∇²p = (ρ/dt) * (∂u*/∂x + ∂v*/∂y)
        
        # Build RHS
        rhs = np.zeros((ny, nx))
        
        # Divergence of tentative velocity at interior nodes
        rhs[1:-1, 1:-1] = (rho / dt) * (
            (u[1:-1, 2:] - u[1:-1, 0:-2]) / (2.0 * dx) +
            (v[2:, 1:-1] - v[0:-2, 1:-1]) / (2.0 * dy)
        )
        
        # BCs for RHS:
        # Top (Dirichlet p=0): rhs = 0
        rhs[-1, :] = 0.0
        # Bottom (Neumann): rhs = 0  (for the equation p[0,i] - p[1,i] = 0)
        rhs[0, :] = 0.0
        # Left (Neumann): rhs = 0
        rhs[:, 0] = 0.0
        # Right (Neumann): rhs = 0
        rhs[:, -1] = 0.0
        
        # Solve
        p_flat = solve_poisson(rhs.flatten())
        p = p_flat.reshape((ny, nx))
        
        # =============================================================
        # Step 3: Corrector — project to divergence-free velocity
        # =============================================================
        # u^{n+1} = u* - (dt/ρ) * ∂p/∂x
        # v^{n+1} = v* - (dt/ρ) * ∂p/∂y
        
        u[1:-1, 1:-1] -= (dt / rho) * (p[1:-1, 2:] - p[1:-1, 0:-2]) / (2.0 * dx)
        v[1:-1, 1:-1] -= (dt / rho) * (p[2:, 1:-1] - p[0:-2, 1:-1]) / (2.0 * dy)
        
        # Re-apply BCs
        u[0, :] = 0.0
        u[-1, :] = U_lid
        u[:, 0] = 0.0
        u[:, -1] = 0.0
        v[0, :] = 0.0
        v[-1, :] = 0.0
        v[:, 0] = 0.0
        v[:, -1] = 0.0
        
        # =============================================================
        # Step 4: Check convergence
        # =============================================================
        du_max = np.max(np.abs(u - un))
        dv_max = np.max(np.abs(v - vn))
        vel_change = max(du_max, dv_max)
        
        if verbose and (iteration % 2000 == 0 or iteration < 10):
            elapsed = time.time() - t0
            div_max = np.max(np.abs(
                (u[1:-1, 2:] - u[1:-1, 0:-2]) / (2*dx) +
                (v[2:, 1:-1] - v[0:-2, 1:-1]) / (2*dy)
            ))
            print(f"  Iter {iteration:6d}: Δvel = {vel_change:.6e}, div = {div_max:.3e}  ({elapsed:.1f}s)")
        
        if vel_change < tol and iteration > 100:
            elapsed = time.time() - t0
            if verbose:
                print(f"  CONVERGED at iter {iteration}, Δvel = {vel_change:.6e} ({elapsed:.1f}s)")
            converged = True
            break
    
    elapsed = time.time() - t0
    
    if not converged and verbose:
        print(f"  NOT converged after {max_iter} iters. Δvel = {vel_change:.6e} ({elapsed:.1f}s)")
    
    return {
        'u': u, 'v': v, 'p': p,
        'x': x, 'y': y,
        'nx': nx, 'ny': ny,
        'Re': Re,
        'converged': converged,
        'iterations': min(iteration + 1, max_iter),
        'elapsed': elapsed,
        'dt': dt
    }


def ghia_data_re100():
    """Ghia, Ghia & Shin (1982) benchmark data for Re=100."""
    # U-velocity along vertical centerline at x=0.5
    ghia_y = np.array([0.0000, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 0.2813,
                        0.4531, 0.5000, 0.6172, 0.7344, 0.8516, 0.9531, 0.9609,
                        0.9688, 1.0000])
    ghia_u = np.array([0.0000, -0.03717, -0.04192, -0.04775, -0.06434, -0.10150, -0.15662,
                        -0.21090, -0.20581, -0.13641, 0.00332, 0.23151, 0.68717, 0.73722,
                        0.78871, 1.0000])

    # V-velocity along horizontal centerline at y=0.5
    ghia_x = np.array([0.0000, 0.0625, 0.0703, 0.0781, 0.0938, 0.1563, 0.2266,
                        0.2344, 0.5000, 0.8047, 0.8594, 0.9063, 0.9453, 0.9531,
                        0.9609, 1.0000])
    ghia_v = np.array([0.0000, 0.09233, 0.10091, 0.10890, 0.12317, 0.16077, 0.17507,
                        0.17527, 0.05454, -0.24533, -0.22445, -0.16914, -0.10313, -0.08864,
                        -0.07391, 0.0000])

    return ghia_y, ghia_u, ghia_x, ghia_v


def compare_with_ghia(result):
    """Compare solver results with Ghia et al. benchmark data."""
    u = result['u']
    v = result['v']
    x = result['x']
    y = result['y']
    nx = result['nx']
    ny = result['ny']

    # Vertical centerline (x = 0.5): extract u[all_j, i_center]
    # Interpolate to exact x=0.5 if grid center doesn't fall exactly there
    i_left = np.searchsorted(x, 0.5) - 1
    i_right = i_left + 1
    if i_left >= 0 and i_right < nx:
        t = (0.5 - x[i_left]) / (x[i_right] - x[i_left])
        u_centerline = (1.0 - t) * u[:, i_left] + t * u[:, i_right]
    else:
        u_centerline = u[:, nx // 2]
    
    # Horizontal centerline (y = 0.5): extract v[j_center, all_i]
    j_left = np.searchsorted(y, 0.5) - 1
    j_right = j_left + 1
    if j_left >= 0 and j_right < ny:
        t = (0.5 - y[j_left]) / (y[j_right] - y[j_left])
        v_centerline = (1.0 - t) * v[j_left, :] + t * v[j_right, :]
    else:
        v_centerline = v[ny // 2, :]

    ghia_y, ghia_u, ghia_x, ghia_v = ghia_data_re100()

    # Interpolate solver results to Ghia data points
    u_interp = np.interp(ghia_y, y, u_centerline)
    v_interp = np.interp(ghia_x, x, v_centerline)

    # Compute errors
    u_errors = np.abs(u_interp - ghia_u)
    v_errors = np.abs(v_interp - ghia_v)

    u_rel_errors = []
    for i in range(len(ghia_u)):
        if abs(ghia_u[i]) > 0.01:
            u_rel_errors.append(abs((u_interp[i] - ghia_u[i]) / ghia_u[i]) * 100)

    v_rel_errors = []
    for i in range(len(ghia_v)):
        if abs(ghia_v[i]) > 0.01:
            v_rel_errors.append(abs((v_interp[i] - ghia_v[i]) / ghia_v[i]) * 100)

    print("\n  U-velocity comparison (vertical centerline, x=0.5):")
    print(f"  {'y':>8s}  {'Ghia':>10s}  {'Solver':>10s}  {'AbsErr':>10s}  {'RelErr':>8s}")
    for i in range(len(ghia_y)):
        ref = ghia_u[i]
        sol = u_interp[i]
        err = abs(sol - ref)
        rel = abs((sol - ref) / ref) * 100 if abs(ref) > 0.01 else 0.0
        marker = " ***" if rel > 5.0 else (" **" if rel > 2.0 else "")
        print(f"  {ghia_y[i]:8.4f}  {ref:10.5f}  {sol:10.5f}  {err:10.6f}  {rel:7.2f}%{marker}")

    print(f"\n  U max absolute error: {np.max(u_errors):.6f}")
    if u_rel_errors:
        print(f"  U max relative error: {max(u_rel_errors):.2f}%")
        print(f"  U mean relative error: {np.mean(u_rel_errors):.2f}%")

    print("\n  V-velocity comparison (horizontal centerline, y=0.5):")
    print(f"  {'x':>8s}  {'Ghia':>10s}  {'Solver':>10s}  {'AbsErr':>10s}  {'RelErr':>8s}")
    for i in range(len(ghia_x)):
        ref = ghia_v[i]
        sol = v_interp[i]
        err = abs(sol - ref)
        rel = abs((sol - ref) / ref) * 100 if abs(ref) > 0.01 else 0.0
        marker = " ***" if rel > 5.0 else (" **" if rel > 2.0 else "")
        print(f"  {ghia_x[i]:8.4f}  {ref:10.5f}  {sol:10.5f}  {err:10.6f}  {rel:7.2f}%{marker}")

    print(f"\n  V max absolute error: {np.max(v_errors):.6f}")
    if v_rel_errors:
        print(f"  V max relative error: {max(v_rel_errors):.2f}%")
        print(f"  V mean relative error: {np.mean(v_rel_errors):.2f}%")

    max_u_rel = max(u_rel_errors) if u_rel_errors else 0.0
    max_v_rel = max(v_rel_errors) if v_rel_errors else 0.0
    overall = max(max_u_rel, max_v_rel)

    print(f"\n  Overall max relative error: {overall:.2f}%")
    print(f"  PASS (<2%): {'YES' if overall < 2.0 else 'NO'}")
    print(f"  PASS (<5%): {'YES' if overall < 5.0 else 'NO'}")

    return {
        'u_interp': u_interp, 'v_interp': v_interp,
        'u_max_rel_error': max_u_rel,
        'v_max_rel_error': max_v_rel,
        'overall_max_error': overall,
        'passed_2pct': overall < 2.0,
        'passed_5pct': overall < 5.0
    }


if __name__ == '__main__':
    print("=" * 70)
    print("Chorin's Projection Cavity Solver (Reference)")
    print("  Method: Fractional step, central differences, scipy sparse LU")
    print("=" * 70)

    results = []
    
    # ---- Test 1: 41x41 for quick check ----
    print("\n--- Test 1: 41x41 grid (Re=100) ---")
    r41 = solve_cavity_chorin(nx=41, ny=41, Re=100.0, 
                               max_iter=50000, tol=1e-7)
    c41 = compare_with_ghia(r41)
    results.append(("41x41", c41))

    # ---- Test 2: 65x65 for better accuracy ----
    print("\n--- Test 2: 65x65 grid (Re=100) ---")
    r65 = solve_cavity_chorin(nx=65, ny=65, Re=100.0,
                               max_iter=80000, tol=1e-7)
    c65 = compare_with_ghia(r65)
    results.append(("65x65", c65))

    # ---- Test 3: 129x129 for high accuracy (match Ghia grid) ----
    print("\n--- Test 3: 129x129 grid (Re=100) ---")
    r129 = solve_cavity_chorin(nx=129, ny=129, Re=100.0,
                                max_iter=200000, tol=1e-7)
    c129 = compare_with_ghia(r129)
    results.append(("129x129", c129))

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    for name, comp in results:
        status = "PASS" if comp['passed_2pct'] else ("MARGINAL" if comp['passed_5pct'] else "FAIL")
        print(f"  {name}: U={comp['u_max_rel_error']:.2f}%, V={comp['v_max_rel_error']:.2f}%, max={comp['overall_max_error']:.2f}% [{status}]")
    print("=" * 70)
