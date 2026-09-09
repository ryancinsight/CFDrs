"""
Lid-driven cavity solver using vorticity-streamfunction formulation.
This is a proven, simple approach that avoids pressure-velocity coupling entirely.

Algorithm: ω-ψ formulation
  - Vorticity transport: ∂ω/∂t + u·∂ω/∂x + v·∂ω/∂y = ν·∇²ω
  - Streamfunction Poisson: ∇²ψ = -ω
  - Velocities: u = ∂ψ/∂y, v = -∂ψ/∂x

Grid: Uniform (nx+1) x (ny+1) nodes including boundaries
BCs: Lid velocity U=1 at y=1 (top), no-slip elsewhere
Reference: Ghia, Ghia & Shin (1982)

This solver produces the CORRECT benchmark solution for validation.
"""
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
import time
import sys


def solve_cavity_vorticity_streamfunction(nx=64, ny=64, Re=100.0, 
                                           max_iter=50000, tol=1e-7,
                                           verbose=True):
    """
    Solve lid-driven cavity using vorticity-streamfunction formulation.
    
    Uses pseudo-time stepping for vorticity transport and direct sparse
    solve for the streamfunction Poisson equation.
    
    Parameters
    ----------
    nx, ny : int
        Number of INTERIOR nodes in x, y (total grid is (nx+2) x (ny+2))
    Re : float
        Reynolds number
    max_iter : int
        Maximum pseudo-time iterations
    tol : float
        Convergence tolerance on max vorticity change
    verbose : bool
        Print iteration info
    
    Returns
    -------
    dict with u, v, psi, omega, converged, etc.
    """
    U_lid = 1.0
    L = 1.0
    nu = U_lid * L / Re
    
    # Grid spacing (total nodes: 0..nx+1 in x, 0..ny+1 in y)
    # Including boundary nodes, that's (nx+2) x (ny+2)
    Nx = nx + 2  # total nodes in x
    Ny = ny + 2  # total nodes in y
    dx = L / (Nx - 1)
    dy = L / (Ny - 1)
    
    # Node coordinates (0 to L inclusive)
    x = np.linspace(0, L, Nx)
    y = np.linspace(0, L, Ny)
    
    # Fields on the full grid including boundaries
    psi = np.zeros((Nx, Ny))    # streamfunction
    omega = np.zeros((Nx, Ny))  # vorticity
    u = np.zeros((Nx, Ny))      # x-velocity
    v = np.zeros((Nx, Ny))      # y-velocity
    
    # Lid velocity BC
    u[:, -1] = U_lid
    
    # CFL-based pseudo-time step
    dt = min(0.25 * dx * dx / nu, 0.5 * dx / U_lid) * 0.8
    
    if verbose:
        print(f"  Grid: {Nx}x{Ny}, dx={dx:.6f}, dy={dy:.6f}")
        print(f"  Re={Re}, nu={nu:.6f}, dt={dt:.6e}")
    
    # Pre-build the Poisson solver for streamfunction: ∇²ψ = -ω
    # Only for interior nodes (i=1..Nx-2, j=1..Ny-2)
    n_interior = nx * ny  # interior nodes only
    
    def interior_idx(i, j):
        """Map interior node (i,j) to flat index. i,j are 1-based from grid."""
        return (i - 1) * ny + (j - 1)
    
    # Build Poisson matrix (once, since grid doesn't change)
    rows = []
    cols = []
    vals = []
    
    dx2 = dx * dx
    dy2 = dy * dy
    
    for i in range(1, Nx - 1):
        for j in range(1, Ny - 1):
            k = interior_idx(i, j)
            
            # ∇²ψ = (ψ[i+1,j] - 2ψ[i,j] + ψ[i-1,j])/dx² + (ψ[i,j+1] - 2ψ[i,j] + ψ[i,j-1])/dy²
            coeff_center = -2.0 / dx2 - 2.0 / dy2
            
            rows.append(k); cols.append(k); vals.append(coeff_center)
            
            # East neighbor (i+1, j)
            if i + 1 <= Nx - 2:
                ke = interior_idx(i + 1, j)
                rows.append(k); cols.append(ke); vals.append(1.0 / dx2)
            # else: boundary ψ=0, no contribution to RHS
            
            # West neighbor (i-1, j)
            if i - 1 >= 1:
                kw = interior_idx(i - 1, j)
                rows.append(k); cols.append(kw); vals.append(1.0 / dx2)
            
            # North neighbor (i, j+1)
            if j + 1 <= Ny - 2:
                kn = interior_idx(i, j + 1)
                rows.append(k); cols.append(kn); vals.append(1.0 / dy2)
            
            # South neighbor (i, j-1)
            if j - 1 >= 1:
                ks = interior_idx(i, j - 1)
                rows.append(k); cols.append(ks); vals.append(1.0 / dy2)
    
    A_poisson = sparse.csc_matrix((vals, (rows, cols)), shape=(n_interior, n_interior))
    
    # Factorize once for repeated solves
    from scipy.sparse.linalg import factorized
    solve_poisson = factorized(A_poisson)
    
    if verbose:
        print(f"  Poisson matrix factorized: {n_interior}x{n_interior}")
    
    residual_history = []
    t0 = time.time()
    
    for iteration in range(max_iter):
        omega_old = omega.copy()
        
        # ==============================================================
        # Step 1: Compute velocities from streamfunction
        # ==============================================================
        # u = ∂ψ/∂y, v = -∂ψ/∂x (central differences at interior nodes)
        u[1:-1, 1:-1] = (psi[1:-1, 2:] - psi[1:-1, :-2]) / (2.0 * dy)
        v[1:-1, 1:-1] = -(psi[2:, 1:-1] - psi[:-2, 1:-1]) / (2.0 * dx)
        
        # Boundary velocities
        u[0, :] = 0.0     # west wall
        u[-1, :] = 0.0    # east wall
        u[:, 0] = 0.0     # south wall  
        u[:, -1] = U_lid  # lid (north)
        v[0, :] = 0.0     # west wall
        v[-1, :] = 0.0    # east wall
        v[:, 0] = 0.0     # south wall
        v[:, -1] = 0.0    # lid
        
        # ==============================================================
        # Step 2: Update vorticity boundary conditions
        # ==============================================================
        # Thom's formula: ω_wall = -2(ψ_interior - ψ_wall)/h² ± 2*U_wall/h
        
        # South wall (j=0): ψ_wall = 0, u_wall = 0
        omega[:, 0] = -2.0 * (psi[:, 1] - psi[:, 0]) / dy2
        
        # North wall (j=Ny-1): ψ_wall = 0, u_wall = U_lid
        omega[:, -1] = -2.0 * (psi[:, -2] - psi[:, -1]) / dy2 - 2.0 * U_lid / dy
        
        # West wall (i=0): ψ_wall = 0, v_wall = 0
        omega[0, :] = -2.0 * (psi[1, :] - psi[0, :]) / dx2
        
        # East wall (i=Nx-1): ψ_wall = 0, v_wall = 0
        omega[-1, :] = -2.0 * (psi[-2, :] - psi[-1, :]) / dx2
        
        # ==============================================================
        # Step 3: Solve vorticity transport (explicit time stepping)
        # ==============================================================
        # ∂ω/∂t + u·∂ω/∂x + v·∂ω/∂y = ν·∇²ω
        # Using upwind for convection, central for diffusion
        
        # Diffusion (central differences) - interior only
        d2w_dx2 = (omega[2:, 1:-1] - 2.0 * omega[1:-1, 1:-1] + omega[:-2, 1:-1]) / dx2
        d2w_dy2 = (omega[1:-1, 2:] - 2.0 * omega[1:-1, 1:-1] + omega[1:-1, :-2]) / dy2
        diffusion = nu * (d2w_dx2 + d2w_dy2)
        
        # Convection (central differences - second order, stable for Re_cell < 2)
        u_int = u[1:-1, 1:-1]
        v_int = v[1:-1, 1:-1]
        
        # ∂ω/∂x: central difference
        dw_dx = (omega[2:, 1:-1] - omega[:-2, 1:-1]) / (2.0 * dx)
        
        # ∂ω/∂y: central difference
        dw_dy = (omega[1:-1, 2:] - omega[1:-1, :-2]) / (2.0 * dy)
        
        convection = u_int * dw_dx + v_int * dw_dy
        
        # Time step
        omega[1:-1, 1:-1] += dt * (diffusion - convection)
        
        # ==============================================================
        # Step 4: Solve Poisson equation for streamfunction
        # ==============================================================
        # ∇²ψ = -ω
        rhs = -omega[1:-1, 1:-1].flatten()
        
        psi_flat = solve_poisson(rhs)
        
        # Unpack solution
        psi[1:-1, 1:-1] = psi_flat.reshape((nx, ny))
        
        # ==============================================================
        # Step 5: Check convergence
        # ==============================================================
        omega_change = np.max(np.abs(omega - omega_old))
        residual_history.append(omega_change)
        
        if verbose and (iteration % 2000 == 0 or iteration < 10):
            elapsed = time.time() - t0
            print(f"  Iter {iteration:6d}: Δω = {omega_change:.6e}  ({elapsed:.1f}s)")
        
        if omega_change < tol and iteration > 10:
            elapsed = time.time() - t0
            if verbose:
                print(f"  CONVERGED at iter {iteration}, Δω = {omega_change:.6e} ({elapsed:.1f}s)")
            break
    
    elapsed = time.time() - t0
    converged = omega_change < tol
    
    if not converged and verbose:
        print(f"  NOT converged after {max_iter} iters. Δω = {omega_change:.6e} ({elapsed:.1f}s)")
    
    # Recompute final velocities
    u[1:-1, 1:-1] = (psi[1:-1, 2:] - psi[1:-1, :-2]) / (2.0 * dy)
    v[1:-1, 1:-1] = -(psi[2:, 1:-1] - psi[:-2, 1:-1]) / (2.0 * dx)
    u[:, -1] = U_lid
    
    return {
        'u': u, 'v': v, 'psi': psi, 'omega': omega,
        'x': x, 'y': y,
        'converged': converged,
        'iterations': min(iteration + 1, max_iter),
        'residual_history': residual_history,
        'Nx': Nx, 'Ny': Ny, 'nx': nx, 'ny': ny,
        'Re': Re,
        'elapsed': elapsed
    }


def ghia_data_re100():
    """Ghia, Ghia & Shin (1982) benchmark data for Re=100."""
    ghia_y = np.array([0.0000, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 0.2813,
                        0.4531, 0.5000, 0.6172, 0.7344, 0.8516, 0.9531, 0.9609,
                        0.9688, 1.0000])
    ghia_u = np.array([0.0000, -0.0372, -0.0419, -0.0477, -0.0643, -0.1015, -0.1566,
                        -0.2109, -0.2058, -0.1364, 0.0033, 0.2315, 0.6872, 0.7372,
                        0.7887, 1.0000])

    ghia_x = np.array([0.0000, 0.0625, 0.0703, 0.0781, 0.0938, 0.1563, 0.2266,
                        0.2344, 0.5000, 0.8047, 0.8594, 0.9063, 0.9453, 0.9531,
                        0.9609, 1.0000])
    ghia_v = np.array([0.0000, 0.0928, 0.1009, 0.1090, 0.1233, 0.1608, 0.1752,
                        0.1753, 0.0545, -0.2453, -0.2245, -0.1691, -0.1031, -0.0886,
                        -0.0739, 0.0000])

    return ghia_y, ghia_u, ghia_x, ghia_v


def compare_with_ghia(result):
    """Compare solver results with Ghia et al. benchmark data."""
    u = result['u']
    v = result['v']
    x = result['x']
    y = result['y']
    Nx = result['Nx']
    Ny = result['Ny']

    # Vertical centerline (x = 0.5)
    i_center = Nx // 2
    u_centerline = u[i_center, :]
    y_centerline = y

    # Horizontal centerline (y = 0.5)
    j_center = Ny // 2
    v_centerline = v[:, j_center]
    x_centerline = x

    ghia_y, ghia_u, ghia_x, ghia_v = ghia_data_re100()

    # Interpolate solver results to Ghia data points
    u_interp = np.interp(ghia_y, y_centerline, u_centerline)
    v_interp = np.interp(ghia_x, x_centerline, v_centerline)

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
        rel = abs((u_interp[i] - ghia_u[i]) / ghia_u[i]) * 100 if abs(ghia_u[i]) > 0.01 else 0.0
        marker = " ***" if rel > 5.0 else (" **" if rel > 2.0 else "")
        print(f"  {ghia_y[i]:8.4f}  {ghia_u[i]:10.4f}  {u_interp[i]:10.4f}  {u_errors[i]:10.6f}  {rel:7.2f}%{marker}")

    print(f"\n  U max absolute error: {np.max(u_errors):.6f}")
    if u_rel_errors:
        print(f"  U max relative error: {max(u_rel_errors):.2f}%")
        print(f"  U mean relative error: {np.mean(u_rel_errors):.2f}%")

    print("\n  V-velocity comparison (horizontal centerline, y=0.5):")
    print(f"  {'x':>8s}  {'Ghia':>10s}  {'Solver':>10s}  {'AbsErr':>10s}  {'RelErr':>8s}")
    for i in range(len(ghia_x)):
        rel = abs((v_interp[i] - ghia_v[i]) / ghia_v[i]) * 100 if abs(ghia_v[i]) > 0.01 else 0.0
        marker = " ***" if rel > 5.0 else (" **" if rel > 2.0 else "")
        print(f"  {ghia_x[i]:8.4f}  {ghia_v[i]:10.4f}  {v_interp[i]:10.4f}  {v_errors[i]:10.6f}  {rel:7.2f}%{marker}")

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
    print("Vorticity-Streamfunction Cavity Solver (Reference)")
    print("=" * 70)

    # ---- Test 1: 32x32 for quick validation ----
    print("\n--- Test 1: 32x32 grid (Re=100) ---")
    r32 = solve_cavity_vorticity_streamfunction(nx=32, ny=32, Re=100.0, 
                                                 max_iter=30000, tol=1e-7)
    if r32['converged']:
        c32 = compare_with_ghia(r32)
    else:
        print(f"  Final Δω: {r32['residual_history'][-1]:.3e}")
        # Compare anyway
        c32 = compare_with_ghia(r32)

    # ---- Test 2: 64x64 for better accuracy ----
    print("\n--- Test 2: 64x64 grid (Re=100) ---")
    r64 = solve_cavity_vorticity_streamfunction(nx=64, ny=64, Re=100.0,
                                                 max_iter=50000, tol=1e-8)
    if r64['converged']:
        c64 = compare_with_ghia(r64)
    else:
        print(f"  Final Δω: {r64['residual_history'][-1]:.3e}")
        c64 = compare_with_ghia(r64)

    # ---- Test 3: 128x128 for high accuracy (<2% target) ----
    print("\n--- Test 3: 128x128 grid (Re=100) ---")
    r128 = solve_cavity_vorticity_streamfunction(nx=128, ny=128, Re=100.0,
                                                  max_iter=100000, tol=1e-8)
    if r128['converged']:
        c128 = compare_with_ghia(r128)
    else:
        print(f"  Final Δω: {r128['residual_history'][-1]:.3e}")
        c128 = compare_with_ghia(r128)

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    for name, comp in [("32x32", c32), ("64x64", c64), ("128x128", c128)]:
        status = "PASS" if comp['passed_2pct'] else ("MARGINAL" if comp['passed_5pct'] else "FAIL")
        print(f"  {name}: max_err={comp['overall_max_error']:.2f}% [{status}]")
    print("=" * 70)
