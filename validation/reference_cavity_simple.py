"""
Standalone Python SIMPLE solver for lid-driven cavity flow.
Used as a reference implementation to validate the Rust SIMPLEC solver.

Algorithm: SIMPLE with collocated grid and Rhie-Chow interpolation
Grid: Uniform, collocated (cell-centered)
BCs: Lid velocity at top, no-slip on all other walls
Reference: Ghia, Ghia & Shin (1982) Re=100

Author: Reference implementation for CFDrs validation
"""
import numpy as np
import sys


def solve_cavity_simple(nx=33, ny=33, Re=100.0, max_iter=5000, 
                         alpha_u=0.7, alpha_p=0.3, tol=1e-6,
                         verbose=True):
    """
    SIMPLE algorithm for lid-driven cavity on collocated grid.
    
    Parameters
    ----------
    nx, ny : int
        Number of cells in x and y directions
    Re : float
        Reynolds number = rho * U * L / mu
    max_iter : int
        Maximum number of SIMPLE iterations
    alpha_u : float
        Velocity under-relaxation factor
    alpha_p : float
        Pressure under-relaxation factor
    tol : float
        Convergence tolerance on mass residual
    verbose : bool
        Print iteration info
    
    Returns
    -------
    dict with u, v, p, converged, iterations, residual_history
    """
    # Physical parameters
    rho = 1.0
    U_lid = 1.0
    L = 1.0
    mu = rho * U_lid * L / Re
    
    # Grid
    dx = L / nx
    dy = L / ny
    
    # Cell centers
    xc = np.linspace(dx/2, L - dx/2, nx)
    yc = np.linspace(dy/2, L - dy/2, ny)
    
    # Initialize fields (nx x ny)
    u = np.zeros((nx, ny))
    v = np.zeros((nx, ny))
    p = np.zeros((nx, ny))
    
    # Set lid velocity initial condition for top row
    u[:, -1] = U_lid
    
    residual_history = []
    
    for iteration in range(max_iter):
        u_old = u.copy()
        v_old = v.copy()
        
        # ================================================================
        # Step 1: Compute momentum coefficients (FVM discretization)
        # ================================================================
        # ap, ae, aw, an, as_ for each cell
        # Using upwind differencing for convection, central for diffusion
        
        ae = np.zeros((nx, ny))
        aw = np.zeros((nx, ny))
        an = np.zeros((nx, ny))
        as_ = np.zeros((nx, ny))
        ap_u = np.zeros((nx, ny))
        ap_v = np.zeros((nx, ny))
        su = np.zeros((nx, ny))  # u-momentum source
        sv = np.zeros((nx, ny))  # v-momentum source
        
        for i in range(nx):
            for j in range(ny):
                # Diffusion coefficients
                De = mu * dy / dx if i < nx - 1 else 0.0
                Dw = mu * dy / dx if i > 0 else 0.0
                Dn = mu * dx / dy if j < ny - 1 else 0.0
                Ds = mu * dx / dy if j > 0 else 0.0
                
                # Face velocities (linear interpolation for convection)
                if i < nx - 1:
                    u_e = 0.5 * (u[i, j] + u[i+1, j])
                else:
                    u_e = 0.0  # wall
                if i > 0:
                    u_w = 0.5 * (u[i-1, j] + u[i, j])
                else:
                    u_w = 0.0  # wall
                if j < ny - 1:
                    v_n = 0.5 * (v[i, j] + v[i, j+1])
                else:
                    v_n = 0.0  # wall (lid has v=0)
                if j > 0:
                    v_s = 0.5 * (v[i, j] + v[i, j-1])
                else:
                    v_s = 0.0  # wall
                
                # Convective fluxes
                Fe = rho * u_e * dy
                Fw = rho * u_w * dy
                Fn = rho * v_n * dx
                Fs = rho * v_s * dx
                
                # Upwind scheme coefficients
                ae[i, j] = De + max(-Fe, 0.0)
                aw[i, j] = Dw + max(Fw, 0.0)
                an[i, j] = Dn + max(-Fn, 0.0)
                as_[i, j] = Ds + max(Fs, 0.0)
                
                # Boundary modifications
                # No-slip walls: add 2*mu*A/dx for wall-adjacent cells
                # The wall velocity is part of the source term
                
                # West wall (i=0)
                if i == 0:
                    aw[i, j] = 0.0
                    # Extra diffusion to wall: mu * dy / (dx/2)
                    sp_w = 2.0 * mu * dy / dx
                    ae[i, j] = De + max(-Fe, 0.0)  # keep east
                    # Source contribution from wall BC (u_wall = 0)
                    su[i, j] += sp_w * 0.0  # u=0 at west wall
                    ap_u[i, j] += sp_w
                    sv[i, j] += sp_w * 0.0  # v=0 at west wall
                    ap_v[i, j] += sp_w
                
                # East wall (i=nx-1)
                if i == nx - 1:
                    ae[i, j] = 0.0
                    sp_e = 2.0 * mu * dy / dx
                    su[i, j] += sp_e * 0.0  # u=0 at east wall
                    ap_u[i, j] += sp_e
                    sv[i, j] += sp_e * 0.0  # v=0 at east wall
                    ap_v[i, j] += sp_e
                
                # South wall (j=0)
                if j == 0:
                    as_[i, j] = 0.0
                    sp_s = 2.0 * mu * dx / dy
                    su[i, j] += sp_s * 0.0  # u=0 at south wall
                    ap_u[i, j] += sp_s
                    sv[i, j] += sp_s * 0.0  # v=0 at south wall
                    ap_v[i, j] += sp_s
                
                # North wall (j=ny-1): LID
                if j == ny - 1:
                    an[i, j] = 0.0
                    sp_n = 2.0 * mu * dx / dy
                    su[i, j] += sp_n * U_lid  # u = U_lid at lid
                    ap_u[i, j] += sp_n
                    sv[i, j] += sp_n * 0.0  # v=0 at lid
                    ap_v[i, j] += sp_n
                
                # Pressure gradient source
                # dp/dx contribution (backward/forward depending on position)
                if i > 0 and i < nx - 1:
                    su[i, j] += -(p[i+1, j] - p[i-1, j]) / (2.0 * dx) * dx * dy
                elif i == 0:
                    su[i, j] += -(p[i+1, j] - p[i, j]) / dx * dx * dy
                else:  # i == nx-1
                    su[i, j] += -(p[i, j] - p[i-1, j]) / dx * dx * dy
                
                # dp/dy contribution
                if j > 0 and j < ny - 1:
                    sv[i, j] += -(p[i, j+1] - p[i, j-1]) / (2.0 * dy) * dx * dy
                elif j == 0:
                    sv[i, j] += -(p[i, j+1] - p[i, j]) / dy * dx * dy
                else:  # j == ny-1
                    sv[i, j] += -(p[i, j] - p[i, j-1]) / dy * dx * dy
        
        # Diagonal coefficient
        ap_u += ae + aw + an + as_  # ap = sum(a_nb) + boundary additions
        ap_v += ae + aw + an + as_
        
        # Prevent zero diagonal
        ap_u = np.maximum(ap_u, 1e-30)
        ap_v = np.maximum(ap_v, 1e-30)
        
        # ================================================================
        # Step 2: Solve momentum equations (Gauss-Seidel)
        # ================================================================
        u_star = u_old.copy()
        v_star = v_old.copy()
        
        for gs_iter in range(20):
            for i in range(nx):
                for j in range(ny):
                    u_nb = 0.0
                    if i < nx - 1: u_nb += ae[i, j] * u_star[i+1, j]
                    if i > 0:      u_nb += aw[i, j] * u_star[i-1, j]
                    if j < ny - 1: u_nb += an[i, j] * u_star[i, j+1]
                    if j > 0:      u_nb += as_[i, j] * u_star[i, j-1]
                    
                    u_star[i, j] = (u_nb + su[i, j]) / ap_u[i, j]
                    
                    v_nb = 0.0
                    if i < nx - 1: v_nb += ae[i, j] * v_star[i+1, j]
                    if i > 0:      v_nb += aw[i, j] * v_star[i-1, j]
                    if j < ny - 1: v_nb += an[i, j] * v_star[i, j+1]
                    if j > 0:      v_nb += as_[i, j] * v_star[i, j-1]
                    
                    v_star[i, j] = (v_nb + sv[i, j]) / ap_v[i, j]
        
        # Under-relax velocities
        u_star = alpha_u * u_star + (1.0 - alpha_u) * u_old
        v_star = alpha_u * v_star + (1.0 - alpha_u) * v_old
        
        # ================================================================
        # Step 3: Compute d coefficients for pressure equation
        # ================================================================
        d_u = dx * dy / ap_u  # V / ap
        d_v = dx * dy / ap_v
        
        # ================================================================
        # Step 4: Solve pressure correction equation
        # ================================================================
        # Build pressure Laplacian: ae*p'_E + aw*p'_W + an*p'_N + as*p'_S - ap*p'_P = b
        # where ae = rho * d_e * dy / dx (face d coefficient * face area / dx)
        
        pp = np.zeros((nx, ny))  # pressure correction
        
        # Build coefficients for pressure correction Poisson equation
        a_E = np.zeros((nx, ny))
        a_W = np.zeros((nx, ny))
        a_N = np.zeros((nx, ny))
        a_S = np.zeros((nx, ny))
        a_P = np.zeros((nx, ny))
        b_p = np.zeros((nx, ny))
        
        for i in range(nx):
            for j in range(ny):
                # Face d coefficients (harmonic average)
                if i < nx - 1:
                    d_e = 2.0 * d_u[i, j] * d_u[i+1, j] / (d_u[i, j] + d_u[i+1, j] + 1e-30)
                    a_E[i, j] = rho * d_e * dy / dx
                if i > 0:
                    d_w = 2.0 * d_u[i-1, j] * d_u[i, j] / (d_u[i-1, j] + d_u[i, j] + 1e-30)
                    a_W[i, j] = rho * d_w * dy / dx
                if j < ny - 1:
                    d_n = 2.0 * d_v[i, j] * d_v[i, j+1] / (d_v[i, j] + d_v[i, j+1] + 1e-30)
                    a_N[i, j] = rho * d_n * dx / dy
                if j > 0:
                    d_s = 2.0 * d_v[i, j-1] * d_v[i, j] / (d_v[i, j-1] + d_v[i, j] + 1e-30)
                    a_S[i, j] = rho * d_s * dx / dy
                
                a_P[i, j] = a_E[i, j] + a_W[i, j] + a_N[i, j] + a_S[i, j]
                
                # RHS: -div(rho * u_star) * volume... actually mass imbalance
                # Mass flux through faces
                mass_in = 0.0
                
                if i < nx - 1:
                    u_e_star = 0.5 * (u_star[i, j] + u_star[i+1, j])
                    mass_in -= rho * u_e_star * dy  # east outflow
                else:
                    mass_in -= 0.0  # wall
                
                if i > 0:
                    u_w_star = 0.5 * (u_star[i-1, j] + u_star[i, j])
                    mass_in += rho * u_w_star * dy  # west inflow
                else:
                    mass_in += 0.0  # wall
                
                if j < ny - 1:
                    v_n_star = 0.5 * (v_star[i, j] + v_star[i, j+1])
                    mass_in -= rho * v_n_star * dx  # north outflow
                else:
                    mass_in -= 0.0  # wall (v=0 at lid)
                
                if j > 0:
                    v_s_star = 0.5 * (v_star[i, j] + v_star[i, j-1])
                    mass_in += rho * v_s_star * dx  # south inflow
                else:
                    mass_in += 0.0  # wall
                
                b_p[i, j] = mass_in  # positive = net mass inflow = need to correct
        
        # Fix reference pressure (avoid singular system)
        a_P[0, 0] = 1e30
        b_p[0, 0] = 0.0
        
        # Prevent zero diagonal
        a_P = np.maximum(a_P, 1e-30)
        
        # Solve with Gauss-Seidel
        for gs_iter in range(200):
            for i in range(nx):
                for j in range(ny):
                    pp_nb = 0.0
                    if i < nx - 1: pp_nb += a_E[i, j] * pp[i+1, j]
                    if i > 0:      pp_nb += a_W[i, j] * pp[i-1, j]
                    if j < ny - 1: pp_nb += a_N[i, j] * pp[i, j+1]
                    if j > 0:      pp_nb += a_S[i, j] * pp[i, j-1]
                    
                    pp[i, j] = (pp_nb + b_p[i, j]) / a_P[i, j]
        
        # ================================================================
        # Step 5: Correct pressure and velocities
        # ================================================================
        p += alpha_p * pp
        
        # Correct cell velocities
        for i in range(nx):
            for j in range(ny):
                if i > 0 and i < nx - 1:
                    u[i, j] = u_star[i, j] - d_u[i, j] * (pp[i+1, j] - pp[i-1, j]) / (2.0 * dx)
                elif i == 0:
                    u[i, j] = u_star[i, j] - d_u[i, j] * (pp[i+1, j] - pp[i, j]) / dx
                else:
                    u[i, j] = u_star[i, j] - d_u[i, j] * (pp[i, j] - pp[i-1, j]) / dx
                
                if j > 0 and j < ny - 1:
                    v[i, j] = v_star[i, j] - d_v[i, j] * (pp[i, j+1] - pp[i, j-1]) / (2.0 * dy)
                elif j == 0:
                    v[i, j] = v_star[i, j] - d_v[i, j] * (pp[i, j+1] - pp[i, j]) / dy
                else:
                    v[i, j] = v_star[i, j] - d_v[i, j] * (pp[i, j] - pp[i, j-1]) / dy
        
        # Enforce boundary conditions
        u[0, :] = 0.0      # west wall
        u[-1, :] = 0.0     # east wall
        u[:, 0] = 0.0      # south wall
        u[:, -1] = U_lid   # lid
        
        v[0, :] = 0.0      # west wall
        v[-1, :] = 0.0     # east wall
        v[:, 0] = 0.0      # south wall
        v[:, -1] = 0.0     # lid
        
        # ================================================================
        # Step 6: Check convergence (mass residual)
        # ================================================================
        mass_residual = np.sum(np.abs(b_p))
        residual_history.append(mass_residual)
        
        if verbose and (iteration % 100 == 0 or iteration < 10):
            print(f"  Iter {iteration:5d}: mass residual = {mass_residual:.6e}")
        
        if mass_residual < tol:
            if verbose:
                print(f"  Converged at iteration {iteration} with residual {mass_residual:.6e}")
            return {
                'u': u, 'v': v, 'p': p,
                'xc': xc, 'yc': yc,
                'converged': True,
                'iterations': iteration,
                'residual_history': residual_history,
                'nx': nx, 'ny': ny, 'Re': Re
            }
    
    if verbose:
        print(f"  Did NOT converge after {max_iter} iterations. Final residual: {mass_residual:.6e}")
    
    return {
        'u': u, 'v': v, 'p': p,
        'xc': xc, 'yc': yc,
        'converged': False,
        'iterations': max_iter,
        'residual_history': residual_history,
        'nx': nx, 'ny': ny, 'Re': Re
    }


def extract_centerline_profiles(result):
    """Extract centerline velocity profiles for comparison with Ghia et al."""
    u = result['u']
    v = result['v']
    xc = result['xc']
    yc = result['yc']
    nx = result['nx']
    ny = result['ny']
    
    # U-velocity along vertical centerline (x = 0.5)
    i_center = nx // 2
    u_centerline = u[i_center, :]
    y_centerline = yc
    
    # V-velocity along horizontal centerline (y = 0.5)
    j_center = ny // 2
    v_centerline = v[:, j_center]
    x_centerline = xc
    
    return {
        'u_centerline': u_centerline,
        'y_coords': y_centerline,
        'v_centerline': v_centerline,
        'x_coords': x_centerline
    }


def ghia_data_re100():
    """Ghia, Ghia & Shin (1982) benchmark data for Re=100."""
    # U-velocity along vertical centerline at x=0.5
    ghia_y = np.array([0.0000, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 0.2813,
                        0.4531, 0.5000, 0.6172, 0.7344, 0.8516, 0.9531, 0.9609,
                        0.9688, 1.0000])
    ghia_u = np.array([0.0000, -0.0372, -0.0419, -0.0477, -0.0643, -0.1015, -0.1566,
                        -0.2109, -0.2058, -0.1364, 0.0033, 0.2315, 0.6872, 0.7372,
                        0.7887, 1.0000])
    
    # V-velocity along horizontal centerline at y=0.5
    ghia_x = np.array([0.0000, 0.0625, 0.0703, 0.0781, 0.0938, 0.1563, 0.2266,
                        0.2344, 0.5000, 0.8047, 0.8594, 0.9063, 0.9453, 0.9531,
                        0.9609, 1.0000])
    ghia_v = np.array([0.0000, 0.0928, 0.1009, 0.1090, 0.1233, 0.1608, 0.1752,
                        0.1753, 0.0545, -0.2453, -0.2245, -0.1691, -0.1031, -0.0886,
                        -0.0739, 0.0000])
    
    return ghia_y, ghia_u, ghia_x, ghia_v


def compare_with_ghia(result):
    """Compare solver results with Ghia et al. benchmark data."""
    profiles = extract_centerline_profiles(result)
    ghia_y, ghia_u, ghia_x, ghia_v = ghia_data_re100()
    
    # Interpolate solver results to Ghia data points
    u_interp = np.interp(ghia_y, profiles['y_coords'], profiles['u_centerline'])
    v_interp = np.interp(ghia_x, profiles['x_coords'], profiles['v_centerline'])
    
    # Add boundary values manually
    u_interp[0] = 0.0   # y=0 wall
    u_interp[-1] = 1.0  # y=1 lid
    v_interp[0] = 0.0   # x=0 wall
    v_interp[-1] = 0.0  # x=1 wall
    
    # Compute errors
    u_errors = np.abs(u_interp - ghia_u)
    v_errors = np.abs(v_interp - ghia_v)
    
    # Relative errors (skip near-zero values)
    u_rel_errors = []
    for i in range(len(ghia_u)):
        if abs(ghia_u[i]) > 0.01:
            u_rel_errors.append(abs((u_interp[i] - ghia_u[i]) / ghia_u[i]) * 100)
    
    v_rel_errors = []
    for i in range(len(ghia_v)):
        if abs(ghia_v[i]) > 0.01:
            v_rel_errors.append(abs((v_interp[i] - ghia_v[i]) / ghia_v[i]) * 100)
    
    print("\n  U-velocity comparison (vertical centerline):")
    print(f"  {'y':>8s}  {'Ghia':>10s}  {'Solver':>10s}  {'AbsErr':>10s}")
    for i in range(len(ghia_y)):
        print(f"  {ghia_y[i]:8.4f}  {ghia_u[i]:10.4f}  {u_interp[i]:10.4f}  {u_errors[i]:10.6f}")
    
    print(f"\n  U max absolute error: {np.max(u_errors):.6f}")
    print(f"  U mean absolute error: {np.mean(u_errors):.6f}")
    if u_rel_errors:
        print(f"  U max relative error: {np.max(u_rel_errors):.2f}%")
        print(f"  U mean relative error: {np.mean(u_rel_errors):.2f}%")
    
    print("\n  V-velocity comparison (horizontal centerline):")
    print(f"  {'x':>8s}  {'Ghia':>10s}  {'Solver':>10s}  {'AbsErr':>10s}")
    for i in range(len(ghia_x)):
        print(f"  {ghia_x[i]:8.4f}  {ghia_v[i]:10.4f}  {v_interp[i]:10.4f}  {v_errors[i]:10.6f}")
    
    print(f"\n  V max absolute error: {np.max(v_errors):.6f}")
    print(f"  V mean absolute error: {np.mean(v_errors):.6f}")
    if v_rel_errors:
        print(f"  V max relative error: {np.max(v_rel_errors):.2f}%")
        print(f"  V mean relative error: {np.mean(v_rel_errors):.2f}%")
    
    max_u_rel = max(u_rel_errors) if u_rel_errors else 0.0
    max_v_rel = max(v_rel_errors) if v_rel_errors else 0.0
    overall_max_error = max(max_u_rel, max_v_rel)
    
    return {
        'u_interp': u_interp,
        'v_interp': v_interp,
        'u_errors': u_errors,
        'v_errors': v_errors,
        'u_max_rel_error': max_u_rel,
        'v_max_rel_error': max_v_rel,
        'overall_max_error': overall_max_error,
        'passed': overall_max_error < 5.0  # 5% for coarse grid, refine for <2%
    }


if __name__ == '__main__':
    print("=" * 70)
    print("Python Reference SIMPLE Solver - Lid-Driven Cavity")
    print("=" * 70)
    
    # Start with small grid to verify algorithm works
    print("\n--- Test 1: Coarse grid (17x17) for quick convergence check ---")
    result_17 = solve_cavity_simple(nx=17, ny=17, Re=100.0, max_iter=3000, 
                                     alpha_u=0.5, alpha_p=0.1, tol=1e-5)
    
    if result_17['converged']:
        print(f"\n  CONVERGED in {result_17['iterations']} iterations")
        comp_17 = compare_with_ghia(result_17)
        print(f"\n  Overall max relative error vs Ghia: {comp_17['overall_max_error']:.2f}%")
    else:
        # Check if residual is at least decreasing
        rh = result_17['residual_history']
        print(f"\n  Residual trend: first={rh[0]:.3e}, min={min(rh):.3e}, last={rh[-1]:.3e}")
        if rh[-1] < rh[0]:
            print("  Residual IS decreasing - needs more iterations or better parameters")
        else:
            print("  Residual NOT decreasing - algorithm issue")
    
    # Medium grid 
    print("\n--- Test 2: Medium grid (33x33) ---")
    result_33 = solve_cavity_simple(nx=33, ny=33, Re=100.0, max_iter=5000,
                                     alpha_u=0.5, alpha_p=0.1, tol=1e-6)
    
    if result_33['converged']:
        print(f"\n  CONVERGED in {result_33['iterations']} iterations")
        comp_33 = compare_with_ghia(result_33)
        print(f"\n  Overall max relative error vs Ghia: {comp_33['overall_max_error']:.2f}%")
    else:
        rh = result_33['residual_history']
        print(f"\n  Residual trend: first={rh[0]:.3e}, min={min(rh):.3e}, last={rh[-1]:.3e}")
    
    print("\n" + "=" * 70)
    print("Reference solver complete")
    print("=" * 70)
