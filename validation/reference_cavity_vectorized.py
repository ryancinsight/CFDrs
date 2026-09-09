"""
Vectorized Python SIMPLE solver for lid-driven cavity flow.
Reference implementation to validate the Rust SIMPLEC solver.

Algorithm: SIMPLE with collocated grid (no Rhie-Chow needed for this formulation)
Uses scipy sparse solver for pressure equation.
Grid: Uniform, collocated (cell-centered)  
BCs: Lid velocity at top, no-slip on all other walls
Reference: Ghia, Ghia & Shin (1982) Re=100
"""
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
import time
import sys


def solve_cavity_simple(nx=33, ny=33, Re=100.0, max_iter=10000,
                         alpha_u=0.5, alpha_p=0.3, tol=1e-6,
                         verbose=True):
    """
    Vectorized SIMPLE algorithm for lid-driven cavity on collocated grid.
    """
    rho = 1.0
    U_lid = 1.0
    L = 1.0
    mu = rho * U_lid * L / Re

    dx = L / nx
    dy = L / ny
    vol = dx * dy

    xc = np.linspace(dx / 2, L - dx / 2, nx)
    yc = np.linspace(dy / 2, L - dy / 2, ny)

    # Fields: u[i, j], v[i, j], p[i, j]
    u = np.zeros((nx, ny))
    v = np.zeros((nx, ny))
    p = np.zeros((nx, ny))

    # Lid BC
    u[:, -1] = U_lid

    residual_history = []
    t0 = time.time()

    for iteration in range(max_iter):
        u_old = u.copy()
        v_old = v.copy()

        # ==============================================================
        # Step 1: Compute FVM coefficients (vectorized)
        # ==============================================================
        # Diffusion coefficients
        D_ew = mu * dy / dx  # east/west diffusion
        D_ns = mu * dx / dy  # north/south diffusion

        # Face velocities for convection (linear interpolation)
        # East faces (i, i+1): u_e[i,j] for i=0..nx-2
        u_e = np.zeros((nx, ny))
        u_e[:-1, :] = 0.5 * (u[:-1, :] + u[1:, :])
        # West faces: u_w[i,j] for i=1..nx-1
        u_w = np.zeros((nx, ny))
        u_w[1:, :] = 0.5 * (u[:-1, :] + u[1:, :])
        # North faces
        v_n = np.zeros((nx, ny))
        v_n[:, :-1] = 0.5 * (v[:, :-1] + v[:, 1:])
        # South faces
        v_s = np.zeros((nx, ny))
        v_s[:, 1:] = 0.5 * (v[:, :-1] + v[:, 1:])

        # Convective fluxes
        Fe = rho * u_e * dy
        Fw = rho * u_w * dy
        Fn = rho * v_n * dx
        Fs = rho * v_s * dx

        # Upwind coefficients
        ae = np.zeros((nx, ny))
        aw = np.zeros((nx, ny))
        an = np.zeros((nx, ny))
        as_ = np.zeros((nx, ny))

        # Interior east/west
        ae[:-1, :] = D_ew + np.maximum(-Fe[:-1, :], 0.0)
        aw[1:, :] = D_ew + np.maximum(Fw[1:, :], 0.0)
        # Interior north/south
        an[:, :-1] = D_ns + np.maximum(-Fn[:, :-1], 0.0)
        as_[:, 1:] = D_ns + np.maximum(Fs[:, 1:], 0.0)

        # Wall boundary modifications
        # West wall (i=0): no neighbor to the west, extra wall diffusion
        sp_w = 2.0 * mu * dy / dx
        aw[0, :] = 0.0

        # East wall (i=nx-1)
        sp_e = 2.0 * mu * dy / dx
        ae[-1, :] = 0.0

        # South wall (j=0)
        sp_s = 2.0 * mu * dx / dy
        as_[:, 0] = 0.0

        # North wall (j=ny-1) = LID
        sp_n = 2.0 * mu * dx / dy
        an[:, -1] = 0.0

        # Source terms for boundaries
        su = np.zeros((nx, ny))
        sv = np.zeros((nx, ny))
        sp = np.zeros((nx, ny))  # diagonal addition from boundaries

        # Wall boundary diagonal additions
        sp[0, :] += sp_w   # west
        sp[-1, :] += sp_e  # east
        sp[:, 0] += sp_s   # south
        sp[:, -1] += sp_n  # north

        # Lid drives u momentum
        su[:, -1] += sp_n * U_lid  # u = U_lid at lid

        # Pressure gradient source terms
        # dp/dx for u-momentum (centered difference where possible)
        # Interior: -(p[i+1,j] - p[i-1,j]) / (2*dx) * vol
        su[1:-1, :] += -(p[2:, :] - p[:-2, :]) / (2.0 * dx) * vol
        # Left boundary: -(p[1,j] - p[0,j]) / dx * vol
        su[0, :] += -(p[1, :] - p[0, :]) / dx * vol
        # Right boundary: -(p[-1,j] - p[-2,j]) / dx * vol
        su[-1, :] += -(p[-1, :] - p[-2, :]) / dx * vol

        # dp/dy for v-momentum
        sv[:, 1:-1] += -(p[:, 2:] - p[:, :-2]) / (2.0 * dy) * vol
        sv[:, 0] += -(p[:, 1] - p[:, 0]) / dy * vol
        sv[:, -1] += -(p[:, -1] - p[:, -2]) / dy * vol

        # Diagonal coefficient
        ap = ae + aw + an + as_ + sp
        ap = np.maximum(ap, 1e-30)

        # ==============================================================
        # Step 2: Solve momentum with Jacobi iterations (vectorized)
        # ==============================================================
        u_star = u_old.copy()
        v_star = v_old.copy()

        for _ in range(30):
            u_nb = np.zeros((nx, ny))
            # East neighbor
            u_nb[:-1, :] += ae[:-1, :] * u_star[1:, :]
            # West neighbor
            u_nb[1:, :] += aw[1:, :] * u_star[:-1, :]
            # North neighbor
            u_nb[:, :-1] += an[:, :-1] * u_star[:, 1:]
            # South neighbor
            u_nb[:, 1:] += as_[:, 1:] * u_star[:, :-1]

            u_star = (u_nb + su) / ap

            v_nb = np.zeros((nx, ny))
            v_nb[:-1, :] += ae[:-1, :] * v_star[1:, :]
            v_nb[1:, :] += aw[1:, :] * v_star[:-1, :]
            v_nb[:, :-1] += an[:, :-1] * v_star[:, 1:]
            v_nb[:, 1:] += as_[:, 1:] * v_star[:, :-1]

            v_star = (v_nb + sv) / ap

        # Under-relax
        u_star = alpha_u * u_star + (1.0 - alpha_u) * u_old
        v_star = alpha_u * v_star + (1.0 - alpha_u) * v_old

        # Enforce BCs on u_star, v_star
        u_star[0, :] = 0.0
        u_star[-1, :] = 0.0
        u_star[:, 0] = 0.0
        u_star[:, -1] = U_lid
        v_star[0, :] = 0.0
        v_star[-1, :] = 0.0
        v_star[:, 0] = 0.0
        v_star[:, -1] = 0.0

        # ==============================================================
        # Step 3: d coefficients
        # ==============================================================
        d = vol / ap  # V / ap for each cell

        # Face d coefficients (harmonic average)
        d_e = np.zeros((nx, ny))
        d_w = np.zeros((nx, ny))
        d_n = np.zeros((nx, ny))
        d_s = np.zeros((nx, ny))

        d_e[:-1, :] = 2.0 * d[:-1, :] * d[1:, :] / (d[:-1, :] + d[1:, :] + 1e-30)
        d_w[1:, :] = d_e[:-1, :]
        d_n[:, :-1] = 2.0 * d[:, :-1] * d[:, 1:] / (d[:, :-1] + d[:, 1:] + 1e-30)
        d_s[:, 1:] = d_n[:, :-1]

        # ==============================================================
        # Step 4: Pressure correction equation (sparse direct solve)
        # ==============================================================
        N = nx * ny

        def idx(i, j):
            return i * ny + j

        # Build sparse pressure Laplacian
        rows = []
        cols = []
        vals = []
        rhs = np.zeros(N)

        for i in range(nx):
            for j in range(ny):
                k = idx(i, j)

                a_e_p = 0.0
                a_w_p = 0.0
                a_n_p = 0.0
                a_s_p = 0.0

                if i < nx - 1:
                    a_e_p = rho * d_e[i, j] * dy / dx
                    rows.append(k)
                    cols.append(idx(i + 1, j))
                    vals.append(a_e_p)
                if i > 0:
                    a_w_p = rho * d_w[i, j] * dy / dx
                    rows.append(k)
                    cols.append(idx(i - 1, j))
                    vals.append(a_w_p)
                if j < ny - 1:
                    a_n_p = rho * d_n[i, j] * dx / dy
                    rows.append(k)
                    cols.append(idx(i, j + 1))
                    vals.append(a_n_p)
                if j > 0:
                    a_s_p = rho * d_s[i, j] * dx / dy
                    rows.append(k)
                    cols.append(idx(i, j - 1))
                    vals.append(a_s_p)

                a_p_p = a_e_p + a_w_p + a_n_p + a_s_p
                rows.append(k)
                cols.append(k)
                vals.append(-a_p_p)

                # RHS: mass imbalance
                # Compute face mass fluxes from u_star, v_star
                mass_out = 0.0

                if i < nx - 1:
                    u_e_f = 0.5 * (u_star[i, j] + u_star[i + 1, j])
                    mass_out += rho * u_e_f * dy
                if i > 0:
                    u_w_f = 0.5 * (u_star[i - 1, j] + u_star[i, j])
                    mass_out -= rho * u_w_f * dy
                if j < ny - 1:
                    v_n_f = 0.5 * (v_star[i, j] + v_star[i, j + 1])
                    mass_out += rho * v_n_f * dx
                if j > 0:
                    v_s_f = 0.5 * (v_star[i, j] + v_star[i, j - 1])
                    mass_out -= rho * v_s_f * dx

                rhs[k] = mass_out  # div(rho * u_star)

        A = sparse.csr_matrix((vals, (rows, cols)), shape=(N, N))

        # Fix reference pressure (cell 0,0)
        ref = idx(0, 0)
        A[ref, :] = 0
        A[ref, ref] = 1.0
        rhs[ref] = 0.0

        pp_flat = spsolve(A, rhs)
        pp = pp_flat.reshape((nx, ny))

        # ==============================================================
        # Step 5: Correct pressure and velocities
        # ==============================================================
        p += alpha_p * pp

        # Correct cell velocities using cell-centered pressure correction gradient
        # Interior cells: centered difference
        u[1:-1, :] = u_star[1:-1, :] - d[1:-1, :] * (pp[2:, :] - pp[:-2, :]) / (2.0 * dx)
        # Boundary cells: one-sided difference
        u[0, :] = u_star[0, :] - d[0, :] * (pp[1, :] - pp[0, :]) / dx
        u[-1, :] = u_star[-1, :] - d[-1, :] * (pp[-1, :] - pp[-2, :]) / dx

        v[:, 1:-1] = v_star[:, 1:-1] - d[:, 1:-1] * (pp[:, 2:] - pp[:, :-2]) / (2.0 * dy)
        v[:, 0] = v_star[:, 0] - d[:, 0] * (pp[:, 1] - pp[:, 0]) / dy
        v[:, -1] = v_star[:, -1] - d[:, -1] * (pp[:, -1] - pp[:, -2]) / dy

        # Enforce BCs
        u[0, :] = 0.0
        u[-1, :] = 0.0
        u[:, 0] = 0.0
        u[:, -1] = U_lid
        v[0, :] = 0.0
        v[-1, :] = 0.0
        v[:, 0] = 0.0
        v[:, -1] = 0.0

        # ==============================================================
        # Step 6: Convergence check
        # ==============================================================
        mass_residual = np.sum(np.abs(rhs))
        residual_history.append(mass_residual)

        if verbose and (iteration % 200 == 0 or iteration < 5):
            elapsed = time.time() - t0
            print(f"  Iter {iteration:5d}: residual = {mass_residual:.6e}  ({elapsed:.1f}s)")

        if mass_residual < tol:
            elapsed = time.time() - t0
            if verbose:
                print(f"  CONVERGED at iter {iteration}, residual {mass_residual:.6e} ({elapsed:.1f}s)")
            return {
                'u': u, 'v': v, 'p': p,
                'xc': xc, 'yc': yc,
                'converged': True,
                'iterations': iteration,
                'residual_history': residual_history,
                'nx': nx, 'ny': ny, 'Re': Re,
                'elapsed': elapsed
            }

    elapsed = time.time() - t0
    if verbose:
        print(f"  NOT converged after {max_iter} iters. Residual: {mass_residual:.6e} ({elapsed:.1f}s)")

    return {
        'u': u, 'v': v, 'p': p,
        'xc': xc, 'yc': yc,
        'converged': False,
        'iterations': max_iter,
        'residual_history': residual_history,
        'nx': nx, 'ny': ny, 'Re': Re,
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
    xc = result['xc']
    yc = result['yc']
    nx = result['nx']
    ny = result['ny']

    i_center = nx // 2
    u_centerline = u[i_center, :]
    j_center = ny // 2
    v_centerline = v[:, j_center]

    ghia_y, ghia_u, ghia_x, ghia_v = ghia_data_re100()

    # Interpolate to Ghia points (add boundary values for interpolation)
    y_with_walls = np.concatenate(([0.0], yc, [1.0]))
    u_with_walls = np.concatenate(([0.0], u_centerline, [1.0]))
    u_interp = np.interp(ghia_y, y_with_walls, u_with_walls)

    x_with_walls = np.concatenate(([0.0], xc, [1.0]))
    v_with_walls = np.concatenate(([0.0], v_centerline, [0.0]))
    v_interp = np.interp(ghia_x, x_with_walls, v_with_walls)

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
    print(f"  {'y':>8s}  {'Ghia':>10s}  {'Solver':>10s}  {'AbsErr':>10s}")
    for i in range(len(ghia_y)):
        marker = " *" if u_errors[i] > 0.02 else ""
        print(f"  {ghia_y[i]:8.4f}  {ghia_u[i]:10.4f}  {u_interp[i]:10.4f}  {u_errors[i]:10.6f}{marker}")

    print(f"\n  U max absolute error: {np.max(u_errors):.6f}")
    print(f"  U mean absolute error: {np.mean(u_errors):.6f}")
    if u_rel_errors:
        print(f"  U max relative error: {max(u_rel_errors):.2f}%")
        print(f"  U mean relative error: {np.mean(u_rel_errors):.2f}%")

    print("\n  V-velocity comparison (horizontal centerline, y=0.5):")
    print(f"  {'x':>8s}  {'Ghia':>10s}  {'Solver':>10s}  {'AbsErr':>10s}")
    for i in range(len(ghia_x)):
        marker = " *" if v_errors[i] > 0.02 else ""
        print(f"  {ghia_x[i]:8.4f}  {ghia_v[i]:10.4f}  {v_interp[i]:10.4f}  {v_errors[i]:10.6f}{marker}")

    print(f"\n  V max absolute error: {np.max(v_errors):.6f}")
    print(f"  V mean absolute error: {np.mean(v_errors):.6f}")
    if v_rel_errors:
        print(f"  V max relative error: {max(v_rel_errors):.2f}%")
        print(f"  V mean relative error: {np.mean(v_rel_errors):.2f}%")

    max_u_rel = max(u_rel_errors) if u_rel_errors else 0.0
    max_v_rel = max(v_rel_errors) if v_rel_errors else 0.0
    overall = max(max_u_rel, max_v_rel)

    return {
        'u_interp': u_interp,
        'v_interp': v_interp,
        'u_max_rel_error': max_u_rel,
        'v_max_rel_error': max_v_rel,
        'overall_max_error': overall,
        'passed_5pct': overall < 5.0,
        'passed_2pct': overall < 2.0
    }


if __name__ == '__main__':
    print("=" * 70)
    print("Python Reference SIMPLE Solver - Lid-Driven Cavity")
    print("=" * 70)

    # ---- Test 1: 17x17 coarse ----
    print("\n--- Test 1: 17x17 coarse grid ---")
    r17 = solve_cavity_simple(nx=17, ny=17, Re=100.0, max_iter=8000,
                               alpha_u=0.5, alpha_p=0.3, tol=1e-5)
    if r17['converged']:
        c17 = compare_with_ghia(r17)
        print(f"\n  Overall max relative error: {c17['overall_max_error']:.2f}%")
        print(f"  PASS (<5%): {c17['passed_5pct']}")
    else:
        rh = r17['residual_history']
        print(f"  Trend: first={rh[0]:.3e} → last={rh[-1]:.3e}")

    # ---- Test 2: 33x33 medium ----
    print("\n--- Test 2: 33x33 medium grid ---")
    r33 = solve_cavity_simple(nx=33, ny=33, Re=100.0, max_iter=10000,
                               alpha_u=0.5, alpha_p=0.3, tol=1e-6)
    if r33['converged']:
        c33 = compare_with_ghia(r33)
        print(f"\n  Overall max relative error: {c33['overall_max_error']:.2f}%")
        print(f"  PASS (<2%): {c33['passed_2pct']}")
    else:
        rh = r33['residual_history']
        print(f"  Trend: first={rh[0]:.3e} → last={rh[-1]:.3e}")

    # ---- Test 3: 65x65 fine (if time permits) ----
    print("\n--- Test 3: 65x65 fine grid ---")
    r65 = solve_cavity_simple(nx=65, ny=65, Re=100.0, max_iter=15000,
                               alpha_u=0.5, alpha_p=0.3, tol=1e-7)
    if r65['converged']:
        c65 = compare_with_ghia(r65)
        print(f"\n  Overall max relative error: {c65['overall_max_error']:.2f}%")
        print(f"  PASS (<2%): {c65['passed_2pct']}")
    else:
        rh = r65['residual_history']
        print(f"  Trend: first={rh[0]:.3e} → last={rh[-1]:.3e}")

    print("\n" + "=" * 70)
    print("Reference solver complete")
    print("=" * 70)
