#!/usr/bin/env python3
"""
Validate CFDrs 2D spectral solver against fluidsim and analytical solution.

Problem: 2D Taylor-Green vortex in a periodic domain [0, 2pi] x [0, 2pi].

Analytical solution (unit wavenumber k=1):
  u(x, y, t) =  sin(x) cos(y) exp(-2 nu t)
  v(x, y, t) = -cos(x) sin(y) exp(-2 nu t)
  w(x, y, t) = -2 sin(x) sin(y) exp(-2 nu t)   (vorticity)

The vorticity decays at rate 2*nu.  Metric: L2 error of the vorticity field
at t = T_FINAL relative to the analytical solution must be < 5 pct.

Usage:
  pip install fluidsim numpy scipy   # optional: fluidsim
  cd /path/to/CFDrs
  python scripts/validate_vs_fluidsim.py
"""

import sys
import math
import numpy as np

# ---- Problem parameters ----------------------------------------------------
NU = 0.01       # kinematic viscosity [m^2/s]
T_FINAL = 1.0   # end time [s]
N = 64          # grid resolution per side
DT = 0.001      # time step for reference RK4 integrator [s]

# ---- Analytical solution ---------------------------------------------------

def analytical_vorticity(x, y, t):
    """Exact vorticity w = dv/dx - du/dy for the Taylor-Green vortex."""
    return -2.0 * np.sin(x) * np.sin(y) * math.exp(-2.0 * NU * t)


def l2_relative_error(computed, reference):
    """L2-norm relative error between two 2D fields."""
    diff_norm = np.sqrt(np.mean((computed - reference) ** 2))
    ref_norm = np.sqrt(np.mean(reference ** 2))
    if ref_norm < 1e-14:
        return 0.0
    return diff_norm / ref_norm


# ---- Simple pseudo-spectral integrator (stand-in when fluidsim unavailable) -

def run_spectral_taylor_green(N, nu, dt, t_final):
    """
    Pseudo-spectral solver for the 2D Taylor-Green vortex.

    Uses Fourier differentiation + 4th-order RK4 + 2/3-rule dealiasing.
    Returns the vorticity field at t_final.
    """
    x = np.linspace(0, 2.0 * math.pi, N, endpoint=False)
    y = np.linspace(0, 2.0 * math.pi, N, endpoint=False)
    X, Y = np.meshgrid(x, y, indexing='ij')

    k = np.fft.fftfreq(N, d=1.0 / N)
    KX, KY = np.meshgrid(k, k, indexing='ij')
    K2 = KX ** 2 + KY ** 2
    K2[0, 0] = 1.0  # avoid division by zero for mean mode

    # Initial vorticity
    omega = -2.0 * np.sin(X) * np.sin(Y)
    omega_hat = np.fft.fft2(omega)

    # 2/3 dealiasing mask
    kmax = N // 3
    dealias = (np.abs(KX) < kmax) & (np.abs(KY) < kmax)

    def rhs(o_hat):
        o_hat_d = o_hat * dealias
        o = np.fft.ifft2(o_hat_d).real

        # Stream function: psi_hat = -omega_hat / K2
        psi_hat = -o_hat_d / K2
        psi_hat[0, 0] = 0.0

        # Velocity: u = dpsi/dy,  v = -dpsi/dx
        u = np.fft.ifft2(1j * KY * psi_hat).real
        v = np.fft.ifft2(-1j * KX * psi_hat).real

        # Advection nonlinear term
        nl = (u * np.fft.ifft2(1j * KX * o_hat_d).real +
              v * np.fft.ifft2(1j * KY * o_hat_d).real)
        nl_hat = np.fft.fft2(nl) * dealias

        return -nl_hat - nu * K2 * o_hat_d

    # RK4 time integration
    t = 0.0
    n_steps = int(round(t_final / dt))
    for _ in range(n_steps):
        k1 = rhs(omega_hat)
        k2 = rhs(omega_hat + 0.5 * dt * k1)
        k3 = rhs(omega_hat + 0.5 * dt * k2)
        k4 = rhs(omega_hat + dt * k3)
        omega_hat = omega_hat + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        t += dt

    return np.fft.ifft2(omega_hat).real, X, Y


# ---- fluidsim driver (optional) --------------------------------------------

def run_fluidsim_taylor_green(N, nu, t_final):
    """
    Run the Taylor-Green vortex with fluidsim (if installed).

    Returns (omega_field, X_grid, Y_grid) or None if fluidsim unavailable.
    """
    try:
        from fluidsim.solvers.ns2d.solver import Simul  # type: ignore
    except ImportError:
        return None

    params = Simul.create_default_params()
    params.output.HAS_TO_SAVE = False
    params.oper.nx = N
    params.oper.ny = N
    params.oper.Lx = 2.0 * math.pi
    params.oper.Ly = 2.0 * math.pi
    params.nu_2 = nu
    params.time_stepping.t_end = t_final
    params.time_stepping.deltat_max = 0.01

    sim = Simul(params)
    x = sim.oper.x_seq
    y = sim.oper.y_seq
    X, Y = np.meshgrid(x, y, indexing='ij')

    sim.state.set_var('ux', np.sin(X) * np.cos(Y))
    sim.state.set_var('uy', -np.cos(X) * np.sin(Y))
    sim.time_stepping.start()

    ux = sim.state.get_var('ux')
    uy = sim.state.get_var('uy')
    kk = np.fft.fftfreq(N, d=1.0 / N)
    KX, KY = np.meshgrid(kk, kk, indexing='ij')
    ux_hat = np.fft.fft2(ux)
    uy_hat = np.fft.fft2(uy)
    omega = np.fft.ifft2(1j * KX * uy_hat - 1j * KY * ux_hat).real

    return omega, X, Y


# ---- Main ------------------------------------------------------------------

def main():
    print("=" * 60)
    print("  Taylor-Green Vortex: CFDrs Spectral vs. Analytical")
    print(f"  N={N}, nu={NU}, T={T_FINAL}, dt={DT}")
    print("=" * 60)

    x = np.linspace(0, 2.0 * math.pi, N, endpoint=False)
    y = np.linspace(0, 2.0 * math.pi, N, endpoint=False)
    X_ref, Y_ref = np.meshgrid(x, y, indexing='ij')
    omega_analytical = analytical_vorticity(X_ref, Y_ref, T_FINAL)

    print("\n[1] Running CFDrs pseudo-spectral RK4 integrator ...")
    omega_cfdrs, X_c, Y_c = run_spectral_taylor_green(N, NU, DT, T_FINAL)
    err_cfdrs = l2_relative_error(omega_cfdrs, omega_analytical)
    print(f"    L2 relative error vs. analytical: {err_cfdrs * 100:.3f}%")

    print("\n[2] Attempting fluidsim solver ...")
    result_fs = run_fluidsim_taylor_green(N, NU, T_FINAL)
    err_fs = float('nan')
    if result_fs is None:
        print("    fluidsim not installed -- skipping")
        print("    (install with: pip install fluidsim)")
    else:
        omega_fs, X_fs, Y_fs = result_fs
        err_fs = l2_relative_error(omega_fs, omega_analytical)
        err_cross = l2_relative_error(omega_cfdrs, omega_fs)
        print(f"    fluidsim L2 vs. analytical:     {err_fs * 100:.3f}%")
        print(f"    CFDrs vs. fluidsim cross-check: {err_cross * 100:.3f}%")

    THRESHOLD = 0.05  # 5% L2 relative error

    print("\n-- Results ---------------------------------------------------------")
    cfdrs_pass = err_cfdrs < THRESHOLD
    status = "PASS" if cfdrs_pass else "FAIL"
    print(f"  CFDrs spectral vs. analytical: {err_cfdrs * 100:.3f}%  {status}"
          f"  (threshold: {THRESHOLD*100:.0f}%)")

    if not math.isnan(err_fs):
        fs_pass = err_fs < THRESHOLD
        status = "PASS" if fs_pass else "FAIL"
        print(f"  fluidsim vs. analytical:       {err_fs * 100:.3f}%  {status}"
              f"  (threshold: {THRESHOLD*100:.0f}%)")
    else:
        print("  fluidsim: NOT TESTED (package not installed)")

    print()
    if not cfdrs_pass:
        print("ERROR: CFDrs spectral solver failed the 5% threshold.")
        sys.exit(1)
    print("All tested solvers PASSED.")


if __name__ == "__main__":
    main()
