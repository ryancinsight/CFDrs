"""Quick grid convergence study for v at (0.5, 0.5)."""
import sys, numpy as np
sys.path.insert(0, 'validation')
from reference_cavity_psiomega import solve_cavity_psi_omega, ghia_data_re100

gy, gu, gx, gv = ghia_data_re100()
i_half_v = np.argmin(np.abs(gx - 0.5))
i_half_u = np.argmin(np.abs(gy - 0.5))

for N in [33, 65, 129]:
    r = solve_cavity_psi_omega(N=N, Re=100.0, max_iter=200000, tol=1e-8,
                                use_woods=True, verbose=False)
    ic = (N - 1) // 2
    u_cl = r['u'][:, ic]
    v_cl = r['v'][ic, :]
    u_i = np.interp(gy, r['y'], u_cl)
    v_i = np.interp(gx, r['x'], v_cl)
    u_err = abs((u_i[i_half_u] - gu[i_half_u]) / gu[i_half_u]) * 100
    v_err = abs((v_i[i_half_v] - gv[i_half_v]) / gv[i_half_v]) * 100
    pmin = r['psi_min']
    iters = r['iterations']
    elapsed = r['elapsed']
    print(f"N={N:4d}: v(0.5,0.5)={v_i[i_half_v]:.5f} err={v_err:.2f}%  "
          f"u(0.5,0.5)={u_i[i_half_u]:.5f} err={u_err:.2f}%  "
          f"psi_min={pmin:.6f}  ({iters} iters, {elapsed:.1f}s)")

print(f"Ghia:   v(0.5,0.5)={gv[i_half_v]:.5f}  "
      f"u(0.5,0.5)={gu[i_half_u]:.5f}  psi_min=-0.103400")
