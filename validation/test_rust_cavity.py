"""Test Rust cavity solver output and compare with Python reference."""
import numpy as np
import pycfdrs
import sys
sys.path.insert(0, 'validation')
from reference_cavity_psiomega import solve_cavity_psi_omega, ghia_data_re100

print("=" * 70)
print("Rust SIMPLEC Cavity Solver Test")
print("=" * 70)

# Test at 33x33 (matching Python reference at N=33)
print("\n--- Rust SIMPLEC solver (33x33) ---")
solver = pycfdrs.CavitySolver2D(reynolds=100.0, nx=33, ny=33)
result = solver.solve()
print(f"  converged={result.converged}, L2_error={result.l2_error:.4e}")
print(f"  u_centerline: len={len(result.u_centerline)}, range=[{min(result.u_centerline):.4f}, {max(result.u_centerline):.4f}]")
print(f"  v_centerline: len={len(result.v_centerline)}, range=[{min(result.v_centerline):.4f}, {max(result.v_centerline):.4f}]")

# Compare with Ghia
ghia_y, ghia_u, ghia_x, ghia_v = ghia_data_re100()
rust_u = np.interp(ghia_y, np.array(result.y_coords), np.array(result.u_centerline))
rust_v = np.interp(ghia_x, np.array(result.x_coords), np.array(result.v_centerline))

print("\n  U along vertical centerline:")
for k in range(len(ghia_y)):
    ref = ghia_u[k]
    sol = rust_u[k]
    rel = abs((sol - ref) / ref) * 100 if abs(ref) > 0.01 else 0
    flag = " ***" if rel > 20 else (" **" if rel > 5 else "")
    print(f"  y={ghia_y[k]:.4f}  Ghia={ref:10.5f}  Rust={sol:10.5f}  err={rel:.1f}%{flag}")

print("\n  V along horizontal centerline:")
for k in range(len(ghia_x)):
    ref = ghia_v[k]
    sol = rust_v[k]
    rel = abs((sol - ref) / ref) * 100 if abs(ref) > 0.01 else 0
    flag = " ***" if rel > 20 else (" **" if rel > 5 else "")
    print(f"  x={ghia_x[k]:.4f}  Ghia={ref:10.5f}  Rust={sol:10.5f}  err={rel:.1f}%{flag}")

# Also get Python reference at N=33 for direct comparison
print("\n--- Python ψ-ω reference (N=33) ---")
ref = solve_cavity_psi_omega(N=33, Re=100.0, max_iter=200000, tol=1e-8, verbose=False)
ref_u = np.interp(ghia_y, ref['y'], ref['u'][:, 16])
ref_v = np.interp(ghia_x, ref['x'], ref['v'][16, :])
print(f"  v at (0.5,0.5): Python={ref_v[8]:.5f}, Rust={rust_v[8]:.5f}, Ghia={ghia_v[8]:.5f}")
