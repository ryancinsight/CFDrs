"""Run MAC solver at multiple resolutions and compute Richardson extrapolation."""
import sys
import numpy as np
sys.path.insert(0, 'validation')
from reference_cavity_mac import solve_cavity_mac, compare_with_ghia, ghia_data_re100

# Test 1: 128x128
print("=== 128x128 quadratic ghost, tol=1e-7 ===")
r128 = solve_cavity_mac(nx=128, ny=128, Re=100.0, max_iter=400000, tol=1e-7)
c128 = compare_with_ghia(r128)

# Test 2: 64x64 (tight tolerance)
print("\n=== 64x64 quadratic ghost, tol=1e-9 ===")
r64 = solve_cavity_mac(nx=64, ny=64, Re=100.0, max_iter=400000, tol=1e-9)
c64 = compare_with_ghia(r64)

# Richardson extrapolation: u_exact ~ (4*u_fine - u_coarse) / 3
gy, gu, gx, gv = ghia_data_re100()
u_rich = (4.0 * c128['u_interp'] - c64['u_interp']) / 3.0
v_rich = (4.0 * c128['v_interp'] - c64['v_interp']) / 3.0

print("\n=== Richardson Extrapolation (64 + 128) ===")
print("  U-velocity at Ghia points:")
for i in range(len(gy)):
    ref = gu[i]
    sol = u_rich[i]
    err = abs(sol - ref)
    rel = abs((sol - ref) / ref) * 100 if abs(ref) > 0.01 else 0.0
    flag = " ***" if rel > 5 else (" **" if rel > 2 else "")
    print(f"  y={gy[i]:.4f}  Ghia={ref:10.5f}  Rich={sol:10.5f}  err={rel:.2f}%{flag}")

print("\n  V-velocity at Ghia points:")
for i in range(len(gx)):
    ref = gv[i]
    sol = v_rich[i]
    err = abs(sol - ref)
    rel = abs((sol - ref) / ref) * 100 if abs(ref) > 0.01 else 0.0
    flag = " ***" if rel > 5 else (" **" if rel > 2 else "")
    print(f"  x={gx[i]:.4f}  Ghia={ref:10.5f}  Rich={sol:10.5f}  err={rel:.2f}%{flag}")

u_rel = [abs((u_rich[i]-gu[i])/gu[i])*100 for i in range(len(gu)) if abs(gu[i]) > 0.01]
v_rel = [abs((v_rich[i]-gv[i])/gv[i])*100 for i in range(len(gv)) if abs(gv[i]) > 0.01]
print(f"\n  Richardson U max rel: {max(u_rel):.3f}%")
print(f"  Richardson V max rel: {max(v_rel):.3f}%")
print(f"  Richardson overall:   {max(max(u_rel), max(v_rel)):.3f}%")

# Convergence order
e64 = c64['overall_max_error']
e128 = c128['overall_max_error']
if e128 > 0 and e64 > 0:
    order = np.log2(e64 / e128)
    print(f"\n  Convergence order (64->128): {order:.2f}")
