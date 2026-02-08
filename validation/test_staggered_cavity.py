"""Test the new staggered-grid Rust cavity solver against Ghia and Python reference."""
import numpy as np
import sys
sys.path.insert(0, 'validation')

print("=" * 70)
print("Rust Staggered-Grid SIMPLE Cavity Solver vs Ghia (1982)")
print("=" * 70)

# Ghia Re=100 benchmark data
ghia_y_u = [0.0000, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719,
            0.2813, 0.4531, 0.5000, 0.6172, 0.7344, 0.8516,
            0.9531, 0.9609, 0.9688, 0.9766, 1.0000]
ghia_u   = [0.0000, -0.03717, -0.04192, -0.04775, -0.06434, -0.10150,
            -0.15662, -0.21090, -0.20581, -0.13641, 0.00332, 0.23151,
            0.68717, 0.73722, 0.78871, 0.84123, 1.00000]

ghia_x_v = [0.0000, 0.0625, 0.0703, 0.0781, 0.0938, 0.1563,
            0.2266, 0.2344, 0.5000, 0.8047, 0.8594, 0.9063,
            0.9453, 0.9531, 0.9609, 0.9688, 1.0000]
ghia_v   = [0.0000, 0.09233, 0.10091, 0.10890, 0.12317, 0.16077,
            0.17507, 0.17527, 0.05454, -0.24533, -0.22445, -0.16914,
            -0.10313, -0.08864, -0.07391, -0.05906, 0.00000]

import pycfdrs

for nx in [32, 64, 128]:
    print(f"\n--- Grid {nx}x{nx} ---")
    solver = pycfdrs.CavitySolver2D(reynolds=100.0, nx=nx, ny=nx)
    result = solver.solve()
    
    y = np.array(result.y_coords)
    u = np.array(result.u_centerline)
    x = np.array(result.x_coords)
    v = np.array(result.v_centerline)
    
    print(f"  converged={result.converged}, L2_error={result.l2_error:.4e}")
    print(f"  u range: [{u.min():.5f}, {u.max():.5f}]")
    print(f"  v range: [{v.min():.5f}, {v.max():.5f}]")
    
    # Interpolate to Ghia points
    u_interp = np.interp(ghia_y_u, y, u)
    v_interp = np.interp(ghia_x_v, x, v)
    
    # Compute errors
    mask_u = np.array([abs(ug) > 0.01 for ug in ghia_u])  
    mask_v = np.array([abs(vg) > 0.01 for vg in ghia_v])
    
    u_err = np.abs((u_interp[mask_u] - np.array(ghia_u)[mask_u]) / np.array(ghia_u)[mask_u]) * 100
    v_err = np.abs((v_interp[mask_v] - np.array(ghia_v)[mask_v]) / np.array(ghia_v)[mask_v]) * 100
    
    print(f"  U max rel error: {u_err.max():.1f}%  mean: {u_err.mean():.1f}%")
    print(f"  V max rel error: {v_err.max():.1f}%  mean: {v_err.mean():.1f}%")
    
    if nx == 64:
        print(f"\n  U along vertical centerline (64x64):")
        for k in range(len(ghia_y_u)):
            ref = ghia_u[k]
            sol = u_interp[k]
            rel = abs((sol - ref) / ref) * 100 if abs(ref) > 0.01 else 0
            flag = " ***" if rel > 20 else (" **" if rel > 5 else "")
            print(f"    y={ghia_y_u[k]:.4f}  Ghia={ref:10.5f}  Rust={sol:10.5f}  err={rel:.1f}%{flag}")
        
        print(f"\n  V along horizontal centerline (64x64):")
        for k in range(len(ghia_x_v)):
            ref = ghia_v[k]
            sol = v_interp[k]
            rel = abs((sol - ref) / ref) * 100 if abs(ref) > 0.01 else 0
            flag = " ***" if rel > 20 else (" **" if rel > 5 else "")
            print(f"    x={ghia_x_v[k]:.4f}  Ghia={ref:10.5f}  Rust={sol:10.5f}  err={rel:.1f}%{flag}")

print("\n" + "=" * 70)
print("DONE")
