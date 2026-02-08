"""Diagnostic: Run cavity solver for minimal steps and check if flow develops."""
import pycfdrs
import time

print("=" * 60)
print("CAVITY SOLVER DIAGNOSTIC")
print("=" * 60)

# Very small grid for fast diagnosis
nx = 9
print(f"\nGrid: {nx}x{nx}, Re=100")
solver = pycfdrs.CavitySolver2D(
    reynolds=100.0,
    nx=nx,
    ny=nx,
    lid_velocity=1.0,
    cavity_size=1.0,
)

t0 = time.time()
result = solver.solve()
elapsed = time.time() - t0

print(f"Time: {elapsed:.2f}s")
print(f"Converged: {result.converged}")
print(f"L2 error vs Ghia: {result.l2_error:.6f}")

# Check u-velocity on vertical centerline
print(f"\nU-velocity on vertical centerline (normalized by U_lid):")
for i, (y, u) in enumerate(zip(result.y_coords, result.u_centerline)):
    print(f"  y={y:.3f}  u={u:+.6f}")

print(f"\nV-velocity on horizontal centerline:")
for i, (x, v) in enumerate(zip(result.x_coords, result.v_centerline)):
    print(f"  x={x:.3f}  v={v:+.6f}")

# Key diagnostic: are velocities non-zero in the interior?
interior_u = [abs(u) for y, u in zip(result.y_coords, result.u_centerline) if 0.1 < y < 0.9]
interior_v = [abs(v) for x, v in zip(result.x_coords, result.v_centerline) if 0.1 < x < 0.9]
max_interior_u = max(interior_u) if interior_u else 0
max_interior_v = max(interior_v) if interior_v else 0

print(f"\nMax interior |u|: {max_interior_u:.6f}")
print(f"Max interior |v|: {max_interior_v:.6f}")

if max_interior_u > 0.01:
    print("PASS: Flow is developing (u-velocity non-trivial)")
else:
    print("FAIL: Flow not developing (u-velocity near zero)")

# Ghia benchmark quick check
u_centerline_max = max(result.u_centerline)
u_centerline_min = min(result.u_centerline)
v_centerline_max = max(result.v_centerline)
v_centerline_min = min(result.v_centerline)

print(f"\nu_max={u_centerline_max:.6f}, u_min={u_centerline_min:.6f}")
print(f"v_max={v_centerline_max:.6f}, v_min={v_centerline_min:.6f}")
print(f"\nGhia Re=100 reference: u_min ~ -0.2109, v_max ~ 0.1753, v_min ~ -0.2453")
