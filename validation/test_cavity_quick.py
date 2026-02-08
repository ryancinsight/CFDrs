"""Quick cavity solver test to check if SIMPLEC now converges."""
import pycfdrs
import time

print("=" * 60)
print("CAVITY SIMPLEC CONVERGENCE TEST")
print("=" * 60)

# Small grid first (17x17) to test convergence quickly
for nx in [17, 33]:
    print(f"\n--- Grid: {nx}x{nx}, Re=100 ---")
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
    
    print(f"  Time: {elapsed:.2f}s")
    print(f"  u_max on centerline: {result.u_centerline_max:.6f}")
    print(f"  v_max on centerline: {result.v_centerline_max:.6f}")
    print(f"  max_pressure: {result.max_pressure:.6e}")
    print(f"  min_pressure: {result.min_pressure:.6e}")
    
    # Ghia et al. (1982) reference for Re=100:
    # u_min on vertical centerline ≈ -0.2109 (at y ≈ 0.453)
    # v_max on horizontal centerline ≈ 0.1753 (at x ≈ 0.234)
    # v_min on horizontal centerline ≈ -0.2453 (at x ≈ 0.805)
    
    # Check if we got meaningful non-zero velocities
    has_flow = abs(result.u_centerline_max) > 0.01 or abs(result.v_centerline_max) > 0.01
    print(f"  Has meaningful flow: {has_flow}")
    
    if has_flow:
        print(f"  PASS: Solver produced non-trivial flow field")
    else:
        print(f"  FAIL: Solver produced near-zero flow (not converging)")

print("\nDone.")
