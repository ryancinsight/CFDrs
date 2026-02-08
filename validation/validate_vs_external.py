
import sys
import os
import numpy as np
import pytest
from datetime import datetime

# Add project root to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import pycfdrs
try:
    import pycfdrs
    print(f"[INFO] Successfully imported pycfdrs version {pycfdrs.__version__}")
except ImportError as e:
    print(f"[ERROR] Failed to import pycfdrs: {e}")
    sys.exit(1)

# Import reference implementation
try:
    from validation.references.cavity_python_cfd import solve_cavity_flow_reference
    print("[INFO] Successfully imported reference solver")
except ImportError as e:
    print(f"[ERROR] Failed to import reference solver: {e}")
    sys.exit(1)

def run_comparison():
    print("\n" + "="*80)
    print("CODE-TO-CODE VALIDATION: Lid-Driven Cavity (Re=100)")
    print("Comparing pycfdrs (Rust) vs. Literature Python Implementation")
    print("="*80)

    # 1. Setup Parameters
    nx = 41
    ny = 41
    reynolds = 100.0
    lid_velocity = 1.0
    cavity_size = 1.0
    
    # 2. Run Reference Solver (Python FDM)
    print(f"\n[REFERENCE] Running Python FDM solver (nx={nx}, ny={ny})...")
    # Reference solver returns full 2D arrays: u, v, p, x, y
    start_time = datetime.now()
    ref_u, ref_v, ref_p, ref_x, ref_y = solve_cavity_flow_reference(
        nx=nx, ny=ny, nt=1000, nit=50, c=lid_velocity, 
        length=cavity_size, rho=1.0, nu=lid_velocity*cavity_size/reynolds
    )
    print(f"[REFERENCE] Finished in {(datetime.now() - start_time).total_seconds():.3f}s")

    # Extract centerlines from Reference
    # For odd grid points (41), center index is (41-1)/2 = 20
    center_idx_x = (nx - 1) // 2
    center_idx_y = (ny - 1) // 2
    
    # U-velocity along vertical centerline (x = 0.5)
    ref_u_centerline = ref_u[:, center_idx_x]
    
    # V-velocity along horizontal centerline (y = 0.5)
    ref_v_centerline = ref_v[center_idx_y, :]

    # 3. Run pycfdrs Solver (Rust FVM/SIMPLEC)
    print(f"\n[PYCFDRS] Running Rust SIMPLEC solver (nx={nx}, ny={ny})...")
    start_time = datetime.now()
    solver = pycfdrs.CavitySolver2D(
        reynolds=reynolds, nx=nx, ny=ny, 
        lid_velocity=lid_velocity, cavity_size=cavity_size
    )
    result = solver.solve()
    print(f"[PYCFDRS] Finished in {(datetime.now() - start_time).total_seconds():.3f}s")
    
    if not result.converged:
        print("[WARNING] pycfdrs did not claim convergence!")

    # 4. Compare Centerlines
    # pycfdrs returns centerlines as 1D arrays
    rust_u_centerline = np.array(result.u_centerline)
    rust_v_centerline = np.array(result.v_centerline)
    
    # Ensure dimensions match
    if len(rust_u_centerline) != len(ref_u_centerline):
        print(f"[ERROR] Dimension mismatch: Rust {len(rust_u_centerline)} vs Ref {len(ref_u_centerline)}")
        # Interpolate if needed (basic check)
        return False

    # Compute Normalized RMS Error
    u_rmse = np.sqrt(np.mean((rust_u_centerline - ref_u_centerline)**2)) / lid_velocity
    v_rmse = np.sqrt(np.mean((rust_v_centerline - ref_v_centerline)**2)) / lid_velocity
    
    print("\n" + "-"*40)
    print("COMPARISON RESULTS")
    print("-"*-40)
    print(f"U-Velocity RMS Error: {u_rmse:.6f} ({u_rmse*100:.3f}%)")
    print(f"V-Velocity RMS Error: {v_rmse:.6f} ({v_rmse*100:.3f}%)")
    
    # 5. Validation Criteria (Allowing < 5% difference between FVM and FDM)
    # Different discretizations (FVM vs FDM) will have some inherent difference.
    # 5% is a reasonable threshold for "code-to-code" agreement of different schemes.
    u_passed = u_rmse < 0.05
    v_passed = v_rmse < 0.05
    
    if u_passed and v_passed:
        print("\n[SUCCESS] Code-to-Code Validation PASSED!")
        return 0
    else:
        print("\n[FAILURE] Validation thresholds exceeded.")
        # Print pointwise comparison for debug
        print("\nDebug U-Velocity (Top 5 deviations):")
        diff = np.abs(rust_u_centerline - ref_u_centerline)
        indices = np.argsort(diff)[-5:]
        for i in indices:
            print(f"  y={ref_y[i]:.3f}: Rust={rust_u_centerline[i]:.3f}, Ref={ref_u_centerline[i]:.3f}, Diff={diff[i]:.3f}")
        return 1

if __name__ == "__main__":
    sys.exit(run_comparison())
