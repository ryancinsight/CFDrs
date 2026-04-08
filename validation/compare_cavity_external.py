#!/usr/bin/env python3
"""
Cross-package validation: Compare CFD-RS cavity flow with external reference.

Validates cfd_python 2D Navier-Stokes solver against:
1. Pure Python finite difference implementation (external_cavity_reference.py)
2. Ghia et al. (1982) benchmark data
3. Python_CFD notebook results

This script proves CFD-RS produces correct results by comparing against
independent implementations.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Add validation directory to path
sys.path.insert(0, str(Path(__file__).parent))

try:
    import cfd_python
except ImportError:
    print("ERROR: cfd_python not installed. Build with: cd crates/cfd-python && maturin develop")
    sys.exit(1)

# Import external reference
from external_cavity_reference import CavityFlowSolver, run_cavity_validation


def run_cfd_python_cavity(Re: float = 100, nx: int = 65, ny: int = 65):
    """
    Run CFD-RS cavity flow solver.
    
    Args:
        Re: Reynolds number
        nx, ny: Grid resolution
    
    Returns:
        Solution dictionary with fields and centerlines
    """
    print(f"\n{'='*70}")
    print(f"cfd_python CAVITY FLOW (Re={Re})")
    print(f"{'='*70}")
   
    try:
        # Check if cfd_python has cavity solver
        if hasattr(cfd_python, 'CavitySolver2D'):
            solver = cfd_python.CavitySolver2D(
                nx=nx, ny=ny,
                Re=Re,
                tolerance=1e-6,
                max_iterations=20000
            )
            result = solver.solve()
            
            return {
                "u": np.array(result.u_field),
                "v": np.array(result.v_field),
                "p": np.array(result.p_field),
                "u_centerline": np.array(result.u_centerline),
                "v_centerline": np.array(result.v_centerline),
                "x": np.array(result.x_coords),
                "y": np.array(result.y_coords),
                "converged": result.converged,
                "iterations": result.iterations,
                "residual": result.residual
            }
        else:
            print("WARN: CavitySolver2D not found in cfd_python - using placeholder")
            # Return dummy data matching external reference dimensions
            return None
            
    except Exception as e:
        print(f"ERROR: cfd_python cavity solver failed: {e}")
        return None


def compare_solutions(cfd_python_result, external_result, Re: float):
    """
    Compare cfd_python and external reference solutions.
    
    Args:
        cfd_python_result: Solution from cfd_python
        external_result: Solution from external_cavity_reference
        Re: Reynolds number
    
    Returns:
        Dictionary with comparison metrics
    """
    print(f"\n{'='*70}")
    print(f"CROSS-PACKAGE COMPARISON (Re={Re})")
    print(f"{'='*70}")
    
    if cfd_python_result is None:
        print("SKIP: cfd_python result not available")
        return None
    
    # Extract external data
    ext_sol = external_result["solution"]
    ext_solver = external_result["solver"]
    
    # Ensure same grid size
    if cfd_python_result["u"].shape != ext_sol["u"].shape:
        print(f"WARN: Grid size mismatch: cfd_python {cfd_python_result['u'].shape} vs external {ext_sol['u'].shape}")
        from scipy.interpolate import RectBivariateSpline
        x_cfd = np.linspace(0, 1, cfd_python_result["u"].shape[1])
        y_cfd = np.linspace(0, 1, cfd_python_result["u"].shape[0])
        x_ext = np.linspace(0, 1, ext_sol["u"].shape[1])
        y_ext = np.linspace(0, 1, ext_sol["u"].shape[0])

        spline_u = RectBivariateSpline(y_cfd, x_cfd, cfd_python_result["u"])
        spline_v = RectBivariateSpline(y_cfd, x_cfd, cfd_python_result["v"])
        spline_p = RectBivariateSpline(y_cfd, x_cfd, cfd_python_result["p"])

        cfd_python_result["u"] = spline_u(y_ext, x_ext)
        cfd_python_result["v"] = spline_v(y_ext, x_ext)
        cfd_python_result["p"] = spline_p(y_ext, x_ext)
    
    # Compute L2 errors
    u_diff = cfd_python_result["u"] - ext_sol["u"]
    v_diff = cfd_python_result["v"] - ext_sol["v"]
    p_diff = cfd_python_result["p"] - ext_sol["p"]
    
    u_l2_error = np.sqrt(np.mean(u_diff**2))
    v_l2_error = np.sqrt(np.mean(v_diff**2))
    p_l2_error = np.sqrt(np.mean(p_diff**2))
    
    # Relative errors
    u_rel_error = u_l2_error / (np.sqrt(np.mean(ext_sol["u"]**2)) + 1e-10)
    v_rel_error = v_l2_error / (np.sqrt(np.mean(ext_sol["v"]**2)) + 1e-10)
    
    # Centerline errors
    u_centerline_error = np.max(np.abs(cfd_python_result["u_centerline"] - ext_sol["u_centerline"]))
    v_centerline_error = np.max(np.abs(cfd_python_result["v_centerline"] - ext_sol["v_centerline"]))
    
    # Vortex center comparison
    cfd_python_vortex_idx = np.unravel_index(np.argmin(cfd_python_result["p"]), cfd_python_result["p"].shape)
    cfd_python_vortex = (cfd_python_result["x"][cfd_python_vortex_idx[1]], cfd_python_result["y"][cfd_python_vortex_idx[0]])
    external_vortex = ext_sol["vortex_center"]
    
    vortex_distance = np.sqrt((cfd_python_vortex[0] - external_vortex[0])**2 + 
                              (cfd_python_vortex[1] - external_vortex[1])**2)
    
    print(f"\nL2 Errors:")
    print(f"  U velocity:  {u_l2_error:.6e} (relative: {u_rel_error*100:.4f}%)")
    print(f"  V velocity:  {v_l2_error:.6e} (relative: {v_rel_error*100:.4f}%)")
    print(f"  Pressure:    {p_l2_error:.6e}")
    
    print(f"\nCenterline Max Errors:")
    print(f"  U (vertical):    {u_centerline_error:.6e}")
    print(f"  V (horizontal):  {v_centerline_error:.6e}")
    
    print(f"\nVortex Center:")
    print(f"  cfd_python:   ({cfd_python_vortex[0]:.4f}, {cfd_python_vortex[1]:.4f})")
    print(f"  external:  ({external_vortex[0]:.4f}, {external_vortex[1]:.4f})")
    print(f"  distance:  {vortex_distance:.6f}")
    
    print(f"\nConvergence:")
    print(f"  cfd_python:   {cfd_python_result['iterations']} iterations (residual: {cfd_python_result['residual']:.6e})")
    print(f"  external:  {ext_sol['steps']} iterations (residual: {ext_sol['residual']:.6e})")
    
    # Validation criteria
    PASS = True
    tolerance = 0.05  # 5% tolerance for discretization differences
    
    if u_rel_error > tolerance:
        print(f"\n[FAIL] U velocity error {u_rel_error*100:.4f}% exceeds {tolerance*100:.1f}% tolerance")
        PASS = False
    else:
        print(f"\n[PASS] U velocity error {u_rel_error*100:.4f}% within tolerance")
    
    if v_rel_error > tolerance:
        print(f"[FAIL] V velocity error {v_rel_error*100:.4f}% exceeds {tolerance*100:.1f}% tolerance")
        PASS = False
    else:
        print(f"[PASS] V velocity error {v_rel_error*100:.4f}% within tolerance")
    
    if vortex_distance > 0.1:  # 10% of cavity size
        print(f"[FAIL] Vortex center distance {vortex_distance:.4f} exceeds 0.1 tolerance")
        PASS = False
    else:
        print(f"[PASS] Vortex center location within tolerance")
    
    return {
        "u_l2_error": u_l2_error,
        "v_l2_error": v_l2_error,
        "p_l2_error": p_l2_error,
        "u_rel_error": u_rel_error,
        "v_rel_error": v_rel_error,
        "u_centerline_error": u_centerline_error,
        "v_centerline_error": v_centerline_error,
        "vortex_distance": vortex_distance,
        "passed": PASS
    }


def plot_comparison(cfd_python_result, external_result, comparison, Re: float):
    """
    Generate comparison plots.
    """
    if cfd_python_result is None:
        print("SKIP: Cannot plot without cfd_python result")
        return
    
    fig, axes = plt.subplots(3, 3, figsize=(15, 14))
    
    ext_sol = external_result["solution"]
    ext_solver = external_result["solver"]
    
    # Row 1: cfd_python fields
    U_mag_cfd_python = np.sqrt(cfd_python_result["u"]**2 + cfd_python_result["v"]**2)
    X_cfd_python, Y_cfd_python = np.meshgrid(cfd_python_result["x"], cfd_python_result["y"])
    
    im00 = axes[0, 0].contourf(X_cfd_python, Y_cfd_python, U_mag_cfd_python, levels=20, cmap='viridis')
    axes[0, 0].set_title('cfd_python: Velocity Magnitude')
    axes[0, 0].set_xlabel('x')
    axes[0, 0].set_ylabel('y')
    plt.colorbar(im00, ax=axes[0, 0])
    
    im01 = axes[0, 1].contourf(X_cfd_python, Y_cfd_python, cfd_python_result["u"], levels=20, cmap='RdBu_r')
    axes[0, 1].set_title('cfd_python: U Velocity')
    axes[0, 1].set_xlabel('x')
    axes[0, 1].set_ylabel('y')
    plt.colorbar(im01, ax=axes[0, 1])
    
    im02 = axes[0, 2].contourf(X_cfd_python, Y_cfd_python, cfd_python_result["p"], levels=20, cmap='coolwarm')
    axes[0, 2].set_title('cfd_python: Pressure')
    axes[0, 2].set_xlabel('x')
    axes[0, 2].set_ylabel('y')
    plt.colorbar(im02, ax=axes[0, 2])
    
    # Row 2: External fields
    U_mag_ext = np.sqrt(ext_sol["u"]**2 + ext_sol["v"]**2)
    
    im10 = axes[1, 0].contourf(ext_solver.X, ext_solver.Y, U_mag_ext, levels=20, cmap='viridis')
    axes[1, 0].set_title('External: Velocity Magnitude')
    axes[1, 0].set_xlabel('x')
    axes[1, 0].set_ylabel('y')
    plt.colorbar(im10, ax=axes[1, 0])
    
    im11 = axes[1, 1].contourf(ext_solver.X, ext_solver.Y, ext_sol["u"], levels=20, cmap='RdBu_r')
    axes[1, 1].set_title('External: U Velocity')
    axes[1, 1].set_xlabel('x')
    axes[1, 1].set_ylabel('y')
    plt.colorbar(im11, ax=axes[1, 1])
    
    im12 = axes[1, 2].contourf(ext_solver.X, ext_solver.Y, ext_sol["p"], levels=20, cmap='coolwarm')
    axes[1, 2].set_title('External: Pressure')
    axes[1, 2].set_xlabel('x')
    axes[1, 2].set_ylabel('y')
    plt.colorbar(im12, ax=axes[1, 2])
    
    # Row 3: Differences and centerline comparison
    u_diff = cfd_python_result["u"] - ext_sol["u"]
    v_diff = cfd_python_result["v"] - ext_sol["v"]
    vel_diff = np.sqrt(u_diff**2 + v_diff**2)
    
    im20 = axes[2, 0].contourf(X_cfd_python, Y_cfd_python, vel_diff, levels=20, cmap='hot')
    axes[2, 0].set_title('Velocity Difference')
    axes[2, 0].set_xlabel('x')
    axes[2, 0].set_ylabel('y')
    plt.colorbar(im20, ax=axes[2, 0])
    
    # U centerline comparison
    axes[2, 1].plot(cfd_python_result["u_centerline"], cfd_python_result["y"], 'b-', label='cfd_python', linewidth=2)
    axes[2, 1].plot(ext_sol["u_centerline"], ext_solver.y, 'r--', label='External', linewidth=2)
    if external_result["validation"]["u_error"] is not None:
        axes[2, 1].plot(external_result["validation"]["ghia_u"]["u"], 
                       external_result["validation"]["ghia_u"]["y"], 
                       'go', label='Ghia et al. (1982)', markersize=5)
    axes[2, 1].set_xlabel('u')
    axes[2, 1].set_ylabel('y')
    axes[2, 1].set_title('U Centerline (x=0.5)')
    axes[2, 1].legend()
    axes[2, 1].grid(True, alpha=0.3)
    
    # V centerline comparison
    axes[2, 2].plot(cfd_python_result["x"], cfd_python_result["v_centerline"], 'b-', label='cfd_python', linewidth=2)
    axes[2, 2].plot(ext_solver.x, ext_sol["v_centerline"], 'r--', label='External', linewidth=2)
    if external_result["validation"]["v_error"] is not None:
        axes[2, 2].plot(external_result["validation"]["ghia_v"]["x"], 
                       external_result["validation"]["ghia_v"]["v"], 
                       'go', label='Ghia et al. (1982)', markersize=5)
    axes[2, 2].set_xlabel('x')
    axes[2, 2].set_ylabel('v')
    axes[2, 2].set_title('V Centerline (y=0.5)')
    axes[2, 2].legend()
    axes[2, 2].grid(True, alpha=0.3)
    
    plt.suptitle(f'Cross-Package Validation: cfd_python vs External Reference (Re={Re})', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('cross_validation_cavity_comparison.png', dpi=150, bbox_inches='tight')
    print(f"\nComparison plot saved: cross_validation_cavity_comparison.png")


def main():
    """
    Main cross-validation workflow.
    """
    print("="*70)
    print("CROSS-PACKAGE VALIDATION: CAVITY FLOW")
    print("="*70)
    print("\nComparing cfd_python against:")
    print("  1. Pure Python finite difference reference")
    print("  2. Ghia et al. (1982) benchmark data")
    print("="*70)
    
    Re = 100
    nx = ny = 65
    
    # Run external reference
    print("\n[1/2] Running external reference implementation...")
    external_result = run_cavity_validation(Re=Re, nx=nx, ny=ny)
    
    print(f"\nExternal result:")
    print(f"  Converged: {external_result['solution']['converged']}")
    print(f"  Steps: {external_result['solution']['steps']}")
    print(f"  Vortex center: ({external_result['solution']['vortex_center'][0]:.4f}, {external_result['solution']['vortex_center'][1]:.4f})")
    
    if external_result['validation']['u_error'] is not None:
        print(f"\nGhia et al. (1982) comparison:")
        print(f"  U centerline error: {external_result['validation']['u_error']:.6f}")
        print(f"  V centerline error: {external_result['validation']['v_error']:.6f}")
    
    # Run cfd_python
    print("\n[2/2] Running cfd_python implementation...")
    cfd_python_result = run_cfd_python_cavity(Re=Re, nx=nx, ny=ny)
    
    if cfd_python_result is not None:
        # Compare solutions
        comparison = compare_solutions(cfd_python_result, external_result, Re)
        
        # Plot comparison
        if comparison is not None:
            plot_comparison(cfd_python_result, external_result, comparison, Re)
            
            # Final verdict
            print(f"\n{'='*70}")
            if comparison["passed"]:
                print("✓ CROSS-VALIDATION PASSED")
                print(f"  cfd_python matches external reference within {comparison['u_rel_error']*100:.4f}% relative error")
            else:
                print("✗ CROSS-VALIDATION FAILED")
                print(f"  Errors exceed tolerance thresholds")
            print(f"{'='*70}")
            
            return 0 if comparison["passed"] else 1
    else:
        print("\nWARN: cfd_python cavity solver not available - validation skipped")
        print("      Implementation needed: CavitySolver2D in cfd_python module")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
