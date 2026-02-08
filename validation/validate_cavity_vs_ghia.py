#!/usr/bin/env python3
"""
Lid-Driven Cavity Validation vs Ghia et al. (1982)

This script validates the pycfdrs CavitySolver2D against the classical
Ghia et al. (1982) benchmark for lid-driven cavity flow at Re=100.

Reference:
    Ghia, U.K.N.G., Ghia, K.N., Shin, C.T. (1982). "High-Re solutions for
    incompressible flow using the Navier-Stokes equations and a multigrid
    method". Journal of Computational Physics, 48(3):387-411.

Acceptance Criteria:
    - L2 error in U-velocity centerline < 1%
    - L2 error in V-velocity centerline < 1%

Usage:
    .\.venv\Scripts\Activate.ps1
    maturin develop --manifest-path crates/pycfdrs/Cargo.toml
    python validation/validate_cavity_vs_ghia.py
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
from dataclasses import dataclass
from typing import Tuple

try:
    import pycfdrs
    from pycfdrs import CavitySolver2D
    HAS_PYCFDRS = True
except ImportError:
    HAS_PYCFDRS = False


@dataclass
class ValidationResult:
    """Result from cavity validation"""
    name: str
    l2_error: float
    tolerance: float
    passed: bool
    details: str


def validate_cavity_re100() -> ValidationResult:
    """
    Validate lid-driven cavity at Re=100 against Ghia et al. (1982).
    
    The Ghia benchmark provides high-Re solutions using a multigrid method
    on a 129×129 uniform grid. We compare centerline velocity profiles.
    
    Physics:
        - Square cavity with moving lid (u = U_lid at top)
        - No-slip walls on other three boundaries
        - Re = U_lid × L / ν = 100
    
    Returns:
        ValidationResult with L2 error and pass/fail status
    """
    print("\n" + "="*70)
    print("LID-DRIVEN CAVITY VALIDATION: Re=100 vs Ghia et al. (1982)")
    print("="*70)
    
    if not HAS_PYCFDRS:
        print("ERROR: pycfdrs not available")
        print("Build with: maturin develop --manifest-path crates/pycfdrs/Cargo.toml")
        return ValidationResult(
            name="Ghia Cavity Re=100",
            l2_error=1.0,
            tolerance=0.01,
            passed=False,
            details="pycfdrs not available"
        )
    
    # Create solver with Ghia benchmark parameters
    # Using 65×65 grid (half of Ghia's 129×129 for faster testing)
    solver = CavitySolver2D(
        reynolds=100.0,
        nx=65,
        ny=65,
        lid_velocity=1.0,
        cavity_size=1.0
    )
    
    print(f"Solver: {solver}")
    print(f"Kinematic viscosity: {solver.viscosity():.6f} m²/s")
    print(f"Grid: {solver.nx}×{solver.ny}")
    
    print("\nRunning SIMPLEC solver...")
    result = solver.solve()
    
    print(f"\nResults:")
    print(f"  L2 Error (vs Ghia): {result.l2_error:.4f}")
    print(f"  Converged: {result.converged}")
    
    # Get Ghia benchmark data
    ghia_u = solver.ghia_u_centerline()
    ghia_v = solver.ghia_v_centerline()
    
    # Create comparison plot
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # U-velocity along vertical centerline
    ax = axes[0]
    ax.plot(ghia_u[:, 1], ghia_u[:, 0], 'ko', markersize=8, label='Ghia et al. (1982)')
    ax.plot(result.u_centerline, result.y_coords, 'b-', linewidth=2, label='pycfdrs (SIMPLEC)')
    ax.set_xlabel('U-velocity (normalized)')
    ax.set_ylabel('Y position')
    ax.set_title(f'U-velocity at x=0.5 (Re={solver.reynolds:.0f})')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # V-velocity along horizontal centerline
    ax = axes[1]
    ax.plot(ghia_v[:, 0], ghia_v[:, 1], 'ko', markersize=8, label='Ghia et al. (1982)')
    ax.plot(result.x_coords, result.v_centerline, 'b-', linewidth=2, label='pycfdrs (SIMPLEC)')
    ax.set_xlabel('X position')
    ax.set_ylabel('V-velocity (normalized)')
    ax.set_title(f'V-velocity at y=0.5 (Re={solver.reynolds:.0f})')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plot_file = 'validation/ghia_cavity_comparison.png'
    plt.savefig(plot_file, dpi=150)
    print(f"\nPlot saved: {plot_file}")
    plt.close()
    
    # Acceptance criteria: L2 error < 1%
    tolerance = 0.01
    passed = result.l2_error < tolerance
    
    print(f"\n{'='*70}")
    print(f"VALIDATION RESULT: {'✓ PASS' if passed else '✗ FAIL'}")
    print(f"  L2 Error: {result.l2_error:.4f}")
    print(f"  Tolerance: {tolerance:.2%}")
    print(f"{'='*70}")
    
    return ValidationResult(
        name="Ghia Cavity Re=100",
        l2_error=result.l2_error,
        tolerance=tolerance,
        passed=passed,
        details=f"Grid: {solver.nx}×{solver.ny}, converged={result.converged}"
    )


def validate_cavity_grid_convergence() -> ValidationResult:
    """
    Verify grid convergence for lid-driven cavity.
    
    As grid resolution increases, L2 error should decrease monotonically.
    This demonstrates numerical consistency.
    """
    print("\n" + "="*70)
    print("GRID CONVERGENCE STUDY: Lid-Driven Cavity Re=100")
    print("="*70)
    
    if not HAS_PYCFDRS:
        return ValidationResult(
            name="Grid Convergence",
            l2_error=1.0,
            tolerance=0.1,
            passed=False,
            details="pycfdrs not available"
        )
    
    grid_sizes = [17, 33, 65]
    errors = []
    
    for n in grid_sizes:
        solver = CavitySolver2D(reynolds=100.0, nx=n, ny=n)
        result = solver.solve()
        errors.append(result.l2_error)
        print(f"  Grid {n}×{n}: L2 error = {result.l2_error:.4f}")
    
    # Check monotonic decrease (allowing 10% tolerance for non-monotonicity)
    is_convergent = all(errors[i] >= errors[i+1] * 0.9 for i in range(len(errors)-1))
    
    print(f"\nGrid convergence: {'✓ Monotonic' if is_convergent else '✗ Non-monotonic'}")
    
    return ValidationResult(
        name="Grid Convergence",
        l2_error=errors[-1],
        tolerance=0.1,
        passed=is_convergent,
        details=f"Errors: {[f'{e:.4f}' for e in errors]}"
    )


def main():
    """Run complete cavity validation suite"""
    print("\n" + "╔" + "═"*68 + "╗")
    print("║" + " "*15 + "GHIA CAVITY BENCHMARK VALIDATION" + " "*20 + "║")
    print("║" + " "*10 + "pycfdrs vs Ghia et al. (1982)" + " "*28 + "║")
    print("╚" + "═"*68 + "╝")
    
    results = []
    
    # Main validation
    results.append(validate_cavity_re100())
    
    # Grid convergence
    results.append(validate_cavity_grid_convergence())
    
    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    
    all_passed = True
    for r in results:
        status = "✓ PASS" if r.passed else "✗ FAIL"
        all_passed = all_passed and r.passed
        print(f"  {r.name:30s}: {status} (L2 error: {r.l2_error:.4f})")
    
    print("="*70)
    if all_passed:
        print("✓ ALL VALIDATIONS PASSED")
        print("  pycfdrs cavity solver matches Ghia et al. (1982) benchmark.")
    else:
        print("✗ SOME VALIDATIONS FAILED")
        print("  Review implementation for failed tests.")
    print("="*70)
    
    return 0 if all_passed else 1


if __name__ == "__main__":
    sys.exit(main())
