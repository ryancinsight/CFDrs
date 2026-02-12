#!/usr/bin/env python3
"""
Lid-Driven Cavity Flow Validation against Ghia et al. (1982)

Compares pycfdrs.CavitySolver2D against benchmark data from:
Ghia, U., Ghia, K. N., & Shin, C. T. (1982). High-Re solutions for
incompressible flow using the Navier-Stokes equations and a multigrid method.
Journal of computational physics, 48(3), 387-411.
"""

import numpy as np
import pycfdrs
import matplotlib.pyplot as plt
from pathlib import Path

def load_ghia_data(re: int) -> tuple:
    """Load Ghia (1982) benchmark data for given Reynolds number."""
    ghia_file = Path(__file__).parent.parent / "external" / "Python_CFD" / "Ghia-1982.txt"
    
    y_coords = []
    u_values = []
    
    with open(ghia_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 8:
                y = float(parts[0])
                re_col_map = {100: 1, 400: 2, 1000: 3, 3200: 4, 
                             5000: 5, 7500: 6, 10000: 7}
                if re in re_col_map:
                    u = float(parts[re_col_map[re]])
                    y_coords.append(y)
                    u_values.append(u)
    
    return np.array(y_coords), np.array(u_values)

def run_cavity_validation(re: int = 100, nx: int = 129, ny: int = 129):
    """Run cavity flow validation at specified Reynolds number."""
    print(f"\n{'='*70}")
    print(f"Lid-Driven Cavity Validation: Re = {re}, Grid = {nx}x{ny}")
    print(f"{'='*70}\n")
    
    # Run pycfdrs solver
    solver = pycfdrs.CavitySolver2D(
        nx=nx,
        ny=ny,
        re=float(re),
        max_iterations=100000,
        tolerance=1e-5
    )
    
    result = solver.solve()
    
    print(f"Solver converged: {result.converged}")
    print(f"Iterations: {result.iterations}")
    print(f"Final residual: {result.residual:.2e}\n")
    
    # Load Ghia benchmark
    y_ghia, u_ghia = load_ghia_data(re)
    
    # Extract u-velocity along vertical centerline
    u_centerline = np.array(result.u_centerline)
    y_centerline = np.array(result.y_coords)
    
    # Interpolate pycfdrs results to Ghia points
    u_interp = np.interp(y_ghia, y_centerline, u_centerline)
    
    # Compute errors
    abs_error = np.abs(u_interp - u_ghia)
    rel_error = abs_error / (np.abs(u_ghia) + 1e-10)
    l2_error = np.sqrt(np.mean((u_interp - u_ghia)**2))
    max_error = np.max(abs_error)
    
    print(f"Error Metrics:")
    print(f"  L2 Error:  {l2_error:.6f}")
    print(f"  Max Error: {max_error:.6f}")
    print(f"  Mean Rel Error: {np.mean(rel_error) * 100:.2f}%\n")
    
    # Validation threshold
    threshold = 0.05  # 5% L2 error
    passed = l2_error < threshold
    
    print(f"Validation Status: {'✓ PASSED' if passed else '✗ FAILED'}")
    print(f"  (Threshold: L2 < {threshold:.2%})\n")
    
    # Plot comparison
    plt.figure(figsize=(10, 6))
    plt.plot(u_centerline, y_centerline, 'b-', linewidth=2, label='pycfdrs')
    plt.plot(u_ghia, y_ghia, 'ro', markersize=8, label='Ghia et al. (1982)')
    plt.xlabel('u-velocity', fontsize=12)
    plt.ylabel('y', fontsize=12)
    plt.title(f'Centerline Velocity Profile: Re = {re}', fontsize=14, fontweight='bold')
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    output_dir = Path(__file__).parent / "results"
    output_dir.mkdir(exist_ok=True)
    plt.savefig(output_dir / f"cavity_re{re}.png", dpi=150)
    print(f"Plot saved: {output_dir / f'cavity_re{re}.png'}\n")
    
    return passed, l2_error

if __name__ == "__main__":
    # Run validation for Re = 100 (coarse grid acceptable)
    passed_100, error_100 = run_cavity_validation(re=100, nx=65, ny=65)
    
    # Run for Re = 400 (requires finer grid)
    # passed_400, error_400 = run_cavity_validation(re=400, nx=129, ny=129)
    
    print(f"\n{'='*70}")
    print(f"Summary:")
    print(f"  Re=100:  {'✓' if passed_100 else '✗'}  (L2 = {error_100:.4f})")
    # print(f"  Re=400:  {'✓' if passed_400 else '✗'}  (L2 = {error_400:.4f})")
    print(f"{'='*70}\n")
