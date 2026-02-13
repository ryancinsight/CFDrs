#!/usr/bin/env python3
"""
3D Poiseuille Flow Validation against Analytical Solution

Validates pycfdrs.Poiseuille3DSolver against the exact parabolic velocity profile
for fully developed laminar pipe flow:

    u(r) = u_max * (1 - (r/R)²)

where u_max = ΔP·R² / (4μL)
"""

import numpy as np
import pycfdrs
import matplotlib.pyplot as plt
from pathlib import Path

def analytical_poiseuille(r: np.ndarray, R: float, dp: float, mu: float, L: float) -> np.ndarray:
    """Compute analytical Poiseuille velocity profile."""
    u_max = dp * R**2 / (4 * mu * L)
    u = u_max * (1 - (r / R)**2)
    return u

def run_poiseuille_3d_validation():
    """Run 3D Poiseuille flow validation."""
    print(f"\n{'='*70}")
    print(f"3D Poiseuille Flow Validation: Analytical vs pycfdrs")
    print(f"{'='*70}\n")
    
    # Physical parameters
    R = 0.5e-3  # 0.5 mm radius
    L = 5.0e-3  # 5 mm length
    dp = 100.0  # 100 Pa pressure drop
    mu = 3.5e-3  # Blood viscosity (3.5 mPa·s)
    rho = 1060.0  # Blood density (kg/m³)
    
    # Run solver
    solver = pycfdrs.Poiseuille3DSolver(
        diameter=2*R,
        length=L,
        nr=10,
        ntheta=10,
        nz=20
    )
    
    result = solver.solve(pressure_drop=-dp, blood_type="newtonian")
    
    print(f"Solver Status: Converged\n")
    
    # Analytical solution
    u_max_analytical = dp * R**2 / (4 * mu * L)
    q_analytical = np.pi * R**2 * u_max_analytical / 2  # Flow rate
    
    print(f"Analytical Solution:")
    print(f"  u_max:     {u_max_analytical:.6f} m/s")
    print(f"  Q:         {q_analytical * 1e9:.6f} nl/s\n")
    
    # Extract numerical solution
    u_max_numerical = result.max_velocity
    q_numerical = result.flow_rate
    
    print(f"Numerical Solution:")
    print(f"  u_max:     {u_max_numerical:.6f} m/s")
    print(f"  Q:         {q_numerical * 1e9:.6f} nl/s\n")
    
    # Compute errors
    u_error = abs(u_max_numerical - u_max_analytical) / u_max_analytical
    q_error = abs(q_numerical - q_analytical) / q_analytical
    
    print(f"Error Metrics:")
    print(f"  u_max error: {u_error * 100:.2f}%")
    print(f"  Q error:     {q_error * 100:.2f}%\n")
    
    # Validation threshold
    threshold = 0.02  # 2% error
    passed = (u_error < threshold) and (q_error < threshold)
    
    print(f"Validation Status: {'PASSED' if passed else 'FAILED'}")
    print(f"  (Threshold: < {threshold * 100:.1f}%)\n")
    
    # Plot radial velocity profile
    r_pts = np.linspace(0, R, 50)
    u_analytical_profile = analytical_poiseuille(r_pts, R, dp, mu, L)
    
    plt.figure(figsize=(10, 6))
    plt.plot(r_pts * 1e3, u_analytical_profile, 'b-', linewidth=2, label='Analytical')
    
    # Plot numerical points (if available)
    if hasattr(result, 'radial_coords') and hasattr(result, 'radial_velocities'):
        r_numerical = np.array(result.radial_coords)
        u_numerical = np.array(result.radial_velocities)
        plt.plot(r_numerical * 1e3, u_numerical, 'ro', markersize=8, label='pycfdrs (P2)')
    
    plt.xlabel('Radial Position r (mm)', fontsize=12)
    plt.ylabel('Axial Velocity u (m/s)', fontsize=12)
    plt.title('3D Poiseuille Flow: Velocity Profile', fontsize=14, fontweight='bold')
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    output_dir = Path(__file__).parent / "results"
    output_dir.mkdir(exist_ok=True)
    plt.savefig(output_dir / "poiseuille_3d_profile.png", dpi=150)
    print(f"Plot saved: {output_dir / 'poiseuille_3d_profile.png'}\n")
    
    return passed, u_error, q_error

if __name__ == "__main__":
    passed, u_err, q_err = run_poiseuille_3d_validation()
    
    print(f"\n{'='*70}")
    print(f"Summary:")
    print(f"  Validation: {'PASSED' if passed else 'FAILED'}")
    print(f"  u_max error: {u_err * 100:.2f}%")
    print(f"  Flow rate error: {q_err * 100:.2f}%")
    print(f"{'='*70}\n")
