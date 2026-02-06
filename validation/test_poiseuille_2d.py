#!/usr/bin/env python3
"""
Test 2D Poiseuille flow solver with non-Newtonian blood

This script validates the 2D Poiseuille solver by:
1. Running the solver with Casson blood model
2. Comparing with analytical Newtonian solution
3. Verifying convergence and physical constraints
"""

import sys

import matplotlib.pyplot as plt
import numpy as np

# Try to import pycfdrs
try:
    import pycfdrs
except ImportError:
    print("ERROR: pycfdrs not installed. Install with:")
    print("  maturin develop --release")
    print("  or")
    print("  pip install target/wheels/pycfdrs-*.whl")
    sys.exit(1)


def analytical_poiseuille(y, H, mu, dP_dx):
    """
    Analytical solution for Newtonian Poiseuille flow

    u(y) = -(dP/dx) * y * (H - y) / (2 * μ)

    Args:
        y: y-coordinates [m]
        H: channel height [m]
        mu: dynamic viscosity [Pa·s]
        dP_dx: pressure gradient [Pa/m]

    Returns:
        velocity profile [m/s]
    """
    return (dP_dx) * y * (H - y) / (2.0 * mu)


def test_poiseuille_2d():
    """Test 2D Poiseuille solver with Casson blood"""

    print("=" * 80)
    print("2D Poiseuille Flow Validation with Non-Newtonian Blood")
    print("=" * 80)

    # Configuration
    height = 0.001  # 1 mm channel
    width = 0.01  # 10 mm width
    length = 0.05  # 50 mm length
    ny = 101  # Grid points
    pressure_gradient = 100000.0  # Pa/m (high shear)

    print(f"\nConfiguration:")
    print(f"  Channel height: {height * 1e3:.2f} mm")
    print(f"  Channel width: {width * 1e3:.2f} mm")
    print(f"  Grid points (y): {ny}")
    print(f"  Pressure gradient: {pressure_gradient:.2e} Pa/m")

    # Create configuration
    config = pycfdrs.PoiseuilleConfig2D(
        height=height,
        width=width,
        length=length,
        ny=ny,
        pressure_gradient=pressure_gradient,
        tolerance=1e-8,
        max_iterations=1000,
        relaxation_factor=0.5,
    )

    print(f"\nSolver config: {config}")

    # Create solver
    solver = pycfdrs.PoiseuilleSolver2D(config)
    print(f"Solver: {solver}")

    # Create Casson blood model
    blood = pycfdrs.CassonBlood()
    print(f"\nBlood model: {blood}")

    # Solve
    print("\nSolving...")
    result = solver.solve(blood)
    print(f"Result: {result}")

    # Extract results
    y_coords = np.array(result.y_coords)
    velocity = np.array(result.velocity)
    shear_rate = np.array(result.shear_rate)
    viscosity = np.array(result.viscosity)

    print(f"\nResults:")
    print(f"  Converged in {result.iterations} iterations")
    print(f"  Flow rate: {result.flow_rate:.6e} m³/s")
    print(f"  Wall shear stress: {result.wall_shear_stress:.3e} Pa")
    print(f"  Max velocity: {velocity.max():.6e} m/s")
    print(f"  Max shear rate: {shear_rate.max():.3e} s⁻¹")
    print(f"  Min viscosity: {viscosity.min():.6e} Pa·s")
    print(f"  Average viscosity: {viscosity.mean():.6e} Pa·s")

    # Compare with analytical solution using minimum viscosity (high-shear)
    mu_eff = viscosity.min()
    u_analytical = analytical_poiseuille(y_coords, height, mu_eff, pressure_gradient)

    # Compute errors (skip boundaries)
    interior_indices = np.arange(1, len(velocity) - 1)
    relative_errors = np.abs(
        (velocity[interior_indices] - u_analytical[interior_indices])
        / u_analytical[interior_indices]
    )
    max_error = relative_errors.max()
    mean_error = relative_errors.mean()

    print(f"\nComparison with Newtonian analytical (μ={mu_eff:.6e} Pa·s):")
    print(f"  Max relative error: {max_error * 100:.2f}%")
    print(f"  Mean relative error: {mean_error * 100:.2f}%")

    # Validation checks
    print("\nValidation checks:")
    checks_passed = 0
    checks_total = 0

    # Check 1: Convergence
    checks_total += 1
    if result.iterations < 1000:
        print(f"  ✓ Solver converged in {result.iterations} iterations")
        checks_passed += 1
    else:
        print(f"  ✗ Solver did not converge (reached max iterations)")

    # Check 2: Positive velocity
    checks_total += 1
    if velocity.min() >= 0:
        print(f"  ✓ All velocities non-negative")
        checks_passed += 1
    else:
        print(f"  ✗ Negative velocities found (min={velocity.min():.3e})")

    # Check 3: Boundary conditions
    checks_total += 1
    if abs(velocity[0]) < 1e-10 and abs(velocity[-1]) < 1e-10:
        print(f"  ✓ No-slip boundary conditions satisfied")
        checks_passed += 1
    else:
        print(
            f"  ✗ Boundary conditions violated (u[0]={velocity[0]:.3e}, u[-1]={velocity[-1]:.3e})"
        )

    # Check 4: Reasonable error
    checks_total += 1
    if max_error < 0.15:  # 15% tolerance for non-Newtonian vs Newtonian
        print(f"  ✓ Error within tolerance ({max_error * 100:.2f}% < 15%)")
        checks_passed += 1
    else:
        print(f"  ✗ Error too large ({max_error * 100:.2f}% >= 15%)")

    # Check 5: Flow rate is positive
    checks_total += 1
    if result.flow_rate > 0:
        print(f"  ✓ Flow rate positive")
        checks_passed += 1
    else:
        print(f"  ✗ Flow rate non-positive ({result.flow_rate:.3e})")

    # Check 6: Shear thinning behavior
    checks_total += 1
    if viscosity.max() > viscosity.min():
        print(
            f"  ✓ Shear-thinning behavior observed (μ_max={viscosity.max():.3e}, μ_min={viscosity.min():.3e})"
        )
        checks_passed += 1
    else:
        print(f"  ✗ No shear-thinning observed")

    print(f"\nValidation summary: {checks_passed}/{checks_total} checks passed")

    # Plot results
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))

    # Plot 1: Velocity profile
    ax1.plot(
        y_coords * 1e3, velocity, "b-", linewidth=2, label="Numerical (Casson blood)"
    )
    ax1.plot(
        y_coords * 1e3,
        u_analytical,
        "r--",
        linewidth=1.5,
        label=f"Analytical (μ={mu_eff:.3e} Pa·s)",
    )
    ax1.set_xlabel("y [mm]")
    ax1.set_ylabel("u [m/s]")
    ax1.set_title("Velocity Profile")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Plot 2: Shear rate profile
    ax2.plot(y_coords * 1e3, shear_rate, "g-", linewidth=2)
    ax2.set_xlabel("y [mm]")
    ax2.set_ylabel("Shear rate [s⁻¹]")
    ax2.set_title("Shear Rate Profile")
    ax2.grid(True, alpha=0.3)

    # Plot 3: Viscosity profile
    ax3.plot(y_coords * 1e3, viscosity * 1e3, "m-", linewidth=2)
    ax3.set_xlabel("y [mm]")
    ax3.set_ylabel("Viscosity [cP]")
    ax3.set_title("Viscosity Profile (Non-Newtonian)")
    ax3.grid(True, alpha=0.3)

    # Plot 4: Relative error
    ax4.plot(y_coords[interior_indices] * 1e3, relative_errors * 100, "k-", linewidth=2)
    ax4.set_xlabel("y [mm]")
    ax4.set_ylabel("Relative Error [%]")
    ax4.set_title("Error vs Analytical Solution")
    ax4.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("poiseuille_2d_validation.png", dpi=150, bbox_inches="tight")
    print(f"\nPlot saved to: poiseuille_2d_validation.png")

    # Return success if all checks passed
    return checks_passed == checks_total


if __name__ == "__main__":
    success = test_poiseuille_2d()
    sys.exit(0 if success else 1)
