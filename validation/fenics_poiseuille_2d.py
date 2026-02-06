#!/usr/bin/env python3
"""
FEniCS validation of 2D Poiseuille flow with non-Newtonian blood

Solves the same problem in FEniCS and compares with pycfdrs to prove correctness.

This implements:
- Steady-state Navier-Stokes with non-Newtonian viscosity
- Iterative solution for shear-dependent viscosity (Casson model)
- Direct comparison of velocity, shear rate, and viscosity profiles
"""

import sys

import matplotlib.pyplot as plt
import numpy as np

# Check for FEniCS
try:
    from dolfin import *

    parameters["linear_algebra_backend"] = "PETSc"
    set_log_level(LogLevel.ERROR)
except ImportError:
    print("ERROR: FEniCS not installed. Install with:")
    print("  conda install -c conda-forge fenics")
    print("  or")
    print("  pip install fenics-dolfinx")
    sys.exit(1)

# Check for pycfdrs
try:
    import pycfdrs
except ImportError:
    print("ERROR: pycfdrs not installed")
    sys.exit(1)


class CassonViscosity:
    """Casson blood model for FEniCS

    μ(γ̇) = (√τ_y / √γ̇ + √μ_∞)²
    """

    def __init__(self, tau_y=0.0056, mu_inf=0.00345):
        self.tau_y = tau_y  # Yield stress [Pa]
        self.mu_inf = mu_inf  # Infinite shear viscosity [Pa·s]

    def viscosity(self, shear_rate):
        """Compute apparent viscosity

        Args:
            shear_rate: Shear rate magnitude [1/s]

        Returns:
            Apparent viscosity [Pa·s]
        """
        # Add small epsilon to avoid division by zero
        gamma = shear_rate + 1e-10
        sqrt_tau_y = np.sqrt(self.tau_y)
        sqrt_mu_inf = np.sqrt(self.mu_inf)
        sqrt_gamma = np.sqrt(gamma)

        mu = (sqrt_tau_y / sqrt_gamma + sqrt_mu_inf) ** 2
        return mu


def solve_fenics_poiseuille(H, ny, dP_dx, tolerance=1e-8, max_iter=1000):
    """
    Solve 2D Poiseuille flow with Casson blood in FEniCS

    Args:
        H: Channel height [m]
        ny: Number of grid points in y-direction
        dP_dx: Pressure gradient [Pa/m]
        tolerance: Convergence tolerance
        max_iter: Maximum iterations

    Returns:
        y_coords, velocity, shear_rate, viscosity, iterations
    """

    print("\n" + "=" * 80)
    print("FEniCS Solution")
    print("=" * 80)

    # Create 1D mesh for fully-developed flow (variation only in y)
    mesh = IntervalMesh(ny - 1, 0, H)

    # Function space for velocity
    V = FunctionSpace(mesh, "CG", 2)  # Continuous Galerkin order 2

    # Boundary conditions (no-slip at walls)
    def boundary(x, on_boundary):
        return on_boundary

    bc = DirichletBC(V, Constant(0.0), boundary)

    # Define variational problem
    u = TrialFunction(V)
    v = TestFunction(V)

    # Initialize viscosity field
    casson = CassonViscosity()
    mu_field = Function(V)
    mu_field.vector()[:] = casson.mu_inf  # Start with constant viscosity

    # Source term from pressure gradient
    f = Constant(dP_dx)

    # Iterative solution for non-Newtonian viscosity
    u_sol = Function(V)
    u_sol.vector()[:] = 0.0

    print(f"Mesh: {ny - 1} elements")
    print(f"DOFs: {V.dim()}")
    print("Iterating for non-Newtonian viscosity...")

    for iteration in range(max_iter):
        # Store old solution
        u_old = u_sol.copy(deepcopy=True)
        mu_old = mu_field.copy(deepcopy=True)

        # Solve linear problem with current viscosity
        # d/dy(μ du/dy) = dP/dx
        a = mu_field * inner(grad(u), grad(v)) * dx
        L = f * v * dx

        solve(a == L, u_sol, bc)

        # Update shear rate: γ̇ = |du/dy|
        # For 1D problem, shear rate is just the gradient magnitude
        u_grad = project(grad(u_sol), VectorFunctionSpace(mesh, "CG", 1))
        shear_rate_vals = np.abs(u_grad.vector().get_local())

        # Update viscosity based on shear rate
        mu_new = np.array([casson.viscosity(gamma) for gamma in shear_rate_vals])
        mu_field.vector()[:] = mu_new

        # Check convergence
        u_diff = u_sol.vector() - u_old.vector()
        mu_diff = mu_field.vector() - mu_old.vector()

        u_norm = np.linalg.norm(u_diff.get_local())
        mu_norm = np.linalg.norm(mu_diff.get_local())
        mu_old_norm = np.linalg.norm(mu_old.vector().get_local())

        residual = mu_norm / (mu_old_norm + 1e-20)

        if iteration % 10 == 0:
            print(f"  Iteration {iteration}: residual = {residual:.6e}")

        if residual < tolerance:
            print(f"Converged in {iteration + 1} iterations (residual={residual:.6e})")
            break
    else:
        print(f"Warning: Did not converge after {max_iter} iterations")
        iteration = max_iter - 1

    # Extract solution at vertices
    y_coords = mesh.coordinates()[:, 0]
    sort_idx = np.argsort(y_coords)
    y_coords = y_coords[sort_idx]

    # Get velocity values
    velocity = u_sol.vector().get_local()[sort_idx]

    # Compute shear rate and viscosity at vertices
    u_grad = project(grad(u_sol), VectorFunctionSpace(mesh, "CG", 1))
    shear_rate = np.abs(u_grad.vector().get_local()[sort_idx])
    viscosity = np.array([casson.viscosity(gamma) for gamma in shear_rate])

    # Compute flow rate (per unit width) using trapezoidal rule
    flow_rate_per_width = np.trapz(velocity, y_coords)

    # Wall shear stress
    wall_shear_stress = viscosity[0] * shear_rate[0]

    print(f"\nFEniCS Results:")
    print(f"  Flow rate (per unit width): {flow_rate_per_width:.6e} m²/s")
    print(f"  Wall shear stress: {wall_shear_stress:.3e} Pa")
    print(f"  Max velocity: {velocity.max():.6e} m/s")
    print(f"  Max shear rate: {shear_rate.max():.3e} s⁻¹")
    print(f"  Min viscosity: {viscosity.min():.6e} Pa·s")

    return y_coords, velocity, shear_rate, viscosity, iteration + 1


def compare_with_pycfdrs():
    """Compare FEniCS solution with pycfdrs"""

    print("\n" + "=" * 80)
    print("FEniCS vs pycfdrs Validation for 2D Poiseuille Flow")
    print("=" * 80)

    # Configuration
    height = 0.001  # 1 mm
    width = 0.01  # 10 mm
    ny = 101
    pressure_gradient = 100000.0  # Pa/m

    print(f"\nTest Configuration:")
    print(f"  Channel height: {height * 1e3:.2f} mm")
    print(f"  Grid points: {ny}")
    print(f"  Pressure gradient: {pressure_gradient:.2e} Pa/m")

    # Solve with FEniCS
    y_fenics, u_fenics, gamma_fenics, mu_fenics, iter_fenics = solve_fenics_poiseuille(
        height, ny, pressure_gradient
    )

    # Solve with pycfdrs
    print("\n" + "=" * 80)
    print("pycfdrs Solution")
    print("=" * 80)

    config = pycfdrs.PoiseuilleConfig2D(
        height=height,
        width=width,
        length=0.05,
        ny=ny,
        pressure_gradient=pressure_gradient,
        tolerance=1e-8,
        max_iterations=1000,
        relaxation_factor=0.5,
    )

    solver = pycfdrs.PoiseuilleSolver2D(config)
    blood = pycfdrs.CassonBlood()
    result = solver.solve(blood)

    y_pycfdrs = np.array(result.y_coords)
    u_pycfdrs = np.array(result.velocity)
    gamma_pycfdrs = np.array(result.shear_rate)
    mu_pycfdrs = np.array(result.viscosity)

    print(f"Converged in {result.iterations} iterations")
    print(f"\npycfdrs Results:")
    print(f"  Flow rate (per unit width): {result.flow_rate / width:.6e} m²/s")
    print(f"  Wall shear stress: {result.wall_shear_stress:.3e} Pa")
    print(f"  Max velocity: {u_pycfdrs.max():.6e} m/s")
    print(f"  Max shear rate: {gamma_pycfdrs.max():.3e} s⁻¹")
    print(f"  Min viscosity: {mu_pycfdrs.min():.6e} Pa·s")

    # Interpolate FEniCS solution to pycfdrs grid for comparison
    u_fenics_interp = np.interp(y_pycfdrs, y_fenics, u_fenics)
    gamma_fenics_interp = np.interp(y_pycfdrs, y_fenics, gamma_fenics)
    mu_fenics_interp = np.interp(y_pycfdrs, y_fenics, mu_fenics)

    # Compute errors (skip boundaries)
    interior = slice(1, -1)

    u_error = np.abs(u_pycfdrs[interior] - u_fenics_interp[interior]) / (
        np.abs(u_fenics_interp[interior]) + 1e-10
    )
    gamma_error = np.abs(gamma_pycfdrs[interior] - gamma_fenics_interp[interior]) / (
        np.abs(gamma_fenics_interp[interior]) + 1e-10
    )
    mu_error = np.abs(mu_pycfdrs[interior] - mu_fenics_interp[interior]) / (
        np.abs(mu_fenics_interp[interior]) + 1e-10
    )

    print("\n" + "=" * 80)
    print("Comparison: pycfdrs vs FEniCS")
    print("=" * 80)
    print(f"Velocity:")
    print(f"  Max relative error: {u_error.max() * 100:.3f}%")
    print(f"  Mean relative error: {u_error.mean() * 100:.3f}%")
    print(f"  RMS error: {np.sqrt(np.mean(u_error**2)) * 100:.3f}%")

    print(f"\nShear rate:")
    print(f"  Max relative error: {gamma_error.max() * 100:.3f}%")
    print(f"  Mean relative error: {gamma_error.mean() * 100:.3f}%")

    print(f"\nViscosity:")
    print(f"  Max relative error: {mu_error.max() * 100:.3f}%")
    print(f"  Mean relative error: {mu_error.mean() * 100:.3f}%")

    # Validation criteria
    print("\n" + "=" * 80)
    print("Validation Checks")
    print("=" * 80)

    checks_passed = 0
    checks_total = 0

    # Check 1: Velocity error < 5%
    checks_total += 1
    if u_error.max() < 0.05:
        print(f"✓ Velocity error acceptable ({u_error.max() * 100:.3f}% < 5%)")
        checks_passed += 1
    else:
        print(f"✗ Velocity error too large ({u_error.max() * 100:.3f}% >= 5%)")

    # Check 2: Both solvers converged
    checks_total += 1
    if result.iterations < 1000 and iter_fenics < 1000:
        print(
            f"✓ Both solvers converged (pycfdrs: {result.iterations}, FEniCS: {iter_fenics})"
        )
        checks_passed += 1
    else:
        print(f"✗ Convergence issues")

    # Check 3: Flow rates match within 5%
    checks_total += 1
    q_fenics = np.trapz(u_fenics, y_fenics)
    q_pycfdrs = result.flow_rate / width
    q_error = abs(q_fenics - q_pycfdrs) / q_fenics
    if q_error < 0.05:
        print(f"✓ Flow rates match ({q_error * 100:.3f}% < 5%)")
        checks_passed += 1
    else:
        print(f"✗ Flow rates differ ({q_error * 100:.3f}% >= 5%)")

    # Check 4: Shear-thinning behavior in both
    checks_total += 1
    if mu_pycfdrs.max() > mu_pycfdrs.min() and mu_fenics.max() > mu_fenics.min():
        print(f"✓ Both show shear-thinning behavior")
        checks_passed += 1
    else:
        print(f"✗ Shear-thinning not consistent")

    print(f"\n{'=' * 80}")
    print(f"VALIDATION RESULT: {checks_passed}/{checks_total} checks passed")
    print(f"{'=' * 80}")

    # Plot comparison
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))

    # Velocity comparison
    ax1.plot(y_fenics * 1e3, u_fenics, "b-", linewidth=2, label="FEniCS")
    ax1.plot(y_pycfdrs * 1e3, u_pycfdrs, "r--", linewidth=1.5, label="pycfdrs")
    ax1.set_xlabel("y [mm]")
    ax1.set_ylabel("Velocity [m/s]")
    ax1.set_title("Velocity Profile Comparison")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Velocity error
    ax2.plot(y_pycfdrs[interior] * 1e3, u_error * 100, "k-", linewidth=2)
    ax2.set_xlabel("y [mm]")
    ax2.set_ylabel("Relative Error [%]")
    ax2.set_title(f"Velocity Error (max={u_error.max() * 100:.3f}%)")
    ax2.grid(True, alpha=0.3)

    # Viscosity comparison
    ax3.plot(y_fenics * 1e3, mu_fenics * 1e3, "b-", linewidth=2, label="FEniCS")
    ax3.plot(y_pycfdrs * 1e3, mu_pycfdrs * 1e3, "r--", linewidth=1.5, label="pycfdrs")
    ax3.set_xlabel("y [mm]")
    ax3.set_ylabel("Viscosity [cP]")
    ax3.set_title("Viscosity Profile Comparison")
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    ax3.set_yscale("log")

    # Shear rate comparison
    ax4.plot(y_fenics * 1e3, gamma_fenics, "b-", linewidth=2, label="FEniCS")
    ax4.plot(y_pycfdrs * 1e3, gamma_pycfdrs, "r--", linewidth=1.5, label="pycfdrs")
    ax4.set_xlabel("y [mm]")
    ax4.set_ylabel("Shear Rate [s⁻¹]")
    ax4.set_title("Shear Rate Comparison")
    ax4.legend()
    ax4.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("fenics_pycfdrs_comparison.png", dpi=150, bbox_inches="tight")
    print(f"\nComparison plot saved: fenics_pycfdrs_comparison.png")

    return checks_passed == checks_total


if __name__ == "__main__":
    success = compare_with_pycfdrs()
    sys.exit(0 if success else 1)
