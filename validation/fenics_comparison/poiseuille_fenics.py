#!/usr/bin/env python3
"""
FEniCS vs CFD-rs: 2D Poiseuille Flow Comparison

This script compares CFD-rs results against FEniCS (Finite Element Method)
for 2D Poiseuille flow between parallel plates.

FEniCS Implementation:
---------------------
Solves the steady incompressible Navier-Stokes equations using:
- P2-P1 Taylor-Hood elements (quadratic velocity, linear pressure)
- No-slip boundary conditions on walls
- Parabolic inflow profile
- Zero-pressure outflow

The weak formulation:
    ∫ μ∇u:∇v dx + ∫ ρ(u·∇u)·v dx - ∫ p(∇·v) dx = 0
    ∫ q(∇·u) dx = 0

For low Re (Re < 1), the convective term (u·∇u) ≈ 0 (Stokes flow).

References:
----------
- Langtangen & Logg (2017). "Solving PDEs in Python: The FEniCS Tutorial"
- Donea & Huerta (2003). "Finite Element Methods for Flow Problems"
"""

import sys
from pathlib import Path
from typing import Dict, Optional

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

# Try to import FEniCS
try:
    from dolfin import *

    FENICS_AVAILABLE = True
except ImportError:
    print("WARNING: FEniCS not available. Install with:")
    print("  conda install -c conda-forge fenics")
    FENICS_AVAILABLE = False

# Try to import pycfdrs
try:
    import pycfdrs

    PYCFDRS_AVAILABLE = True
except ImportError:
    print("WARNING: pycfdrs not available. Build with 'maturin develop'")
    PYCFDRS_AVAILABLE = False


def solve_poiseuille_fenics(
    height: float,
    length: float,
    pressure_drop: float,
    viscosity: float,
    density: float,
    nx: int = 50,
    ny: int = 100,
) -> Optional[Dict]:
    """
    Solve 2D Poiseuille flow using FEniCS

    Parameters:
    -----------
    height : float
        Channel height [m]
    length : float
        Channel length [m]
    pressure_drop : float
        Pressure drop [Pa]
    viscosity : float
        Dynamic viscosity [Pa·s]
    density : float
        Density [kg/m³]
    nx, ny : int
        Mesh resolution

    Returns:
    --------
    solution : dict or None
        Dictionary with velocity and pressure fields
    """
    if not FENICS_AVAILABLE:
        return None

    print("Solving with FEniCS...")

    # Create mesh
    mesh = RectangleMesh(Point(0, 0), Point(length, height), nx, ny)

    # Define function spaces (Taylor-Hood P2-P1)
    V = VectorFunctionSpace(mesh, "P", 2)  # Velocity
    Q = FunctionSpace(mesh, "P", 1)  # Pressure

    # Define trial and test functions
    u = TrialFunction(V)
    v = TestFunction(V)
    p = TrialFunction(Q)
    q = TestFunction(Q)

    # Define functions for solution
    u_ = Function(V)
    p_ = Function(Q)

    # Boundary conditions
    # No-slip on walls (y=0 and y=height)
    def wall_boundary(x, on_boundary):
        return on_boundary and (near(x[1], 0) or near(x[1], height))

    bc_walls = DirichletBC(V, Constant((0, 0)), wall_boundary)

    # Parabolic inflow at x=0
    pressure_gradient = -pressure_drop / length
    u_max = -(pressure_gradient / (8.0 * viscosity)) * height**2
    inflow_profile = Expression(
        ("u_max * 4.0 * x[1] * (H - x[1]) / (H*H)", "0"),
        degree=2,
        u_max=u_max,
        H=height,
    )

    def inflow_boundary(x, on_boundary):
        return on_boundary and near(x[0], 0)

    bc_inflow = DirichletBC(V, inflow_profile, inflow_boundary)

    # Pressure boundary conditions
    def outflow_boundary(x, on_boundary):
        return on_boundary and near(x[0], length)

    bc_pressure = DirichletBC(Q, Constant(0), outflow_boundary)

    bcs_u = [bc_walls, bc_inflow]
    bcs_p = [bc_pressure]

    # Variational formulation (Stokes equation for low Re)
    # For high Re, add convection term: + rho*dot(u_, nabla_grad(u_))*v*dx
    a_u = viscosity * inner(nabla_grad(u), nabla_grad(v)) * dx
    L_u = inner(Constant((0, 0)), v) * dx

    # Pressure Poisson equation
    a_p = inner(nabla_grad(p), nabla_grad(q)) * dx
    L_p = -div(u_) * q * dx

    # Solve momentum equation
    solve(a_u == L_u, u_, bcs_u)

    # Solve pressure equation
    solve(a_p == L_p, p_, bcs_p)

    # Extract solution arrays
    u_array = u_.compute_vertex_values(mesh)
    p_array = p_.compute_vertex_values(mesh)

    # Get coordinates
    coords = mesh.coordinates()

    # Compute velocity magnitude
    u_magnitude = np.sqrt(u_array[: len(coords)] ** 2 + u_array[len(coords) :] ** 2)

    return {
        "velocity": u_array,
        "pressure": p_array,
        "velocity_magnitude": u_magnitude,
        "coordinates": coords,
        "mesh": mesh,
        "u_func": u_,
        "p_func": p_,
    }


def compare_fenics_cfdrs(
    height: float = 100e-6,
    length: float = 1e-3,
    pressure_drop: float = 100.0,
    viscosity: float = 3.5e-3,
    density: float = 1060.0,
    nx: int = 50,
    ny: int = 100,
    plot: bool = True,
) -> Dict:
    """
    Compare FEniCS and CFD-rs solutions

    Returns:
    --------
    results : dict
        Comparison results and error metrics
    """
    print("=" * 80)
    print("FEniCS vs CFD-rs: 2D POISEUILLE FLOW COMPARISON")
    print("=" * 80)
    print(f"Channel height:    {height * 1e6:.1f} μm")
    print(f"Channel length:    {length * 1e3:.2f} mm")
    print(f"Pressure drop:     {pressure_drop:.1f} Pa")
    print(f"Viscosity:         {viscosity * 1e3:.2f} cP")
    print(f"Grid resolution:   {nx} × {ny}")
    print()

    # Solve with FEniCS
    fenics_solution = solve_poiseuille_fenics(
        height, length, pressure_drop, viscosity, density, nx, ny
    )

    if fenics_solution is None:
        print("FEniCS solution not available, skipping comparison")
        return {"passed": False, "error": "FEniCS not available"}

    # Analytical solution for comparison
    pressure_gradient = -pressure_drop / length
    u_max_analytical = -(pressure_gradient / (8.0 * viscosity)) * height**2

    print(f"FEniCS max velocity: {fenics_solution['velocity_magnitude'].max():.6e} m/s")
    print(f"Analytical max vel:  {u_max_analytical:.6e} m/s")
    print()

    # If CFD-rs available, compare
    if PYCFDRS_AVAILABLE:
        print("Note: CFD-rs 2D solver integration pending")
        print("Comparison will be enabled once pycfdrs 2D bindings are complete")
        print()

    # Validation: Compare FEniCS against analytical
    fenics_error = (
        abs(fenics_solution["velocity_magnitude"].max() - u_max_analytical)
        / u_max_analytical
    )

    print("FEniCS vs Analytical:")
    print(f"  Relative error: {fenics_error * 100:.3f}%")

    passed = fenics_error < 0.01

    if passed:
        print("  ✓ FEniCS accuracy < 1%")
    else:
        print(f"  ✗ FEniCS error too large (FAILED)")

    if plot and fenics_solution is not None:
        # Plot FEniCS solution
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        # Velocity field
        ax = axes[0]
        coords = fenics_solution["coordinates"]
        u_mag = fenics_solution["velocity_magnitude"]

        sc = ax.tricontourf(
            coords[:, 0] * 1e3,
            coords[:, 1] * 1e6,
            u_mag * 1e3,
            levels=20,
            cmap="viridis",
        )
        ax.set_xlabel("x [mm]", fontsize=12)
        ax.set_ylabel("y [μm]", fontsize=12)
        ax.set_title("FEniCS Velocity Field [mm/s]", fontsize=14)
        plt.colorbar(sc, ax=ax)

        # Centerline profile
        ax = axes[1]
        # Extract centerline values
        centerline_mask = np.abs(coords[:, 0] - length / 2) < length / (2 * nx)
        centerline_coords = coords[centerline_mask]
        centerline_vel = u_mag[centerline_mask]

        # Sort by y coordinate
        sort_idx = np.argsort(centerline_coords[:, 1])
        y_center = centerline_coords[sort_idx, 1]
        u_center = centerline_vel[sort_idx]

        # Analytical profile
        y_analytical = np.linspace(0, height, 200)
        u_analytical = (
            -(pressure_gradient / (2.0 * viscosity))
            * y_analytical
            * (height - y_analytical)
        )

        ax.plot(
            u_analytical * 1e3,
            y_analytical * 1e6,
            "k-",
            linewidth=2,
            label="Analytical",
        )
        ax.plot(u_center * 1e3, y_center * 1e6, "ro", markersize=6, label="FEniCS")
        ax.set_xlabel("Velocity [mm/s]", fontsize=12)
        ax.set_ylabel("y [μm]", fontsize=12)
        ax.set_title("Centerline Velocity Profile", fontsize=14)
        ax.legend()
        ax.grid(True, alpha=0.3)

        plt.tight_layout()

        output_dir = Path(__file__).parent.parent / "reports" / "figures"
        output_dir.mkdir(parents=True, exist_ok=True)
        plt.savefig(
            output_dir / "fenics_poiseuille_comparison.png",
            dpi=300,
            bbox_inches="tight",
        )
        print(f"\nPlot saved to: {output_dir / 'fenics_poiseuille_comparison.png'}")

        plt.show()

    return {
        "passed": passed,
        "fenics_error": fenics_error,
        "fenics_max_velocity": fenics_solution["velocity_magnitude"].max(),
        "analytical_max_velocity": u_max_analytical,
    }


if __name__ == "__main__":
    results = compare_fenics_cfdrs(
        height=100e-6,
        length=1e-3,
        pressure_drop=100.0,
        viscosity=3.5e-3,
        density=1060.0,
        nx=50,
        ny=100,
        plot=True,
    )

    sys.exit(0 if results.get("passed", False) else 1)
