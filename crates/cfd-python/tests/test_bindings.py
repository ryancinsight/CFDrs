"""Value-semantic smoke tests for the installed cfd-python wheel."""

import importlib.metadata
import math

import cfd_python


def test_module_version_matches_wheel_metadata() -> None:
    """The Python module reports the version encoded in its wheel."""
    assert cfd_python.__version__ == importlib.metadata.version("cfd-python")


def test_casson_viscosity_matches_its_constitutive_equation() -> None:
    """The binding preserves the input-sensitive Casson constitutive law."""
    blood = cfd_python.CassonBlood()
    shear_rate = 25.0
    expected = (
        math.sqrt(blood.yield_stress() / shear_rate)
        + math.sqrt(blood.viscosity_high_shear())
    ) ** 2

    assert math.isclose(
        blood.apparent_viscosity(shear_rate),
        expected,
        rel_tol=2.0e-15,
        abs_tol=0.0,
    )
    assert blood.apparent_viscosity(1.0) > blood.apparent_viscosity(100.0)


def test_womersley_number_matches_dimensionless_definition() -> None:
    """The binding computes alpha = R * sqrt(omega * rho / mu)."""
    radius = 5.0e-3
    omega = 2.0 * math.pi
    density = 1060.0
    viscosity = 4.0e-3
    womersley = cfd_python.WomersleyNumber(radius, omega, density, viscosity)
    expected = radius * math.sqrt(omega * density / viscosity)

    assert math.isclose(womersley.value(), expected, rel_tol=2.0e-15, abs_tol=0.0)
    assert womersley.flow_regime()


def test_poiseuille_solver_returns_analytical_profile_values() -> None:
    """The installed binding returns the analytical channel-flow values."""
    height = 2.0e-3
    width = 1.0e-2
    length = 5.0e-2
    pressure_drop = 100.0
    viscosity = 3.5e-3
    solver = cfd_python.Poiseuille2DSolver(height, width, length, 9, 5)
    result = solver.solve(pressure_drop, "newtonian")
    expected_max = height**2 * (pressure_drop / length) / (8.0 * viscosity)
    expected_flow_rate = (2.0 / 3.0) * expected_max * height * width

    assert math.isclose(result.max_velocity, expected_max, rel_tol=2.0e-15)
    assert math.isclose(result.flow_rate, expected_flow_rate, rel_tol=2.0e-15)
    assert result.u_centerline[0] == 0.0
    assert math.isclose(result.u_centerline[2], result.max_velocity, rel_tol=2.0e-15)
    assert result.u_centerline[-1] == 0.0
