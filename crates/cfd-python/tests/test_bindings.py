"""Value-semantic smoke tests for the installed cfd-python wheel."""

import importlib.metadata
import math
import threading
from typing import Dict

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
    assert womersley.flow_regime() == "Inertial"


def test_poiseuille_solver_returns_analytical_profile_values() -> None:
    """The installed binding returns the analytical channel-flow values."""
    height = 2.0e-3
    width = 1.0e-2
    length = 5.0e-2
    pressure_drop = 100.0
    # The Newtonian branch uses cfd_core's canonical infinite-shear blood
    # viscosity (0.00345 Pa·s; Cho & Kensey, 1991).
    viscosity = 3.45e-3
    solver = cfd_python.Poiseuille2DSolver(height, width, length, 9, 5)
    result = solver.solve(pressure_drop, "newtonian")
    expected_max = height**2 * (pressure_drop / length) / (8.0 * viscosity)
    expected_flow_rate = (2.0 / 3.0) * expected_max * height * width

    assert math.isclose(result.max_velocity, expected_max, rel_tol=2.0e-15)
    assert math.isclose(result.flow_rate, expected_flow_rate, rel_tol=2.0e-15)
    assert result.u_centerline[0] == 0.0
    assert math.isclose(result.u_centerline[2], result.max_velocity, rel_tol=2.0e-15)
    assert result.u_centerline[-1] == 0.0


def test_solver_releases_gil_for_another_python_thread() -> None:
    """A Rust solver does not block another Python thread from progressing."""
    solver = cfd_python.Poiseuille2DSolver(2.0e-3, 1.0e-2, 5.0e-2, 1024, 512)
    solver_started = threading.Event()
    solver_finished = threading.Event()
    progress_finished = threading.Event()
    outcome: Dict[str, object] = {}

    def run_solver() -> None:
        solver_started.set()
        outcome["result"] = solver.solve(100.0, "newtonian")
        solver_finished.set()

    def run_python_progress() -> None:
        assert solver_started.wait(timeout=10.0)
        outcome["progress_sum"] = sum(range(100_000))
        outcome["progress_before_solver_finished"] = not solver_finished.is_set()
        progress_finished.set()

    solver_thread = threading.Thread(target=run_solver, daemon=True)
    progress_thread = threading.Thread(target=run_python_progress, daemon=True)
    solver_thread.start()
    assert solver_started.wait(timeout=10.0)
    progress_thread.start()

    assert progress_finished.wait(timeout=30.0)
    solver_thread.join(timeout=30.0)
    progress_thread.join(timeout=30.0)

    assert not solver_thread.is_alive()
    assert not progress_thread.is_alive()
    assert outcome["progress_sum"] == 4_999_950_000
    assert outcome["progress_before_solver_finished"] is True

    result = outcome["result"]
    assert isinstance(result, cfd_python.Poiseuille2DResult)
    assert math.isfinite(result.max_velocity)
    assert result.u_centerline[0] == 0.0
    assert result.u_centerline[-1] == 0.0
    assert result.u_centerline[256] > 0.0
