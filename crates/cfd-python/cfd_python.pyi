"""Typed interface for the ``cfd_python`` PyO3 extension module."""

from types import ModuleType
from typing import Sequence, Tuple, Union

import numpy as np
from numpy.typing import NDArray

__version__: str


class CassonBlood:
    def __init__(self) -> None: ...
    def apparent_viscosity(self, gamma: float) -> float: ...
    def viscosity(self, gamma: float) -> float: ...
    def yield_stress(self) -> float: ...
    def density(self) -> float: ...
    def viscosity_high_shear(self) -> float: ...


class CarreauYasudaBlood:
    def __init__(self) -> None: ...
    def apparent_viscosity(self, gamma: float) -> float: ...
    def viscosity(self, gamma: float) -> float: ...
    def density(self) -> float: ...
    def viscosity_zero_shear(self) -> float: ...
    def viscosity_high_shear(self) -> float: ...


class CrossBlood:
    def __init__(self) -> None: ...
    def apparent_viscosity(self, gamma: float) -> float: ...
    def viscosity(self, gamma: float) -> float: ...
    def density(self) -> float: ...
    def viscosity_zero_shear(self) -> float: ...
    def viscosity_high_shear(self) -> float: ...
    def time_constant(self) -> float: ...
    def rate_index(self) -> float: ...


class FahraeuasLindqvist:
    def __init__(self, diameter: float, hematocrit: float) -> None: ...
    def is_significant(self) -> bool: ...
    def relative_viscosity(self) -> float: ...
    def apparent_viscosity(self) -> float: ...
    def tube_hematocrit(self) -> float: ...
    def plasma_viscosity(self) -> float: ...


BloodModel = Union[CassonBlood, CarreauYasudaBlood, CrossBlood]


class BifurcationResult:
    q_parent: float
    q_1: float
    q_2: float
    p_parent: float
    p_1: float
    p_2: float
    dp_1: float
    dp_2: float
    gamma_1: float
    gamma_2: float
    mu_1: float
    mu_2: float
    wss_1: float
    wss_2: float
    mass_conservation_error: float
    pressure_continuity_error: float
    def flow_split_ratio(self) -> float: ...
    def is_valid(self, tolerance: float) -> bool: ...


class BifurcationSolver:
    d_parent: float
    d_daughter1: float
    d_daughter2: float
    length: float
    flow_split_ratio: float
    def __init__(
        self,
        d_parent: float,
        d_daughter1: float,
        d_daughter2: float,
        length: float = ...,
        flow_split_ratio: float = ...,
    ) -> None: ...
    def solve(self, flow_rate: float, pressure: float, blood_type: str) -> BifurcationResult: ...
    def murray_law_deviation(self) -> float: ...
    def areas(self) -> Tuple[float, float, float]: ...


class TrifurcationResult:
    q_parent: float
    q_daughters: Sequence[float]
    p_parent: float
    p_daughters: Sequence[float]
    mass_conservation_error: float


class TrifurcationSolver:
    def __init__(self, *args: float, **kwargs: float) -> None: ...
    def solve(self, flow_rate: float, pressure: float, blood_type: str) -> TrifurcationResult: ...


class PoiseuilleConfig2D:
    height: float
    width: float
    length: float
    ny: int
    pressure_gradient: float
    tolerance: float
    max_iterations: int
    relaxation_factor: float
    def __init__(
        self,
        height: float,
        width: float,
        length: float,
        ny: int,
        pressure_gradient: float,
        tolerance: float = ...,
        max_iterations: int = ...,
        relaxation_factor: float = ...,
    ) -> None: ...


class Poiseuille2DResult:
    max_velocity: float
    flow_rate: float
    reynolds_number: float
    pressure_drop: float
    wall_shear_stress: float
    u_centerline: Sequence[float]
    y_coords: Sequence[float]


class PoiseuilleResult2D:
    y_coords: Sequence[float]
    velocity: Sequence[float]
    shear_rate: Sequence[float]
    viscosity: Sequence[float]
    flow_rate: float
    wall_shear_stress: float
    iterations: int


class Poiseuille2DSolver:
    def __init__(self, height: float, width: float, length: float, nx: int, ny: int) -> None: ...
    def solve(self, pressure_drop: float, blood_type: str) -> Poiseuille2DResult: ...


class PoiseuilleSolver2D_Legacy:
    def __init__(self, config: PoiseuilleConfig2D) -> None: ...
    def solve(self, blood: BloodModel) -> PoiseuilleResult2D: ...


class WomersleyNumber:
    def __init__(self, radius: float, omega: float, density: float, viscosity: float) -> None: ...
    @classmethod
    def from_heart_rate(cls, diameter: float, heart_rate_hz: float, density: float, viscosity: float) -> "WomersleyNumber": ...
    @classmethod
    def human_aorta(cls) -> "WomersleyNumber": ...
    @classmethod
    def human_femoral(cls) -> "WomersleyNumber": ...
    def value(self) -> float: ...
    def stokes_layer_thickness(self) -> float: ...
    def radius_to_stokes_ratio(self) -> float: ...
    def flow_regime(self) -> str: ...


class WomersleyProfile:
    def __init__(self, radius: float, omega: float, density: float, viscosity: float, pressure_amplitude: float) -> None: ...
    def velocity(self, xi: float, t: float) -> float: ...
    def centerline_velocity(self, t: float) -> float: ...
    def wall_shear_stress(self, t: float) -> float: ...
    def flow_rate(self, t: float) -> float: ...


class WomersleyFlow:
    def __init__(self, *args: float, **kwargs: float) -> None: ...
    def womersley_number(self) -> float: ...
    def velocity(self, xi: float, t: float) -> float: ...
    def impedance_magnitude(self) -> float: ...


class TrifurcationSolver2D:
    def __init__(self, width: float, length: float, angle: float, nx: int) -> None: ...
    def solve(self, flow_rate: float, blood_type: str) -> "TrifurcationResult2D": ...


class TrifurcationResult2D:
    execution_time: float
    mass_conservation_error: float


class BifurcationSolver2D:
    def __init__(self, *args: float, **kwargs: float) -> None: ...
    def murray_law_deviation(self) -> float: ...
    def solve(self, inlet_velocity: float, blood_type: str) -> "BifurcationResult2D": ...


class BifurcationResult2D:
    q_parent: float
    q_daughter1: float
    q_daughter2: float
    mass_balance_error: float
    flow_split_ratio: float


class CavitySolver2D:
    def __init__(self, reynolds: float, nx: int, ny: int, lid_velocity: float, cavity_size: float) -> None: ...
    def viscosity(self) -> float: ...
    def ghia_u_centerline(self) -> NDArray[np.float64]: ...
    def ghia_v_centerline(self) -> NDArray[np.float64]: ...
    def solve(self) -> "CavityResult2D": ...


class CavityResult2D:
    l2_error: float
    u_centerline: Sequence[float]
    v_centerline: Sequence[float]
    y_coords: Sequence[float]
    x_coords: Sequence[float]
    converged: bool


class SerpentineSolver1D:
    def __init__(self, *args: float, **kwargs: float) -> None: ...
    def solve(self, velocity: float, blood_type: str) -> "SerpentineResult1D": ...


class SerpentineResult1D:
    pressure_drop: float
    resistance: float
    dean_number: float
    reynolds_number: float
    apparent_viscosity: float


class VenturiSolver2D:
    def __init__(self, *args: float, **kwargs: float) -> None: ...
    @classmethod
    def iso_5167_standard(cls, nx: int, ny: int) -> "VenturiSolver2D": ...
    def pressure_coefficient_analytical(self) -> float: ...
    def area_ratio(self) -> float: ...
    def total_length(self) -> float: ...
    def solve(self, inlet_velocity: float, blood_type: str) -> "VenturiResult2D": ...


class VenturiResult2D:
    cp_throat: float
    pressure_recovery: float
    velocity_ratio: float
    mass_conservation_error: float


class VenturiSolver1D:
    def __init__(self, *args: float, **kwargs: float) -> None: ...
    def beta(self) -> float: ...
    def solve(self, velocity: float, blood_type: str) -> "VenturiResult1D": ...


class VenturiResult1D:
    pressure_drop: float
    resistance: float
    dp_bernoulli: float
    reynolds_number: float
    beta: float
    apparent_viscosity: float


class Bifurcation3DSolver:
    def __init__(self, *args: float, **kwargs: float) -> None: ...
    def murray_law_deviation(self) -> float: ...
    def grid_size(self) -> Tuple[int, int, int]: ...
    def solve(self, flow_rate: float, blood_type: str) -> "Bifurcation3DResult": ...


class Bifurcation3DResult:
    max_wss: float
    min_wss: float
    mean_wss: float
    wss_ratio: float
    flow_split_ratio: float
    mass_conservation_error: float


class Trifurcation3DSolver:
    def __init__(self, d_parent: float, d_daughter: float, length: float) -> None: ...
    def solve(self, flow_rate: float, blood_type: str) -> "Trifurcation3DResult": ...


class Trifurcation3DResult:
    max_wss: float
    min_wss: float
    flow_rates: Sequence[float]
    mass_conservation_error: float


class Poiseuille3DSolver:
    def __init__(self, diameter: float, length: float, nr: int, ntheta: int, nz: int) -> None: ...
    def analytical_max_velocity(self, pressure_gradient: float, viscosity: float) -> float: ...
    def analytical_flow_rate(self, pressure_gradient: float, viscosity: float) -> float: ...
    def solve(self, pressure_drop: float, blood_type: str) -> "Poiseuille3DResult": ...


class Poiseuille3DResult:
    max_velocity: float
    flow_rate: float
    reynolds_number: float
    wall_shear_stress: float


class Venturi3DSolver:
    def __init__(self, *args: float, **kwargs: Union[float, bool]) -> None: ...
    def solve(self, flow_rate: float, blood_type: str) -> "Venturi3DResult": ...


class Venturi3DResult:
    u_inlet: float
    u_throat: float
    p_inlet: float
    p_throat: float
    p_outlet: float
    dp_throat: float
    dp_recovery: float
    cp_throat: float
    cp_recovery: float


class Serpentine3DSolver:
    def __init__(self, diameter: float, wavelength: float, amplitude: float, cycles: int, circular: bool) -> None: ...
    def solve(self, flow_rate: float, blood_type: str) -> "Serpentine3DResult": ...


class Serpentine3DResult:
    u_inlet: float
    p_inlet: float
    dp_total: float
    dean_number: float


cfd_python: ModuleType
validation: ModuleType
