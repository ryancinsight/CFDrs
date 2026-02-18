//! Womersley pulsatile flow PyO3 wrappers

use cfd_1d::vascular::womersley::{
    WomersleyFlow as RustWomersleyFlow, WomersleyNumber as RustWomersleyNumber,
    WomersleyProfile as RustWomersleyProfile,
};
use pyo3::prelude::*;

/// Womersley number calculator and flow regime classifier
///
/// alpha = R * sqrt(omega * rho / mu)
///
/// Characterises the ratio of unsteady inertial to viscous forces
/// in pulsatile flow.
///
/// # Flow Regimes
/// - alpha < 1: quasi-steady (Poiseuille-like)
/// - 1 <= alpha < 3: transitional
/// - 3 <= alpha < 10: inertia-dominated
/// - alpha >= 10: plug flow with thin Stokes layer
#[pyclass(name = "WomersleyNumber")]
pub struct PyWomersleyNumber {
    inner: RustWomersleyNumber<f64>,
}

#[pymethods]
impl PyWomersleyNumber {
    /// Create Womersley number calculator
    ///
    /// # Arguments
    /// - `radius`: Vessel radius [m]
    /// - `omega`: Angular frequency [rad/s]
    /// - `density`: Fluid density [kg/m^3]
    /// - `viscosity`: Dynamic viscosity [Pa.s]
    #[new]
    fn new(radius: f64, omega: f64, density: f64, viscosity: f64) -> Self {
        PyWomersleyNumber {
            inner: RustWomersleyNumber::<f64>::new(radius, omega, density, viscosity),
        }
    }

    /// Create from vessel diameter and heart rate (Hz)
    #[staticmethod]
    fn from_heart_rate(diameter: f64, heart_rate_hz: f64, density: f64, viscosity: f64) -> Self {
        PyWomersleyNumber {
            inner: RustWomersleyNumber::<f64>::from_heart_rate(
                diameter,
                heart_rate_hz,
                density,
                viscosity,
            ),
        }
    }

    /// Human aorta at 72 bpm (alpha ~ 18)
    #[staticmethod]
    fn human_aorta() -> Self {
        PyWomersleyNumber {
            inner: RustWomersleyNumber::<f64>::human_aorta(),
        }
    }

    /// Human femoral artery at 72 bpm (alpha ~ 3.3)
    #[staticmethod]
    fn human_femoral() -> Self {
        PyWomersleyNumber {
            inner: RustWomersleyNumber::<f64>::human_femoral(),
        }
    }

    /// Calculate the Womersley number alpha [-]
    fn value(&self) -> f64 {
        self.inner.value()
    }

    /// Stokes layer thickness delta [m]
    fn stokes_layer_thickness(&self) -> f64 {
        self.inner.stokes_layer_thickness()
    }

    /// Ratio R / delta [-]
    fn radius_to_stokes_ratio(&self) -> f64 {
        self.inner.radius_to_stokes_ratio()
    }

    /// Classify flow regime as string
    fn flow_regime(&self) -> String {
        format!("{:?}", self.inner.flow_regime())
    }

    fn __str__(&self) -> String {
        format!(
            "WomersleyNumber(alpha={:.4}, regime={})",
            self.value(),
            self.flow_regime()
        )
    }

    fn __repr__(&self) -> String {
        self.__str__()
    }
}

/// Womersley velocity profile calculator
///
/// Computes the unsteady velocity field u(r/R, t) for pulsatile pipe flow
/// using asymptotic approximations valid for low, intermediate, and high
/// Womersley numbers.
#[pyclass(name = "WomersleyProfile")]
pub struct PyWomersleyProfile {
    inner: RustWomersleyProfile<f64>,
}

#[pymethods]
impl PyWomersleyProfile {
    /// Create velocity profile calculator
    ///
    /// # Arguments
    /// - `radius`: Vessel radius [m]
    /// - `omega`: Angular frequency [rad/s]
    /// - `density`: Fluid density [kg/m^3]
    /// - `viscosity`: Dynamic viscosity [Pa.s]
    /// - `pressure_amplitude`: Oscillatory pressure gradient amplitude [Pa/m]
    #[new]
    fn new(
        radius: f64,
        omega: f64,
        density: f64,
        viscosity: f64,
        pressure_amplitude: f64,
    ) -> Self {
        let wom = RustWomersleyNumber::<f64>::new(radius, omega, density, viscosity);
        PyWomersleyProfile {
            inner: RustWomersleyProfile::<f64>::new(wom, pressure_amplitude),
        }
    }

    /// Calculate velocity at radial position xi = r/R and time t
    ///
    /// # Arguments
    /// - `xi`: Dimensionless radial position r/R [0, 1]
    /// - `t`: Time [s]
    ///
    /// # Returns
    /// Axial velocity u [m/s]
    fn velocity(&self, xi: f64, t: f64) -> f64 {
        self.inner.velocity(xi, t)
    }

    /// Centerline velocity at time t [m/s]
    fn centerline_velocity(&self, t: f64) -> f64 {
        self.inner.centerline_velocity(t)
    }

    /// Wall shear stress at time t [Pa]
    fn wall_shear_stress(&self, t: f64) -> f64 {
        self.inner.wall_shear_stress(t)
    }

    /// Volumetric flow rate at time t [m^3/s]
    fn flow_rate(&self, t: f64) -> f64 {
        self.inner.flow_rate(t)
    }

    fn __str__(&self) -> String {
        format!(
            "WomersleyProfile(alpha={:.4})",
            self.inner.womersley.value()
        )
    }

    fn __repr__(&self) -> String {
        self.__str__()
    }
}

/// Complete Womersley flow solver for arterial segments
///
/// Provides time-varying flow solutions combining mean (steady Poiseuille)
/// and pulsatile (Womersley) components.
#[pyclass(name = "WomersleyFlow")]
pub struct PyWomersleyFlow {
    inner: RustWomersleyFlow<f64>,
}

#[pymethods]
impl PyWomersleyFlow {
    /// Create Womersley flow solver
    ///
    /// # Arguments
    /// - `radius`: Vessel radius [m]
    /// - `length`: Vessel length [m]
    /// - `density`: Fluid density [kg/m^3]
    /// - `viscosity`: Dynamic viscosity [Pa.s]
    /// - `omega`: Angular frequency [rad/s]
    /// - `inlet_pressure_amplitude`: Pulsatile inlet pressure amplitude [Pa]
    /// - `mean_pressure_gradient`: Steady mean pressure gradient [Pa/m]
    #[new]
    fn new(
        radius: f64,
        length: f64,
        density: f64,
        viscosity: f64,
        omega: f64,
        inlet_pressure_amplitude: f64,
        mean_pressure_gradient: f64,
    ) -> Self {
        PyWomersleyFlow {
            inner: RustWomersleyFlow::<f64>::new(
                radius,
                length,
                density,
                viscosity,
                omega,
                inlet_pressure_amplitude,
                mean_pressure_gradient,
            ),
        }
    }

    /// Womersley number alpha [-]
    fn womersley_number(&self) -> f64 {
        self.inner.womersley_number().value()
    }

    /// Total velocity (mean + pulsatile) at xi = r/R and time t [m/s]
    fn velocity(&self, xi: f64, t: f64) -> f64 {
        self.inner.velocity(xi, t)
    }

    /// Impedance magnitude |Z| [Pa.s/m^3]
    fn impedance_magnitude(&self) -> f64 {
        self.inner.impedance_magnitude()
    }

    fn __str__(&self) -> String {
        format!(
            "WomersleyFlow(alpha={:.4}, R={:.1} um, L={:.1} mm)",
            self.womersley_number(),
            self.inner.radius * 1e6,
            self.inner.length * 1e3,
        )
    }

    fn __repr__(&self) -> String {
        self.__str__()
    }
}
