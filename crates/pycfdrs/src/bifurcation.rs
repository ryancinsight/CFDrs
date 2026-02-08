//! Bifurcation solver PyO3 wrapper
//!
//! Provides Python interface to the 1D bifurcation solver with non-Newtonian blood flow.

use cfd_1d::bifurcation::{BifurcationJunction, TrifurcationJunction};
use cfd_1d::channel::{Channel, ChannelGeometry};
use cfd_core::physics::fluid::blood::{CarreauYasudaBlood, CassonBlood};
use pyo3::exceptions::{PyRuntimeError, PyTypeError};
use pyo3::prelude::*;

use crate::result_types::PyBifurcationResult;

/// 1D bifurcation flow solver with non-Newtonian blood
///
/// This solver computes pressure and flow distribution in bifurcating vessels using
/// the Hagen-Poiseuille equation with shear-rate dependent viscosity.
///
/// # Physics
///
/// The bifurcation satisfies:
/// - Mass conservation: Q_parent = Q_1 + Q_2
/// - Pressure-flow relationship: ΔP = (128μQL)/(πD⁴)
/// - Non-Newtonian viscosity: μ(γ̇) from blood model
/// - Wall shear rate: γ̇ = 32Q/(πD³)
///
/// # Example
///
/// ```python
/// import pycfdrs
///
/// # Create symmetric bifurcation
/// bifurc = pycfdrs.BifurcationSolver(
///     d_parent=100e-6,  # 100 μm
///     d_daughter1=80e-6,
///     d_daughter2=80e-6,
///     length=1e-3  # 1 mm
/// )
///
/// # Create blood model
/// blood = pycfdrs.CassonBlood()
///
/// # Solve
/// result = bifurc.solve(
///     flow_rate=3e-8,  # 30 nL/s
///     pressure=40.0,   # 40 Pa
///     blood=blood
/// )
///
/// print(f"Flow split: {result.flow_split_ratio():.3f}")
/// print(f"Pressure drop D1: {result.dp_1:.3f} Pa")
/// ```
#[pyclass(name = "BifurcationSolver")]
pub struct PyBifurcationSolver {
    #[pyo3(get)]
    d_parent: f64,
    #[pyo3(get)]
    d_daughter1: f64,
    #[pyo3(get)]
    d_daughter2: f64,
    #[pyo3(get)]
    length: f64,
    #[pyo3(get)]
    flow_split_ratio: f64,
}

#[pymethods]
impl PyBifurcationSolver {
    /// Create new bifurcation solver
    ///
    /// # Arguments
    /// - `d_parent`: Parent vessel diameter [m]
    /// - `d_daughter1`: Daughter 1 vessel diameter [m]
    /// - `d_daughter2`: Daughter 2 vessel diameter [m]
    /// - `length`: Vessel length [m] (default: 1e-3)
    /// - `flow_split_ratio`: Flow distribution ratio Q_1/(Q_1+Q_2) (default: 0.5)
    #[new]
    #[pyo3(signature = (d_parent, d_daughter1, d_daughter2, length=1e-3, flow_split_ratio=0.5))]
    fn new(
        d_parent: f64,
        d_daughter1: f64,
        d_daughter2: f64,
        length: f64,
        flow_split_ratio: f64,
    ) -> Self {
        PyBifurcationSolver {
            d_parent,
            d_daughter1,
            d_daughter2,
            length,
            flow_split_ratio,
        }
    }

    /// Solve bifurcation with Casson blood model
    ///
    /// # Arguments
    /// - `flow_rate`: Parent inlet flow rate [m³/s]
    /// - `pressure`: Parent inlet pressure [Pa]
    /// - `blood`: Blood model (CassonBlood or CarreauYasudaBlood)
    ///
    /// # Returns
    /// - BifurcationResult containing all flow and pressure data
    ///
    /// # Raises
    /// - PyTypeError if blood model is not recognized
    /// - RuntimeError if solver fails to converge
    fn solve(
        &self,
        flow_rate: f64,
        pressure: f64,
        blood_type: &str,
    ) -> PyResult<PyBifurcationResult> {
        let roughness = 0.0; // Smooth walls for microfluidic channels

        let parent_geom = ChannelGeometry::circular(self.length, self.d_parent, roughness);
        let d1_geom = ChannelGeometry::circular(self.length, self.d_daughter1, roughness);
        let d2_geom = ChannelGeometry::circular(self.length, self.d_daughter2, roughness);

        let parent = Channel::new(parent_geom);
        let d1 = Channel::new(d1_geom);
        let d2 = Channel::new(d2_geom);

        let bifurc = BifurcationJunction::new(parent, d1, d2, self.flow_split_ratio);

        match blood_type {
            "casson" => {
                let blood = CassonBlood::<f64>::normal_blood();
                // Clone to satisfy Copy bound - blood models can't be Copy due to String field
                let blood_copy = blood.clone();
                match bifurc.solve(blood_copy, flow_rate, pressure) {
                    Ok(solution) => {
                        let wss_1 = (4.0 * solution.mu_1 * solution.q_1)
                            / (std::f64::consts::PI * self.d_daughter1.powi(3));
                        let wss_2 = (4.0 * solution.mu_2 * solution.q_2)
                            / (std::f64::consts::PI * self.d_daughter2.powi(3));

                        Ok(PyBifurcationResult::new(
                            solution.q_parent,
                            solution.q_1,
                            solution.q_2,
                            solution.p_parent,
                            solution.p_1,
                            solution.p_2,
                            solution.dp_1,
                            solution.dp_2,
                            solution.gamma_1,
                            solution.gamma_2,
                            solution.mu_1,
                            solution.mu_2,
                            wss_1,
                            wss_2,
                            solution.mass_conservation_error,
                            solution.junction_pressure_error,
                        ))
                    }
                    Err(e) => Err(PyTypeError::new_err(format!("Solver failed: {}", e))),
                }
            }
            "carreau_yasuda" => {
                let blood = CarreauYasudaBlood::<f64>::normal_blood();
                match bifurc.solve(blood, flow_rate, pressure) {
                    Ok(solution) => {
                        let wss_1 = (4.0 * solution.mu_1 * solution.q_1)
                            / (std::f64::consts::PI * self.d_daughter1.powi(3));
                        let wss_2 = (4.0 * solution.mu_2 * solution.q_2)
                            / (std::f64::consts::PI * self.d_daughter2.powi(3));

                        Ok(PyBifurcationResult::new(
                            solution.q_parent,
                            solution.q_1,
                            solution.q_2,
                            solution.p_parent,
                            solution.p_1,
                            solution.p_2,
                            solution.dp_1,
                            solution.dp_2,
                            solution.gamma_1,
                            solution.gamma_2,
                            solution.mu_1,
                            solution.mu_2,
                            wss_1,
                            wss_2,
                            solution.mass_conservation_error,
                            solution.junction_pressure_error,
                        ))
                    }
                    Err(e) => Err(PyTypeError::new_err(format!("Solver failed: {}", e))),
                }
            }
            _ => Err(PyTypeError::new_err(
                "blood_type must be 'casson' or 'carreau_yasuda'",
            )),
        }
    }

    /// Calculate Murray's law deviation
    ///
    /// Returns: |D_parent³ - (D_1³ + D_2³)| / D_parent³
    fn murray_law_deviation(&self) -> f64 {
        let d0_cubed = self.d_parent.powi(3);
        let daughters_cubed = self.d_daughter1.powi(3) + self.d_daughter2.powi(3);
        (d0_cubed - daughters_cubed).abs() / d0_cubed
    }

    /// Get vessel area ratios
    ///
    /// Returns: (A_parent, A_daughter1, A_daughter2)
    fn areas(&self) -> (f64, f64, f64) {
        let pi = std::f64::consts::PI;
        let a_parent = pi * (self.d_parent / 2.0).powi(2);
        let a_d1 = pi * (self.d_daughter1 / 2.0).powi(2);
        let a_d2 = pi * (self.d_daughter2 / 2.0).powi(2);
        (a_parent, a_d1, a_d2)
    }

    fn __str__(&self) -> String {
        format!(
            "BifurcationSolver(d_p={:.1} μm, d_1={:.1} μm, d_2={:.1} μm, L={:.2} mm)",
            self.d_parent * 1e6,
            self.d_daughter1 * 1e6,
            self.d_daughter2 * 1e6,
            self.length * 1e3
        )
    }

    fn __repr__(&self) -> String {
        self.__str__()
    }
}
/// 1D trifurcation flow solver
#[pyclass(name = "TrifurcationSolver")]
pub struct PyTrifurcationSolver {
    #[pyo3(get)]
    d_parent: f64,
    #[pyo3(get)]
    d_daughter1: f64,
    #[pyo3(get)]
    d_daughter2: f64,
    #[pyo3(get)]
    d_daughter3: f64,
    #[pyo3(get)]
    length: f64,
}

#[pymethods]
impl PyTrifurcationSolver {
    #[new]
    #[pyo3(signature = (d_parent, d_daughter1, d_daughter2, d_daughter3, length=1e-3))]
    fn new(
        d_parent: f64,
        d_daughter1: f64,
        d_daughter2: f64,
        d_daughter3: f64,
        length: f64,
    ) -> Self {
        PyTrifurcationSolver {
            d_parent,
            d_daughter1,
            d_daughter2,
            d_daughter3,
            length,
        }
    }

    /// Solve trifurcation with Casson blood
    fn solve(
        &self,
        flow_rate: f64,
        pressure: f64,
        _blood_type: &str,
    ) -> PyResult<PyTrifurcationResult> {
        let roughness = 0.0; // Smooth walls for microfluidic channels

        let parent_geom = ChannelGeometry::circular(self.length, self.d_parent, roughness);
        let d1_geom = ChannelGeometry::circular(self.length, self.d_daughter1, roughness);
        let d2_geom = ChannelGeometry::circular(self.length, self.d_daughter2, roughness);
        let d3_geom = ChannelGeometry::circular(self.length, self.d_daughter3, roughness);

        let parent = Channel::new(parent_geom);
        let d1 = Channel::new(d1_geom);
        let d2 = Channel::new(d2_geom);
        let d3 = Channel::new(d3_geom);

        let trifurc =
            TrifurcationJunction::new(parent, d1, d2, d3, (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0));
        let blood = CassonBlood::<f64>::normal_blood();

        match trifurc.solve(blood, flow_rate, pressure) {
            Ok(solution) => Ok(PyTrifurcationResult {
                q_parent: solution.q_parent,
                q_daughters: [solution.q_1, solution.q_2, solution.q_3],
                p_parent: solution.p_parent,
                p_daughters: [solution.p_1, solution.p_2, solution.p_3],
                mass_conservation_error: solution.mass_conservation_error,
            }),
            Err(e) => Err(pyo3::exceptions::PyRuntimeError::new_err(format!(
                "Solver failed: {}",
                e
            ))),
        }
    }

    fn __str__(&self) -> String {
        format!(
            "TrifurcationSolver(d_p={:.1} μm, daughters=[{:.1}, {:.1}, {:.1}] μm, L={:.2} mm)",
            self.d_parent * 1e6,
            self.d_daughter1 * 1e6,
            self.d_daughter2 * 1e6,
            self.d_daughter3 * 1e6,
            self.length * 1e3
        )
    }
}

#[pyclass(name = "TrifurcationResult")]
pub struct PyTrifurcationResult {
    #[pyo3(get)]
    pub q_parent: f64,
    #[pyo3(get)]
    pub q_daughters: [f64; 3],
    #[pyo3(get)]
    pub p_parent: f64,
    #[pyo3(get)]
    pub p_daughters: [f64; 3],
    #[pyo3(get)]
    pub mass_conservation_error: f64,
}

#[pymethods]
impl PyTrifurcationResult {
    fn __str__(&self) -> String {
        format!(
            "TrifurcationResult(Q_daughters=[{:.2e}, {:.2e}, {:.2e}] m³/s)",
            self.q_daughters[0], self.q_daughters[1], self.q_daughters[2]
        )
    }
}
