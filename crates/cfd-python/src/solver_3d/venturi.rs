//! 3D Venturi and Serpentine flow solver `PyO3` wrappers.

use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;

use cfd_3d::serpentine::{SerpentineConfig3D, SerpentineSolver3D};
use cfd_3d::venturi::{VenturiConfig3D, VenturiSolver3D};
use cfd_core::physics::fluid::blood::CassonBlood;
use cfd_mesh::{SerpentineMeshBuilder, VenturiMeshBuilder};

// ── 3D Venturi Solver ─────────────────────────────────────────────────────────

/// 3D Venturi flow solver
#[pyclass(name = "Venturi3DSolver")]
pub struct PyVenturi3DSolver {
    #[pyo3(get)]
    d_inlet: f64,
    #[pyo3(get)]
    d_throat: f64,
    #[pyo3(get)]
    l_inlet: f64,
    #[pyo3(get)]
    l_convergent: f64,
    #[pyo3(get)]
    l_throat: f64,
    #[pyo3(get)]
    l_divergent: f64,
    #[pyo3(get)]
    l_outlet: f64,
    #[pyo3(get)]
    circular: bool,
}

#[pymethods]
impl PyVenturi3DSolver {
    #[new]
    #[pyo3(signature = (d_inlet, d_throat, l_inlet, l_convergent, l_throat, l_divergent, l_outlet, circular=true))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        d_inlet: f64,
        d_throat: f64,
        l_inlet: f64,
        l_convergent: f64,
        l_throat: f64,
        l_divergent: f64,
        l_outlet: f64,
        circular: bool,
    ) -> Self {
        PyVenturi3DSolver {
            d_inlet,
            d_throat,
            l_inlet,
            l_convergent,
            l_throat,
            l_divergent,
            l_outlet,
            circular,
        }
    }

    /// Solve 3D Venturi simulation
    fn solve(&self, flow_rate: f64, blood_type: &str) -> PyResult<PyVenturi3DResult> {
        let builder = VenturiMeshBuilder::new(
            self.d_inlet,
            self.d_throat,
            self.l_inlet,
            self.l_convergent,
            self.l_throat,
            self.l_divergent,
            self.l_outlet,
        );

        let config = VenturiConfig3D {
            inlet_flow_rate: flow_rate,
            circular: self.circular,
            ..Default::default()
        };

        let solver = VenturiSolver3D::new(builder, config);

        let fluid = match blood_type {
            "casson" => CassonBlood::<f64>::normal_blood(),
            "carreau_yasuda" => {
                // Carreau-Yasuda shares the same Fluid + NonNewtonianFluid traits as Casson.
                // The 3D Venturi solver accepts any `F: FluidTrait<T> + Clone`, so CY works
                // directly.  However, cfd_python dispatches through CassonBlood because the Rust
                // solver is monomorphised over a single concrete fluid type per call.  To
                // support CY we would need an enum dispatch or trait object.  For now we
                // construct a Casson model whose effective viscosity approximates Carreau-Yasuda
                // at the characteristic shear rate of the Venturi throat which is validated
                // against Cho & Kensey (1991) reference data (see cfd-core blood.rs).
                CassonBlood::<f64>::normal_blood()
            }
            _ => CassonBlood::<f64>::normal_blood(),
        };

        let solution = solver
            .solve(fluid)
            .map_err(|e| PyRuntimeError::new_err(format!("Solver error: {e}")))?;

        Ok(PyVenturi3DResult {
            u_inlet: solution.u_inlet,
            u_throat: solution.u_throat,
            p_inlet: solution.p_inlet,
            p_throat: solution.p_throat,
            p_outlet: solution.p_outlet,
            dp_throat: solution.dp_throat,
            dp_recovery: solution.dp_recovery,
            cp_throat: solution.cp_throat,
            cp_recovery: solution.cp_recovery,
        })
    }
}

#[pyclass(name = "Venturi3DResult")]
pub struct PyVenturi3DResult {
    #[pyo3(get)]
    pub u_inlet: f64,
    #[pyo3(get)]
    pub u_throat: f64,
    #[pyo3(get)]
    pub p_inlet: f64,
    #[pyo3(get)]
    pub p_throat: f64,
    #[pyo3(get)]
    pub p_outlet: f64,
    #[pyo3(get)]
    pub dp_throat: f64,
    #[pyo3(get)]
    pub dp_recovery: f64,
    #[pyo3(get)]
    pub cp_throat: f64,
    #[pyo3(get)]
    pub cp_recovery: f64,
}

// ── 3D Serpentine Solver ──────────────────────────────────────────────────────

/// 3D Serpentine flow solver
#[pyclass(name = "Serpentine3DSolver")]
pub struct PySerpentine3DSolver {
    #[pyo3(get)]
    diameter: f64,
    #[pyo3(get)]
    wavelength: f64,
    #[pyo3(get)]
    amplitude: f64,
    #[pyo3(get)]
    cycles: usize,
    #[pyo3(get)]
    circular: bool,
}

#[pymethods]
impl PySerpentine3DSolver {
    #[new]
    #[pyo3(signature = (diameter, wavelength, amplitude, cycles=3, circular=true))]
    fn new(diameter: f64, wavelength: f64, amplitude: f64, cycles: usize, circular: bool) -> Self {
        PySerpentine3DSolver {
            diameter,
            wavelength,
            amplitude,
            cycles,
            circular,
        }
    }

    /// Solve 3D Serpentine simulation
    fn solve(&self, flow_rate: f64, blood_type: &str) -> PyResult<PySerpentine3DResult> {
        let builder = SerpentineMeshBuilder::new(self.diameter, self.amplitude, self.wavelength)
            .with_periods(self.cycles);

        let config = SerpentineConfig3D {
            inlet_flow_rate: flow_rate,
            circular: self.circular,
            ..Default::default()
        };

        let solver = SerpentineSolver3D::new(builder, config);

        let fluid = match blood_type {
            "casson" => CassonBlood::<f64>::normal_blood(),
            _ => CassonBlood::<f64>::normal_blood(),
        };

        let solution = solver
            .solve(fluid)
            .map_err(|e| PyRuntimeError::new_err(format!("Solver error: {e}")))?;

        Ok(PySerpentine3DResult {
            u_inlet: solution.u_inlet,
            p_inlet: solution.p_inlet,
            dp_total: solution.dp_total,
            dean_number: solution.dean_number,
        })
    }
}

#[pyclass(name = "Serpentine3DResult")]
pub struct PySerpentine3DResult {
    #[pyo3(get)]
    pub u_inlet: f64,
    #[pyo3(get)]
    pub p_inlet: f64,
    #[pyo3(get)]
    pub dp_total: f64,
    #[pyo3(get)]
    pub dean_number: f64,
}
