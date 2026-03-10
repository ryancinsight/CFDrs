//! 1D Serpentine resistance solver `PyO3` wrapper.

use pyo3::prelude::*;
use pyo3::exceptions::PyRuntimeError;

// ── Solver ───────────────────────────────────────────────────────────────────

/// 1D serpentine channel resistance model.
///
/// Computes total pressure drop through a serpentine channel accounting for:
/// - Straight section friction (Shah-London for rectangular, Hagen-Poiseuille for circular)
/// - Dean flow enhancement in curved sections (White 1929, Ito 1959)
/// - Bend minor losses (Idelchik 2007)
///
/// # Physics
///
/// Dean number: De = Re * `sqrt(D_h` / `2R_c`)
/// Curved friction enhancement: `f_curved/f_straight` depends on De regime
/// Bend loss: `K_bend` = C1 + C2/Re (per 180° turn)
/// Total: ΔP = f*(`L/D_h`)*(ρV²/2) + N*K_bend*(ρV²/2)
#[pyclass(name = "SerpentineSolver1D")]
pub struct PySerpentineSolver1D {
    #[pyo3(get)]
    width: f64,
    #[pyo3(get)]
    height: f64,
    #[pyo3(get)]
    straight_length: f64,
    #[pyo3(get)]
    num_segments: usize,
    #[pyo3(get)]
    bend_radius: f64,
}

#[pymethods]
impl PySerpentineSolver1D {
    #[new]
    #[pyo3(signature = (width, height, straight_length, num_segments=10, bend_radius=0.002))]
    fn new(width: f64, height: f64, straight_length: f64, num_segments: usize, bend_radius: f64) -> Self {
        PySerpentineSolver1D { width, height, straight_length, num_segments, bend_radius }
    }

    /// Solve serpentine resistance for given flow conditions.
    fn solve(&self, velocity: f64, blood_type: &str) -> PyResult<PySerpentineResult1D> {
        use cfd_1d::physics::resistance::models::{
            FlowConditions, ResistanceModel, SerpentineCrossSection, SerpentineModel,
        };
        use cfd_core::physics::fluid::blood::CassonBlood as RustCasson;
        use cfd_core::physics::fluid::blood::CarreauYasudaBlood as RustCY;

        let cross_section = SerpentineCrossSection::Rectangular {
            width: self.width,
            height: self.height,
        };
        let total_straight = self.straight_length * self.num_segments as f64;
        let model = SerpentineModel::new(
            total_straight, self.num_segments, cross_section, self.bend_radius,
        );

        let dh = cross_section.hydraulic_diameter();
        let density = 1060.0;

        let (mu, re) = match blood_type {
            "casson" => {
                let blood = RustCasson::<f64>::normal_blood();
                let mu = blood.apparent_viscosity(velocity / dh * 8.0);
                let re = density * velocity * dh / mu;
                (mu, re)
            }
            "carreau_yasuda" => {
                let blood = RustCY::<f64>::normal_blood();
                let mu = blood.apparent_viscosity(velocity / dh * 8.0);
                let re = density * velocity * dh / mu;
                (mu, re)
            }
            _ => {
                let mu = 0.0035;
                let re = density * velocity * dh / mu;
                (mu, re)
            }
        };

        let mut conditions = FlowConditions::new(velocity);
        conditions.reynolds_number = Some(re);

        let fluid = cfd_core::physics::fluid::ConstantPropertyFluid::new(
            "blood".to_string(), density, mu, 3617.0, 0.52, 1570.0,
        );

        let resistance = model.calculate_resistance(&fluid, &conditions)
            .map_err(|e| PyRuntimeError::new_err(format!("Serpentine solver error: {e}")))?;

        let area = self.width * self.height;
        let flow_rate = velocity * area;
        let pressure_drop = resistance * flow_rate;
        let dean_number = re * (dh / (2.0 * self.bend_radius)).sqrt();

        Ok(PySerpentineResult1D {
            pressure_drop,
            resistance,
            dean_number,
            reynolds_number: re,
            apparent_viscosity: mu,
        })
    }

    fn __str__(&self) -> String {
        format!(
            "SerpentineSolver1D(w={:.0} μm, h={:.0} μm, {} segments, R_bend={:.1} mm)",
            self.width * 1e6,
            self.height * 1e6,
            self.num_segments,
            self.bend_radius * 1e3
        )
    }
}

// ── Result ───────────────────────────────────────────────────────────────────

/// Result from 1D serpentine resistance calculation
#[pyclass(name = "SerpentineResult1D")]
#[derive(Debug, Clone, Copy)]
pub struct PySerpentineResult1D {
    #[pyo3(get)]
    pub pressure_drop: f64,
    #[pyo3(get)]
    pub resistance: f64,
    #[pyo3(get)]
    pub dean_number: f64,
    #[pyo3(get)]
    pub reynolds_number: f64,
    #[pyo3(get)]
    pub apparent_viscosity: f64,
}

#[pymethods]
impl PySerpentineResult1D {
    fn __str__(&self) -> String {
        format!(
            "SerpentineResult1D(ΔP={:.2} Pa, De={:.1}, Re={:.1})",
            self.pressure_drop, self.dean_number, self.reynolds_number
        )
    }
}
