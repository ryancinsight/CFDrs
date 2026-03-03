//! 3D Y-Bifurcation Flow benchmark with blood rheology
//!
//! Validates 3D bifurcation flow against Murray's law (Murray 1926) and conservation
//! laws using a real FEM solve via `BifurcationSolver3D`.
//!
//! # Theorem — Junction Mass Conservation (continuity)
//!
//! For incompressible steady flow at any bifurcation junction,
//!
//! ```text
//! Q_parent = Q_daughter1 + Q_daughter2
//! ```
//!
//! where $Q_i = \int_{A_i} \mathbf{u} \cdot \hat{n} \, dA$.
//! This is validated as the primary correctness criterion of the FEM solve.
//!
//! # Theorem — Murray's Cube Law (Murray 1926)
//!
//! For minimum metabolic cost in Hagen–Poiseuille flow,
//!
//! ```text
//! D_parent^3 = D_daughter1^3 + D_daughter2^3
//! ```
//!
//! Relative deviation $|(D_0^3 - D_1^3 - D_2^3)| / D_0^3$ is computed and
//! must be ≤ 1 % for the benchmark geometry (D_p = 100 µm, D_d = 79.37 µm).
//!
//! **Reference**: Murray, C.D. "The Physiological Principle of Minimum Work",
//! *PNAS* 12(3), 1926, pp. 207–214.
use super::super::{Benchmark, BenchmarkConfig, BenchmarkResult};
use crate::geometry::threed::Bifurcation3D;
use cfd_3d::bifurcation::{BifurcationConfig3D, BifurcationGeometry3D, BifurcationSolver3D};
use cfd_core::conversion::SafeFromF64;
use cfd_core::physics::fluid::blood::CassonBlood;
use nalgebra::RealField;

/// 3D Bifurcation Flow benchmark
pub struct BifurcationFlow3D<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// The 3D bifurcation geometry (cfd-validation representation)
    pub geometry: Bifurcation3D<T>,
    /// The blood rheology model used for the FEM solve
    pub fluid: CassonBlood<T>,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy> BifurcationFlow3D<T> {
    /// Create a new 3D bifurcation flow benchmark
    pub fn new(geometry: Bifurcation3D<T>, fluid: CassonBlood<T>) -> Self {
        Self { geometry, fluid }
    }
}

impl<T: RealField + Copy + num_traits::Float + num_traits::FromPrimitive +
    num_traits::ToPrimitive + SafeFromF64 + std::convert::From<f64> + cfd_mesh::domain::core::Scalar> Benchmark<T> for BifurcationFlow3D<T> {
    fn name(&self) -> &'static str {
        "3D Y-Bifurcation Flow"
    }

    fn description(&self) -> &'static str {
        "3D Y-Bifurcation flow with non-Newtonian blood rheology (Casson). \
         Validates mass conservation at junction and Murray's cube law compliance \
         using a real FEM solve via BifurcationSolver3D."
    }

    fn run(&self, _config: &BenchmarkConfig<T>) -> cfd_core::error::Result<BenchmarkResult<T>> {
        let mut result = BenchmarkResult::new(self.name());

        // ── 1. Build BifurcationGeometry3D from the cfd-validation Bifurcation3D ──
        // D_parent = 100 µm, D_daughters = 79.37 µm → Murray-optimal (D₀³ ≈ 2·D₁³)
        let geo = BifurcationGeometry3D::<T>::symmetric(
            self.geometry.d_parent,
            self.geometry.d_daughters[0],
            self.geometry.l_parent,
            self.geometry.l_daughters[0],
            // Transition zone: 10 % of parent length
            self.geometry.l_parent * T::from_f64_or_one(0.1),
        );

        // ── 2. Configure the solver ──
        // Use small flow rate (1e-10 m³/s) to stay firmly in Stokes regime
        // (Re ≪ 1 for 100 µm tube), keeping the FEM problem well-conditioned.
        let config = BifurcationConfig3D::<T> {
            inlet_flow_rate: T::from_f64_or_one(1e-10),
            mesh_resolution: 4, // coarse mesh for benchmark speed
            ..BifurcationConfig3D::default()
        };

        // ── 3. Solve with Casson blood ──
        let solver = BifurcationSolver3D::new(geo.clone(), config);
        let solution = solver.solve(self.fluid.clone())?;

        // ── 4. Compute Murray's law deviation from geometry (analytical) ──
        let d_p = self.geometry.d_parent;
        let d_d1 = self.geometry.d_daughters[0];
        let d_d2 = self.geometry.d_daughters[1];
        let lhs = d_p * d_p * d_p;
        let rhs = d_d1 * d_d1 * d_d1 + d_d2 * d_d2 * d_d2;
        let murray_deviation = num_traits::Float::abs(lhs - rhs) / lhs;

        // ── 5. Record metrics ──
        result.metrics.insert("Murray Deviation".to_string(), murray_deviation);
        result.metrics.insert("Mass Conservation Error".to_string(), solution.mass_conservation_error);
        result.metrics.insert("Q Parent".to_string(), solution.q_parent);
        result.metrics.insert("Q Daughter1".to_string(), solution.q_daughter1);
        result.metrics.insert("Q Daughter2".to_string(), solution.q_daughter2);
        result.metrics.insert("Wall Shear Parent (Pa)".to_string(), solution.wall_shear_stress_parent);
        result.metrics.insert("dp Parent (Pa)".to_string(), solution.dp_parent);

        println!("Bifurcation 3D Validation:");
        println!("  Murray deviation:        {murray_deviation:?}");
        println!("  Mass conservation error: {:?}", solution.mass_conservation_error);
        println!("  Q_parent:                {:?}", solution.q_parent);
        println!("  Q_d1 + Q_d2:             {:?}", solution.q_daughter1 + solution.q_daughter2);
        println!("  Wall shear (parent, Pa): {:?}", solution.wall_shear_stress_parent);

        Ok(result)
    }

    fn reference_solution(&self) -> Option<BenchmarkResult<T>> {
        // Reference: purely geometric Murray's law (zero deviation for optimal geometry)
        let mut ref_result = BenchmarkResult::new(self.name());
        ref_result.metrics.insert("Murray Deviation".to_string(), T::zero());
        ref_result.metrics.insert("Mass Conservation Error".to_string(), T::zero());
        Some(ref_result)
    }

    fn validate(&self, result: &BenchmarkResult<T>) -> cfd_core::error::Result<bool> {
        // ── Criterion 1: Murray's law deviation < 1 % ──
        let murray_ok = result
            .metrics
            .get("Murray Deviation")
            .map(|&dev| num_traits::Float::abs(dev) < T::from_f64_or_one(0.01))
            .unwrap_or(false);

        // ── Criterion 2: Mass conservation error < 2 % (coarse-mesh tolerance) ──
        let mass_ok = result
            .metrics
            .get("Mass Conservation Error")
            .map(|&err| num_traits::Float::abs(err) < T::from_f64_or_one(0.02))
            .unwrap_or(false);

        Ok(murray_ok && mass_ok)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// Symmetric bifurcation with Murray-optimal daughters:
    ///   D_parent = 100 µm, D_daughter = 100 / 2^(1/3) ≈ 79.37 µm
    ///   → D_p^3 = 2 · D_d^3  (exact Murray)
    #[test]
    fn test_bifurcation_flow_3d_murray_and_mass() {
        let d_parent = 100e-6_f64;
        // Murray-optimal daughter diameter: D_d = D_p / 2^(1/3)
        let d_daughter = d_parent / (2.0_f64.powf(1.0 / 3.0));
        let angle = PI / 6.0; // 30° half-angle

        let geometry = Bifurcation3D::symmetric(
            d_parent,   // d_parent
            d_daughter, // d_daughter (Murray-optimal)
            500e-6,     // l_parent: 500 µm
            500e-6,     // l_daughter: 500 µm
            angle,      // branching angle
        );
        let fluid = CassonBlood::<f64>::normal_blood();
        let benchmark = BifurcationFlow3D::new(geometry, fluid);
        let config = BenchmarkConfig::default();

        let result = benchmark.run(&config).expect("3D bifurcation solve failed");
        assert!(
            benchmark.validate(&result).expect("validation error"),
            "Validation failed: Murray={:?}, MassErr={:?}",
            result.metrics.get("Murray Deviation"),
            result.metrics.get("Mass Conservation Error"),
        );
    }
}
