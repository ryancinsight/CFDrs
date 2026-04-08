//! Venturi tube resistance model for converging-diverging channels.
//!
//! ## Venturi Effect
//!
//! **Theorem**: In a converging-diverging (Venturi) tube, the pressure drop
//! between the upstream section (1) and the throat section (2) is given by:
//!
//! ΔP₁₂ = (ρ V₂²/2)(1 − β⁴) / C_d²
//!
//! where β = D₂/D₁ is the diameter ratio.
//!
//! **Proof sketch**: continuity gives $V_2 = V_1 / \beta^2$ for incompressible
//! flow, and Bernoulli between the upstream section and the throat gives the
//! ideal acceleration pressure drop $\frac{1}{2}\rho V_2^2 (1 - \beta^4)$. The
//! discharge coefficient $C_d$ accounts for contraction losses and viscous
//! non-idealities, yielding the corrected denominator $C_d^2$ used by the
//! reduced-order resistance model.
//!
//! ## Theorem - Diffuser Loss Monotonicity
//!
//! For the calibrated gradual-expansion map in this module, the diffuser loss
//! coefficient $K_{exp}(\theta)$ is monotone nondecreasing with half-angle
//! $\theta$, and the idealized recovery efficiency $\eta_r = 1 - K_{exp}$ is
//! therefore monotone nonincreasing.
//!
//! **Proof sketch**: the implementation is a piecewise-constant calibration
//! table ordered by increasing angle bands. Each successive band assigns an
//! equal or larger loss coefficient than the previous band, so $K_{exp}$ cannot
//! decrease as the diffuser becomes more abrupt.
//!
//! ## Module structure
//!
//! | Module | Contents |
//! |--------|----------|
//! | [`model`] | `VenturiModel<T>` struct, `ResistanceModel<T>` impl |
//! | [`analysis`] | `VenturiAnalysis<T>`, `analyze()` |
//!
//! ## References
//!
//! - ISO 5167-4:2003. Venturi tubes.
//! - Idelchik, I. E. (2007). *Handbook of Hydraulic Resistance* (4th ed.), §5.1-5.8.
//! - Reader-Harris, M. (2015). *Orifice Plates and Venturi Tubes.* Springer.

pub mod analysis;
pub mod model;

// Re-import traits needed by child modules
use super::traits;
use serde::{Deserialize, Serialize};

// ── Constants ───────────────────────────────────────────────────────────────

pub(crate) const LAMINAR_LIMIT_RE: f64 = 2300.0;
pub(crate) const LAMINAR_FRICTION_COEFF: f64 = 64.0;
pub(crate) const BLASIUS_COEFF: f64 = 0.3164;
pub(crate) const BLASIUS_EXP: f64 = 0.25;

// ── Enums ───────────────────────────────────────────────────────────────────

/// Venturi tube geometry specification
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum VenturiGeometry {
    /// Machined convergent section (ISO 5167-4, C_d ≈ 0.995)
    MachinedConvergent,
    /// Rough-cast convergent section (ISO 5167-4, C_d ≈ 0.984)
    RoughCastConvergent,
    /// Rough-welded convergent section (ISO 5167-4, C_d ≈ 0.985)
    RoughWeldedConvergent,
    /// Custom discharge coefficient
    Custom {
        /// Discharge coefficient C_d
        discharge_coefficient: f64,
    },
}

impl VenturiGeometry {
    /// Get the nominal discharge coefficient for this geometry type
    #[must_use]
    pub fn discharge_coefficient(&self) -> f64 {
        match self {
            Self::MachinedConvergent => 0.995,
            Self::RoughCastConvergent => 0.984,
            Self::RoughWeldedConvergent => 0.985,
            Self::Custom {
                discharge_coefficient,
            } => *discharge_coefficient,
        }
    }
}

/// Expansion geometry type for the diverging section
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum ExpansionType {
    /// Sudden expansion (Borda-Carnot, K_exp = 1.0)
    Sudden,
    /// Gradual expansion with specified half-angle [degrees]
    Gradual {
        /// Diffuser half-angle in degrees
        half_angle_deg: f64,
    },
}

impl ExpansionType {
    /// Get the expansion loss coefficient K_exp
    ///
    /// Based on Idelchik (2007) Table 5-1 for conical diffusers
    #[must_use]
    pub fn loss_coefficient(&self) -> f64 {
        match self {
            Self::Sudden => 1.0,
            Self::Gradual { half_angle_deg } => {
                let angle = *half_angle_deg;
                if angle <= 3.0 {
                    0.10 // Optimal diffuser
                } else if angle <= 5.0 {
                    0.14
                } else if angle <= 7.0 {
                    0.20
                } else if angle <= 10.0 {
                    0.30
                } else if angle <= 15.0 {
                    0.45
                } else if angle <= 20.0 {
                    0.55
                } else if angle <= 30.0 {
                    0.70
                } else if angle <= 45.0 {
                    0.85
                } else {
                    1.0 // Approaches sudden expansion
                }
            }
        }
    }

    /// Get pressure recovery efficiency η_r = 1 - K_exp
    #[must_use]
    pub fn recovery_efficiency(&self) -> f64 {
        1.0 - self.loss_coefficient()
    }
}

// ── Re-exports ──────────────────────────────────────────────────────────────

pub use analysis::VenturiAnalysis;
pub use model::VenturiModel;

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use proptest::prelude::*;
    use traits::{FlowConditions, ResistanceModel};

    proptest! {
        #[test]
        fn gradual_expansion_loss_is_monotone_nondecreasing(
            small_half_angle_deg in 0.0_f64..89.0,
            large_half_angle_deg in 0.0_f64..89.0,
        ) {
            prop_assume!(small_half_angle_deg <= large_half_angle_deg);

            let small = ExpansionType::Gradual {
                half_angle_deg: small_half_angle_deg,
            }
            .loss_coefficient();
            let large = ExpansionType::Gradual {
                half_angle_deg: large_half_angle_deg,
            }
            .loss_coefficient();

            prop_assert!(small <= large);
            prop_assert!((0.10..=1.0).contains(&small));
            prop_assert!((0.10..=1.0).contains(&large));
        }
    }

    proptest! {
        #[test]
        fn gradual_expansion_recovery_is_monotone_nonincreasing(
            small_half_angle_deg in 0.0_f64..89.0,
            large_half_angle_deg in 0.0_f64..89.0,
        ) {
            prop_assume!(small_half_angle_deg <= large_half_angle_deg);

            let small = ExpansionType::Gradual {
                half_angle_deg: small_half_angle_deg,
            }
            .recovery_efficiency();
            let large = ExpansionType::Gradual {
                half_angle_deg: large_half_angle_deg,
            }
            .recovery_efficiency();

            prop_assert!(small >= large);
            prop_assert!((0.0..=0.90).contains(&small));
            prop_assert!((0.0..=0.90).contains(&large));
        }
    }

    #[test]
    fn test_venturi_bernoulli_limit() -> cfd_core::error::Result<()> {
        let model = VenturiModel::symmetric(
            0.01,  // 10mm inlet
            0.005, // 5mm throat (β = 0.5)
            0.01,  // 10mm throat length
            0.05,  // 50mm total
        )
        .with_geometry(VenturiGeometry::Custom {
            discharge_coefficient: 1.0,
        })
        .with_expansion(ExpansionType::Gradual {
            half_angle_deg: 3.0,
        });

        let fluid = cfd_core::physics::fluid::database::water_20c::<f64>()?;

        let velocity = 2.0; // m/s
        let mut conditions = FlowConditions::new(velocity);
        conditions.reynolds_number = Some(998.2 * velocity * 0.01 / 1.002e-3);

        let analysis = model.analyze(&fluid, &conditions)?;

        // V_throat = V_inlet * (D_inlet/D_throat)² = 2 * 4 = 8 m/s
        assert_relative_eq!(analysis.throat_velocity, 8.0, epsilon = 0.01);

        // β⁴ = 0.0625 → (1 − β⁴) = 0.9375
        // Ideal Bernoulli: ΔP = ½ * 998.2 * 64 * 0.9375 ≈ 29,946 Pa
        let dp_ideal = 0.5 * 998.2 * 64.0 * 0.9375;
        assert_relative_eq!(analysis.dp_contraction, dp_ideal, epsilon = dp_ideal * 0.05);

        Ok(())
    }

    #[test]
    fn test_venturi_millifluidic_blood() -> cfd_core::error::Result<()> {
        let model = VenturiModel::millifluidic(
            0.001,  // 1mm inlet
            0.0005, // 0.5mm throat
            0.002,  // 2mm throat length
        );

        let fluid = cfd_core::physics::fluid::database::water_20c::<f64>()?;

        let mut conditions = FlowConditions::new(0.01); // 10 mm/s
        let density = 998.2;
        let viscosity = 1.002e-3;
        let re = density * 0.01 * 0.001 / viscosity;
        conditions.reynolds_number = Some(re);

        let analysis = model.analyze(&fluid, &conditions)?;

        assert!(analysis.dp_friction > 0.0);
        assert!(analysis.dp_contraction > 0.0);
        assert!(analysis.dp_total > 0.0);

        // V_throat = V_inlet * (A_inlet/A_throat) = 0.01 * 4 = 0.04 m/s
        assert_relative_eq!(analysis.throat_velocity, 0.04, epsilon = 0.001);

        Ok(())
    }

    #[test]
    fn test_venturi_area_ratio() {
        let model = VenturiModel::<f64>::symmetric(0.01, 0.005, 0.01, 0.05);

        // β² = (D_t/D_i)² = (0.005/0.01)² = 0.25
        assert_relative_eq!(model.beta_squared(), 0.25, epsilon = 1e-10);

        let a_inlet = std::f64::consts::PI * 0.01_f64.powi(2) / 4.0;
        assert_relative_eq!(model.inlet_area(), a_inlet, epsilon = 1e-10);

        let a_throat = std::f64::consts::PI * 0.005_f64.powi(2) / 4.0;
        assert_relative_eq!(model.throat_area(), a_throat, epsilon = 1e-10);
    }

    #[test]
    fn test_venturi_friction_factor() {
        let model = VenturiModel::<f64>::symmetric(0.01, 0.005, 0.01, 0.05);

        // Laminar: f = 64/Re
        let f_lam = model.throat_friction_factor(100.0);
        assert_relative_eq!(f_lam, 0.64, epsilon = 1e-6);

        // Blasius: f = 0.3164 / Re^0.25
        let f_turb = model.throat_friction_factor(10000.0);
        let expected = 0.3164 / 10000.0_f64.powf(0.25);
        assert_relative_eq!(f_turb, expected, epsilon = 1e-6);
    }

    #[test]
    fn test_venturi_recovery_reynolds_correction() {
        let model = VenturiModel::<f64>::symmetric(0.01, 0.005, 0.01, 0.05)
            .with_expansion(ExpansionType::Gradual { half_angle_deg: 7.0 });

        assert_relative_eq!(
            VenturiModel::<f64>::diffuser_recovery_reynolds_correction(5000.0),
            1.0,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            VenturiModel::<f64>::diffuser_recovery_reynolds_correction(200.0),
            0.6,
            epsilon = 1e-12
        );
        assert!(model.effective_recovery_efficiency(150.0) < model.effective_recovery_efficiency(5000.0));
    }

    #[test]
    fn test_venturi_expansion_coefficients() {
        assert_relative_eq!(
            ExpansionType::Sudden.loss_coefficient(),
            1.0,
            epsilon = 1e-10
        );
        assert_relative_eq!(
            ExpansionType::Gradual {
                half_angle_deg: 3.0
            }
            .loss_coefficient(),
            0.10,
            epsilon = 1e-10
        );
        assert_relative_eq!(
            ExpansionType::Gradual {
                half_angle_deg: 7.0
            }
            .loss_coefficient(),
            0.20,
            epsilon = 1e-10
        );
    }

    #[test]
    fn test_venturi_validate_invariants() {
        let model = VenturiModel::<f64>::symmetric(0.01, 0.005, 0.01, 0.05);
        let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().unwrap();
        let conditions = FlowConditions::new(0.1);

        assert!(model.validate_invariants(&fluid, &conditions).is_ok());

        // Invalid: throat >= inlet
        let bad_model = VenturiModel::<f64>::symmetric(0.005, 0.01, 0.01, 0.05);
        assert!(bad_model.validate_invariants(&fluid, &conditions).is_err());
    }

    #[test]
    fn test_venturi_resistance_model_trait() -> cfd_core::error::Result<()> {
        let model = VenturiModel::symmetric(0.01_f64, 0.005, 0.01, 0.05);
        let fluid = cfd_core::physics::fluid::database::water_20c::<f64>()?;

        let mut conditions = FlowConditions::new(0.1);
        conditions.reynolds_number = Some(500.0);

        let resistance = model.calculate_resistance(&fluid, &conditions)?;
        assert!(resistance > 0.0);

        let (r, k) = model.calculate_coefficients(&fluid, &conditions)?;
        assert!(r >= 0.0 || k >= 0.0); // At least one should be positive

        Ok(())
    }

    #[test]
    fn test_iso5167_discharge_coefficients() {
        assert_relative_eq!(
            VenturiGeometry::MachinedConvergent.discharge_coefficient(),
            0.995,
            epsilon = 1e-10
        );
        assert_relative_eq!(
            VenturiGeometry::RoughCastConvergent.discharge_coefficient(),
            0.984,
            epsilon = 1e-10
        );
        assert_relative_eq!(
            VenturiGeometry::RoughWeldedConvergent.discharge_coefficient(),
            0.985,
            epsilon = 1e-10
        );
    }
}
