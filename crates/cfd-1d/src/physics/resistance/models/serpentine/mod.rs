//! Serpentine channel resistance model with Dean flow corrections.
//!
//! ## Dean Flow in Curved Channels
//!
//! **Dean Number**: When fluid flows through a curved channel, centrifugal forces
//! create secondary flow vortices (Dean vortices). The Dean number quantifies
//! the ratio of centrifugal to viscous forces:
//!
//! De = Re √(D_h / 2R_c)
//!
//! ### Friction Factor Enhancement
//!
//! The friction factor in curved channels is higher than in straight channels
//! due to the secondary flow. The enhancement factor depends on the Dean number:
//!
//! **Laminar regime** (De < 17): Negligible curvature effect
//! **Low Dean number** (17 ≤ De < 370): Ito (1959) or Dean (1928) perturbation
//! **High Dean number** (De ≥ 370): Ito (1959) boundary layer limit
//!
//! ### Bend Loss Coefficient
//!
//! Each 180° bend introduces an additional minor loss:
//! K_bend = C₁ + C₂/Re (Idelchik 2007, §6.2)
//!
//! ## Module structure
//!
//! | Module | Contents |
//! |--------|----------|
//! | [`model`] | `SerpentineModel<T>` struct, `ResistanceModel<T>` impl |
//! | [`analysis`] | `SerpentineAnalysis<T>`, `bayat_rezai_enhancement()` |
//!
//! ## References
//!
//! - Dean, W. R. (1928). *Proc. R. Soc. Lond. A*, 121(787), 402-420.
//! - Ito, H. (1959). *ASME J. Basic Eng.*, 81(2), 123-134.
//! - Shah & London (1978). *Laminar Flow Forced Convection in Ducts*.
//! - Bayat, P. & Rezai, P. (2017). *Sci. Rep.* 7:13655.
//! - Idelchik, I. E. (2007). *Handbook of Hydraulic Resistance* (4th ed.), §6.1-6.4.

pub mod analysis;
pub mod model;

// Re-import traits needed by this module's types
use super::traits;
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;
use serde::{Deserialize, Serialize};

// ── Enums ───────────────────────────────────────────────────────────────────

/// Bend type in the serpentine channel
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum BendType {
    /// Sharp 180° U-turn (no fillet radius)
    Sharp,
    /// Smooth 180° bend with specified R/D ratio
    Smooth {
        /// Ratio of bend radius to hydraulic diameter (R/D_h)
        radius_to_dh_ratio: f64,
    },
}

impl BendType {
    /// Minor loss coefficient for a single bend: K = C1 + C2/Re
    ///
    /// Returns (C1, C2) constants based on Idelchik (2007) §6.2
    #[must_use]
    pub fn loss_constants(&self) -> (f64, f64) {
        match self {
            Self::Sharp => (2.2, 250.0),
            Self::Smooth { radius_to_dh_ratio } => {
                let r = *radius_to_dh_ratio;
                if r <= 1.5 {
                    (0.9, 200.0)
                } else if r <= 3.0 {
                    (0.4, 100.0)
                } else if r <= 5.0 {
                    (0.3, 75.0)
                } else {
                    (0.2, 50.0)
                }
            }
        }
    }

    /// Compute bend loss coefficient K at the given Reynolds number
    pub fn loss_coefficient<T: RealField + Copy + FromPrimitive>(&self, reynolds: T) -> T {
        let (c1, c2) = self.loss_constants();
        let c1_t = T::from_f64(c1).expect("Mathematical constant conversion compromised");
        let c2_t = T::from_f64(c2).expect("Mathematical constant conversion compromised");
        c1_t + c2_t / reynolds
    }
}

/// Cross-section type for the serpentine channel
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum SerpentineCrossSection {
    /// Circular cross-section with given diameter [m]
    Circular {
        /// Channel diameter [m]
        diameter: f64,
    },
    /// Rectangular cross-section with width × height
    Rectangular {
        /// Channel width [m]
        width: f64,
        /// Channel height (depth) [m]
        height: f64,
    },
}

impl SerpentineCrossSection {
    /// Hydraulic diameter D_h
    #[must_use]
    pub fn hydraulic_diameter(&self) -> f64 {
        match self {
            Self::Circular { diameter } => *diameter,
            Self::Rectangular { width, height } => 2.0 * width * height / (width + height),
        }
    }

    /// Cross-sectional area [m²]
    #[must_use]
    pub fn area(&self) -> f64 {
        match self {
            Self::Circular { diameter } => std::f64::consts::PI * diameter * diameter / 4.0,
            Self::Rectangular { width, height } => width * height,
        }
    }

    /// Aspect ratio for rectangular channels (max/min dimension)
    #[must_use]
    pub fn aspect_ratio(&self) -> f64 {
        match self {
            Self::Circular { .. } => 1.0,
            Self::Rectangular { width, height } => {
                let (a, b) = if width > height {
                    (*width, *height)
                } else {
                    (*height, *width)
                };
                a / b
            }
        }
    }

    /// Shah-London correction factor for rectangular channels
    ///
    /// f·Re product for fully developed laminar flow in rectangular ducts
    /// compared to circular pipes:
    ///
    /// f·Re = C(α) where α is the aspect ratio
    ///
    /// Reference: Shah & London (1978), Table 43
    #[must_use]
    pub fn shah_london_fre_factor(&self) -> f64 {
        match self {
            Self::Circular { .. } => 1.0,
            Self::Rectangular { width, height } => {
                let alpha = {
                    let (a, b) = if width > height {
                        (*width, *height)
                    } else {
                        (*height, *width)
                    };
                    // Use the smaller/larger ratio so alpha ∈ (0, 1]
                    b / a
                };
                // Shah-London polynomial fit for fRe_rect / fRe_circular
                // fRe_circular = 64, fRe_rect(α) = 96(1 - 1.3553α + 1.9467α² - 1.7012α³ + 0.9564α⁴ - 0.2537α⁵)
                // Correction factor = fRe_rect / 64
                let fre_rect = 96.0
                    * (1.0 - 1.3553 * alpha + 1.9467 * alpha.powi(2) - 1.7012 * alpha.powi(3)
                        + 0.9564 * alpha.powi(4)
                        - 0.2537 * alpha.powi(5));
                fre_rect / 64.0
            }
        }
    }
}

// ── Re-exports ──────────────────────────────────────────────────────────────

pub use analysis::{bayat_rezai_enhancement, SerpentineAnalysis};
pub use model::SerpentineModel;

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use cfd_core::physics::fluid::FluidTrait;
    use traits::{FlowConditions, ResistanceModel};

    fn water() -> impl FluidTrait<f64> {
        cfd_core::physics::fluid::database::water_20c::<f64>().unwrap()
    }

    #[test]
    fn test_dean_number_calculation() {
        let model = SerpentineModel {
            straight_length: 0.02,
            num_segments: 5,
            cross_section: SerpentineCrossSection::Circular { diameter: 0.001 },
            bend_radius: 0.005,
            bend_type: BendType::Smooth {
                radius_to_dh_ratio: 5.0,
            },
        };

        // De = Re × √(D_h / (2 R_c)) = 100 × √(0.001 / 0.01) = 100 × √0.1 ≈ 31.6
        let dean = model.dean_number(100.0);
        let expected = 100.0 * (0.001 / (2.0 * 0.005_f64)).sqrt();
        assert_relative_eq!(dean, expected, epsilon = 1e-6);
    }

    #[test]
    fn test_curvature_enhancement_low_de() {
        let model = SerpentineModel::<f64> {
            straight_length: 0.02,
            num_segments: 5,
            cross_section: SerpentineCrossSection::Circular { diameter: 0.001 },
            bend_radius: 0.005,
            bend_type: BendType::Smooth {
                radius_to_dh_ratio: 5.0,
            },
        };

        // Exact Dean (1928) limit for De=5:
        // Enhancement should be effectively 1.0 at very low De.
        let enhancement = model.curvature_enhancement(5.0);
        assert!(
            (enhancement - 1.0).abs() < 0.005,
            "Enhancement at De=5 should be virtually 1.0, got {enhancement}"
        );
    }

    /// Validate curvature_enhancement against the exact Boundary Layer and Dean behaviors
    #[test]
    fn test_curvature_enhancement_exact_behaviors() {
        let model = SerpentineModel::<f64> {
            straight_length: 0.02,
            num_segments: 5,
            cross_section: SerpentineCrossSection::Circular { diameter: 0.001 },
            bend_radius: 0.005,
            bend_type: BendType::Smooth {
                radius_to_dh_ratio: 5.0,
            },
        };

        // De = 0 or negative → enhancement = 1.0 (no curvature)
        assert_relative_eq!(model.curvature_enhancement(0.0_f64), 1.0, epsilon = 1e-9);

        // De = 100: Transition to boundary layer exact limit mapping  = 0.1033 * 100.0^(1/2) = 1.033
        let e100 = model.curvature_enhancement(100.0_f64);
        assert_relative_eq!(e100, 1.033, epsilon = 0.002);

        // De = 370: Boundary layer limit
        let e370 = model.curvature_enhancement(370.0_f64);
        assert_relative_eq!(e370, 0.1033 * 370.0_f64.sqrt(), epsilon = 0.002);
    }

    #[test]
    fn test_curvature_enhancement_moderate_de() {
        let model = SerpentineModel::<f64> {
            straight_length: 0.02,
            num_segments: 5,
            cross_section: SerpentineCrossSection::Circular { diameter: 0.001 },
            bend_radius: 0.005,
            bend_type: BendType::Smooth {
                radius_to_dh_ratio: 5.0,
            },
        };

        // At De = 100, enhancement should be > 1.0 (secondary flow increases friction)
        let enhancement = model.curvature_enhancement(100.0);
        assert!(enhancement > 1.0);
        assert!(enhancement < 5.0);
    }

    #[test]
    fn test_shah_london_square() {
        // For square channel (α = 1): fRe ≈ 56.9, factor = 56.9/64 ≈ 0.889
        let cs = SerpentineCrossSection::Rectangular {
            width: 0.001,
            height: 0.001,
        };
        let factor = cs.shah_london_fre_factor();
        // Shah-London: 96(1 - 1.3553 + 1.9467 - 1.7012 + 0.9564 - 0.2537) / 64
        // = 96 × 0.5929 / 64 = 56.9184 / 64 ≈ 0.8894
        assert_relative_eq!(factor, 0.8894, epsilon = 0.01);
    }

    #[test]
    fn test_bend_loss_coefficient() {
        // Sharp bend at Re = 100: K = 2.2 + 250/100 = 4.7
        let k = BendType::Sharp.loss_coefficient(100.0_f64);
        assert_relative_eq!(k, 4.7, epsilon = 1e-6);

        // Smooth bend (R/D=5) at Re = 1000: K = 0.3 + 75/1000 = 0.375
        let k = BendType::Smooth {
            radius_to_dh_ratio: 5.0,
        }
        .loss_coefficient(1000.0_f64);
        assert_relative_eq!(k, 0.375, epsilon = 1e-6);
    }

    #[test]
    fn test_serpentine_resistance_positive() -> cfd_core::error::Result<()> {
        let model = SerpentineModel {
            straight_length: 0.02_f64,
            num_segments: 5,
            cross_section: SerpentineCrossSection::Circular { diameter: 0.001 },
            bend_radius: 0.005,
            bend_type: BendType::Smooth {
                radius_to_dh_ratio: 5.0,
            },
        };

        let fluid = water();
        let conditions = FlowConditions::new(0.1); // 100 mm/s

        let resistance = model.calculate_resistance(&fluid, &conditions)?;
        assert!(resistance > 0.0, "Resistance must be positive");

        let (r, k) = model.calculate_coefficients(&fluid, &conditions)?;
        assert!(r >= 0.0, "Linear coefficient must be non-negative");
        assert!(k >= 0.0, "Quadratic coefficient must be non-negative");

        Ok(())
    }

    #[test]
    fn test_serpentine_analysis() -> cfd_core::error::Result<()> {
        let model = SerpentineModel {
            straight_length: 0.02_f64,
            num_segments: 5,
            cross_section: SerpentineCrossSection::Circular { diameter: 0.001 },
            bend_radius: 0.005,
            bend_type: BendType::Smooth {
                radius_to_dh_ratio: 5.0,
            },
        };

        let fluid = water();
        let mut conditions = FlowConditions::new(0.1);
        conditions.reynolds_number = Some(100.0);

        let analysis = model.analyze(&fluid, &conditions)?;

        assert!(analysis.reynolds > 0.0);
        assert!(analysis.dean_number > 0.0);
        assert!(analysis.curvature_enhancement >= 1.0);
        assert!(analysis.dp_total > 0.0);
        assert_eq!(analysis.num_bends, 4); // 5 segments → 4 bends

        Ok(())
    }

    #[test]
    fn test_serpentine_more_bends_more_resistance() -> cfd_core::error::Result<()> {
        let fluid = water();
        let conditions = FlowConditions::new(0.1);

        let model_3 = SerpentineModel {
            straight_length: 0.02_f64,
            num_segments: 3,
            cross_section: SerpentineCrossSection::Circular { diameter: 0.001 },
            bend_radius: 0.005,
            bend_type: BendType::Smooth {
                radius_to_dh_ratio: 5.0,
            },
        };

        let model_10 = SerpentineModel {
            straight_length: 0.02_f64,
            num_segments: 10,
            cross_section: SerpentineCrossSection::Circular { diameter: 0.001 },
            bend_radius: 0.005,
            bend_type: BendType::Smooth {
                radius_to_dh_ratio: 5.0,
            },
        };

        let r3 = model_3.calculate_resistance(&fluid, &conditions)?;
        let r10 = model_10.calculate_resistance(&fluid, &conditions)?;

        // More bends → more resistance (same straight length, but more minor losses)
        assert!(r10 > r3, "More bends should increase resistance");

        Ok(())
    }

    #[test]
    fn test_serpentine_validate_invariants() {
        let fluid = water();
        let conditions = FlowConditions::new(0.1);

        let good_model = SerpentineModel {
            straight_length: 0.02_f64,
            num_segments: 5,
            cross_section: SerpentineCrossSection::Circular { diameter: 0.001 },
            bend_radius: 0.005,
            bend_type: BendType::Smooth {
                radius_to_dh_ratio: 5.0,
            },
        };
        assert!(good_model.validate_invariants(&fluid, &conditions).is_ok());

        let bad_model = SerpentineModel {
            straight_length: -0.01_f64,
            num_segments: 5,
            cross_section: SerpentineCrossSection::Circular { diameter: 0.001 },
            bend_radius: 0.005,
            bend_type: BendType::Smooth {
                radius_to_dh_ratio: 5.0,
            },
        };
        assert!(bad_model.validate_invariants(&fluid, &conditions).is_err());
    }

    // ─── Bayat & Rezai (2017) enhancement tests ─────────────────────────

    /// De = 0 gives enhancement = 1.0 (no curvature effect).
    #[test]
    fn test_bayat_rezai_no_curvature() {
        let e = bayat_rezai_enhancement(0.0);
        assert_relative_eq!(e, 1.0, epsilon = 1e-12);
    }

    /// De = 50 gives enhancement ≈ 1 + 0.085·50^0.48 ≈ 1.55.
    #[test]
    fn test_bayat_rezai_moderate_dean() {
        let de = 50.0_f64;
        let e = bayat_rezai_enhancement(de);
        let expected = 1.0 + 0.085 * de.powf(0.48);
        assert_relative_eq!(e, expected, max_relative = 1e-10);
        // Sanity: should be around 1.55
        assert!(e > 1.4 && e < 1.7, "Enhancement at De=50 should be ~1.55, got {e}");
    }

    /// Compare Bayat & Rezai (2017) vs Ito (1959) at De = 20.
    ///
    /// For rectangular millifluidic channels, Bayat gives lower enhancement
    /// than the classical Ito circular-tube formula.
    #[test]
    fn test_bayat_rezai_vs_ito() {
        let de = 20.0_f64;
        let bayat = bayat_rezai_enhancement(de);
        let model = SerpentineModel::<f64> {
            straight_length: 0.02,
            num_segments: 5,
            cross_section: SerpentineCrossSection::Circular { diameter: 0.001 },
            bend_radius: 0.005,
            bend_type: BendType::Smooth {
                radius_to_dh_ratio: 5.0,
            },
        };
        let ito = model.curvature_enhancement(de);

        assert!(
            bayat > 1.0,
            "Bayat enhancement at De=20 should be > 1.0, got {bayat}"
        );
        assert!(
            (bayat - ito).abs() > 0.01,
            "Bayat ({bayat}) and Ito ({ito}) should differ at De=20"
        );
    }

    /// Verify that the millifluidic wrapper method gives the same result as the standalone function.
    #[test]
    fn test_millifluidic_enhancement_consistent_with_standalone() {
        let model = SerpentineModel::<f64> {
            straight_length: 0.02,
            num_segments: 5,
            cross_section: SerpentineCrossSection::Circular { diameter: 0.001 },
            bend_radius: 0.005,
            bend_type: BendType::Smooth {
                radius_to_dh_ratio: 5.0,
            },
        };

        for &de in &[0.0, 1.0, 10.0, 50.0, 100.0] {
            let from_method = model.curvature_enhancement_millifluidic(de);
            let from_standalone = bayat_rezai_enhancement(de);
            assert_relative_eq!(
                from_method, from_standalone,
                epsilon = 1e-15,
            );
        }
    }
}
