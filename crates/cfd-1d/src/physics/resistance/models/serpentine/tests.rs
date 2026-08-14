#![expect(
    clippy::unwrap_used,
    reason = "ratchet CFDRS-UNWRAP-1: pre-existing debt"
)]

use super::traits::{FlowConditions, ResistanceModel};
use super::*;
use cfd_core::physics::fluid::FluidTrait;
use eunomia::assert_relative_eq;

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

    let enhancement = model.curvature_enhancement(5.0);
    assert!(
        (enhancement - 1.0).abs() < 0.005,
        "Enhancement at De=5 should be virtually 1.0, got {enhancement}"
    );
}

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

    assert_relative_eq!(model.curvature_enhancement(0.0_f64), 1.0, epsilon = 1e-9);

    let e100 = model.curvature_enhancement(100.0_f64);
    assert_relative_eq!(e100, 1.033, epsilon = 0.002);

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

    let enhancement = model.curvature_enhancement(100.0);
    assert!(enhancement > 1.0);
    assert!(enhancement < 5.0);
}

#[test]
fn test_shah_london_square() {
    let cs = SerpentineCrossSection::Rectangular {
        width: 0.001,
        height: 0.001,
    };
    let factor = cs.shah_london_fre_factor();
    assert_relative_eq!(factor, 0.8894, epsilon = 0.01);
}

#[test]
fn test_bend_loss_coefficient() {
    let k = BendType::Sharp.loss_coefficient(100.0_f64);
    assert_relative_eq!(k, 4.7, epsilon = 1e-6);

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
    let conditions = FlowConditions::new(0.1);

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
    assert_eq!(analysis.num_bends, 4);

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

    assert!(r10 > r3, "More bends should increase resistance");

    Ok(())
}

#[test]
fn serpentine_coefficients_are_invariant_under_reverse_flow() -> cfd_core::error::Result<()> {
    let fluid = water();
    let model = SerpentineModel {
        straight_length: 0.02_f64,
        num_segments: 6,
        cross_section: SerpentineCrossSection::Rectangular {
            width: 1.2e-3,
            height: 0.45e-3,
        },
        bend_radius: 1.5e-3,
        bend_type: BendType::Smooth {
            radius_to_dh_ratio: 3.0,
        },
    };

    let forward = FlowConditions::new(0.08);
    let reverse = FlowConditions::new(-0.08);
    let (r_forward, k_forward) = model.calculate_coefficients(&fluid, &forward)?;
    let (r_reverse, k_reverse) = model.calculate_coefficients(&fluid, &reverse)?;
    let resistance_forward = model.calculate_resistance(&fluid, &forward)?;
    let resistance_reverse = model.calculate_resistance(&fluid, &reverse)?;

    assert_relative_eq!(r_forward, r_reverse, max_relative = 1e-12);
    assert_relative_eq!(k_forward, k_reverse, max_relative = 1e-12);
    assert_relative_eq!(resistance_forward, resistance_reverse, max_relative = 1e-12);

    Ok(())
}

#[test]
fn serpentine_analysis_uses_flow_magnitude_for_scalar_losses() -> cfd_core::error::Result<()> {
    let fluid = water();
    let model = SerpentineModel {
        straight_length: 0.02_f64,
        num_segments: 6,
        cross_section: SerpentineCrossSection::Circular { diameter: 0.8e-3 },
        bend_radius: 2.0e-3,
        bend_type: BendType::Smooth {
            radius_to_dh_ratio: 2.5,
        },
    };

    let forward = model.analyze(&fluid, &FlowConditions::new(0.12))?;
    let reverse = model.analyze(&fluid, &FlowConditions::new(-0.12))?;

    assert_relative_eq!(forward.reynolds, reverse.reynolds, max_relative = 1e-12);
    assert_relative_eq!(
        forward.wall_shear_rate,
        reverse.wall_shear_rate,
        max_relative = 1e-12
    );
    assert_relative_eq!(
        forward.dean_number,
        reverse.dean_number,
        max_relative = 1e-12
    );
    assert_relative_eq!(forward.dp_total, reverse.dp_total, max_relative = 1e-12);
    assert!(reverse.reynolds > 0.0);
    assert!(reverse.wall_shear_rate > 0.0);
    assert!(reverse.dp_total > 0.0);

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

#[test]
fn test_bayat_rezai_no_curvature() {
    let e = bayat_rezai_enhancement(0.0);
    assert_relative_eq!(e, 1.0, epsilon = 1e-12);
}

#[test]
fn test_bayat_rezai_moderate_dean() {
    let de = 50.0_f64;
    let e = bayat_rezai_enhancement(de);
    let expected = 1.0 + 0.085 * de.powf(0.48);
    assert_relative_eq!(e, expected, max_relative = 1e-10);
    assert!(
        e > 1.4 && e < 1.7,
        "Enhancement at De=50 should be ~1.55, got {e}"
    );
}

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
        assert_relative_eq!(from_method, from_standalone, epsilon = 1e-15);
    }
}
