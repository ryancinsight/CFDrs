use super::super::traits::LESTurbulenceModel;
use leto::Array2;

use super::*;

fn create_test_fields(nx: usize, ny: usize) -> (Array2<f64>, Array2<f64>, Array2<f64>) {
    let mut velocity_u = Array2::zeros([nx, ny]);
    let velocity_v = Array2::zeros([nx, ny]);
    let pressure = Array2::zeros([nx, ny]);

    // Simple shear flow
    for i in 0..nx {
        for j in 0..ny {
            velocity_u[[i, j]] = (j as f64) * 0.1;
        }
    }

    (velocity_u, velocity_v, pressure)
}

#[test]
fn test_des_creation() {
    let config = DESConfig::default();
    let des = DetachedEddySimulation::new(10, 10, 0.1, 0.1, config, &[]);

    assert_eq!(des.sgs_viscosity.shape()[0], 10);
    assert_eq!(des.sgs_viscosity.shape()[1], 10);
    assert_eq!(des.des_length_scale.shape()[0], 10);
    assert_eq!(des.des_length_scale.shape()[1], 10);
}

#[test]
fn test_des_variants() {
    let variants = vec![DESVariant::DES97, DESVariant::DDES, DESVariant::IDDES];

    for variant in variants {
        let config = DESConfig {
            variant,
            ..Default::default()
        };
        let des = DetachedEddySimulation::new(10, 10, 0.1, 0.1, config, &[]);

        let name = des.get_model_name();
        assert!(name.contains("DES"));
    }
}

#[test]
fn test_des_length_scale_computation() {
    let config = DESConfig {
        variant: DESVariant::DES97, // Use simpler DES97 for testing
        ..DESConfig::default()
    };
    let des = DetachedEddySimulation::new(10, 10, 0.1, 0.1, config, &[]);

    // Create dummy velocity fields
    let velocity_u = Array2::from_elem([10, 10], 1.0);
    let velocity_v = Array2::from_elem([10, 10], 0.5);

    let length_scale = des.compute_des_length_scale(&velocity_u, &velocity_v, 0.1, 0.1);

    // Check dimensions
    assert_eq!(length_scale.shape()[0], 10);
    assert_eq!(length_scale.shape()[1], 10);

    // Length scale should be positive and reasonable
    assert!(length_scale.iter().all(|&l| l > 0.0 && l.is_finite()));
}

#[test]
fn test_des_model_update() {
    let config = DESConfig::default();
    let mut des = DetachedEddySimulation::new(10, 10, 0.1, 0.1, config, &[]);
    let (velocity_u, velocity_v, pressure) = create_test_fields(10, 10);

    let result = des.update(
        &velocity_u,
        &velocity_v,
        &pressure,
        1.0,
        0.01,
        0.001,
        0.1,
        0.1,
    );
    result
        .as_ref()
        .expect("the DES eddy-viscosity update must succeed");
}

#[test]
fn test_des_mode_detection() {
    let config = DESConfig::default();
    let des = DetachedEddySimulation::new(10, 10, 0.1, 0.1, config, &[]);

    // Test a point - should work without panicking
    let is_les = des.is_les_mode(5, 5);
    // Result depends on computed fields, just check it doesn't crash
    let _ = is_les;
}

#[test]
fn test_des_model_constants() {
    let config = DESConfig {
        des_constant: 0.7,
        max_sgs_ratio: 0.6,
        ..Default::default()
    };
    let des = DetachedEddySimulation::new(10, 10, 0.1, 0.1, config, &[]);

    let constants = des.get_model_constants();
    assert!(!constants.is_empty());

    // Should include DES-specific constants
    assert!(constants.iter().any(|(name, _)| *name == "DES Constant"));
    assert!(constants.iter().any(|(name, _)| *name == "Max SGS Ratio"));
}

#[test]
fn test_iddes_length_scale_computation() {
    let config = DESConfig {
        variant: DESVariant::IDDES,
        ..Default::default()
    };
    let des = DetachedEddySimulation::new(10, 10, 0.1, 0.1, config, &[]);

    let velocity_u = Array2::from_elem([10, 10], 1.0);
    let velocity_v = Array2::from_elem([10, 10], 0.5);

    let length_scale = des.compute_des_length_scale(&velocity_u, &velocity_v, 0.1, 0.1);

    assert_eq!(length_scale.shape()[0], 10);
    assert_eq!(length_scale.shape()[1], 10);
    assert!(length_scale.iter().all(|&l| l > 0.0 && l.is_finite()));
}

#[test]
fn test_ddes_shielding_behavior() {
    // Setup DDES
    let config = DESConfig {
        variant: DESVariant::DDES,
        des_constant: 0.65,
        ..Default::default()
    };

    let mut des = DetachedEddySimulation::new(10, 10, 1.0, 1.0, config, &[]);

    // Manually overwrite wall distance for a specific point (5, 5)
    // Set d_w = 1.0.
    des.wall_distance[[5, 5]] = 1.0;

    // l_RANS = 1.0.
    // l_LES = 0.65 * 1.0 = 0.65.
    // l_RANS > l_LES, so shielding logic enters.

    // Case 1: High Eddy Viscosity (Shielding Active -> RANS mode)
    // nu_tilde needs to be high enough.
    des.nu_tilde[[5, 5]] = 10.0;

    // Create velocity field that gives S=1.0 at (5,5).
    // Simple shear: u = y. du/dy = 1.
    let mut velocity_u = Array2::zeros([10, 10]);
    for j in 0..10 {
        for i in 0..10 {
            velocity_u[[i, j]] = j as f64;
        }
    }
    let velocity_v = Array2::zeros([10, 10]);

    let length_scale = des.compute_des_length_scale(&velocity_u, &velocity_v, 1.0, 1.0);

    let l = length_scale[[5, 5]];
    // Expect RANS length scale (1.0) because of shielding
    assert!(
        (l - 1.0).abs() < 1e-3,
        "With high viscosity, DDES should shield and return RANS length (1.0), got {l}"
    );

    // Case 2: Zero Eddy Viscosity (Shielding Inactive -> LES mode)
    des.nu_tilde[[5, 5]] = 0.0;
    let length_scale_les = des.compute_des_length_scale(&velocity_u, &velocity_v, 1.0, 1.0);
    let l_les = length_scale_les[[5, 5]];

    // Expect LES length scale (0.65)
    assert!(
        (l_les - 0.65).abs() < 1e-3,
        "With zero viscosity, DDES should return LES length (0.65), got {l_les}"
    );
}

/// **Positive**: `try_new` accepts valid arguments.
#[test]
fn des_try_new_accepts_valid_arguments() {
    let des = DetachedEddySimulation::try_new(8, 12, 1e-3, 1e-3, DESConfig::default(), &[])
        .expect("valid must succeed");
    assert_eq!(des.dx, 1e-3);
    assert_eq!(des.dy, 1e-3);
}

/// **Adversarial**: zero `dx` is rejected.
#[test]
fn des_try_new_rejects_zero_dx() {
    match DetachedEddySimulation::try_new(8, 12, 0.0, 1e-3, DESConfig::default(), &[]) {
        Err(e) => assert!(e.to_string().contains("dx"), "error must mention dx: {e}"),
        Ok(_) => panic!("zero dx must be rejected"),
    }
}

/// **Adversarial**: NaN `dy` is rejected.
#[test]
fn des_try_new_rejects_nan_dy() {
    match DetachedEddySimulation::try_new(8, 12, 1e-3, f64::NAN, DESConfig::default(), &[]) {
        Err(e) => assert!(e.to_string().contains("dy"), "error must mention dy: {e}"),
        Ok(_) => panic!("NaN dy must be rejected"),
    }
}

/// **Adversarial**: zero `des_constant` is rejected.
#[test]
fn des_try_new_rejects_zero_des_constant() {
    let cfg = DESConfig {
        des_constant: 0.0,
        ..DESConfig::default()
    };
    match DetachedEddySimulation::try_new(8, 12, 1e-3, 1e-3, cfg, &[]) {
        Err(e) => assert!(
            e.to_string().contains("des_constant"),
            "error must mention des_constant: {e}"
        ),
        Ok(_) => panic!("zero des_constant must be rejected"),
    }
}

/// **Boundary**: `new` panics on invalid `dx` (thin wrapper contract).
#[test]
#[should_panic(expected = "dx")]
fn des_new_panics_on_invalid_dx() {
    let _ = DetachedEddySimulation::new(8, 12, 0.0, 1e-3, DESConfig::default(), &[]);
}
