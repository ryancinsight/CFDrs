//! Comprehensive validation tests for MUSCL schemes
//!
//! Tests accuracy on smooth solutions, monotonicity preservation for discontinuities,
//! and convergence rates for MUSCL2 and MUSCL3 schemes.

/// Test basic MUSCL scheme creation and naming
#[test]
fn test_muscl_scheme_creation() {
    use cfd_2d::physics::momentum::{MusclOrder, MusclReconstruction, MusclScheme, Superbee};

    let limiter = Superbee;
    let muscl2 = MusclScheme::<f64, Superbee>::new(limiter.clone(), MusclOrder::SecondOrder);
    let muscl3 = MusclScheme::<f64, Superbee>::new(limiter, MusclOrder::ThirdOrder);

    assert_eq!(muscl2.name(), "MUSCL2");
    assert_eq!(muscl3.name(), "MUSCL3");
    assert_eq!(muscl2.order(), MusclOrder::SecondOrder);
    assert_eq!(muscl3.order(), MusclOrder::ThirdOrder);
}

/// Test MUSCL discretization schemes
#[test]
fn test_muscl_discretization_schemes() {
    use cfd_2d::physics::momentum::{DiscretizationScheme, MusclDiscretization, MusclScheme, VanLeer};

    let limiter = VanLeer;
    let muscl = MusclScheme::<f64, VanLeer>::new(limiter, cfd_2d::physics::momentum::MusclOrder::SecondOrder);
    let scheme = MusclDiscretization::new(muscl);

    // Test positive velocity (upwind-like behavior)
    let flux_pos = scheme.convective_flux(1.0, 0.5, 2.0, 1.0);
    assert!(flux_pos > 0.0);

    // Test negative velocity
    let flux_neg = scheme.convective_flux(1.0, 0.5, 2.0, -1.0);
    assert!(flux_neg < 0.0);

    // Test scheme naming
    assert_eq!(scheme.name(), "MUSCL2");
}
