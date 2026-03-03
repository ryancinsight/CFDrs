//! Integration tests for `CascadeSolver3D` — validates multi-channel CIF
//! network 3D FEM solves with both Newtonian and non-Newtonian fluids.

use cfd_3d::cascade::{CascadeChannelSpec, CascadeConfig3D, CascadeSolver3D};
use cfd_core::physics::fluid::blood::CassonBlood;
use cfd_core::physics::fluid::ConstantPropertyFluid;

/// Helper: creates a CIF-like asymmetric bifurcation with a wide outer arm
/// and a narrow central arm that has a Venturi throat.
fn cif_two_channel_specs() -> Vec<CascadeChannelSpec> {
    vec![
        // Outer (wide) arm — no venturi, higher flow fraction
        CascadeChannelSpec {
            id: "outer".into(),
            length: 10.0e-3,
            width: 2.0e-3,
            height: 1.0e-3,
            flow_rate_m3_s: 5e-8,
            is_venturi_throat: false,
            throat_width: None,
            local_hematocrit: None,
        },
        // Central (narrow) arm — venturi throat for SDT
        CascadeChannelSpec {
            id: "center".into(),
            length: 10.0e-3,
            width: 1.0e-3,
            height: 1.0e-3,
            flow_rate_m3_s: 2e-8,
            is_venturi_throat: true,
            throat_width: Some(0.5e-3),
            local_hematocrit: None,
        },
    ]
}

/// Low-resolution config for fast CI tests.
fn ci_config() -> CascadeConfig3D {
    CascadeConfig3D {
        outlet_pressure: 0.0,
        resolution: (20, 5, 5),
        max_picard_iterations: 5,
        picard_tolerance: 1e-2,
    }
}

// ── Newtonian water tests ─────────────────────────────────────────────

#[test]
fn cascade_water_solves_two_channels() {
    let water = ConstantPropertyFluid::<f64>::water_20c().unwrap();
    let solver = CascadeSolver3D::new(ci_config(), water);
    let result = solver.solve(&cif_two_channel_specs()).unwrap();

    assert_eq!(result.channel_results.len(), 2);
    for cr in &result.channel_results {
        assert!(
            cr.wall_shear_mean_pa >= 0.0,
            "Channel {} has negative mean shear: {}",
            cr.channel_id,
            cr.wall_shear_mean_pa
        );
        assert!(
            cr.max_velocity > 0.0,
            "Channel {} has zero max velocity",
            cr.channel_id
        );
    }
}

#[test]
fn cascade_channels_produce_distinct_shear() {
    // At coarse resolution, both channels should produce distinct non-zero
    // wall shear stress values reflecting their different geometries.
    let water = ConstantPropertyFluid::<f64>::water_20c().unwrap();
    let solver = CascadeSolver3D::new(ci_config(), water);
    let result = solver.solve(&cif_two_channel_specs()).unwrap();

    let outer = result
        .channel_results
        .iter()
        .find(|r| r.channel_id == "outer")
        .unwrap();
    let center = result
        .channel_results
        .iter()
        .find(|r| r.channel_id == "center")
        .unwrap();

    // Both channels must produce non-zero wall shear.
    assert!(
        outer.wall_shear_max_pa > 0.0,
        "Outer channel has zero max shear"
    );
    assert!(
        center.wall_shear_max_pa > 0.0,
        "Center channel has zero max shear"
    );
    // The shear values should differ (different geometry).
    assert!(
        (outer.wall_shear_max_pa - center.wall_shear_max_pa).abs() > 1e-10,
        "Both channels report identical shear — suspicious"
    );
    // The max_shear_channel_id must reference a valid channel.
    assert!(
        result.max_shear_channel_id == "outer"
            || result.max_shear_channel_id == "center",
        "max_shear_channel_id is '{}'",
        result.max_shear_channel_id
    );
}

// ── Non-Newtonian blood tests ─────────────────────────────────────────

#[test]
fn cascade_blood_nonnewtonian_converges() {
    let blood = CassonBlood::normal_blood();
    let solver = CascadeSolver3D::new(ci_config(), blood);
    let result = solver.solve(&cif_two_channel_specs()).unwrap();

    assert_eq!(result.channel_results.len(), 2);
    for cr in &result.channel_results {
        assert!(
            cr.converged,
            "Channel {} did not converge after {} Picard iterations",
            cr.channel_id,
            cr.picard_iterations
        );
    }
}

// ── Edge cases ────────────────────────────────────────────────────────

#[test]
fn cascade_empty_channels_returns_error() {
    let water = ConstantPropertyFluid::<f64>::water_20c().unwrap();
    let solver = CascadeSolver3D::new(ci_config(), water);
    assert!(solver.solve(&[]).is_err());
}

#[test]
fn cascade_single_straight_channel() {
    let water = ConstantPropertyFluid::<f64>::water_20c().unwrap();
    let solver = CascadeSolver3D::new(ci_config(), water);
    let specs = vec![CascadeChannelSpec {
        id: "straight".into(),
        length: 5.0e-3,
        width: 1.5e-3,
        height: 1.0e-3,
        flow_rate_m3_s: 3e-8,
        is_venturi_throat: false,
        throat_width: None,
        local_hematocrit: None,
    }];
    let result = solver.solve(&specs).unwrap();
    assert_eq!(result.channel_results.len(), 1);
    assert!(result.channel_results[0].max_velocity > 0.0);
}

// ── Hematocrit coupling (D2) ──────────────────────────────────────────

#[test]
fn cascade_hematocrit_sets_local_field() {
    // Verify that `local_hematocrit` is correctly propagated to results.
    let water = ConstantPropertyFluid::<f64>::water_20c().unwrap();
    let solver = CascadeSolver3D::new(ci_config(), water);

    let specs = vec![
        CascadeChannelSpec {
            id: "bypass".into(),
            length: 10.0e-3,
            width: 1.5e-3,
            height: 1.0e-3,
            flow_rate_m3_s: 3e-8,
            is_venturi_throat: false,
            throat_width: None,
            local_hematocrit: Some(0.20),
        },
        CascadeChannelSpec {
            id: "center".into(),
            length: 10.0e-3,
            width: 1.5e-3,
            height: 1.0e-3,
            flow_rate_m3_s: 3e-8,
            is_venturi_throat: false,
            throat_width: None,
            local_hematocrit: Some(0.55),
        },
        CascadeChannelSpec {
            id: "default".into(),
            length: 10.0e-3,
            width: 1.5e-3,
            height: 1.0e-3,
            flow_rate_m3_s: 3e-8,
            is_venturi_throat: false,
            throat_width: None,
            local_hematocrit: None,
        },
    ];
    let result = solver.solve(&specs).unwrap();

    let bypass = result.channel_results.iter().find(|r| r.channel_id == "bypass").unwrap();
    let center = result.channel_results.iter().find(|r| r.channel_id == "center").unwrap();
    let default_ch = result.channel_results.iter().find(|r| r.channel_id == "default").unwrap();

    assert_eq!(bypass.local_hematocrit, 0.20);
    assert_eq!(center.local_hematocrit, 0.55);
    assert_eq!(default_ch.local_hematocrit, 0.45); // feed default
}
