use super::*;
use crate::physics::resistance::models::FlowConditions;
use aequitas::systems::si::quantities::{ReciprocalTime, Velocity, VolumetricFlowRate};
use cfd_core::physics::fluid::CassonBlood;
use eunomia::assert_relative_eq;

fn make_model(jtype: JunctionType, direction: JunctionFlowDirection) -> JunctionLossModel {
    JunctionLossModel {
        junction_type: jtype,
        flow_direction: direction,
        branch_diameter: Length::from_base(500e-6),
        branch_area: Area::from_base(std::f64::consts::PI * (250e-6_f64).powi(2)),
        junction_length: Length::from_base(500e-6),
        density: MassDensity::from_base(1000.0),
        run_area: None,
    }
}

#[test]
fn tee_dividing_k_factor() {
    let m = make_model(JunctionType::Tee, JunctionFlowDirection::Dividing);
    assert!((m.loss_coefficient().into_base() - 1.0).abs() < 1e-10);
}

#[test]
fn cross_combining_k_factor() {
    let m = make_model(JunctionType::Cross, JunctionFlowDirection::Combining);
    assert!((m.loss_coefficient().into_base() - 1.2).abs() < 1e-10);
}

#[test]
fn linearised_resistance_positive() {
    let m = make_model(JunctionType::Cross, JunctionFlowDirection::Dividing);
    let r = m
        .linearised_resistance(VolumetricFlowRate::from_base(1e-6))
        .into_base();
    assert!(r > 0.0, "linearised resistance must be positive");
}

#[test]
fn resistance_model_trait_works() {
    let m = make_model(JunctionType::Tee, JunctionFlowDirection::Combining);
    let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().expect("expected value");
    let cond = FlowConditions::new(Velocity::from_base(0.01));
    let r = m
        .calculate_resistance(&fluid, &cond)
        .expect("expected value")
        .into_base();
    assert!(r > 0.0, "resistance must be positive: got {r}");
}

#[test]
fn coefficients_have_quadratic_term() {
    let m = make_model(JunctionType::Cross, JunctionFlowDirection::Dividing);
    let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().expect("expected value");
    let cond = FlowConditions::new(Velocity::from_base(0.01));
    let (r, k) = m
        .calculate_coefficients(&fluid, &cond)
        .expect("expected value");
    let (r, k) = (r.into_base(), k.into_base());
    assert!(r > 0.0, "linear resistance must be positive");
    assert!(k > 0.0, "quadratic loss coefficient must be positive");
}

#[test]
fn non_newtonian_resistance_uses_explicit_shear_rate() {
    let m = make_model(JunctionType::Tee, JunctionFlowDirection::Combining);
    let blood = CassonBlood::<f64>::normal_blood();
    let mut low_shear = FlowConditions::new(Velocity::from_base(0.01));
    low_shear.shear_rate = Some(ReciprocalTime::from_base(10.0));
    let mut high_shear = FlowConditions::new(Velocity::from_base(0.01));
    high_shear.shear_rate = Some(ReciprocalTime::from_base(1000.0));

    let (low_r, low_k) = m
        .calculate_coefficients(&blood, &low_shear)
        .expect("positive explicit shear rate must compute junction loss");
    let (low_r, low_k) = (low_r.into_base(), low_k.into_base());
    let (high_r, high_k) = m
        .calculate_coefficients(&blood, &high_shear)
        .expect("positive explicit shear rate must compute junction loss");
    let (high_r, high_k) = (high_r.into_base(), high_k.into_base());

    assert!(
        low_r > high_r,
        "Casson shear-thinning requires higher junction viscous resistance at 10 1/s than at 1000 1/s: {low_r} <= {high_r}"
    );
    assert_relative_eq!(low_k, high_k, max_relative = 1e-12);
    assert!(low_k > 0.0);
}

#[test]
fn negative_explicit_shear_rate_rejected() {
    let m = make_model(JunctionType::Cross, JunctionFlowDirection::Dividing);
    let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().expect("expected value");
    let mut cond = FlowConditions::new(Velocity::from_base(0.01));
    cond.shear_rate = Some(ReciprocalTime::from_base(-1.0));

    let err = m
        .calculate_coefficients(&fluid, &cond)
        .expect_err("negative explicit wall shear rate violates junction rheology");
    assert!(
        err.to_string().contains("shear rate"),
        "error must identify the rejected physical invariant: {err}"
    );
}

#[test]
fn non_newtonian_reverse_flow_uses_velocity_magnitude() {
    let m = make_model(JunctionType::Tee, JunctionFlowDirection::Combining);
    let blood = CassonBlood::<f64>::normal_blood();
    let forward = FlowConditions::new(Velocity::from_base(0.01));
    let reverse = FlowConditions::new(Velocity::from_base(-0.01));

    let (forward_r, forward_k) = m
        .calculate_coefficients(&blood, &forward)
        .expect("forward non-Newtonian junction loss must compute");
    let (forward_r, forward_k) = (forward_r.into_base(), forward_k.into_base());
    let (reverse_r, reverse_k) = m
        .calculate_coefficients(&blood, &reverse)
        .expect("reverse non-Newtonian junction loss must compute");
    let (reverse_r, reverse_k) = (reverse_r.into_base(), reverse_k.into_base());

    assert_relative_eq!(forward_r, reverse_r, max_relative = 1e-12);
    assert_relative_eq!(forward_k, reverse_k, max_relative = 1e-12);
    assert!(forward_k > 0.0);
}

#[test]
fn combining_higher_loss_than_dividing_for_tee() {
    let div = make_model(JunctionType::Tee, JunctionFlowDirection::Dividing);
    let comb = make_model(JunctionType::Tee, JunctionFlowDirection::Combining);
    assert!(
        comb.loss_coefficient().into_base() > div.loss_coefficient().into_base(),
        "combining flow typically has higher losses than dividing"
    );
}

#[test]
fn diameter_ratio_correction_narrows_branch_reduces_k() {
    let equal = make_model(JunctionType::Cross, JunctionFlowDirection::Dividing);
    let k_equal = equal.loss_coefficient().into_base();
    let narrow_branch = JunctionLossModel {
        run_area: Some(equal.branch_area * 4.0), // run 4× larger than branch
        ..equal
    };
    let k_narrow = narrow_branch.loss_coefficient().into_base();
    assert!(
        k_narrow < k_equal,
        "narrow branch into wide run should reduce K: {k_narrow} vs {k_equal}"
    );
    // area ratio = 0.25, sqrt = 0.5 → K_narrow = K_equal × 0.5
    assert!((k_narrow - k_equal * 0.5).abs() < 1e-10);
}

#[test]
fn diameter_ratio_none_gives_equal_area_k() {
    let m = make_model(JunctionType::Tee, JunctionFlowDirection::Dividing);
    assert!(m.run_area.is_none());
    assert!((m.loss_coefficient().into_base() - 1.0).abs() < 1e-10);
}

#[test]
fn nonfinite_geometry_is_rejected() {
    let subscriber = tracing_subscriber::fmt()
        .with_test_writer()
        .with_span_events(tracing_subscriber::fmt::format::FmtSpan::FULL)
        .finish();
    tracing::subscriber::with_default(subscriber, || {
        let (mut junction, fluid, conditions) = tracing::info_span!("arrange").in_scope(|| {
            (
                make_model(JunctionType::Tee, JunctionFlowDirection::Dividing),
                cfd_core::physics::fluid::database::water_20c::<f64>().expect("water"),
                FlowConditions::new(Velocity::from_base(0.01)),
            )
        });
        let result = tracing::info_span!("act").in_scope(|| {
            junction.branch_area = Area::from_base(f64::NAN);
            junction.validate_invariants(&fluid, &conditions)
        });
        tracing::info_span!("assert").in_scope(|| {
            assert!(
                result
                    .expect_err("NaN area is invalid")
                    .to_string()
                    .contains("branch area")
            );
        });
    });
}

#[test]
fn invalid_geometry_rejected() {
    let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().expect("expected value");
    let cond = FlowConditions::new(Velocity::from_base(0.01));

    let mut m = make_model(JunctionType::Tee, JunctionFlowDirection::Combining);
    m.branch_area = Area::from_base(-1e-6);
    assert!(m.validate_invariants(&fluid, &cond).is_err());

    let mut m2 = make_model(JunctionType::Tee, JunctionFlowDirection::Combining);
    m2.run_area = Some(Area::from_base(-1e-6));
    assert!(m2.validate_invariants(&fluid, &cond).is_err());
}

#[test]
fn junction_reference_flow_has_the_signed_loss_derivative_and_existing_wire_shape() {
    let subscriber = tracing_subscriber::fmt()
        .with_test_writer()
        .with_span_events(tracing_subscriber::fmt::format::FmtSpan::FULL)
        .finish();
    tracing::subscriber::with_default(subscriber, || {
        let input = tracing::info_span!("arrange").in_scope(|| {
            serde_json::json!({
                "junction_type": "Tee",
                "flow_direction": "Dividing",
                "branch_diameter_m": 0.5,
                "branch_area_m2": 0.25,
                "junction_length_m": 0.125,
                "density_kg_m3": 1024.0,
                "run_area_m2": 1.0,
            })
        });
        let (model, slope, reverse, quadratic) = tracing::info_span!("act").in_scope(|| {
            let model: JunctionLossModel = serde_json::from_value(input.clone()).expect("junction");
            let flow = VolumetricFlowRate::from_base(0.125);
            (
                serde_json::to_value(&model).expect("junction JSON"),
                model.linearised_resistance(flow),
                model.linearised_resistance(-flow),
                minor_loss_coefficient(model.loss_coefficient(), model.density, model.branch_area),
            )
        });
        tracing::info_span!("assert").in_scope(|| {
            assert_eq!(input, model);
            // K = sqrt(1/4) = 1/2. C = K * 1024 / (2 * (1/4)^2) = 4096.
            // The derivative of C Q |Q| is 2 C |Q|, including for reverse flow.
            assert_eq!(quadratic.into_base(), 4096.0);
            assert_eq!(slope, HydraulicResistance::from_base(1024.0));
            assert_eq!(reverse, slope);
        });
    });
}
