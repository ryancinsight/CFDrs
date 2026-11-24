use super::*;
use crate::physics::resistance::models::DarcyWeisbachModel;
use aequitas::systems::si::quantities::{
    DynamicViscosity, MassDensity, SpecificHeatCapacity, ThermalConductivity,
};
use cfd_core::physics::fluid::ConstantPropertyFluid;

#[test]
fn mach_limit_uses_velocity_magnitude_and_rejects_nonfinite_samples() {
    let subscriber = tracing_subscriber::fmt()
        .with_test_writer()
        .with_span_events(tracing_subscriber::fmt::format::FmtSpan::FULL)
        .finish();
    tracing::subscriber::with_default(subscriber, || {
        let (fluid, model, just_above) = tracing::info_span!("arrange").in_scope(|| {
            (
                ConstantPropertyFluid::new(
                    "constant liquid".into(),
                    MassDensity::from_base(1000.0),
                    DynamicViscosity::from_base(0.001),
                    SpecificHeatCapacity::from_base(4000.0),
                    ThermalConductivity::from_base(0.5),
                    Velocity::from_base(1000.0),
                ),
                DarcyWeisbachModel::circular(0.01, 0.1, 0.0),
                300.0_f64.next_up(),
            )
        });
        let passes = tracing::info_span!("act").in_scope(|| {
            [300.0, -300.0, just_above, -just_above, f64::NAN].map(|speed| {
                let conditions = FlowConditions::new(Velocity::from_base(speed));
                model.validate_mach_number(&fluid, &conditions).is_ok()
            })
        });
        tracing::info_span!("assert").in_scope(|| {
            // Ma = |v|/c, so both signs at 300 m/s are on the limit.
            assert_eq!(passes, [true, true, false, false, false]);
        });
    });
}

#[test]
fn model_returns_dimensioned_resistance_and_signed_pressure_drop() -> Result<()> {
    use aequitas::systems::si::quantities::{
        HydraulicResistance, Pressure, QuadraticHydraulicResistance, VolumetricFlowRate,
    };
    use cfd_core::physics::fluid::database::water_20c;

    let _trace = tracing::subscriber::set_default(
        tracing_subscriber::fmt()
            .with_test_writer()
            .with_span_events(tracing_subscriber::fmt::format::FmtSpan::FULL)
            .finish(),
    );
    let span = tracing::info_span!("arrange");
    let (model, fluid, flow) = span.in_scope(|| -> Result<_> {
        Ok((
            DarcyWeisbachModel::circular(0.01, 0.1, 0.0),
            water_20c::<f64>()?,
            VolumetricFlowRate::from_base(-1.0e-4),
        ))
    })?;
    let (resistance, pressure) = tracing::info_span!("act").in_scope(|| -> Result<_> {
        let conditions = FlowConditions::from_flow_rate(flow);
        let (r, k): (HydraulicResistance<f64>, QuadraticHydraulicResistance<f64>) =
            model.calculate_coefficients(&fluid, &conditions)?;
        let effective: HydraulicResistance<f64> =
            model.calculate_resistance(&fluid, &conditions)?;
        let pressure: Pressure<f64> = r * flow + k * flow * flow.abs();
        Ok((effective, pressure))
    })?;
    tracing::info_span!("assert").in_scope(|| {
        assert!(resistance.into_base() > 0.0);
        assert!(pressure.into_base() < 0.0);
        eunomia::assert_relative_eq!(
            pressure.into_base(),
            (resistance * flow).into_base(),
            max_relative = 2.0 * f64::EPSILON
        );
    });
    Ok(())
}
