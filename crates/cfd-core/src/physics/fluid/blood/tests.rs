use super::*;
use crate::physics::fluid::traits::{Fluid as FluidTrait, NonNewtonianFluid};

#[test]
fn test_fluid_trait_implementations() {
    let casson = CassonBlood::<f64>::normal_blood();
    let carreau = CarreauYasudaBlood::<f64>::normal_blood();
    let cross_blood = CrossBlood::<f64>::normal_blood();

    // All should implement FluidTrait
    let casson_state = casson
        .properties_at(310.0, 101325.0)
        .expect("expected value");
    let carreau_state = carreau
        .properties_at(310.0, 101325.0)
        .expect("expected value");
    let cross_state = cross_blood
        .properties_at(310.0, 101325.0)
        .expect("expected value");

    assert_eq!(casson_state.density.into_base(), constants::BLOOD_DENSITY);
    assert_eq!(carreau_state.density.into_base(), constants::BLOOD_DENSITY);
    assert_eq!(cross_state.density.into_base(), constants::BLOOD_DENSITY);

    // All have positive viscosity
    assert!(casson_state.dynamic_viscosity.into_base() > 0.0);
    assert!(carreau_state.dynamic_viscosity.into_base() > 0.0);
    assert!(cross_state.dynamic_viscosity.into_base() > 0.0);
}

#[test]
fn test_non_newtonian_trait() {
    let casson = CassonBlood::<f64>::normal_blood();
    let carreau = CarreauYasudaBlood::<f64>::normal_blood();

    // Casson has yield stress
    assert!(casson.has_yield_stress());
    assert!(casson.yield_stress().is_some());

    // Carreau-Yasuda does not have yield stress
    assert!(!carreau.has_yield_stress());
    assert!(carreau.yield_stress().is_none());
}

#[test]
fn test_model_comparison_at_intermediate_shear() {
    let casson = CassonBlood::<f64>::normal_blood();
    let carreau = CarreauYasudaBlood::<f64>::normal_blood();
    let cross_blood = CrossBlood::<f64>::normal_blood();

    // At γ̇ = 100 s⁻¹ (typical arterial condition), all models should agree within factor of 2
    let gamma = 100.0;
    let mu_casson = casson.apparent_viscosity(gamma);
    let mu_carreau = carreau.apparent_viscosity(gamma);
    let mu_cross = cross_blood.apparent_viscosity(gamma);

    // All should be in reasonable range (3-10 mPa·s)
    for (name, mu) in [
        ("Casson", mu_casson),
        ("Carreau", mu_carreau),
        ("Cross", mu_cross),
    ] {
        assert!(
            mu > 0.003 && mu < 0.010,
            "{name} viscosity at 100 s⁻¹ = {mu} Pa·s out of expected range"
        );
    }
}
