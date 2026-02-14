use cfd_1d::components::{Component, ComponentFactory, OrganCompartment, PorousMembrane};
use cfd_1d::resistance::{FlowConditions, ResistanceCalculator};
use std::collections::HashMap;

#[test]
fn porous_membrane_component_returns_positive_resistance() {
    let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().expect("water");
    let membrane = PorousMembrane::new(10e-6, 1e-3, 1e-3, 0.5e-6, 0.2);
    let r = membrane.resistance(&fluid);
    assert!(r > 0.0);
}

#[test]
fn resistance_calculator_supports_membrane_model() {
    let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().expect("water");
    let calc = ResistanceCalculator::<f64>::new();
    let conditions = FlowConditions::new(0.0);

    let r = calc
        .calculate_membrane_porous(10e-6, 1e-3, 1e-3, 0.5e-6, 0.2, &fluid, &conditions)
        .expect("membrane resistance");

    assert!(r > 0.0);
}

#[test]
fn component_factory_creates_membrane_and_organ_components() {
    let mut membrane_params = HashMap::new();
    membrane_params.insert("thickness".to_string(), 10e-6);
    membrane_params.insert("width".to_string(), 1e-3);
    membrane_params.insert("height".to_string(), 1e-3);
    membrane_params.insert("pore_radius".to_string(), 0.5e-6);
    membrane_params.insert("porosity".to_string(), 0.2);

    let membrane = ComponentFactory::create::<f64>("PorousMembrane", &membrane_params)
        .expect("factory membrane");
    assert_eq!(membrane.component_type(), "PorousMembrane");

    let mut organ_params = HashMap::new();
    organ_params.insert("length".to_string(), 1e-3);
    organ_params.insert("width".to_string(), 1e-3);
    organ_params.insert("height".to_string(), 0.5e-3);
    organ_params.insert("hydraulic_resistance".to_string(), 5e11);

    let organ =
        ComponentFactory::create::<f64>("OrganCompartment", &organ_params).expect("factory organ");
    assert_eq!(organ.component_type(), "OrganCompartment");
}

#[test]
fn organ_compartment_exposes_volume() {
    let organ = OrganCompartment::new(2e-3, 1e-3, 0.5e-3, 1e12);
    let volume = organ.volume().expect("volume");
    assert!(volume > 0.0);
}
