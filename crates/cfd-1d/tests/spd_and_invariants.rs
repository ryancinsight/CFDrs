use cfd_1d::domain::network::{
    Edge, EdgeProperties, EdgeType, Network, NetworkBuilder, ResistanceUpdatePolicy,
};
use cfd_1d::solver::core::{LinearSolverMethod, NetworkProblem, NetworkSolver};
use cfd_1d::{
    durst_resistance_multiplier, ChannelGeometry, CrossSection, FlowConditions,
    ResistanceCalculator, ResistanceChannelGeometry, SurfaceProperties, Wettability,
};
use cfd_core::compute::solver::Solver;
use cfd_core::error::Result;
use cfd_core::physics::fluid::blood::{CarreauYasudaBlood, FahraeuasLindqvist};
use cfd_core::physics::fluid::{ConstantPropertyFluid, FluidTrait};

fn network_two_node<T: nalgebra::RealField + Copy + num_traits::FromPrimitive>() -> (
    Network<T>,
    petgraph::graph::EdgeIndex,
    petgraph::graph::NodeIndex,
    petgraph::graph::NodeIndex,
) {
    let mut builder = NetworkBuilder::<T>::new();
    let inlet = builder.add_inlet("in".to_string());
    let outlet = builder.add_outlet("out".to_string());
    let edge = builder.add_edge(inlet, outlet, Edge::new("e".to_string(), EdgeType::Pipe));
    let graph = builder.build().expect("test invariant");
    let fluid = ConstantPropertyFluid::<T>::water_20c().expect("test invariant");
    let net = Network::new(graph, fluid);
    (net, edge, inlet, outlet)
}

#[test]
fn spd_heuristic_selects_cg() -> Result<()> {
    type F = f64;
    let (mut net, edge, inlet, _outlet) = network_two_node::<F>();
    if let Some(e) = net.graph.edge_weight_mut(edge) {
        e.resistance = 4.0;
        e.quad_coeff = 0.0;
    }
    net.set_pressure(inlet, 10.0);
    let problem = NetworkProblem::new(net.clone());
    let solver = NetworkSolver::<F>::new();
    let solved = solver.solve(&problem)?;
    match solved.last_solver_method.expect("test invariant") {
        LinearSolverMethod::ConjugateGradient => (),
        _ => panic!("expected CG method selection"),
    }
    Ok(())
}

#[test]
fn negative_coefficients_rejected_in_update() {
    type F = f64;
    let (mut net, edge, inlet, outlet) = network_two_node::<F>();
    if let Some(e) = net.graph.edge_weight_mut(edge) {
        e.resistance = -1.0;
        e.quad_coeff = 0.0;
    }
    let mut x = nalgebra::DVector::<F>::zeros(net.node_count());
    x[inlet.index()] = 5.0;
    x[outlet.index()] = 0.0;
    let err = net.update_from_solution(&x).unwrap_err();
    match err {
        cfd_core::error::Error::InvalidConfiguration(_) => (),
        _ => panic!("unexpected error type"),
    }
}

#[test]
fn flow_invariant_edges_skip_resistance_recomputation() -> Result<()> {
    type F = f64;
    let (mut net, edge, _inlet, _outlet) = network_two_node::<F>();
    if let Some(e) = net.graph.edge_weight_mut(edge) {
        e.resistance = 7.5;
        e.quad_coeff = 0.0;
    }
    net.set_flow_rate(edge, 3.0e-9);
    net.add_edge_properties(
        edge,
        EdgeProperties {
            id: "e".to_string(),
            component_type: cfd_1d::ComponentType::Pipe,
            length: 1.0e-2,
            area: 1.0e-6,
            hydraulic_diameter: Some(1.0e-3),
            resistance: 7.5,
            geometry: Some(ChannelGeometry {
                channel_type: cfd_1d::ChannelType::Straight,
                length: 1.0e-2,
                cross_section: CrossSection::Circular { diameter: 1.0e-3 },
                surface: SurfaceProperties {
                    roughness: 0.0,
                    contact_angle: None,
                    surface_energy: None,
                    wettability: Wettability::Hydrophilic,
                },
                variations: Vec::new(),
            }),
            resistance_update_policy: ResistanceUpdatePolicy::FlowInvariant,
            properties: std::collections::HashMap::new(),
        },
    );

    net.update_resistances()?;
    let edge_ref = net.graph.edge_weight(edge).expect("edge must exist");
    assert!((edge_ref.resistance - 7.5).abs() < 1.0e-12);
    assert!(edge_ref.quad_coeff.abs() < 1.0e-12);
    Ok(())
}

#[test]
fn flow_dependent_short_channel_reapplies_durst_correction() -> Result<()> {
    type F = f64;
    let (mut net, edge, _inlet, _outlet) = network_two_node::<F>();
    let diameter = 100.0e-6;
    let length = 1.0e-3;
    let area = std::f64::consts::PI * diameter * diameter / 4.0;
    let flow_rate = 2.0e-9;

    net.set_flow_rate(edge, flow_rate);
    net.add_edge_properties(
        edge,
        EdgeProperties {
            id: "e".to_string(),
            component_type: cfd_1d::ComponentType::Pipe,
            length,
            area,
            hydraulic_diameter: Some(diameter),
            resistance: 0.0,
            geometry: Some(ChannelGeometry {
                channel_type: cfd_1d::ChannelType::Straight,
                length,
                cross_section: CrossSection::Circular { diameter },
                surface: SurfaceProperties {
                    roughness: 0.0,
                    contact_angle: None,
                    surface_energy: None,
                    wettability: Wettability::Hydrophilic,
                },
                variations: Vec::new(),
            }),
            resistance_update_policy: ResistanceUpdatePolicy::FlowDependent,
            properties: std::collections::HashMap::new(),
        },
    );

    net.update_resistances()?;

    let mut conditions = FlowConditions::new(0.0);
    conditions.velocity = None;
    conditions.flow_rate = Some(flow_rate);
    conditions.temperature =
        cfd_core::physics::constants::physics::thermo::T_STANDARD;

    let calculator = ResistanceCalculator::<F>::new();
    let (base_resistance, _) = calculator.calculate_coefficients_auto(
        &ResistanceChannelGeometry::Circular { diameter, length },
        net.fluid(),
        &conditions,
    )?;

    let fluid_state = net
        .fluid()
        .properties_at(conditions.temperature, conditions.pressure)?;
    let velocity = flow_rate / area;
    let reynolds = fluid_state.density * velocity * diameter / fluid_state.dynamic_viscosity;
    let expected_multiplier = durst_resistance_multiplier(reynolds, length / diameter);

    let edge_ref = net.graph.edge_weight(edge).expect("edge must exist");
    assert!((edge_ref.resistance - base_resistance * expected_multiplier).abs() < 1.0e-9);
    assert!(edge_ref.resistance > base_resistance);
    assert!(edge_ref.quad_coeff.abs() < 1.0e-12);
    Ok(())
}

#[test]
fn flow_dependent_blood_microchannel_reapplies_fahraeus_lindqvist_reduction() -> Result<()> {
    type F = f64;
    let mut builder = NetworkBuilder::<F>::new();
    let inlet = builder.add_inlet("in".to_string());
    let outlet = builder.add_outlet("out".to_string());
    let edge = builder.add_edge(inlet, outlet, Edge::new("e".to_string(), EdgeType::Pipe));
    let graph = builder.build().expect("test invariant");
    let fluid = CarreauYasudaBlood::<F>::normal_blood();
    let mut net = Network::new(graph, fluid);

    let diameter = 100.0e-6;
    let length = 10.0e-3;
    let area = std::f64::consts::PI * diameter * diameter / 4.0;
    let flow_rate = 2.0e-9;

    net.set_flow_rate(edge, flow_rate);
    net.add_edge_properties(
        edge,
        EdgeProperties {
            id: "e".to_string(),
            component_type: cfd_1d::ComponentType::Pipe,
            length,
            area,
            hydraulic_diameter: Some(diameter),
            resistance: 0.0,
            geometry: Some(ChannelGeometry {
                channel_type: cfd_1d::ChannelType::Straight,
                length,
                cross_section: CrossSection::Circular { diameter },
                surface: SurfaceProperties {
                    roughness: 0.0,
                    contact_angle: None,
                    surface_energy: None,
                    wettability: Wettability::Hydrophilic,
                },
                variations: Vec::new(),
            }),
            resistance_update_policy: ResistanceUpdatePolicy::FlowDependent,
            properties: std::collections::HashMap::new(),
        },
    );

    net.update_resistances()?;

    let mut conditions = FlowConditions::new(0.0);
    conditions.velocity = None;
    conditions.flow_rate = Some(flow_rate);
    conditions.temperature =
        cfd_core::physics::constants::physics::thermo::T_STANDARD;

    let calculator = ResistanceCalculator::<F>::new();
    let (base_resistance, _) = calculator.calculate_coefficients_auto(
        &ResistanceChannelGeometry::Circular { diameter, length },
        net.fluid(),
        &conditions,
    )?;

    let fl = FahraeuasLindqvist::<F>::new(diameter, 0.45);
    let expected_reduction = (fl.relative_viscosity() / 3.2).clamp(0.5, 1.0);

    let edge_ref = net.graph.edge_weight(edge).expect("edge must exist");
    assert!((edge_ref.resistance - base_resistance * expected_reduction).abs() < 1.0e-9);
    assert!(edge_ref.resistance < base_resistance);
    assert!(edge_ref.quad_coeff.abs() < 1.0e-12);
    Ok(())
}

#[test]
fn flow_dependent_recompute_uses_shear_aware_reynolds_for_blood() -> Result<()> {
    type F = f64;
    let mut builder = NetworkBuilder::<F>::new();
    let inlet = builder.add_inlet("in".to_string());
    let outlet = builder.add_outlet("out".to_string());
    let edge = builder.add_edge(inlet, outlet, Edge::new("e".to_string(), EdgeType::Pipe));
    let graph = builder.build().expect("test invariant");
    let fluid = CarreauYasudaBlood::<F>::normal_blood();
    let mut net = Network::new(graph, fluid);

    let diameter = 0.04;
    let length = 2.5;
    let area = std::f64::consts::PI * diameter * diameter / 4.0;
    let fluid_state = net
        .fluid()
        .properties_at(
            cfd_core::physics::constants::physics::thermo::T_STANDARD,
            cfd_core::physics::constants::physics::thermo::P_ATM,
        )?;
    let target_reference_re = 2400.0;
    let flow_rate = target_reference_re * std::f64::consts::PI * diameter
        * fluid_state.dynamic_viscosity
        / (4.0 * fluid_state.density);
    let velocity = flow_rate / area;
    let shear_rate = 8.0 * velocity / diameter;
    let apparent_viscosity = net.fluid().viscosity_at_shear(
        shear_rate,
        cfd_core::physics::constants::physics::thermo::T_STANDARD,
        cfd_core::physics::constants::physics::thermo::P_ATM,
    )?;
    let shear_aware_re = fluid_state.density * velocity * diameter / apparent_viscosity;

    assert!(target_reference_re > 2300.0);
    assert!(shear_aware_re < 2300.0, "expected shear-aware Reynolds to remain laminar, got {shear_aware_re}");

    net.set_flow_rate(edge, flow_rate);
    net.add_edge_properties(
        edge,
        EdgeProperties {
            id: "e".to_string(),
            component_type: cfd_1d::ComponentType::Pipe,
            length,
            area,
            hydraulic_diameter: Some(diameter),
            resistance: 0.0,
            geometry: Some(ChannelGeometry {
                channel_type: cfd_1d::ChannelType::Straight,
                length,
                cross_section: CrossSection::Circular { diameter },
                surface: SurfaceProperties {
                    roughness: 0.0,
                    contact_angle: None,
                    surface_energy: None,
                    wettability: Wettability::Hydrophilic,
                },
                variations: Vec::new(),
            }),
            resistance_update_policy: ResistanceUpdatePolicy::FlowDependent,
            properties: std::collections::HashMap::new(),
        },
    );

    net.update_resistances()?;

    let mut conditions = FlowConditions::new(0.0);
    conditions.velocity = None;
    conditions.flow_rate = Some(flow_rate);
    conditions.temperature =
        cfd_core::physics::constants::physics::thermo::T_STANDARD;
    conditions.shear_rate = Some(shear_rate);
    conditions.reynolds_number = Some(shear_aware_re);

    let calculator = ResistanceCalculator::<F>::new();
    let (expected_resistance, expected_k) = calculator.calculate_hagen_poiseuille_coefficients(
        diameter,
        length,
        net.fluid(),
        &conditions,
    )?;

    let edge_ref = net.graph.edge_weight(edge).expect("edge must exist");
    assert!((edge_ref.resistance - expected_resistance).abs() < 1.0e-9);
    assert!((edge_ref.quad_coeff - expected_k).abs() < 1.0e-12);
    assert!(edge_ref.quad_coeff.abs() < 1.0e-12);
    Ok(())
}
