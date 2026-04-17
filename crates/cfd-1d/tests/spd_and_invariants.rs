use cfd_1d::domain::network::{
    Edge, EdgeProperties, EdgeType, Network, NetworkBuilder, ResistanceUpdatePolicy,
    EDGE_PROPERTY_HEMATOCRIT, EDGE_PROPERTY_LOCAL_APPARENT_VISCOSITY_PA_S,
    EDGE_PROPERTY_LOCAL_HEMATOCRIT, EDGE_PROPERTY_PLASMA_VISCOSITY_PA_S,
};
use cfd_1d::solver::core::{LinearSolverMethod, NetworkProblem, NetworkSolver};
use cfd_1d::{
    blood_microchannel_apparent_viscosity, durst_resistance_multiplier, pries_phase_separation,
    ChannelGeometry, CrossSection, FlowConditions, ResistanceCalculator, ResistanceChannelGeometry,
    SurfaceProperties, Wettability,
};
use cfd_core::compute::solver::Solver;
use cfd_core::error::Result;
use cfd_core::physics::fluid::blood::CarreauYasudaBlood;
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
    conditions.temperature = cfd_core::physics::constants::physics::thermo::T_STANDARD;

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
fn flow_dependent_blood_microchannel_uses_hematocrit_aware_apparent_viscosity() -> Result<()> {
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
    conditions.temperature = cfd_core::physics::constants::physics::thermo::T_STANDARD;

    let calculator = ResistanceCalculator::<F>::new();
    let (base_resistance, _) = calculator.calculate_coefficients_auto(
        &ResistanceChannelGeometry::Circular { diameter, length },
        net.fluid(),
        &conditions,
    )?;

    let fluid_state = net
        .fluid()
        .properties_at(conditions.temperature, conditions.pressure)?;
    let target_mu = blood_microchannel_apparent_viscosity(
        diameter,
        flow_rate,
        area,
        0.45,
        fluid_state.dynamic_viscosity / 3.2,
    )
    .expect("blood apparent viscosity estimate must exist");
    let shear_rate = 8.0 * (flow_rate / area) / diameter;
    let base_mu =
        net.fluid()
            .viscosity_at_shear(shear_rate, conditions.temperature, conditions.pressure)?;
    let expected_factor = target_mu / base_mu;

    let edge_ref = net.graph.edge_weight(edge).expect("edge must exist");
    assert!((edge_ref.resistance - base_resistance * expected_factor).abs() < 1.0e-9);
    assert!(edge_ref.quad_coeff.abs() < 1.0e-12);
    Ok(())
}

#[test]
fn blood_microchannel_override_properties_change_resistance() -> Result<()> {
    type F = f64;
    let mut builder = NetworkBuilder::<F>::new();
    let inlet = builder.add_inlet("in".to_string());
    let outlet = builder.add_outlet("out".to_string());
    let edge = builder.add_edge(inlet, outlet, Edge::new("e".to_string(), EdgeType::Pipe));
    let graph = builder.build().expect("test invariant");
    let fluid = CarreauYasudaBlood::<F>::normal_blood();
    let mut net = Network::new(graph, fluid);

    let diameter = 80.0e-6;
    let length = 8.0e-3;
    let area = std::f64::consts::PI * diameter * diameter / 4.0;
    let flow_rate = 1.5e-9;

    net.set_flow_rate(edge, flow_rate);
    let mut properties = std::collections::HashMap::new();
    properties.insert(EDGE_PROPERTY_HEMATOCRIT.to_string(), 0.25);
    properties.insert(EDGE_PROPERTY_PLASMA_VISCOSITY_PA_S.to_string(), 1.05e-3);
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
            properties,
        },
    );

    net.update_resistances()?;
    let overridden_resistance = net
        .graph
        .edge_weight(edge)
        .expect("edge must exist")
        .resistance;

    let mut net_default = net.clone();
    let mut default_props = net_default
        .properties
        .get_mut(&edge)
        .expect("edge properties must exist")
        .properties
        .clone();
    default_props.remove(EDGE_PROPERTY_HEMATOCRIT);
    default_props.remove(EDGE_PROPERTY_PLASMA_VISCOSITY_PA_S);
    net_default
        .properties
        .get_mut(&edge)
        .expect("edge properties must exist")
        .properties = default_props;
    net_default.update_resistances()?;
    let default_resistance = net_default
        .graph
        .edge_weight(edge)
        .expect("edge must exist")
        .resistance;

    assert!(
        (overridden_resistance - default_resistance).abs() > 1.0e-12,
        "hematocrit/plasma-viscosity override should change the edge resistance"
    );
    Ok(())
}

#[test]
fn local_apparent_viscosity_override_changes_blood_edge_resistance() -> Result<()> {
    type F = f64;
    let mut builder = NetworkBuilder::<F>::new();
    let inlet = builder.add_inlet("in".to_string());
    let outlet = builder.add_outlet("out".to_string());
    let edge = builder.add_edge(inlet, outlet, Edge::new("e".to_string(), EdgeType::Pipe));
    let graph = builder.build().expect("test invariant");
    let fluid = CarreauYasudaBlood::<F>::normal_blood();
    let mut net = Network::new(graph, fluid);

    let diameter = 80.0e-6;
    let length = 8.0e-3;
    let area = std::f64::consts::PI * diameter * diameter / 4.0;
    let flow_rate = 1.5e-9;

    net.set_flow_rate(edge, flow_rate);
    let mut properties = std::collections::HashMap::new();
    properties.insert(EDGE_PROPERTY_HEMATOCRIT.to_string(), 0.25);
    properties.insert(EDGE_PROPERTY_PLASMA_VISCOSITY_PA_S.to_string(), 1.05e-3);
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
            properties,
        },
    );

    net.update_resistances()?;
    let baseline = net
        .graph
        .edge_weight(edge)
        .expect("edge must exist")
        .resistance;

    net.properties
        .get_mut(&edge)
        .expect("edge properties must exist")
        .properties
        .insert(
            EDGE_PROPERTY_LOCAL_APPARENT_VISCOSITY_PA_S.to_string(),
            1.0e-2,
        );
    net.update_resistances()?;
    let overridden = net
        .graph
        .edge_weight(edge)
        .expect("edge must exist")
        .resistance;

    assert!((overridden - baseline).abs() > 1.0e-12);
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
    let fluid_state = net.fluid().properties_at(
        cfd_core::physics::constants::physics::thermo::T_STANDARD,
        cfd_core::physics::constants::physics::thermo::P_ATM,
    )?;
    let target_reference_re = 2400.0;
    let flow_rate =
        target_reference_re * std::f64::consts::PI * diameter * fluid_state.dynamic_viscosity
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
    assert!(
        shear_aware_re < 2300.0,
        "expected shear-aware Reynolds to remain laminar, got {shear_aware_re}"
    );

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
    conditions.temperature = cfd_core::physics::constants::physics::thermo::T_STANDARD;
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

#[test]
fn bifurcation_transport_propagates_phase_separated_hematocrit() -> Result<()> {
    type F = f64;
    let mut builder = NetworkBuilder::<F>::new();
    let inlet = builder.add_inlet("in".to_string());
    let junction = builder.add_junction("j".to_string());
    let outlet_a = builder.add_outlet("out_a".to_string());
    let outlet_b = builder.add_outlet("out_b".to_string());
    let parent = builder.add_edge(
        inlet,
        junction,
        Edge::new("parent".to_string(), EdgeType::Pipe),
    );
    let daughter_a = builder.add_edge(
        junction,
        outlet_a,
        Edge::new("a".to_string(), EdgeType::Pipe),
    );
    let daughter_b = builder.add_edge(
        junction,
        outlet_b,
        Edge::new("b".to_string(), EdgeType::Pipe),
    );
    let graph = builder.build().expect("test invariant");
    let fluid = CarreauYasudaBlood::<F>::normal_blood();
    let mut net = Network::new(graph, fluid);

    let parent_d = 120.0e-6;
    let daughter_a_d = 90.0e-6;
    let daughter_b_d = 60.0e-6;
    let length = 6.0e-3;
    let feed_hct = 0.45;
    let q_parent = 3.0e-9;
    let q_a = 2.0e-9;
    let q_b = 1.0e-9;

    net.set_flow_rate(parent, q_parent);
    net.set_flow_rate(daughter_a, q_a);
    net.set_flow_rate(daughter_b, q_b);

    let make_props =
        |id: &str, diameter: f64, properties: std::collections::HashMap<String, f64>| {
            EdgeProperties {
                id: id.to_string(),
                component_type: cfd_1d::ComponentType::Pipe,
                length,
                area: std::f64::consts::PI * diameter * diameter / 4.0,
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
                properties,
            }
        };

    let mut parent_props = std::collections::HashMap::new();
    parent_props.insert(EDGE_PROPERTY_HEMATOCRIT.to_string(), feed_hct);
    net.add_edge_properties(parent, make_props("parent", parent_d, parent_props));
    net.add_edge_properties(
        daughter_a,
        make_props("a", daughter_a_d, std::collections::HashMap::new()),
    );
    net.add_edge_properties(
        daughter_b,
        make_props("b", daughter_b_d, std::collections::HashMap::new()),
    );

    net.update_resistances()?;

    let expected_a = pries_phase_separation(
        feed_hct,
        q_a / (q_a + q_b),
        daughter_a_d * 1.0e6,
        daughter_b_d * 1.0e6,
        parent_d * 1.0e6,
    )?;
    let expected_b =
        ((1.0 - expected_a.cell_fraction) * feed_hct / (q_b / (q_a + q_b))).clamp(0.0, 1.0);

    let hct_a = net
        .properties
        .get(&daughter_a)
        .and_then(|props| props.properties.get(EDGE_PROPERTY_LOCAL_HEMATOCRIT))
        .copied()
        .expect("daughter a should receive transported hematocrit");
    let hct_b = net
        .properties
        .get(&daughter_b)
        .and_then(|props| props.properties.get(EDGE_PROPERTY_LOCAL_HEMATOCRIT))
        .copied()
        .expect("daughter b should receive transported hematocrit");

    assert!((hct_a - expected_a.daughter_hematocrit).abs() < 1.0e-9);
    assert!((hct_b - expected_b).abs() < 1.0e-9);
    assert!(
        hct_a > hct_b,
        "larger daughter should carry higher hematocrit"
    );
    Ok(())
}

#[test]
fn cascade_transport_uses_parent_daughter_hematocrit_as_next_feed() -> Result<()> {
    type F = f64;
    let mut builder = NetworkBuilder::<F>::new();
    let inlet = builder.add_inlet("in".to_string());
    let j1 = builder.add_junction("j1".to_string());
    let j2 = builder.add_junction("j2".to_string());
    let outlet_b = builder.add_outlet("out_b".to_string());
    let outlet_c = builder.add_outlet("out_c".to_string());
    let outlet_d = builder.add_outlet("out_d".to_string());
    let parent = builder.add_edge(inlet, j1, Edge::new("parent".to_string(), EdgeType::Pipe));
    let branch_b = builder.add_edge(j1, outlet_b, Edge::new("b".to_string(), EdgeType::Pipe));
    let trunk = builder.add_edge(j1, j2, Edge::new("trunk".to_string(), EdgeType::Pipe));
    let branch_c = builder.add_edge(j2, outlet_c, Edge::new("c".to_string(), EdgeType::Pipe));
    let branch_d = builder.add_edge(j2, outlet_d, Edge::new("d".to_string(), EdgeType::Pipe));
    let graph = builder.build().expect("test invariant");
    let fluid = CarreauYasudaBlood::<F>::normal_blood();
    let mut net = Network::new(graph, fluid);

    let parent_d = 120.0e-6;
    let branch_b_d = 70.0e-6;
    let trunk_d = 100.0e-6;
    let branch_c_d = 80.0e-6;
    let branch_d_d = 55.0e-6;
    let length = 6.0e-3;
    let feed_hct = 0.45;
    let q_parent = 4.0e-9;
    let q_b = 1.0e-9;
    let q_trunk = 3.0e-9;
    let q_c = 2.0e-9;
    let q_d = 1.0e-9;

    net.set_flow_rate(parent, q_parent);
    net.set_flow_rate(branch_b, q_b);
    net.set_flow_rate(trunk, q_trunk);
    net.set_flow_rate(branch_c, q_c);
    net.set_flow_rate(branch_d, q_d);

    let make_props =
        |id: &str, diameter: f64, properties: std::collections::HashMap<String, f64>| {
            EdgeProperties {
                id: id.to_string(),
                component_type: cfd_1d::ComponentType::Pipe,
                length,
                area: std::f64::consts::PI * diameter * diameter / 4.0,
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
                properties,
            }
        };

    let mut parent_props = std::collections::HashMap::new();
    parent_props.insert(EDGE_PROPERTY_HEMATOCRIT.to_string(), feed_hct);
    net.add_edge_properties(parent, make_props("parent", parent_d, parent_props));
    net.add_edge_properties(
        branch_b,
        make_props("b", branch_b_d, std::collections::HashMap::new()),
    );
    net.add_edge_properties(
        trunk,
        make_props("trunk", trunk_d, std::collections::HashMap::new()),
    );
    net.add_edge_properties(
        branch_c,
        make_props("c", branch_c_d, std::collections::HashMap::new()),
    );
    net.add_edge_properties(
        branch_d,
        make_props("d", branch_d_d, std::collections::HashMap::new()),
    );

    net.update_resistances()?;

    let first_split = pries_phase_separation(
        feed_hct,
        q_b / (q_b + q_trunk),
        branch_b_d * 1.0e6,
        trunk_d * 1.0e6,
        parent_d * 1.0e6,
    )?;
    let trunk_hct = ((1.0 - first_split.cell_fraction) * feed_hct / (q_trunk / (q_b + q_trunk)))
        .clamp(0.0, 1.0);
    let second_split = pries_phase_separation(
        trunk_hct,
        q_c / (q_c + q_d),
        branch_c_d * 1.0e6,
        branch_d_d * 1.0e6,
        trunk_d * 1.0e6,
    )?;
    let expected_d =
        ((1.0 - second_split.cell_fraction) * trunk_hct / (q_d / (q_c + q_d))).clamp(0.0, 1.0);

    let local_trunk = net
        .properties
        .get(&trunk)
        .and_then(|props| props.properties.get(EDGE_PROPERTY_LOCAL_HEMATOCRIT))
        .copied()
        .expect("trunk should receive transported hematocrit");
    let local_c = net
        .properties
        .get(&branch_c)
        .and_then(|props| props.properties.get(EDGE_PROPERTY_LOCAL_HEMATOCRIT))
        .copied()
        .expect("branch c should receive transported hematocrit");
    let local_d = net
        .properties
        .get(&branch_d)
        .and_then(|props| props.properties.get(EDGE_PROPERTY_LOCAL_HEMATOCRIT))
        .copied()
        .expect("branch d should receive transported hematocrit");

    assert!((local_trunk - trunk_hct).abs() < 1.0e-9);
    assert!((local_c - second_split.daughter_hematocrit).abs() < 1.0e-9);
    assert!((local_d - expected_d).abs() < 1.0e-9);
    Ok(())
}

#[test]
fn reconverging_transport_mixes_incoming_rbc_fluxes() -> Result<()> {
    type F = f64;
    let mut builder = NetworkBuilder::<F>::new();
    let inlet_a = builder.add_inlet("in_a".to_string());
    let inlet_b = builder.add_inlet("in_b".to_string());
    let merge = builder.add_junction("merge".to_string());
    let outlet = builder.add_outlet("out".to_string());
    let edge_a = builder.add_edge(inlet_a, merge, Edge::new("a".to_string(), EdgeType::Pipe));
    let edge_b = builder.add_edge(inlet_b, merge, Edge::new("b".to_string(), EdgeType::Pipe));
    let outlet_edge = builder.add_edge(merge, outlet, Edge::new("out".to_string(), EdgeType::Pipe));
    let graph = builder.build().expect("test invariant");
    let fluid = CarreauYasudaBlood::<F>::normal_blood();
    let mut net = Network::new(graph, fluid);

    let diameter = 90.0e-6;
    let length = 5.0e-3;
    let q_a = 1.0e-9;
    let q_b = 3.0e-9;
    let hct_a = 0.20;
    let hct_b = 0.45;

    net.set_flow_rate(edge_a, q_a);
    net.set_flow_rate(edge_b, q_b);
    net.set_flow_rate(outlet_edge, q_a + q_b);

    let make_props =
        |id: &str, properties: std::collections::HashMap<String, f64>| EdgeProperties {
            id: id.to_string(),
            component_type: cfd_1d::ComponentType::Pipe,
            length,
            area: std::f64::consts::PI * diameter * diameter / 4.0,
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
            properties,
        };

    let mut props_a = std::collections::HashMap::new();
    props_a.insert(EDGE_PROPERTY_HEMATOCRIT.to_string(), hct_a);
    let mut props_b = std::collections::HashMap::new();
    props_b.insert(EDGE_PROPERTY_HEMATOCRIT.to_string(), hct_b);
    net.add_edge_properties(edge_a, make_props("a", props_a));
    net.add_edge_properties(edge_b, make_props("b", props_b));
    net.add_edge_properties(
        outlet_edge,
        make_props("out", std::collections::HashMap::new()),
    );

    net.update_resistances()?;

    let expected_mixed = (q_a * hct_a + q_b * hct_b) / (q_a + q_b);
    let local_out = net
        .properties
        .get(&outlet_edge)
        .and_then(|props| props.properties.get(EDGE_PROPERTY_LOCAL_HEMATOCRIT))
        .copied()
        .expect("outlet edge should receive mixed hematocrit");

    assert!((local_out - expected_mixed).abs() < 1.0e-9);
    Ok(())
}
