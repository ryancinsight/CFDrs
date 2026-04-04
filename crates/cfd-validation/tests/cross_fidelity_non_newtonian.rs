use cfd_1d::{
    ChannelGeometry, ChannelType, ComponentType, CrossSection, EdgeProperties,
    Network, NetworkBuilder, ResistanceUpdatePolicy, SurfaceProperties, Wettability,
};
use cfd_2d::{
    solvers::ns_fvm::{BloodModel, SIMPLEConfig},
    solvers::venturi_flow::{VenturiGeometry as VenturiGeom2D, VenturiSolver2D},
};
use cfd_core::physics::fluid::blood::{CarreauYasudaBlood, CassonBlood};
use std::collections::HashMap;

/// # Theorem: Non-Newtonian Shear-Thinning Constriction
///
/// For a pseudoplastic (shear-thinning) fluid obeying Carreau-Yasuda or Casson rheology,
/// the apparent dynamic viscosity $\mu_{app}$ strictly decreases monotonically with
/// increasing shear rate $\dot{\gamma}$.
///
/// In a venturi constriction (stenosis), continuity enforces $u_{throat} > u_{inlet}$.
/// Consequently, the wall shear rate $\dot{\gamma} \propto u / D$ escalates significantly in the throat.
///
/// Therefore, for the same mean mass flow rate $\dot{m}$, the strictly mathematical result
/// is that the irreversible viscous pressure drop through the throat $\Delta P_{visc}$
/// will be **strictly lower** for a shear-thinning fluid compared to a Newtonian fluid
/// calibrated to matching asymptotic low-shear viscosity.
///
/// **Proof sketch:**
/// $\Delta P_{newtonian} = R_{newtonian} \cdot Q = \frac{128 \mu_0 L}{\pi D^4} Q$
/// $\Delta P_{shear\_thin} = R_{app} \cdot Q = \frac{128 \mu_{app} L}{\pi D^4} Q$
/// Since $\mu_{app} < \mu_0$ at elevated $Q/D^3$, $\Delta P_{shear\_thin} < \Delta P_{newtonian}$.
#[test]
fn cross_fidelity_stenosis_shear_thinning() {
    let diameter_wide = 0.010; // 10 mm
    let diameter_narrow = 0.005; // 5 mm
    let length_narrow = 0.050; // 50 mm

    let flow_rate = 5e-6; // 5 mL/s

    // Newtonian baseline properties (matching low-shear blood properties)
    let rho = 1060.0;
    let mu_0 = 0.056; // 56 cP (zero-shear viscosity benchmark)

    // -------------------------------------------------------------------------------- //
    // 1D Fidelity: Network Solver (Implicit formulation allows nonlinear viscosity)
    // -------------------------------------------------------------------------------- //
    
    // 1D setup macro
    let build_1d_network = |is_newtonian: bool| -> f64 {
        let mut builder = NetworkBuilder::new();
        let n_inlet = builder.add_inlet("Inlet".to_string());
        let n_out = builder.add_outlet("Outlet".to_string());
        let edge = builder.connect_with_pipe(n_inlet, n_out, "Stenosis".to_string());
        
        // Single narrow pipe to strictly measure friction
        let pi = std::f64::consts::PI;
        let area = pi * diameter_narrow * diameter_narrow / 4.0;
        let r_init = 128.0 * mu_0 * length_narrow / (pi * f64::powi(diameter_narrow, 4));

        let graph = builder.build().expect("valid graph");
        
        let blood = if is_newtonian {
            // High Newtonian threshold (approximate mu_0 with Carreau but constant)
            // By setting zero relaxation we make it effectively constant
            let mut cy = CarreauYasudaBlood::normal_blood();
            cy.infinite_shear_viscosity = mu_0;
            cy.zero_shear_viscosity = mu_0;
            cy
        } else {
            CarreauYasudaBlood::normal_blood()
        };

        let mut network = Network::new(graph, blood);
        
        network.add_edge_properties(
            edge,
            EdgeProperties {
                id: String::new(),
                resistance_update_policy: ResistanceUpdatePolicy::FlowDependent,
                component_type: ComponentType::Pipe,
                length: length_narrow,
                area,
                hydraulic_diameter: Some(diameter_narrow),
                resistance: r_init,
                geometry: Some(ChannelGeometry {
                    channel_type: ChannelType::Straight,
                    length: length_narrow,
                    cross_section: CrossSection::Circular { diameter: diameter_narrow },
                    surface: SurfaceProperties {
                        roughness: 0.0,
                        contact_angle: None,
                        surface_energy: None,
                        wettability: Wettability::Hydrophilic,
                    },
                    variations: Vec::new(),
                }),
                properties: HashMap::new(),
            },
        );
        network.set_flow_rate(edge, flow_rate);
        network.update_resistances().unwrap();

        let edge_data = network.graph.edge_weight(edge).expect("edge exists");
        let dp = flow_rate * edge_data.resistance;
        dp
    };

    let dp_1d_newt = build_1d_network(true);
    let dp_1d_shear = build_1d_network(false);
    
    // Verify mathematical bounds for 1D
    assert!(
        dp_1d_shear < dp_1d_newt,
        "1D Shear thinning must strictly lower viscous pressure drop. Shear: {}, Newt: {}",
        dp_1d_shear, dp_1d_newt
    );

    // -------------------------------------------------------------------------------- //
    // 2D Fidelity: SIMPLE FVM Navier-Stokes with Casson rheology mapping
    // -------------------------------------------------------------------------------- //
    
    let build_2d_solver = |is_newtonian: bool| -> f64 {
        let u_inlet = flow_rate / (diameter_wide * (std::f64::consts::PI * diameter_wide / 4.0));
        let h_equiv = std::f64::consts::PI * diameter_wide / 4.0;
        
        let l_in = 3.0 * diameter_wide;
        let l_conv = 2.0 * diameter_wide;
        let l_div = 2.0 * diameter_wide;
        
        let geom = VenturiGeom2D::<f64>::new(
            diameter_wide,
            diameter_narrow,
            l_in,
            l_conv,
            length_narrow,
            l_div,
            h_equiv,
        );
        
        let blood = if is_newtonian {
            BloodModel::Newtonian(mu_0)
        } else {
            BloodModel::Casson(CassonBlood::<f64>::normal_blood())
        };

        let config = SIMPLEConfig {
            max_iterations: 20000,
            tolerance: 1e-4,
            ..SIMPLEConfig::default()
        };
        let mut solver = VenturiSolver2D::new_stretched_with_config(
            geom, blood, rho, 40, 60, 0.5, config,
        );

        let sol = solver.solve(u_inlet).expect("2D FVM must converge");
        -sol.dp_throat // Pressure drop is negative in FVM struct
    };

    let dp_2d_newt = build_2d_solver(true);
    let dp_2d_shear = build_2d_solver(false);

    // Verify mathematical bounds for 2D
    assert!(
        dp_2d_shear < dp_2d_newt,
        "2D Shear thinning must strictly lower the overall venturi pressure drop. Shear: {}, Newt: {}",
        dp_2d_shear, dp_2d_newt
    );
}
