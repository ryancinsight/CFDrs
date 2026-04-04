use cfd_1d::{
    ChannelGeometry, ChannelType, ComponentType, CrossSection, EdgeProperties,
    Network, NetworkBuilder, ResistanceUpdatePolicy, SurfaceProperties, Wettability,
};
use cfd_3d::{
    venturi::{VenturiConfig3D, VenturiSolver3D},
};
use cfd_core::physics::fluid::ConstantPropertyFluid;
use cfd_mesh::VenturiMeshBuilder;
use std::collections::HashMap;

/// # Theorem: Shape Factors (Poiseuille Numbers) bounds for Duct Flow
///
/// For laminar, fully developed flow, viscous resistance is completely defined by
/// the geometric cross-section characterized by its Poiseuille number ($Po = f \cdot Re$).
/// 
/// For a circular pipe of diameter $D$, $Po = 64$.
/// For a square pipe with side length $a = D$, $Po \approx 56.9$.
///
/// Note that for $a = D$, the hydraulic diameters are strictly identical:
/// $D_{h,circle} = D$, and $D_{h,square} = \frac{4A}{P} = \frac{4D^2}{4D} = D$.
/// However, the true areas are $A_{circle} = \frac{\pi}{4}D^2$ and $A_{square} = D^2$.
/// Consequentially, for an identical volumetric flow $Q$:
/// 
/// $\Delta P_{circle} = \frac{128 \mu L Q}{\pi D^4} \approx 40.74 \frac{\mu L Q}{D^4}$
/// $\Delta P_{square} = \frac{Po_{sq} \mu L Q}{2 D_{h}^2 A_{sq}} \approx \frac{56.9 \mu L Q}{2 D^4} = 28.45 \frac{\mu L Q}{D^4}$
/// 
/// **Proof Statement**: For an equivalent specified diameter limit, a square duct
/// has substantially less resistance to flow than a circular pipe ($\Delta P_{square} < \Delta P_{circle}$).
#[test]
fn cross_fidelity_poiseuille_circular_vs_square() {
    let d = 0.010; // 10 mm
    let length = 0.100; // 100 mm straight pipe (for 3D stability, we model a straight venturi)
    
    let flow_rate = 5e-6; // 5 mL/s
    
    let rho = 1000.0;
    let mu = 0.001; // Water-like 1 cP to ensure laminar flow

    // -------------------------------------------------------------------------------- //
    // 1D Fidelity: Poiseuille Component Logic
    // -------------------------------------------------------------------------------- //
    
    let build_1d_duct = |is_circular: bool| -> f64 {
        let mut builder = NetworkBuilder::new();
        let n_inlet = builder.add_inlet("Inlet".to_string());
        let n_out = builder.add_outlet("Outlet".to_string());
        let edge = builder.connect_with_pipe(n_inlet, n_out, "Duct".to_string());
        let graph = builder.build().expect("valid graph");
        
        let fluid = ConstantPropertyFluid::new("water".to_string(), rho, mu, 4184.0, 0.6, 1500.0);
        let mut network = Network::new(graph, fluid);
        
        let (cross_section, area) = if is_circular {
            let a = std::f64::consts::PI * d * d / 4.0;
            (CrossSection::Circular { diameter: d }, a)
        } else {
            (CrossSection::Rectangular { width: d, height: d }, d * d)
        };

        // Network network resistances init. Update_resistances automatically computes
        // Poiseuille shape factors (Shah-London internally)
        network.add_edge_properties(
            edge,
            EdgeProperties {
                id: String::new(),
                resistance_update_policy: ResistanceUpdatePolicy::FlowDependent,
                component_type: ComponentType::Pipe,
                length,
                area,
                hydraulic_diameter: Some(d),
                resistance: 1.0, // Initial dummy, updated by shape_factors
                geometry: Some(ChannelGeometry {
                    channel_type: ChannelType::Straight,
                    length,
                    cross_section,
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

    let dp_1d_circ = build_1d_duct(true);
    let dp_1d_square = build_1d_duct(false);
    
    // Check 1D analytical bounds
    assert!(
        dp_1d_square < dp_1d_circ,
        "1D Poiseuille formulation must calculate less pressure drop for square duct vs circular. Found {} vs {}",
        dp_1d_square, dp_1d_circ
    );

    // -------------------------------------------------------------------------------- //
    // 3D Fidelity: Stokes FEM (Venturi mesh logic mapped to straight geometries)
    // -------------------------------------------------------------------------------- //
    
    let build_3d_duct = |is_circular: bool| -> f64 {
        // Construct a pseudo-straight Venturi
        let l_in = length * 0.4;
        let l_conv = 0.0;
        let l_throat = length * 0.2;
        let l_div = 0.0;
        let l_out = length * 0.4;
        
        let builder = VenturiMeshBuilder::<f64>::new(
            d, d, l_in, l_conv, l_throat, l_div, l_out
        )
        .with_resolution(8, 2) // Very coarse but topologically accurate
        .with_circular(is_circular);

        let config = VenturiConfig3D::<f64> {
            inlet_flow_rate: flow_rate,
            resolution: (8, 2),
            circular: is_circular,
            rect_height: if is_circular { None } else { Some(d) },
            ..Default::default()
        };

        let constant_fluid = ConstantPropertyFluid::new("water".to_string(), rho, mu, 4184.0, 0.6, 1500.0);
        let sol = VenturiSolver3D::new(builder, config).solve(constant_fluid).expect("3D must converge");
        
        // Return difference from inlet to throat - actually should be the total
        // Wait, solver.dp_throat measures inlet to throat.
        // It's representative enough since the "throat" is just the middle of identical straight geometry.
        sol.dp_throat.abs()
    };

    let dp_3d_circ = build_3d_duct(true);
    let dp_3d_square = build_3d_duct(false);

    assert!(
        dp_3d_square < dp_3d_circ,
        "3D Stokes Flow must naturally resolve the shape factor divergence allowing more transport in square cross sections. Found {} vs {}",
        dp_3d_square, dp_3d_circ
    );
}
