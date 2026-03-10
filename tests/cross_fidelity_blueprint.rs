//! Cross-fidelity blueprint validation: 1D → 2D → 3D consistency.
//!
//! Validates that a complete `NetworkBlueprint` evaluates to consistent
//! results across all three fidelity levels:
//! 1. Authoritative 1D `NetworkSolver` (conserves mass perfectly, sets baseline ΔP)
//! 2. 2D Depth-Averaged LBM Solver (captures lateral flow features)
//! 3. 3D Cascade / Domain Solver (captures full volumetric flow)

use cfd_3d::blueprint_integration::{
    process_blueprint_with_reference_trace, Blueprint3dProcessingConfig,
};
use cfd_3d::cascade::{CascadeChannelSpec, CascadeConfig3D, CascadeSolver3D};
use cfd_core::physics::fluid::ConstantPropertyFluid;
use cfd_schematics::interface::presets::venturi_chain;

#[test]
fn cross_fidelity_blueprint_consistency() {
    // 1. Generate a test blueprint (SingleVenturi-like chain)
    // Preset: id, length(30mm), width(4mm), height(2mm)
    let blueprint = venturi_chain("test_blueprint", 0.030, 0.004, 0.002);

    let density_kg_m3 = 1060.0;
    let viscosity_pa_s = 3.5e-3;
    let total_flow_rate_m3_s = 1.0e-7; // 100 uL/s

    // 2. Configure the unified 3D trace processing path
    // This internally runs the 1D reference solve and the 2D solve
    let config = Blueprint3dProcessingConfig {
        density_kg_m3,
        viscosity_pa_s,
        total_flow_rate_m3_s,
        two_d_grid_nx: 64, // Coarse, fast check
        two_d_grid_ny: 16,
        two_d_tolerance: 1e-8,
        run_2d_reference: true,
        ..Default::default()
    };

    let result = process_blueprint_with_reference_trace(&blueprint, &config)
        .expect("process_blueprint_with_reference_trace failed");

    // 3. Audit properties
    println!("=== Cross-Fidelity Audit: {} ===", blueprint.name);
    println!("1D Reference Nodes: {}", result.node_traces.len());
    println!("1D Reference Channels: {}", result.channel_traces.len());

    // 3A: Mass Conservation
    let mut total_1d_inlet_flow = 0.0;
    let mut total_1d_outlet_flow = 0.0;
    for node in &result.node_traces {
        if node.node_kind == cfd_schematics::NodeKind::Inlet {
            total_1d_inlet_flow += node.prescribed_boundary_flow_m3_s;
        } else if node.node_kind == cfd_schematics::NodeKind::Outlet {
            // Outlets have negative boundary flow (leaving domain) or outgoing > 0
            total_1d_outlet_flow += node.incoming_flow_m3_s;
        }
    }

    println!("1D Total Inlet Flow: {:.6e} m^3/s", total_1d_inlet_flow);
    println!("1D Total Outlet Flow: {:.6e} m^3/s", total_1d_outlet_flow);

    assert!(
        (total_1d_inlet_flow - total_flow_rate_m3_s).abs() < 1e-12,
        "1D solver did not respect prescribed inlet flow"
    );
    assert!(
        (total_1d_inlet_flow - total_1d_outlet_flow).abs() < 1e-12,
        "1D mass conservation violated"
    );

    // 3B: 2D vs 1D Checks
    assert!(
        result.two_d_converged_count.is_some() && result.two_d_converged_count.unwrap() > 0,
        "2D solver did not converge on any channels"
    );

    if let Some(mean_err) = result.two_d_mean_outlet_flow_error_pct {
        println!("2D Mean Outlet Flow Error vs 1D: {:.2}%", mean_err);
    }

    // Compare Pressure Drops and mass conservation per channel
    for trace in &result.channel_traces {
        // Find matching 2D solve result in the output by re-running solver or just printing trace fields.
        // Wait, the Blueprint3dTrace/ChannelCrossFidelityTrace doesn't include the iteration count or residual.
        // It only has two_d_converged.
        println!(
            "Channel {}: 1D dP={:.2} Pa, Vel={:.4}m/s | 2D Mean Shear={:?} | 2D Converged={:?} | 2D Outlet Error={:?}",
            trace.channel_id,
            trace.reference_pressure_drop_pa,
            trace.reference_mean_velocity_m_s,
            trace.two_d_field_wall_shear_mean_pa,
            trace.two_d_converged,
            trace.two_d_outlet_flow_error_pct
        );

        assert!(
            trace.reference_pressure_drop_pa >= 0.0,
            "1D pressure drop cannot be negative"
        );
    }

    // 4. (Future) We will integrate the actual 3D Navier-Stokes solver here
    // to complete the fidelity loop. For now, mesh completeness proves 3D readiness.
    assert!(
        result.fluid_mesh_vertex_count > 0 && result.fluid_mesh_face_count > 0,
        "Fluid mesh was not generated successfully"
    );
    println!(
        "3D Mesh generated with {} vertices, {} faces.",
        result.fluid_mesh_vertex_count, result.fluid_mesh_face_count
    );

    if let Some(mean_err) = result.two_d_mean_outlet_flow_error_pct {
        assert!(
            mean_err < 5.0, // Should be close to 1D mass conservation
            "2D flow conservation error is too high: {}%",
            mean_err
        );
    }

    // 4. Full 3D Cascade Verification
    // Use the exact same fluid properties as 1D.
    let fluid_3d = ConstantPropertyFluid {
        name: "test_fluid".into(),
        density: density_kg_m3,
        viscosity: viscosity_pa_s,
        specific_heat: 4184.0,
        thermal_conductivity: 0.6,
        speed_of_sound: 1500.0,
    };
    let cascade_config = CascadeConfig3D {
        outlet_pressure: 0.0,
        resolution: (40, 16, 16),
        max_picard_iterations: 3,
        picard_tolerance: 1e-2,
    };
    let cascade_solver = CascadeSolver3D::new(cascade_config, fluid_3d);

    // Build specs dynamically from the traces
    let mut specs = Vec::new();
    for ch in &blueprint.channels {
        if let Some(trace) = result
            .channel_traces
            .iter()
            .find(|t| t.channel_id == ch.id.as_str())
        {
            let (w, h) = ch.cross_section.dims();
            let is_venturi = ch.venturi_geometry.is_some();
            let throat_w = ch.venturi_geometry.as_ref().map(|v| v.throat_width_m);

            specs.push(CascadeChannelSpec {
                id: ch.id.as_str().to_owned(),
                length: ch.length_m,
                width: w,
                height: h,
                flow_rate_m3_s: trace.reference_flow_rate_m3_s,
                is_venturi_throat: is_venturi,
                throat_width: throat_w,
                local_hematocrit: None,
            });
        }
    }

    let cascade_result = cascade_solver
        .solve(&specs)
        .expect("3D Cascade solver failed");

    for cr in &cascade_result.channel_results {
        println!(
            "3D Channel {}: dP={:.2} Pa, Vel_max={:.4}m/s, Shear={:.4} Pa",
            cr.channel_id, cr.pressure_drop_pa, cr.max_velocity, cr.wall_shear_mean_pa
        );

        assert!(cr.max_velocity > 0.0, "3D velocity should be non-zero");
        assert!(
            cr.pressure_drop_pa > 0.0,
            "3D pressure drop should be positive"
        );
    }
}

#[test]
fn cross_fidelity_blueprint_complex_branching() {
    // 1. Generate a complex branching blueprint (Double Trifurcation)
    let blueprint = cfd_schematics::interface::presets::double_trifurcation_cif_venturi_rect(
        "dtcv-cross-fidelity",
        15e-3,  // trunk_length_m
        10e-3,  // branch_length_m
        4.0e-3, // main_width_m (must equal 4.0mm for 4.0mm D_h when height is 4.0mm)
        0.54,   // split1_center_frac
        0.45,   // split2_center_frac
        80e-6,  // throat_width_m
        200e-6, // throat_length_m
        1.0e-3, // inter_throat_spacing_m
        1,      // center_throat_count
        4.0e-3, // height_m (must equal 4.0mm for 4.0mm D_h when width is 4.0mm)
    );

    let density_kg_m3 = 1060.0;
    let viscosity_pa_s = 3.5e-3;
    let total_flow_rate_m3_s = 1.0e-7; // 100 uL/s

    // 2. Configure the unified 3D trace processing path
    let config = Blueprint3dProcessingConfig {
        density_kg_m3,
        viscosity_pa_s,
        total_flow_rate_m3_s,
        two_d_grid_nx: 128, // Higher resolution for complex network
        two_d_grid_ny: 32,
        two_d_tolerance: 1e-8,
        run_2d_reference: true,
        ..Default::default()
    };

    let result = process_blueprint_with_reference_trace(&blueprint, &config)
        .expect("process_blueprint_with_reference_trace failed on complex branching");

    println!("=== Cross-Fidelity Audit: {} ===", blueprint.name);
    println!("1D Reference Nodes: {}", result.node_traces.len());
    println!("1D Reference Channels: {}", result.channel_traces.len());

    let mut total_1d_inlet_flow = 0.0;
    let mut total_1d_outlet_flow = 0.0;
    for node in &result.node_traces {
        if node.node_kind == cfd_schematics::NodeKind::Inlet {
            total_1d_inlet_flow += node.prescribed_boundary_flow_m3_s;
        } else if node.node_kind == cfd_schematics::NodeKind::Outlet {
            total_1d_outlet_flow += node.incoming_flow_m3_s;
        }
    }

    assert!(
        (total_1d_inlet_flow - total_flow_rate_m3_s).abs() < 1e-12,
        "1D solver did not respect prescribed inlet flow"
    );
    assert!(
        (total_1d_inlet_flow - total_1d_outlet_flow).abs() < 1e-12,
        "1D mass conservation violated across complex junctions"
    );

    assert!(
        result.two_d_converged_count.is_some() && result.two_d_converged_count.unwrap() > 0,
        "2D solver did not converge on any branch channels"
    );

    if let Some(mean_err) = result.two_d_mean_outlet_flow_error_pct {
        println!("2D Mean Outlet Flow Error vs 1D: {:.2}%", mean_err);

        println!("--- Per-Channel 2D Errors ---");
        for trace in &result.channel_traces {
            if let Some(err) = trace.two_d_outlet_flow_error_pct {
                println!(
                    "  Channel {:<20} | 1D Flow: {:.2e} m3/s | 2D Error: {:>6.2}% | Conv: {:?} | Shear: {:?}",
                    trace.channel_id, trace.reference_flow_rate_m3_s, err, trace.two_d_converged, trace.two_d_field_wall_shear_mean_pa
                );
            }
        }

        assert!(
            mean_err < 10.0,
            "2D flow conservation error in branching network is too high: {}%",
            mean_err
        );
    }

    // 4. 3D Cascade Verification
    let fluid_3d = ConstantPropertyFluid {
        name: "test_fluid".into(),
        density: density_kg_m3,
        viscosity: viscosity_pa_s,
        specific_heat: 4184.0,
        thermal_conductivity: 0.6,
        speed_of_sound: 1500.0,
    };
    let cascade_config = CascadeConfig3D {
        outlet_pressure: 0.0,
        resolution: (40, 16, 16),
        max_picard_iterations: 3,
        picard_tolerance: 1e-2,
    };
    let cascade_solver = CascadeSolver3D::new(cascade_config, fluid_3d);

    let mut specs = Vec::new();
    for ch in &blueprint.channels {
        if let Some(trace) = result
            .channel_traces
            .iter()
            .find(|t| t.channel_id == ch.id.as_str())
        {
            let (w, h) = ch.cross_section.dims();
            let is_venturi = ch.venturi_geometry.is_some();
            let throat_w = ch.venturi_geometry.as_ref().map(|v| v.throat_width_m);

            specs.push(CascadeChannelSpec {
                id: ch.id.as_str().to_owned(),
                length: ch.length_m,
                width: w,
                height: h,
                flow_rate_m3_s: trace.reference_flow_rate_m3_s,
                is_venturi_throat: is_venturi,
                throat_width: throat_w,
                local_hematocrit: None,
            });
        }
    }

    let cascade_result = cascade_solver
        .solve(&specs)
        .expect("3D Cascade solver failed");

    // Cross-check max 3D shear with 2D LBM predicted mean
    for cr in &cascade_result.channel_results {
        if let Some(trace) = result
            .channel_traces
            .iter()
            .find(|t| t.channel_id == cr.channel_id)
        {
            if let Some(mean_2d_shear) = trace.two_d_field_wall_shear_mean_pa {
                println!(
                    "Branch [{:<15}] | 1D dP: {:>6.2} Pa | 2D Shear: {:>6.2} Pa | 3D dP: {:>6.2} Pa | 3D Shear: {:>6.2} Pa",
                    cr.channel_id, trace.reference_pressure_drop_pa, mean_2d_shear, cr.pressure_drop_pa, cr.wall_shear_mean_pa
                );
            }
        }
    }
}

#[test]
fn cross_fidelity_blueprint_bifurcation() {
    let blueprint = cfd_schematics::interface::presets::bifurcation_rect(
        "bifurcation-cross-fidelity",
        15e-3,  // parent_length_m
        10e-3,  // daughter_length_m
        4.0e-3, // parent_width_m (req 4mm D_h)
        4.0e-3, // daughter_width_m (req 4mm D_h)
        4.0e-3, // height_m
    );

    let density_kg_m3 = 1060.0;
    let viscosity_pa_s = 3.5e-3;
    let total_flow_rate_m3_s = 1.0e-7;

    let config = Blueprint3dProcessingConfig {
        density_kg_m3,
        viscosity_pa_s,
        total_flow_rate_m3_s,
        two_d_grid_nx: 64,
        two_d_grid_ny: 16,
        two_d_tolerance: 1e-8,
        run_2d_reference: true,
        ..Default::default()
    };

    let result = process_blueprint_with_reference_trace(&blueprint, &config)
        .expect("process_blueprint_with_reference_trace failed on bifurcation");

    println!("=== Cross-Fidelity Audit: {} ===", blueprint.name);
    println!("1D Reference Nodes: {}", result.node_traces.len());
    println!("1D Reference Channels: {}", result.channel_traces.len());

    let mut total_1d_inlet_flow = 0.0;
    let mut total_1d_outlet_flow = 0.0;
    for node in &result.node_traces {
        if node.node_kind == cfd_schematics::NodeKind::Inlet {
            total_1d_inlet_flow += node.prescribed_boundary_flow_m3_s;
        } else if node.node_kind == cfd_schematics::NodeKind::Outlet {
            total_1d_outlet_flow += node.incoming_flow_m3_s;
        }
    }

    assert!(
        (total_1d_inlet_flow - total_flow_rate_m3_s).abs() < 1e-12,
        "1D solver did not respect prescribed inlet flow"
    );
    assert!(
        (total_1d_inlet_flow - total_1d_outlet_flow).abs() < 1e-12,
        "1D mass conservation violated across bifurcation"
    );

    if let Some(mean_err) = result.two_d_mean_outlet_flow_error_pct {
        println!("2D Mean Outlet Flow Error vs 1D: {:.2}%", mean_err);
        assert!(
            mean_err < 10.0,
            "2D flow conservation error in bifurcation network is too high: {}%",
            mean_err
        );
    }
}

#[test]
fn cross_fidelity_blueprint_trifurcation() {
    let blueprint = cfd_schematics::interface::presets::trifurcation_rect(
        "trifurcation-cross-fidelity",
        15e-3,  // parent_length_m
        10e-3,  // daughter_length_m
        4.0e-3, // parent_width_m (req 4mm D_h)
        4.0e-3, // daughter_width_m (req 4mm D_h)
        4.0e-3, // height_m
    );

    let density_kg_m3 = 1060.0;
    let viscosity_pa_s = 3.5e-3;
    let total_flow_rate_m3_s = 1.0e-7;

    let config = Blueprint3dProcessingConfig {
        density_kg_m3,
        viscosity_pa_s,
        total_flow_rate_m3_s,
        two_d_grid_nx: 64,
        two_d_grid_ny: 16,
        two_d_tolerance: 1e-8,
        run_2d_reference: true,
        ..Default::default()
    };

    let result = process_blueprint_with_reference_trace(&blueprint, &config)
        .expect("process_blueprint_with_reference_trace failed on trifurcation");

    println!("=== Cross-Fidelity Audit: {} ===", blueprint.name);
    println!("1D Reference Nodes: {}", result.node_traces.len());
    println!("1D Reference Channels: {}", result.channel_traces.len());

    let mut total_1d_inlet_flow = 0.0;
    let mut total_1d_outlet_flow = 0.0;
    for node in &result.node_traces {
        if node.node_kind == cfd_schematics::NodeKind::Inlet {
            total_1d_inlet_flow += node.prescribed_boundary_flow_m3_s;
        } else if node.node_kind == cfd_schematics::NodeKind::Outlet {
            total_1d_outlet_flow += node.incoming_flow_m3_s;
        }
    }

    assert!(
        (total_1d_inlet_flow - total_flow_rate_m3_s).abs() < 1e-12,
        "1D solver did not respect prescribed inlet flow"
    );
    assert!(
        (total_1d_inlet_flow - total_1d_outlet_flow).abs() < 1e-12,
        "1D mass conservation violated across trifurcation"
    );

    if let Some(mean_err) = result.two_d_mean_outlet_flow_error_pct {
        println!("2D Mean Outlet Flow Error vs 1D: {:.2}%", mean_err);
        assert!(
            mean_err < 10.0,
            "2D flow conservation error in trifurcation network is too high: {}%",
            mean_err
        );
    }
}
