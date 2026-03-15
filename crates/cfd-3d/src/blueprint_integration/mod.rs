//! Canonical blueprint-driven 3D preprocessing and cross-fidelity tracing.

use std::collections::HashMap;

use cfd_2d::network::{solve_reference_trace, Network2dBuilderSink};
use cfd_core::error::{Error, Result as CfdResult};
use cfd_core::physics::fluid::BloodModel;
use cfd_mesh::application::pipeline::TopologyClass;
use cfd_mesh::application::pipeline::{
    BlueprintMeshPipeline, PipelineConfig, PipelineVolumeTrace, SegmentCenterline,
};
use cfd_schematics::application::ports::GraphSink;
use cfd_schematics::{NetworkBlueprint, NodeKind};

/// Configuration for the canonical `NetworkBlueprint -> cfd-mesh -> cfd-3d` trace path.
#[derive(Debug, Clone)]
pub struct Blueprint3dProcessingConfig {
    /// Mesh-pipeline configuration used for `cfd-mesh` conversion.
    pub mesh: PipelineConfig,
    /// Fluid density used for cfd-1d/cfd-2d reference solves [kg/m^3].
    pub density_kg_m3: f64,
    /// Dynamic viscosity used for cfd-1d/cfd-2d reference solves [Pa·s].
    pub viscosity_pa_s: f64,
    /// Target total inlet flow rate for scaled cross-fidelity references [m^3/s].
    pub total_flow_rate_m3_s: f64,
    /// Structured-grid X resolution for optional cfd-2d comparison.
    pub two_d_grid_nx: usize,
    /// Structured-grid Y resolution for optional cfd-2d comparison.
    pub two_d_grid_ny: usize,
    /// Solver tolerance for the optional cfd-2d comparison.
    pub two_d_tolerance: f64,
    /// Whether to run the full cfd-2d comparison in addition to the cfd-1d trace.
    pub run_2d_reference: bool,
}

impl Default for Blueprint3dProcessingConfig {
    fn default() -> Self {
        Self {
            mesh: PipelineConfig::default(),
            density_kg_m3: 1060.0,
            viscosity_pa_s: 3.5e-3,
            total_flow_rate_m3_s: 1.0e-9,
            two_d_grid_nx: 32,
            two_d_grid_ny: 12,
            two_d_tolerance: 1e-6,
            run_2d_reference: true,
        }
    }
}

/// Per-channel trace spanning schematic, mesh, cfd-1d, and optional cfd-2d diagnostics.
#[derive(Debug, Clone)]
pub struct ChannelCrossFidelityTrace {
    /// Blueprint channel identifier.
    pub channel_id: String,
    /// Upstream blueprint node identifier.
    pub from_node_id: String,
    /// Downstream blueprint node identifier.
    pub to_node_id: String,
    /// Schematic volume carried by this channel [mm^3].
    pub schematic_volume_mm3: f64,
    /// Meshed 3D path volume for this channel [mm^3].
    pub meshed_volume_mm3: f64,
    /// Relative 3D meshed volume error against the schematic contract [%].
    pub mesh_volume_error_pct: f64,
    /// Authoritative cfd-1d reference flow rate [m^3/s].
    pub reference_flow_rate_m3_s: f64,
    /// Authoritative cfd-1d pressure drop [Pa].
    pub reference_pressure_drop_pa: f64,
    /// Authoritative cfd-1d mean velocity [m/s].
    pub reference_mean_velocity_m_s: f64,
    /// cfd-2d outlet-flow error against the cfd-1d reference [%], when computed.
    pub two_d_outlet_flow_error_pct: Option<f64>,
    /// Mean wall shear extracted from the cfd-2d field [Pa], when computed.
    pub two_d_field_wall_shear_mean_pa: Option<f64>,
    /// Whether the cfd-2d solve converged for this channel, when computed.
    pub two_d_converged: Option<bool>,
}

/// Per-node trace derived from the authoritative cfd-1d solve.
#[derive(Debug, Clone)]
pub struct NodeCrossFidelityTrace {
    /// Blueprint node identifier.
    pub node_id: String,
    /// Blueprint node classification.
    pub node_kind: NodeKind,
    /// Reference nodal pressure [Pa].
    pub pressure_pa: f64,
    /// Sum of incoming flow rates [m^3/s].
    pub incoming_flow_m3_s: f64,
    /// Sum of outgoing flow rates [m^3/s].
    pub outgoing_flow_m3_s: f64,
    /// Prescribed boundary source or sink [m^3/s].
    pub prescribed_boundary_flow_m3_s: f64,
    /// Continuity residual [m^3/s].
    pub continuity_residual_m3_s: f64,
}

/// Full blueprint-processing trace for the 3D pipeline entry point.
#[derive(Debug, Clone)]
pub struct Blueprint3dTrace {
    /// Topology class reported by `cfd-mesh`.
    pub topology_class: TopologyClass,
    /// Blueprint channel count.
    pub segment_count: usize,
    /// Synthesized layout segments projected onto the chip plane.
    pub layout_segments: Vec<SegmentCenterline>,
    /// Mesh/schematic volume diagnostics.
    pub volume_trace: PipelineVolumeTrace,
    /// Per-channel cross-fidelity traces.
    pub channel_traces: Vec<ChannelCrossFidelityTrace>,
    /// Per-node continuity and pressure traces.
    pub node_traces: Vec<NodeCrossFidelityTrace>,
    /// Fluid-mesh vertex count.
    pub fluid_mesh_vertex_count: usize,
    /// Fluid-mesh face count.
    pub fluid_mesh_face_count: usize,
    /// Number of converged cfd-2d channels, when that comparison is enabled.
    pub two_d_converged_count: Option<usize>,
    /// Maximum cfd-2d outlet-flow error against the cfd-1d reference [%].
    pub two_d_max_outlet_flow_error_pct: Option<f64>,
    /// Mean cfd-2d outlet-flow error against the cfd-1d reference [%].
    pub two_d_mean_outlet_flow_error_pct: Option<f64>,
}

/// Process a schematic blueprint through the canonical mesh path and attach cross-fidelity traces.
pub fn process_blueprint_with_reference_trace(
    blueprint: &NetworkBlueprint,
    config: &Blueprint3dProcessingConfig,
) -> CfdResult<Blueprint3dTrace> {
    let pipeline_output = BlueprintMeshPipeline::run(blueprint, &config.mesh).map_err(|err| {
        Error::InvalidInput(format!(
            "Blueprint3dTrace failed during cfd-mesh blueprint conversion: {err}"
        ))
    })?;

    let reference_trace = solve_reference_trace::<f64>(
        blueprint,
        config.density_kg_m3,
        config.viscosity_pa_s,
        config.total_flow_rate_m3_s,
    )?;
    let reference_by_channel_id: HashMap<&str, _> = reference_trace
        .channel_traces
        .iter()
        .map(|trace| (trace.channel_id.as_str(), trace))
        .collect();

    let two_d_result = if config.run_2d_reference {
        let builder = Network2dBuilderSink::new(
            BloodModel::Newtonian(config.viscosity_pa_s),
            config.density_kg_m3,
            config.total_flow_rate_m3_s,
            config.two_d_grid_nx,
            config.two_d_grid_ny,
        );
        let mut solver = builder.build(blueprint)?;
        Some(solver.solve_all(config.two_d_tolerance)?)
    } else {
        None
    };
    let two_d_by_channel_id: HashMap<&str, _> = two_d_result
        .as_ref()
        .map(|result| {
            result
                .channels
                .iter()
                .map(|channel| (channel.channel_id.as_str(), channel))
                .collect()
        })
        .unwrap_or_default();

    let channel_traces = pipeline_output
        .volume_trace
        .channel_traces
        .iter()
        .map(|mesh_trace| {
            let reference = reference_by_channel_id
                .get(mesh_trace.channel_id.as_str())
                .copied()
                .ok_or_else(|| {
                    Error::InvalidInput(format!(
                        "Blueprint3dTrace missing cfd-1d reference for channel '{}'",
                        mesh_trace.channel_id
                    ))
                })?;
            let two_d = two_d_by_channel_id
                .get(mesh_trace.channel_id.as_str())
                .copied();

            Ok(ChannelCrossFidelityTrace {
                channel_id: mesh_trace.channel_id.clone(),
                from_node_id: mesh_trace.from_node_id.clone(),
                to_node_id: mesh_trace.to_node_id.clone(),
                schematic_volume_mm3: mesh_trace.schematic_volume_mm3,
                meshed_volume_mm3: mesh_trace.meshed_volume_mm3,
                mesh_volume_error_pct: mesh_trace.volume_error_pct,
                reference_flow_rate_m3_s: reference.flow_rate_m3_s,
                reference_pressure_drop_pa: reference.pressure_drop_pa,
                reference_mean_velocity_m_s: reference.mean_velocity_m_s,
                two_d_outlet_flow_error_pct: two_d
                    .map(|channel| channel.field_outlet_flow_error_pct),
                two_d_field_wall_shear_mean_pa: two_d
                    .map(|channel| channel.field_wall_shear_mean_pa),
                two_d_converged: two_d.map(|channel| channel.solve_result.converged),
            })
        })
        .collect::<CfdResult<Vec<_>>>()?;

    let node_traces = reference_trace
        .node_traces
        .iter()
        .map(|trace| NodeCrossFidelityTrace {
            node_id: trace.node_id.clone(),
            node_kind: trace.node_kind,
            pressure_pa: trace.pressure_pa,
            incoming_flow_m3_s: trace.incoming_flow_m3_s,
            outgoing_flow_m3_s: trace.outgoing_flow_m3_s,
            prescribed_boundary_flow_m3_s: trace.prescribed_boundary_flow_m3_s,
            continuity_residual_m3_s: trace.continuity_residual_m3_s,
        })
        .collect();

    Ok(Blueprint3dTrace {
        topology_class: pipeline_output.topology_class,
        segment_count: pipeline_output.segment_count,
        layout_segments: pipeline_output.layout_segments,
        volume_trace: pipeline_output.volume_trace,
        channel_traces,
        node_traces,
        fluid_mesh_vertex_count: pipeline_output.fluid_mesh.vertex_count(),
        fluid_mesh_face_count: pipeline_output.fluid_mesh.face_count(),
        two_d_converged_count: two_d_result.as_ref().map(|result| result.converged_count),
        two_d_max_outlet_flow_error_pct: two_d_result
            .as_ref()
            .map(|result| result.max_field_outlet_flow_error_pct),
        two_d_mean_outlet_flow_error_pct: two_d_result
            .as_ref()
            .map(|result| result.mean_field_outlet_flow_error_pct),
    })
}
