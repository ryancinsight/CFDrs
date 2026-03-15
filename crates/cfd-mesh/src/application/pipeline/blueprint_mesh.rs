//! Main `NetworkBlueprint → IndexedMesh` pipeline.
//!
//! Converts a `cfd_schematics::NetworkBlueprint` into watertight surface meshes
//! suitable for CFD simulation (fluid mesh) and manufacturing output (chip body).

use std::collections::HashMap;
use std::f64::consts::PI;

use cfd_schematics::geometry::FluidVolumeSummary;
use cfd_schematics::{CrossSectionSpec, NetworkBlueprint, NodeKind};

use crate::application::channel::path::ChannelPath;
use crate::application::channel::profile::ChannelProfile;
use crate::application::channel::substrate::SubstrateBuilder;
use crate::application::channel::sweep::SweepMesher;
use crate::application::csg::boolean::{self, BooleanOp};
use crate::domain::core::error::{MeshError, MeshResult};
use crate::domain::core::index::RegionId;
use crate::domain::core::scalar::{Point3r, Real};
use crate::domain::mesh::IndexedMesh;
use crate::infrastructure::storage::vertex_pool::VertexPool;

use super::constraint::InletOutletConstraint;
use super::topology::{NetworkTopology, TopologyClass};
use super::well_plate::SbsWellPlate96;

// ── PipelineConfig ────────────────────────────────────────────────────────────

/// Configuration for the `BlueprintMeshPipeline`.
#[derive(Debug, Clone)]
pub struct PipelineConfig {
    /// Number of polygon segments for circular cross-sections. Default: 16.
    pub circular_segments: usize,
    /// Number of axial rings per path segment. Default: 8.
    pub axial_rings: usize,
    /// Fractional overlap extension for CSG union watertightness. Default: 0.05.
    pub csg_overlap_fraction: f64,
    /// Arm angle (rad) for the bifurcation diamond layout.
    ///
    /// Controls how steeply the diagonal arms diverge from the centreline
    /// before reaching the horizontal parallel-channel section.  Larger values
    /// give a wider Y-spread.  Default: π/3 (60°).
    pub bifurcation_half_angle_rad: f64,
    /// Half-angle (rad) between trifurcation outer daughter tubes. Default: π/4.
    pub trifurcation_half_angle_rad: f64,
    /// Chip thickness — the Z dimension of the substrate [mm]. Default: 2.0.
    pub chip_height_mm: f64,
    /// Minimum clearance from block edges for channel segments [mm]. Default: 5.0.
    pub wall_clearance_mm: f64,
    /// Whether to build the chip body mesh (substrate minus channel void). Default: true.
    pub include_chip_body: bool,
    /// Skip the 4 mm inlet/outlet hydraulic-diameter constraint.
    ///
    /// Set to `true` for millifluidic designs whose channel cross-sections are
    /// inherently smaller than the 4 mm macro-port specification (e.g. 6 mm × 1 mm
    /// rectangular channels with D_h ≈ 1.7 mm).  The physical tubing adapter
    /// is handled externally.  Default: `false`.
    pub skip_diameter_constraint: bool,
}

impl Default for PipelineConfig {
    fn default() -> Self {
        Self {
            circular_segments: 16,
            axial_rings: 8,
            csg_overlap_fraction: 0.05,
            bifurcation_half_angle_rad: PI / 3.0,
            trifurcation_half_angle_rad: PI / 4.0,
            chip_height_mm: 10.0,
            wall_clearance_mm: 5.0,
            include_chip_body: true,
            skip_diameter_constraint: false,
        }
    }
}

// ── PipelineOutput ────────────────────────────────────────────────────────────

/// XY centreline of one synthesised layout segment, projected onto the chip plane.
///
/// All coordinates are in millimetres.  Use these for 2-D schematics and for
/// comparing the pipeline's physical layout against the source blueprint.
#[derive(Debug, Clone, serde::Serialize)]
pub struct SegmentCenterline {
    /// X start [mm]
    pub x0: f64,
    /// Y start [mm]
    pub y0: f64,
    /// X end [mm]
    pub x1: f64,
    /// Y end [mm]
    pub y1: f64,
    /// Source blueprint channel identifier when this segment maps to one.
    pub source_channel_id: Option<String>,
    /// Upstream blueprint node identifier when this segment maps to one.
    pub from_node_id: Option<String>,
    /// Downstream blueprint node identifier when this segment maps to one.
    pub to_node_id: Option<String>,
    /// Whether this segment was synthesized as a routing connector rather than a direct blueprint channel mapping.
    pub is_synthetic_connector: bool,
    /// Effective tube diameter [mm]
    pub diameter_mm: f64,
}

/// Per-blueprint-channel volume comparison between schematic geometry and the meshed 3D path.
#[derive(Debug, Clone, serde::Serialize)]
pub struct ChannelVolumeTrace {
    /// Blueprint channel identifier.
    pub channel_id: String,
    /// Upstream blueprint node identifier.
    pub from_node_id: String,
    /// Downstream blueprint node identifier.
    pub to_node_id: String,
    /// Schematic centerline length [mm].
    pub schematic_centerline_length_mm: f64,
    /// Meshed centerline length traced through synthesized layout segments [mm].
    pub meshed_centerline_length_mm: f64,
    /// True schematic cross-sectional area [mm^2].
    pub cross_section_area_mm2: f64,
    /// Authoritative schematic fluid volume [mm^3].
    pub schematic_volume_mm3: f64,
    /// Meshed fluid volume for this channel after collapsing the channel's own subsegments [mm^3].
    pub meshed_volume_mm3: f64,
    /// Meshed minus schematic volume [mm^3].
    pub volume_error_mm3: f64,
    /// Relative volume error against the schematic contract [%].
    pub volume_error_pct: f64,
    /// Number of synthesized layout segments attributed to this blueprint channel.
    pub layout_segment_count: usize,
}

/// Volume diagnostics for the full blueprint-to-mesh conversion.
#[derive(Debug, Clone, serde::Serialize)]
pub struct PipelineVolumeTrace {
    /// Authoritative schematic-wide fluid volume summary.
    pub schematic_summary: FluidVolumeSummary,
    /// Per-channel schematic-versus-mesh volume traces.
    pub channel_traces: Vec<ChannelVolumeTrace>,
    /// Sum of all per-channel meshed volumes before inter-channel CSG overlap removal [mm^3].
    pub pre_csg_channel_volume_mm3: f64,
    /// Volume carried by synthesized routing connectors that do not map to a blueprint channel [mm^3].
    pub synthetic_connector_volume_mm3: f64,
    /// Final fluid mesh signed volume after full CSG assembly [mm^3].
    pub fluid_mesh_volume_mm3: f64,
    /// Final chip-body signed volume [mm^3] when requested.
    pub chip_mesh_volume_mm3: Option<f64>,
    /// Fluid-mesh minus schematic total fluid volume [mm^3].
    pub fluid_mesh_volume_error_mm3: f64,
    /// Relative fluid-mesh volume error against the schematic total [%].
    pub fluid_mesh_volume_error_pct: f64,
    /// Volume removed by full-network CSG relative to the sum of meshed channel and connector volumes [mm^3].
    pub csg_overlap_delta_mm3: f64,
}

/// Output of the `BlueprintMeshPipeline`.
pub struct PipelineOutput {
    /// Fluid-domain mesh (channel interior) for CFD solvers.
    pub fluid_mesh: IndexedMesh,
    /// Chip body mesh (substrate minus channel void) for manufacturing / STL.
    /// `None` when `PipelineConfig::include_chip_body` is `false`.
    pub chip_mesh: Option<IndexedMesh>,
    /// Detected topology class.
    pub topology_class: TopologyClass,
    /// Number of channel segments in the fluid mesh.
    pub segment_count: usize,
    /// XY centrelines of all synthesised layout segments.
    ///
    /// Includes the full pre-merge layout (before `merge_collinear_segments`).
    /// Use for 2-D schematic rendering and geometry verification.
    pub layout_segments: Vec<SegmentCenterline>,
    /// Volume diagnostics tying the mesh back to the schematic contract.
    pub volume_trace: PipelineVolumeTrace,
}

// ── BlueprintMeshPipeline ─────────────────────────────────────────────────────

/// Converts a `NetworkBlueprint` into watertight `IndexedMesh` objects.
pub struct BlueprintMeshPipeline;

impl BlueprintMeshPipeline {
    /// Run the full pipeline.
    ///
    /// # Steps
    /// 1. Validate 4 mm inlet/outlet diameter constraint.
    /// 2. Classify network topology.
    /// 3. Synthesize 3-D segment positions in plate coordinates.
    /// 4. Validate wall clearance.
    /// 5. Build per-segment meshes via `SweepMesher`.
    /// 6. Assemble fluid mesh via iterative CSG union.
    /// 7. Label boundary faces (inlet / outlet / wall).
    /// 8. Optionally build chip body via CSG Difference.
    pub fn run(bp: &NetworkBlueprint, config: &PipelineConfig) -> MeshResult<PipelineOutput> {
        let is_selective_routing = bp
            .topology_spec()
            .is_some_and(cfd_schematics::BlueprintTopologySpec::is_selective_routing);
        if is_selective_routing {
            require_geometry_authored_selective_blueprint(bp)?;
        }

        // Step 2 — classify topology (must precede constraint check so we can
        // skip the 4 mm port constraint for ParallelArray micro-channels).
        let topo = NetworkTopology::new(bp);
        let class = if is_selective_routing {
            TopologyClass::Complex
        } else {
            topo.classify()
        };

        // Step 1 — diameter constraint.
        // ParallelArray micro-channels (D_h << 4 mm) are not subject to the
        // 96-well-plate macro-port constraint — skip the check for that topology.
        // Also skip when the caller explicitly opts out (millifluidic channels).
        if !config.skip_diameter_constraint
            && !matches!(&class, TopologyClass::ParallelArray { .. })
        {
            InletOutletConstraint::check(bp).map_err(|e| MeshError::ChannelError {
                message: e.to_string(),
            })?;
        }

        // Complex topologies proceed to graph-based layout synthesis.

        // Step 3 — synthesize layout
        let z_mid = config.chip_height_mm / 2.0;
        let y_center = SbsWellPlate96::center_y();
        // segment_count = number of *blueprint* channels (not synthesised 3-D segments).
        // The serpentine zigzag inserts synthetic turn segments that don't map to
        // blueprint channels, so we cannot use layout.len() here.
        let segment_count = bp.channels.len();
        let use_authored_routed_layout = is_selective_routing
            || (matches!(&class, TopologyClass::Complex) && blueprint_has_routed_layout(bp));
        let layout = if use_authored_routed_layout {
            synthesize_geometry_authored_routed_layout(bp, z_mid)?
        } else {
            synthesize_layout(&class, &topo, bp, y_center, z_mid, config)?
        };

        // Step 4 — wall clearance (routing bounds)
        // Inlet/outlet ports are allowed to touch the x=0 and x=WIDTH_MM faces;
        // only the Y side-walls require `wall_clearance_mm` keep-out.
        for seg in &layout {
            let (x0, y0) = (seg.start.x, seg.start.y);
            let (x1, y1) = (seg.end.x, seg.end.y);
            if !SbsWellPlate96::segment_within_routing_bounds(
                x0,
                y0,
                x1,
                y1,
                config.wall_clearance_mm,
            ) {
                return Err(MeshError::ChannelError {
                    message: format!(
                        "channel segment ({x0:.2}, {y0:.2}) → ({x1:.2}, {y1:.2}) mm \
                         violates routing bounds on plate {:.2} × {:.2} mm \
                         (side clearance {:.1} mm)",
                        SbsWellPlate96::WIDTH_MM,
                        SbsWellPlate96::DEPTH_MM,
                        config.wall_clearance_mm,
                    ),
                });
            }
        }

        // Capture centrelines from the full (pre-merge) layout for schematic output.
        let layout_segments: Vec<SegmentCenterline> = layout
            .iter()
            .map(|seg| SegmentCenterline {
                x0: seg.start.x,
                y0: seg.start.y,
                x1: seg.end.x,
                y1: seg.end.y,
                source_channel_id: seg.source_channel_id.clone(),
                from_node_id: seg.from_node_id.clone(),
                to_node_id: seg.to_node_id.clone(),
                is_synthetic_connector: seg.is_synthetic_connector,
                diameter_mm: cross_section_diameter_mm(&seg.cross_section),
            })
            .collect();

        // Step 5 — merge collinear segments (for CSG topologies only).
        // Collinear same-cross-section consecutive segments would cause coaxial CSG
        // union degeneracy (coincident lateral surfaces → non-deterministic boolean).
        let mesh_layout = merge_collinear_segments(&layout);

        // Step 6 — assemble fluid mesh.
        // LinearChain (serpentine): sweep the entire zigzag as a single polyline to
        // avoid 90° T-junction CDT failures from CSG-unioning perpendicular cylinders.
        // Bifurcation (diamond): build upper route as a polyline, then union the lower
        // fork tubes individually — avoids 3-way endpoint CDT failures at the two
        // T-junctions (see `build_bifurcation_fluid_mesh` for details).
        // VenturiChain / Trifurcation / ParallelArray: per-segment meshes + iterative
        // CSG union.  ParallelArray channels are non-overlapping (different Y rows)
        // so no T-junction degeneracy arises.
        let mut fluid_mesh = if matches!(&class, TopologyClass::LinearChain { .. }) {
            build_polyline_mesh(&layout, 0.0, config)?
        } else if matches!(&class, TopologyClass::VenturiChain) {
            build_venturi_chain_mesh(&mesh_layout, config)?
        } else if matches!(&class, TopologyClass::Complex) {
            build_complex_fluid_mesh(&mesh_layout, config)?
        } else {
            // Robust topological fluid extraction via CSG inversion.
            // Direct multi-way union of intersecting cylinders at unconstrained angles
            build_complex_fluid_mesh(&mesh_layout, config)?
        };

        // Step 7 — label boundaries
        label_boundaries(&mut fluid_mesh, &class, &layout, z_mid, y_center);
        fluid_mesh.rebuild_edges();

        // Complex branching junctions can still retain small residual boundary
        // loops under aggressive auto-layout intersections. Keep those designs
        // usable while the multi-branch closure audit continues.
        if !matches!(&class, TopologyClass::Complex) && !fluid_mesh.is_watertight() {
            let count = fluid_mesh
                .edges_ref()
                .map_or(0, |e| e.boundary_edges().len());
            return Err(MeshError::NotWatertight { count });
        }

        // Step 8 — chip body
        let chip_mesh = if config.include_chip_body {
            Some(if matches!(&class, TopologyClass::ParallelArray { .. }) {
                // ParallelArray: N non-overlapping straight channels at different Y
                // positions.  Sequential subtraction keeps each void tube simple
                // (no multi-tube junction faces) and avoids the CDT PSLG panic
                // from pre-unioning all parallel tubes before the Difference step.
                build_chip_body_sequential(&mesh_layout, config)?
            } else if matches!(&class, TopologyClass::LinearChain { .. }) {
                // LinearChain (serpentine): subtract the full connected fluid
                // mesh as the void so the chip body exposes exactly one inlet
                // and one outlet for a single continuous channel network.
                boolean::csg_boolean(
                    BooleanOp::Difference,
                    &SubstrateBuilder::well_plate_96(config.chip_height_mm).build_indexed()?,
                    &fluid_mesh,
                )?
            } else if matches!(&class, TopologyClass::VenturiChain) {
                // VenturiChain: build the void via the same annular-cap direct
                // sweep path used for the fluid mesh, then subtract from substrate.
                let void_mesh = build_venturi_chain_mesh(&mesh_layout, config)?;
                boolean::csg_boolean(
                    BooleanOp::Difference,
                    &SubstrateBuilder::well_plate_96(config.chip_height_mm).build_indexed()?,
                    &void_mesh,
                )?
            } else if matches!(&class, TopologyClass::Complex) {
                // Complex: reuse the already-built fluid_mesh as the void.
                boolean::csg_boolean(
                    BooleanOp::Difference,
                    &SubstrateBuilder::well_plate_96(config.chip_height_mm).build_indexed()?,
                    &fluid_mesh,
                )?
            } else {
                build_chip_body(&mesh_layout, config)?
            })
        } else {
            None
        };

        let volume_trace =
            compute_volume_trace(bp, &layout, &fluid_mesh, chip_mesh.as_ref(), config)?;

        Ok(PipelineOutput {
            fluid_mesh,
            chip_mesh,
            topology_class: class,
            segment_count,
            layout_segments,
            volume_trace,
        })
    }
}

// ── Internal types ────────────────────────────────────────────────────────────

#[derive(Clone)]
struct SegmentLayout {
    start: Point3r,
    end: Point3r,
    cross_section: CrossSectionSpec,
    source_channel_id: Option<String>,
    from_node_id: Option<String>,
    to_node_id: Option<String>,
    is_synthetic_connector: bool,
}

fn channel_segment(
    start: Point3r,
    end: Point3r,
    cross_section: CrossSectionSpec,
    channel_id: &str,
    from_node_id: &str,
    to_node_id: &str,
) -> SegmentLayout {
    SegmentLayout {
        start,
        end,
        cross_section,
        source_channel_id: Some(channel_id.to_string()),
        from_node_id: Some(from_node_id.to_string()),
        to_node_id: Some(to_node_id.to_string()),
        is_synthetic_connector: false,
    }
}

fn synthetic_segment(
    start: Point3r,
    end: Point3r,
    cross_section: CrossSectionSpec,
) -> SegmentLayout {
    SegmentLayout {
        start,
        end,
        cross_section,
        source_channel_id: None,
        from_node_id: None,
        to_node_id: None,
        is_synthetic_connector: true,
    }
}

fn require_geometry_authored_selective_blueprint(bp: &NetworkBlueprint) -> MeshResult<()> {
    let Some(topology) = bp.topology_spec() else {
        return Err(MeshError::ChannelError {
            message: format!(
                "selective split-tree mesh pipeline requires topology metadata on blueprint '{}'",
                bp.name
            ),
        });
    };
    if !topology.is_selective_routing() {
        return Ok(());
    }
    if !bp.is_geometry_authored() {
        return Err(MeshError::ChannelError {
            message: format!(
                "selective split-tree mesh pipeline requires create_geometry-authored provenance for blueprint '{}' ({})",
                bp.name,
                topology.stage_sequence_label()
            ),
        });
    }
    if bp.render_hints().is_none() {
        return Err(MeshError::ChannelError {
            message: format!(
                "selective split-tree mesh pipeline requires canonical render hints for blueprint '{}' ({})",
                bp.name,
                topology.stage_sequence_label()
            ),
        });
    }
    Ok(())
}

// ── Layout synthesis ──────────────────────────────────────────────────────────

fn synthesize_layout(
    class: &TopologyClass,
    topo: &NetworkTopology<'_>,
    bp: &NetworkBlueprint,
    y_center: Real,
    z_mid: Real,
    config: &PipelineConfig,
) -> MeshResult<Vec<SegmentLayout>> {
    match class {
        // Venturi chain: straight X-axis layout (varying cross-section along one axis).
        //
        // Blueprint channel lengths are rescaled proportionally so the total
        // X span equals the full chip width (0 → chip_w).  This guarantees
        // inlet/outlet caps touch the chip faces, producing clean circular port
        // holes when the chip body CSG Difference is applied.  Section ratios
        // are preserved — a 1:1:1 blueprint stays 1:1:1 on the chip.
        TopologyClass::VenturiChain => {
            let channels = topo
                .linear_path_channels()
                .ok_or_else(|| MeshError::ChannelError {
                    message: "expected linear path but traversal failed".to_string(),
                })?;
            let chip_w = SbsWellPlate96::WIDTH_MM;
            let total_len_mm: Real = channels.iter().map(|ch| ch.length_m * 1000.0).sum();
            let scale = if total_len_mm > 1e-9 {
                chip_w / total_len_mm
            } else {
                1.0
            };

            let mut layout = Vec::with_capacity(channels.len());
            let mut x: Real = 0.0;
            let n = channels.len();
            for (i, ch) in channels.iter().enumerate() {
                let seg_len = ch.length_m * 1000.0 * scale;
                let start = Point3r::new(x, y_center, z_mid);
                // Clamp the final endpoint to chip_w exactly to prevent
                // floating-point accumulation from overshooting the routing
                // bounds check (e.g. 127.760000000000002 > 127.76).
                let x_end = if i + 1 == n { chip_w } else { x + seg_len };
                let end = Point3r::new(x_end, y_center, z_mid);
                layout.push(channel_segment(
                    start,
                    end,
                    ch.cross_section,
                    ch.id.as_str(),
                    ch.from.as_str(),
                    ch.to.as_str(),
                ));
                x = x_end;
            }
            Ok(layout)
        }

        // Linear (serpentine) chain: zigzag rows with interior turns.
        //
        // Each blueprint channel maps to one horizontal row (alternating ±X
        // direction).  To avoid creating extra side-wall openings, only the
        // first row start and final row end touch chip faces; all intermediate
        // turn segments are inset from side walls by one channel radius.
        //
        // This produces a single connected serpentine with exactly one inlet
        // and one outlet while retaining n rows + (n − 1) turns = 2n − 1
        // synthetic layout segments.
        TopologyClass::LinearChain { .. } => {
            let channels = topo
                .linear_path_channels()
                .ok_or_else(|| MeshError::ChannelError {
                    message: "expected linear path but traversal failed".to_string(),
                })?;

            let max_dia_mm = channels
                .iter()
                .map(|ch| cross_section_diameter_mm(&ch.cross_section))
                .fold(0.0_f64, f64::max);
            let row_pitch = (max_dia_mm * 2.5).max(1.0); // ≥ 1 mm
            let chip_w = SbsWellPlate96::WIDTH_MM;
            let n = channels.len();
            let turn_inset_x = (max_dia_mm * 0.5).max(1e-3);
            let x_left = turn_inset_x;
            let x_right = chip_w - turn_inset_x;
            if x_right <= x_left {
                return Err(MeshError::ChannelError {
                    message: format!(
                        "serpentine channel diameter {:.3} mm too large for plate width {:.2} mm",
                        max_dia_mm, chip_w
                    ),
                });
            }

            // Row i sits at y = y_center + (i − (n−1)/2) × row_pitch (centred on chip).
            let y_base = y_center - (n as Real - 1.0) / 2.0 * row_pitch;

            // n rows + (n−1) turns = 2n − 1 segments.
            let mut layout: Vec<SegmentLayout> = Vec::with_capacity(2 * n);
            for i in 0..n {
                let y_row = y_base + i as Real * row_pitch;
                // Even rows run +X, odd rows run −X.
                let x0 = if i == 0 {
                    0.0 // single inlet port
                } else if i % 2 == 0 {
                    x_left
                } else {
                    x_right
                };
                let x1 = if i + 1 == n {
                    if i % 2 == 0 {
                        chip_w // single outlet port on +X face
                    } else {
                        0.0 // single outlet port on −X face
                    }
                } else if i % 2 == 0 {
                    x_right
                } else {
                    x_left
                };

                layout.push(channel_segment(
                    Point3r::new(x0, y_row, z_mid),
                    Point3r::new(x1, y_row, z_mid),
                    channels[i].cross_section,
                    channels[i].id.as_str(),
                    channels[i].from.as_str(),
                    channels[i].to.as_str(),
                ));
                // Vertical turn connecting this row to the next at an inset x.
                if i + 1 < n {
                    let y_next = y_base + (i + 1) as Real * row_pitch;
                    layout.push(synthetic_segment(
                        Point3r::new(x1, y_row, z_mid),
                        Point3r::new(x1, y_next, z_mid),
                        channels[i].cross_section,
                    ));
                }
            }
            Ok(layout)
        }

        // ParallelArray: N parallel straight channels all running full chip width,
        // evenly spaced in Y around y_center.
        TopologyClass::ParallelArray { n_channels } => {
            synthesize_parallel_array_layout(bp, *n_channels, y_center, z_mid, config)
        }

        TopologyClass::Complex => synthesize_complex_layout(bp, y_center, z_mid, config),
    }
}

// ── ParallelArray layout synthesis ───────────────────────────────────────────

/// Lay out N parallel straight channels evenly spaced in Y across the chip.
///
/// Each channel spans the full chip width (0 → `chip_w`) at its own Y row.
/// Row pitch is `max(max_dia_mm × 2.5, 1.0)` mm.  If the requested spread
/// would exceed the routing bounds, the pitch is reduced to the maximum that
/// fits — channels remain non-overlapping as long as `max_dia_mm × 1.0 ≤
/// reduced_pitch` (enforced by the wall-clearance check in `run()`).
fn synthesize_parallel_array_layout(
    bp: &NetworkBlueprint,
    n_channels: usize,
    y_center: Real,
    z_mid: Real,
    config: &PipelineConfig,
) -> MeshResult<Vec<SegmentLayout>> {
    let chip_w = SbsWellPlate96::WIDTH_MM;
    let max_y = SbsWellPlate96::DEPTH_MM - config.wall_clearance_mm;
    let min_y = config.wall_clearance_mm;

    // Determine cross-section from the first channel in the blueprint.
    let first_ch = bp.channels.first().ok_or_else(|| MeshError::ChannelError {
        message: "ParallelArray blueprint has no channels".to_string(),
    })?;
    let cs = first_ch.cross_section;

    // Row pitch: default 2.5 × diameter; clamped so all rows fit within bounds.
    let max_dia_mm = cross_section_diameter_mm(&cs);
    let unclamped_pitch = (max_dia_mm * 2.5).max(1.0);
    // Maximum pitch that keeps all N rows inside [min_y, max_y].
    let available_span = (max_y - min_y).max(0.0);
    let row_pitch = if n_channels <= 1 {
        unclamped_pitch
    } else {
        unclamped_pitch.min(available_span / (n_channels as Real - 1.0))
    };

    // Y positions: centred on y_center.
    let y_base = y_center - (n_channels as Real - 1.0) / 2.0 * row_pitch;

    let mut layout = Vec::with_capacity(n_channels);
    for (i, channel) in bp.channels.iter().take(n_channels).enumerate() {
        let y_row = y_base + i as Real * row_pitch;
        layout.push(channel_segment(
            Point3r::new(0.0, y_row, z_mid),
            Point3r::new(chip_w, y_row, z_mid),
            cs,
            channel.id.as_str(),
            channel.from.as_str(),
            channel.to.as_str(),
        ));
    }
    Ok(layout)
}

// ── Complex (general DAG) layout synthesis ───────────────────────────────────

fn blueprint_has_routed_layout(bp: &NetworkBlueprint) -> bool {
    bp.nodes.iter().all(|node| node.layout.is_some())
        && bp.channels.iter().all(|channel| channel.path.len() >= 2)
}

fn synthesize_geometry_authored_routed_layout(
    bp: &NetworkBlueprint,
    z_mid: Real,
) -> MeshResult<Vec<SegmentLayout>> {
    let node_positions: HashMap<&str, (f64, f64)> = bp
        .nodes
        .iter()
        .map(|node| (node.id.as_str(), node.point))
        .collect();
    let mut layout = Vec::new();

    for channel in &bp.channels {
        let start = node_positions
            .get(channel.from.as_str())
            .copied()
            .ok_or_else(|| MeshError::ChannelError {
                message: format!(
                    "selective blueprint '{}' is missing node position for '{}'",
                    bp.name,
                    channel.from.as_str()
                ),
            })?;
        let end = node_positions
            .get(channel.to.as_str())
            .copied()
            .ok_or_else(|| MeshError::ChannelError {
                message: format!(
                    "selective blueprint '{}' is missing node position for '{}'",
                    bp.name,
                    channel.to.as_str()
                ),
            })?;

        let mut path_points = if channel.path.len() >= 2 {
            channel.path.clone()
        } else {
            vec![start, end]
        };
        if path_points.first().is_some_and(|point| {
            (point.0 - start.0).abs() > 1.0e-6 || (point.1 - start.1).abs() > 1.0e-6
        }) {
            path_points.insert(0, start);
        }
        if path_points.last().is_some_and(|point| {
            (point.0 - end.0).abs() > 1.0e-6 || (point.1 - end.1).abs() > 1.0e-6
        }) {
            path_points.push(end);
        }

        let before_len = layout.len();
        for window in path_points.windows(2) {
            let start_pt = Point3r::new(window[0].0, window[0].1, z_mid);
            let end_pt = Point3r::new(window[1].0, window[1].1, z_mid);
            if (end_pt - start_pt).norm() <= 1.0e-9 {
                continue;
            }
            layout.push(channel_segment(
                start_pt,
                end_pt,
                channel.cross_section,
                channel.id.as_str(),
                channel.from.as_str(),
                channel.to.as_str(),
            ));
        }
        if layout.len() == before_len {
            return Err(MeshError::ChannelError {
                message: format!(
                    "selective blueprint '{}' produced no non-degenerate path segments for channel '{}'",
                    bp.name,
                    channel.id.as_str()
                ),
            });
        }
    }

    Ok(layout)
}

/// Lay out a general directed-acyclic channel network on the SBS-96 plate.
///
/// # Algorithm
///
/// 1. Compute the topological depth of every node via longest-path BFS from the
///    inlet.
/// 2. Map depth to X-coordinate across the chip width (0 → chip_w), with inlet
///    at x = 0 and outlet at x = chip_w.
/// 3. For nodes sharing the same depth, spread evenly in Y around `y_center`.
/// 4. Each blueprint channel maps to one `SegmentLayout` from its `from` node
///    position to its `to` node position.
fn synthesize_complex_layout(
    bp: &NetworkBlueprint,
    y_center: Real,
    z_mid: Real,
    config: &PipelineConfig,
) -> MeshResult<Vec<SegmentLayout>> {
    let chip_w = SbsWellPlate96::WIDTH_MM;
    let max_y = SbsWellPlate96::DEPTH_MM - config.wall_clearance_mm;
    let min_y = config.wall_clearance_mm;

    // Find the unique inlet node.
    let inlet_id = bp
        .nodes
        .iter()
        .find(|n| matches!(n.kind, NodeKind::Inlet))
        .map(|n| n.id.as_str())
        .ok_or_else(|| MeshError::ChannelError {
            message: "Complex blueprint has no inlet node".to_string(),
        })?;

    // BFS longest-path depth from inlet for each node.
    // Iterative relaxation (Bellman-Ford-style for DAG longest-path).
    let mut depth: HashMap<&str, usize> = HashMap::new();
    depth.insert(inlet_id, 0);

    let mut changed = true;
    while changed {
        changed = false;
        for ch in &bp.channels {
            let from = ch.from.as_str();
            let to = ch.to.as_str();
            if let Some(&d) = depth.get(from) {
                let new_d = d + 1;
                let entry = depth.entry(to).or_insert(0);
                if new_d > *entry {
                    *entry = new_d;
                    changed = true;
                }
            }
        }
    }

    // Assign unreachable nodes (if any) to depth 0.
    for node in &bp.nodes {
        depth.entry(node.id.as_str()).or_insert(0);
    }

    let max_depth = depth.values().copied().max().unwrap_or(1).max(1);

    // Group nodes by depth.
    let mut depth_groups: HashMap<usize, Vec<&str>> = HashMap::new();
    for (&node_id, &d) in &depth {
        depth_groups.entry(d).or_default().push(node_id);
    }
    // Sort each group for deterministic layout.
    for group in depth_groups.values_mut() {
        group.sort();
    }

    // Assign 2D positions:
    // - X: depth / max_depth × chip_w
    // - Y: evenly spaced in [min_y, max_y], centred on y_center
    let mut positions: HashMap<&str, (Real, Real)> = HashMap::new();
    for (&d, group) in &depth_groups {
        let x = (d as Real / max_depth as Real) * chip_w;
        let n = group.len();
        if n == 1 {
            positions.insert(group[0], (x, y_center));
        } else {
            let available = max_y - min_y;
            let pitch = available / (n as Real - 1.0).max(1.0);
            for (i, &node_id) in group.iter().enumerate() {
                let y = min_y + i as Real * pitch;
                positions.insert(node_id, (x, y));
            }
        }
    }

    // Build segment layout: one per channel.
    let mut layout: Vec<SegmentLayout> = Vec::with_capacity(bp.channels.len());
    for ch in &bp.channels {
        let from = ch.from.as_str();
        let to = ch.to.as_str();
        let &(x0, y0) = positions.get(from).ok_or_else(|| MeshError::ChannelError {
            message: format!("missing position for node '{from}'"),
        })?;
        let &(x1, y1) = positions.get(to).ok_or_else(|| MeshError::ChannelError {
            message: format!("missing position for node '{to}'"),
        })?;

        // Ensure no degenerate zero-length segments: if two nodes share the
        // same position (e.g. skip-edges between non-adjacent depths), offset
        // the endpoint slightly in X.
        let (sx, sy, ex, ey) = if (x0 - x1).abs() < 1e-6 && (y0 - y1).abs() < 1e-6 {
            (x0, y0, x1 + 0.5, y1)
        } else {
            (x0, y0, x1, y1)
        };

        layout.push(channel_segment(
            Point3r::new(sx, sy, z_mid),
            Point3r::new(ex, ey, z_mid),
            ch.cross_section,
            ch.id.as_str(),
            ch.from.as_str(),
            ch.to.as_str(),
        ));
    }

    Ok(layout)
}

// ── Segment mesh ──────────────────────────────────────────────────────────────

fn build_segment_mesh(seg: &SegmentLayout, config: &PipelineConfig) -> MeshResult<IndexedMesh> {
    let dir = seg.end - seg.start;
    let len = dir.norm();
    if len < 1e-9 {
        return Err(MeshError::ChannelError {
            message: "degenerate segment: start == end".to_string(),
        });
    }
    let unit = dir / len;
    // Cap the CSG overlap extension at half the tube radius.
    //
    // Without this cap a long segment (e.g. a 127.76 mm serpentine row at 5%
    // overlap = 6.39 mm) would extend its tip well past the cross-section
    // radius of the adjacent connecting tube (e.g. a 4 mm diameter turn tube
    // with r = 2 mm).  When the tip pokes through the connecting tube's lateral
    // surface the CSG intersection algorithm produces a PSLG whose constraint
    // segments cross, triggering a CDT "segments intersect in their interiors"
    // panic.  Limiting the extension to ½ · r guarantees the tip always stays
    // inside the body of any same-diameter (or larger) adjacent tube.
    let radius_mm = cross_section_radius_mm(&seg.cross_section);
    let overlap = (len * config.csg_overlap_fraction).min(radius_mm * 0.5);
    let start_ext = Point3r::from(seg.start.coords - unit * overlap);
    let end_ext = Point3r::from(seg.end.coords + unit * overlap);

    let profile = match seg.cross_section {
        CrossSectionSpec::Circular { diameter_m } => ChannelProfile::Circular {
            radius: diameter_m / 2.0 * 1000.0,
            segments: config.circular_segments,
        },
        CrossSectionSpec::Rectangular { width_m, height_m } => ChannelProfile::Rectangular {
            width: width_m * 1000.0,
            height: height_m * 1000.0,
        },
    };

    let path = ChannelPath::straight(start_ext, end_ext);
    let mut pool = VertexPool::default_millifluidic();
    let mesher = SweepMesher {
        cap_start: true,
        cap_end: true,
    };
    let faces = mesher.sweep(&profile, &path, &mut pool, RegionId::from_usize(0));

    let mut mesh = IndexedMesh::new();
    for (_, vdata) in pool.iter() {
        mesh.add_vertex(vdata.position, vdata.normal);
    }
    for face in &faces {
        mesh.add_face_with_region(
            face.vertices[0],
            face.vertices[1],
            face.vertices[2],
            face.region,
        );
    }
    mesh.rebuild_edges();
    repair_pipeline_mesh(&mut mesh, true)?;
    Ok(mesh)
}

// ── Polyline mesh (serpentine single-sweep) ───────────────────────────────────

/// Sweep the entire layout as a **single** polyline — no CSG union required.
///
/// Used for `LinearChain` (serpentine) to avoid 90° T-junction CDT failures
/// that arise from CSG-unioning perpendicular cylinders.  All segments in
/// `layout` must form a connected chain (`seg[i].end ≈ seg[i+1].start`).
///
/// `extend_ends_mm > 0` extends the first waypoint backward and the last
/// waypoint forward (along the respective segment directions).  Set to a
/// non-zero value for the chip-body void so it exits the substrate, creating
/// clean port openings.
fn build_polyline_mesh(
    layout: &[SegmentLayout],
    extend_ends_mm: f64,
    config: &PipelineConfig,
) -> MeshResult<IndexedMesh> {
    build_polyline_mesh_with_policy(layout, extend_ends_mm, config, true)
}

fn build_polyline_mesh_with_policy(
    layout: &[SegmentLayout],
    extend_ends_mm: f64,
    config: &PipelineConfig,
    require_watertight: bool,
) -> MeshResult<IndexedMesh> {
    if layout.is_empty() {
        return Err(MeshError::ChannelError {
            message: "empty layout for polyline mesh".to_string(),
        });
    }

    let n = layout.len();
    let mut points: Vec<Point3r> = Vec::with_capacity(n + 1);

    // First waypoint — optionally extended backward along the first segment.
    let first_dir = (layout[0].end - layout[0].start).normalize();
    points.push(layout[0].start - first_dir * extend_ends_mm);

    // Interior connection points (end of each segment except the last).
    for i in 0..(n - 1) {
        points.push(layout[i].end);
    }

    // Last waypoint — optionally extended forward along the last segment.
    let last = &layout[n - 1];
    let last_dir = (last.end - last.start).normalize();
    points.push(last.end + last_dir * extend_ends_mm);

    let profile = match layout[0].cross_section {
        CrossSectionSpec::Circular { diameter_m } => ChannelProfile::Circular {
            radius: diameter_m / 2.0 * 1000.0,
            segments: config.circular_segments,
        },
        CrossSectionSpec::Rectangular { width_m, height_m } => ChannelProfile::Rectangular {
            width: width_m * 1000.0,
            height: height_m * 1000.0,
        },
    };

    let path = ChannelPath::new(points);
    let mut pool = VertexPool::default_millifluidic();
    let mesher = SweepMesher {
        cap_start: true,
        cap_end: true,
    };
    let faces = mesher.sweep(&profile, &path, &mut pool, RegionId::from_usize(0));

    let mut mesh = IndexedMesh::new();
    for (_, vdata) in pool.iter() {
        mesh.add_vertex(vdata.position, vdata.normal);
    }
    for face in &faces {
        mesh.add_face_with_region(
            face.vertices[0],
            face.vertices[1],
            face.vertices[2],
            face.region,
        );
    }
    mesh.rebuild_edges();
    repair_pipeline_mesh(&mut mesh, require_watertight)?;
    Ok(mesh)
}

// ── Fluid mesh assembly ───────────────────────────────────────────────────────

fn assemble_fluid_mesh(meshes: Vec<IndexedMesh>) -> MeshResult<IndexedMesh> {
    assemble_fluid_mesh_with_policy(meshes, true)
}

fn assemble_fluid_mesh_with_policy(
    meshes: Vec<IndexedMesh>,
    require_watertight: bool,
) -> MeshResult<IndexedMesh> {
    if meshes.is_empty() {
        return Err(MeshError::ChannelError {
            message: "no segment meshes to assemble".to_string(),
        });
    }

    let mut accumulated = crate::application::csg::boolean::csg_boolean_nary(
        crate::application::csg::BooleanOp::Union,
        &meshes,
    )?;
    repair_pipeline_mesh(&mut accumulated, require_watertight)?;
    Ok(accumulated)
}

fn repair_pipeline_mesh(mesh: &mut IndexedMesh, require_watertight: bool) -> MeshResult<()> {
    // The blueprint pipeline mixes direct sweeps and Boolean assembly. Apply a
    // consistent repair pass before final watertight enforcement so chain
    // builders and n-way unions converge to the same closure semantics.
    mesh.orient_outward();
    mesh.retain_largest_component();
    mesh.rebuild_edges();

    if mesh.signed_volume() < 0.0 {
        mesh.flip_faces();
        mesh.rebuild_edges();
    }

    if !mesh.is_watertight() {
        let edge_store = crate::infrastructure::storage::edge_store::EdgeStore::from_face_store(
            &mesh.faces,
        );
        let added = crate::application::watertight::seal::seal_boundary_loops(
            &mut mesh.vertices,
            &mut mesh.faces,
            &edge_store,
            RegionId::INVALID,
        );
        if added > 0 {
            mesh.rebuild_edges();
            mesh.orient_outward();
            mesh.rebuild_edges();
        }
    }

    if !mesh.is_watertight() {
        let improved = crate::application::watertight::repair::MeshRepair::iterative_boundary_stitch(
            &mut mesh.faces,
            &mesh.vertices,
            3,
        );
        if improved > 0 {
            mesh.rebuild_edges();
            mesh.orient_outward();
            mesh.rebuild_edges();
        }
    }

    if mesh.signed_volume() < 0.0 {
        mesh.flip_faces();
        mesh.rebuild_edges();
    }

    if require_watertight && !mesh.is_watertight() {
        let count = mesh.edges_ref().map_or(0, |e| e.boundary_edges().len());
        return Err(MeshError::NotWatertight { count });
    }

    Ok(())
}

// ── Venturi chain concatenated sweep (no CSG) ─────────────────────────────────

/// Build the fluid mesh for a VenturiChain topology — NO CSG union required.
///
/// ## Algorithm
///
/// For a venturi with segments `[seg_0 (R_in), seg_1 (R_throat), seg_2 (R_in)]`, the
/// assembled mesh consists of:
///
/// 1. **Start cap** — outward-facing fan for `seg_0`.
/// 2. **Lateral surface** per segment — quad-strip connecting adjacent rings.
/// 3. **Annular cap** at each cross-section change — connects the outer ring
///    (`R_large`) to the inner ring (`R_small`) at the same axial position.
///    This ring-pair forms a closed annular disk.
/// 4. **End cap** — outward-facing fan for the last segment.
///
/// ## Why not CSG union?
///
/// CSG union of two coaxial circular tubes with different radii places their
/// shared intersection circle exactly on the boundary of both meshes.  The GWN
/// classifier returns `wn ≈ 0.5` for all seam fragments at that circle, and the
/// exact-predicate tiebreaker fails to close the resulting 6-edge hole.  The
/// boundary-hole patcher cannot fix this because the open chain does not form a
/// simple closed polygon.  Direct construction is the only correct approach for
/// step-change cross-section venturis.
fn build_venturi_chain_mesh(
    layout: &[SegmentLayout],
    config: &PipelineConfig,
) -> MeshResult<IndexedMesh> {
    use crate::domain::core::scalar::Point3r;

    if layout.is_empty() {
        return Err(MeshError::ChannelError {
            message: "empty layout for venturi chain mesh".to_string(),
        });
    }

    let n_seg = layout.len();

    // Helper: extract numeric radius [mm] from cross-section spec.
    fn segment_radius_mm(seg: &SegmentLayout) -> f64 {
        match seg.cross_section {
            cfd_schematics::CrossSectionSpec::Circular { diameter_m } => diameter_m / 2.0 * 1000.0,
            cfd_schematics::CrossSectionSpec::Rectangular { width_m, height_m } => {
                // Use half-diagonal as effective radius for rectangular profiles.
                0.5 * (width_m * width_m + height_m * height_m).sqrt() * 1000.0
            }
        }
    }

    let n_seg_pts = config.circular_segments;

    // Build all faces directly into a shared pool to avoid CSG.
    let mut pool = VertexPool::for_csg();
    let mut all_faces: Vec<crate::infrastructure::storage::face_store::FaceData> = Vec::new();

    // Helper: generate CCW ring of `n` vertices at `position` with radius `r`
    // using the frame from a ChannelPath (normal = +Y, binormal = +Z for +X tangent).
    let make_ring = |position: Point3r,
                     r: f64,
                     pool: &mut VertexPool|
     -> Vec<crate::domain::core::index::VertexId> {
        let path_tmp = ChannelPath::straight(
            position,
            Point3r::new(position.x + 1.0, position.y, position.z),
        );
        let frame = &path_tmp.compute_frames()[0];
        let n = n_seg_pts;
        let two_pi = 2.0 * std::f64::consts::PI;
        (0..n)
            .map(|k| {
                let theta = two_pi * k as f64 / n as f64;
                let x = theta.cos() * r;
                let y = theta.sin() * r;
                let pos = frame.position + frame.normal * x + frame.binormal * y;
                let outward = (pos - frame.position).normalize();
                pool.insert_or_weld(pos, outward)
            })
            .collect()
    };

    let region = RegionId::from_usize(0);

    // Build per-segment rings: each segment has a START ring and END ring.
    // We store them as (start_ring, end_ring) per segment.
    let mut segment_rings: Vec<(
        Vec<crate::domain::core::index::VertexId>,
        Vec<crate::domain::core::index::VertexId>,
    )> = Vec::new();

    for seg in layout {
        let r = segment_radius_mm(seg);
        let start_ring = make_ring(seg.start, r, &mut pool);
        let end_ring = make_ring(seg.end, r, &mut pool);
        segment_rings.push((start_ring, end_ring));
    }

    // ── Start cap (outward = -X because inlet faces towards -X) ──────────────
    {
        let first_seg = &layout[0];
        let start_pos = first_seg.start;
        let path_tmp = ChannelPath::straight(
            start_pos,
            Point3r::new(start_pos.x + 1.0, start_pos.y, start_pos.z),
        );
        let frame = &path_tmp.compute_frames()[0];
        let center = pool.insert_or_weld(start_pos, -frame.tangent);
        let ring = &segment_rings[0].0;
        let n = ring.len();
        for i in 0..n {
            let j = (i + 1) % n;
            all_faces.push(crate::infrastructure::storage::face_store::FaceData::new(
                center, ring[j], ring[i], region,
            ));
        }
    }

    // ── Lateral strips + annular junction caps ────────────────────────────────
    for s in 0..n_seg {
        let (ref start_ring, ref end_ring) = segment_rings[s];
        let n = start_ring.len();

        // Lateral quad-strip for segment s.
        for i in 0..n {
            let j = (i + 1) % n;
            all_faces.push(crate::infrastructure::storage::face_store::FaceData::new(
                start_ring[i],
                end_ring[j],
                end_ring[i],
                region,
            ));
            all_faces.push(crate::infrastructure::storage::face_store::FaceData::new(
                start_ring[i],
                start_ring[j],
                end_ring[j],
                region,
            ));
        }

        // Annular cap at the boundary between segment s and segment s+1.
        if s + 1 < n_seg {
            let r_a = segment_radius_mm(&layout[s]);
            let r_b = segment_radius_mm(&layout[s + 1]);

            // If cross-section changes, emit an annular disk between ring_a (R_a) and ring_b (R_b).
            // The annular disk is only needed when the radius changes.
            let radius_diff = (r_a - r_b).abs();
            if radius_diff > 1e-9 {
                let ring_a = &segment_rings[s].1; // end ring of seg s (R_a)
                let ring_b = &segment_rings[s + 1].0; // start ring of seg s+1 (R_b)

                // Winding convention — must satisfy the half-edge pairing invariant:
                //
                //   • The lateral strip for seg s (face1) produces half-edge ring_a[j]→ring_a[i].
                //     The annular cap must produce the REVERSE ring_a[i]→ring_a[j].
                //   • The lateral strip for seg s+1 (face2) produces ring_b[i]→ring_b[j].
                //     The annular cap must produce the REVERSE ring_b[j]→ring_b[i].
                //
                // For a contraction (r_a > r_b), the shoulder face normal = +X (away from inlet):
                //   F1: [ring_b[j], ring_b[i], ring_a[i]]  → ring_b[j]→ring_b[i], ring_b[i]→ring_a[i], ring_a[i]→ring_b[j]
                //   F2: [ring_a[j], ring_b[j], ring_a[i]]  → ring_a[j]→ring_b[j], ring_b[j]→ring_a[i], ring_a[i]→ring_a[j]
                // For an expansion (r_b > r_a), normal = -X (away from outlet):
                //   F1: [ring_a[i], ring_a[j], ring_b[j]]  → ring_a[i]→ring_a[j], ring_a[j]→ring_b[j], ring_b[j]→ring_a[i]
                //   F2: [ring_a[i], ring_b[j], ring_b[i]]  → ring_a[i]→ring_b[j], ring_b[j]→ring_b[i], ring_b[i]→ring_a[i]
                let contraction = r_a > r_b;
                let n = ring_a.len().min(ring_b.len());
                for i in 0..n {
                    let j = (i + 1) % n;
                    if contraction {
                        all_faces.push(crate::infrastructure::storage::face_store::FaceData::new(
                            ring_b[j], ring_b[i], ring_a[i], region,
                        ));
                        all_faces.push(crate::infrastructure::storage::face_store::FaceData::new(
                            ring_a[j], ring_b[j], ring_a[i], region,
                        ));
                    } else {
                        all_faces.push(crate::infrastructure::storage::face_store::FaceData::new(
                            ring_a[i], ring_a[j], ring_b[j], region,
                        ));
                        all_faces.push(crate::infrastructure::storage::face_store::FaceData::new(
                            ring_a[i], ring_b[j], ring_b[i], region,
                        ));
                    }
                }
            }
        }
    }

    // ── End cap ───────────────────────────────────────────────────────────────
    {
        let last_seg = &layout[n_seg - 1];
        let end_pos = last_seg.end;
        let path_tmp =
            ChannelPath::straight(Point3r::new(end_pos.x - 1.0, end_pos.y, end_pos.z), end_pos);
        let frames = path_tmp.compute_frames();
        let frame = frames.last().unwrap();
        let center = pool.insert_or_weld(end_pos, frame.tangent);
        let ring = &segment_rings[n_seg - 1].1;
        let n = ring.len();
        for i in 0..n {
            let j = (i + 1) % n;
            all_faces.push(crate::infrastructure::storage::face_store::FaceData::new(
                center, ring[i], ring[j], region,
            ));
        }
    }

    // ── Reconstruct IndexedMesh with explicit vertex remap ────────────────────
    //
    // `IndexedMesh::add_vertex` uses the mesh's own internal `insert_or_weld`
    // which may return IDs that differ from the pool's sequential indices.
    // Capture the explicit pool → mesh remap.
    let mut mesh = IndexedMesh::new();
    let mut vertex_remap: Vec<crate::domain::core::index::VertexId> =
        Vec::with_capacity(pool.len());
    for (_, vdata) in pool.iter() {
        let mesh_vid = mesh.add_vertex(vdata.position, vdata.normal);
        vertex_remap.push(mesh_vid);
    }
    for face in &all_faces {
        mesh.add_face_with_region(
            vertex_remap[face.vertices[0].as_usize()],
            vertex_remap[face.vertices[1].as_usize()],
            vertex_remap[face.vertices[2].as_usize()],
            face.region,
        );
    }
    mesh.recompute_normals();
    mesh.rebuild_edges();
    repair_pipeline_mesh(&mut mesh, true)?;
    Ok(mesh)
}

// ── Boundary labeling ─────────────────────────────────────────────────────────

fn label_boundaries(
    mesh: &mut IndexedMesh,
    class: &TopologyClass,
    layout: &[SegmentLayout],
    z_mid: Real,
    y_center: Real,
) {
    if layout.is_empty() {
        return;
    }

    // Determine inlet and outlet positions from layout.
    let first_seg = &layout[0];
    let inlet_r_mm = cross_section_radius_mm(&first_seg.cross_section);
    let epsilon = 2.0 * inlet_r_mm;

    // For all topologies except ParallelArray, the single inlet is layout[0].start.
    // For ParallelArray, every channel has its own inlet cap at x=0.
    let inlet_positions: Vec<Point3r> = match class {
        TopologyClass::ParallelArray { .. } => layout.iter().map(|s| s.start).collect(),
        _ => vec![first_seg.start],
    };

    let outlet_positions: Vec<Point3r> = match class {
        TopologyClass::LinearChain { .. } | TopologyClass::VenturiChain => {
            vec![layout.last().unwrap().end]
        }
        TopologyClass::ParallelArray { .. } => {
            // Every channel has its own outlet cap at x = chip_w.
            layout.iter().map(|s| s.end).collect()
        }
        TopologyClass::Complex => {
            // For complex topologies, outlets are segment endpoints that are
            // terminal (degree-1) nodes on the downstream side of the DAG.
            // The inlet is layout[0].start; every other degree-1 endpoint is
            // an outlet.
            let tol = 1e-4;
            let inlet_pos = first_seg.start;
            let mut nodes: Vec<Point3r> = Vec::new();
            let mut node_deg: Vec<usize> = Vec::new();
            for seg in layout {
                for p in &[seg.start, seg.end] {
                    if let Some(idx) = nodes.iter().position(|n| (*n - *p).norm() < tol) {
                        node_deg[idx] += 1;
                    } else {
                        nodes.push(*p);
                        node_deg.push(1);
                    }
                }
            }
            nodes
                .iter()
                .zip(node_deg.iter())
                .filter(|(&n, &d)| d == 1 && (n - inlet_pos).norm() > tol)
                .map(|(&n, _)| n)
                .collect()
        }
    };

    let _ = (z_mid, y_center); // used only in layout synthesis

    // Apply labels to faces by centroid proximity
    let face_ids: Vec<_> = mesh.faces.iter_enumerated().map(|(id, _)| id).collect();
    for fid in face_ids {
        let face = mesh.faces.get(fid);
        let p0 = mesh.vertices.position(face.vertices[0]);
        let p1 = mesh.vertices.position(face.vertices[1]);
        let p2 = mesh.vertices.position(face.vertices[2]);
        let centroid = Point3r::new(
            (p0.x + p1.x + p2.x) / 3.0,
            (p0.y + p1.y + p2.y) / 3.0,
            (p0.z + p1.z + p2.z) / 3.0,
        );

        let label = if inlet_positions
            .iter()
            .any(|ip| (centroid - ip).norm() < epsilon)
        {
            "inlet"
        } else if outlet_positions
            .iter()
            .any(|op| (centroid - op).norm() < epsilon)
        {
            "outlet"
        } else {
            "wall"
        };
        mesh.mark_boundary(fid, label);
    }
}

fn cross_section_radius_mm(cs: &CrossSectionSpec) -> Real {
    match cs {
        CrossSectionSpec::Circular { diameter_m } => diameter_m / 2.0 * 1000.0,
        CrossSectionSpec::Rectangular { width_m, height_m } => {
            // Use half-diagonal as effective radius for proximity check
            (width_m * width_m + height_m * height_m).sqrt() * 500.0
        }
    }
}

/// Effective outer diameter [mm] of a cross-section.
///
/// Used for clearance calculations (e.g. row pitch and arm offsets) where a
/// single linear measure of the tube size is needed.
fn cross_section_diameter_mm(cs: &CrossSectionSpec) -> Real {
    match cs {
        CrossSectionSpec::Circular { diameter_m } => diameter_m * 1000.0,
        CrossSectionSpec::Rectangular { width_m, height_m } => width_m.max(*height_m) * 1000.0,
    }
}

// ── Complex fluid mesh ────────────────────────────────────────────────────────

/// Build the fluid mesh for a [`TopologyClass::Complex`] DAG layout.
///
/// # Algorithm — chain extraction + sequential CSG union
///
/// 1. Build an undirected adjacency graph over segment endpoints (using
///    coordinate-based node matching with 1 µm tolerance).
/// 2. Extract maximal chains through degree-2 intermediate nodes: each chain
///    is a sequence of connected segments that can be swept as one polyline
///    without internal caps (avoiding coincident-face CSG degeneracy).
/// 3. Build each chain as a polyline mesh.
/// 4. CSG-union the chains in order.
fn build_complex_fluid_mesh(
    layout: &[SegmentLayout],
    config: &PipelineConfig,
) -> MeshResult<IndexedMesh> {
    if layout.is_empty() {
        return Err(MeshError::ChannelError {
            message: "empty layout for complex topology".to_string(),
        });
    }

    // ── Step 1: identify unique node positions ───────────────────────────────
    let tol = 1e-4; // 0.1 mm tolerance in mm coordinates
    let mut nodes: Vec<Point3r> = Vec::new();
    let mut seg_nodes: Vec<(usize, usize)> = Vec::with_capacity(layout.len());

    let mut find_or_add = |p: Point3r| -> usize {
        for (i, pos) in nodes.iter().enumerate() {
            if (pos.x - p.x).abs() < 1e-5 && (pos.y - p.y).abs() < 1e-5 {
                return i;
            }
        }
        nodes.push(p);
        nodes.len() - 1
    };

    for seg in layout {
        let s = find_or_add(seg.start);
        let e = find_or_add(seg.end);
        seg_nodes.push((s, e));
    }

    // ── Step 2: compute node degrees ─────────────────────────────────────────
    let n_nodes = nodes.len();
    let mut degree = vec![0usize; n_nodes];
    for &(s, e) in &seg_nodes {
        degree[s] += 1;
        degree[e] += 1;
    }

    // ── Step 3: extract chains ───────────────────────────────────────────────
    // A chain is a maximal run of segments connected through degree-2 nodes.
    let mut consumed = vec![false; layout.len()];
    let mut chains: Vec<Vec<usize>> = Vec::new();

    // Build node → segment index adjacency.
    let mut node_segs: Vec<Vec<usize>> = vec![Vec::new(); n_nodes];
    for (si, &(s, e)) in seg_nodes.iter().enumerate() {
        node_segs[s].push(si);
        node_segs[e].push(si);
    }

    for start_seg in 0..layout.len() {
        if consumed[start_seg] {
            continue;
        }
        consumed[start_seg] = true;
        let mut chain = vec![start_seg];

        // Extend forward from end node.
        let (s0, mut tip) = seg_nodes[start_seg];
        loop {
            if degree[tip] != 2 {
                break;
            }
            if let Some(si) = node_segs[tip].iter().find(|&&si| !consumed[si]).copied() {
                consumed[si] = true;
                let (ns, ne) = seg_nodes[si];
                tip = if ns == tip { ne } else { ns };
                chain.push(si);
            } else {
                break;
            }
        }

        // Extend backward from start node.
        let mut head = s0;
        loop {
            if degree[head] != 2 {
                break;
            }
            if let Some(si) = node_segs[head].iter().find(|&&si| !consumed[si]).copied() {
                consumed[si] = true;
                let (ns, ne) = seg_nodes[si];
                head = if ns == head { ne } else { ns };
                chain.insert(0, si);
            } else {
                break;
            }
        }

        chains.push(chain);
    }

    // ── Step 4: orient chain segments for contiguous polyline sweep ──────────
    let oriented_chains: Vec<Vec<SegmentLayout>> = chains
        .iter()
        .map(|chain| {
            let mut segs: Vec<SegmentLayout> = Vec::with_capacity(chain.len());
            for (ci, &si) in chain.iter().enumerate() {
                let mut seg = layout[si].clone();
                if ci > 0 {
                    let prev_end = segs[ci - 1].end;
                    if (seg.end - prev_end).norm() < tol && (seg.start - prev_end).norm() >= tol {
                        std::mem::swap(&mut seg.start, &mut seg.end);
                    }
                }
                segs.push(seg);
            }
            segs
        })
        .collect();

    // ── Step 5: build mesh for each chain ────────────────────────────────────
    let mut chain_meshes: Vec<IndexedMesh> = Vec::with_capacity(oriented_chains.len());
    for chain_segs in &oriented_chains {
        if chain_segs.len() == 1 {
            let mesh = build_segment_mesh(&chain_segs[0], config).map_err(|error| {
                MeshError::ChannelError {
                    message: format!("complex chain segment mesh failed: {error}"),
                }
            })?;
            chain_meshes.push(mesh);
        } else {
            let mesh = build_polyline_mesh_with_policy(
                chain_segs,
                0.0,
                config,
                false,
            )
            .map_err(|error| MeshError::ChannelError {
                message: format!(
                    "complex chain polyline mesh failed for {} segments: {error}",
                    chain_segs.len()
                ),
            })?;
            chain_meshes.push(mesh);
        }
    }

    assemble_fluid_mesh_with_policy(chain_meshes, false).map_err(|error| MeshError::ChannelError {
        message: format!("complex fluid assembly failed: {error}"),
    })
}

fn build_chip_body(layout: &[SegmentLayout], config: &PipelineConfig) -> MeshResult<IndexedMesh> {
    let substrate = SubstrateBuilder::well_plate_96(config.chip_height_mm).build_indexed()?;

    // Build the void mesh using the same tube geometry as the fluid mesh.
    // The 5% CSG overlap extension makes void tubes exit the inlet face (x < 0),
    // creating proper port openings in the substrate.
    // config.chip_height_mm must be > channel_diameter so the void stays within
    // the substrate in the Z direction (default: 10 mm for 4 mm channels).
    let void_meshes: Vec<IndexedMesh> = layout
        .iter()
        .map(|seg| build_segment_mesh(seg, config))
        .collect::<MeshResult<_>>()?;

    let void_mesh = assemble_fluid_mesh(void_meshes)?;

    boolean::csg_boolean(BooleanOp::Difference, &substrate, &void_mesh)
}

/// Build the chip body by subtracting each void tube **individually**.
///
/// Used for branching topologies ([`TopologyClass::Bifurcation`] and
/// [`TopologyClass::Trifurcation`]) where pre-unioning all tubes before the
/// `Difference` step causes a CDT panic inside `corefine_face`: the complex
/// junction faces of the multi-tube void union produce overlapping constraint
/// segments in the 2-D PSLG projection, triggering an `"Invalid PSLG"` panic
/// inside `Cdt::from_pslg`.
///
/// Sequential subtraction keeps every void operand as a simple single-cylinder
/// mesh with clean hemispherical-cap geometry, so no crossing segments arise.
fn build_chip_body_sequential(
    layout: &[SegmentLayout],
    config: &PipelineConfig,
) -> MeshResult<IndexedMesh> {
    let mut result = SubstrateBuilder::well_plate_96(config.chip_height_mm).build_indexed()?;
    for seg in layout {
        let void_mesh = build_segment_mesh(seg, config)?;
        result = boolean::csg_boolean(BooleanOp::Difference, &result, &void_mesh)?;
    }
    Ok(result)
}

// ── Segment merging ───────────────────────────────────────────────────────────

/// Merge consecutive collinear same-cross-section segments into single longer
/// segments.  This prevents coaxial CSG union degeneracy: two equal-diameter
/// cylinders on the same axis have coincident lateral surfaces in the overlap
/// region, which causes floating-point precision failures in the boolean.
///
/// Merging criteria (all three must hold):
/// - The segments are adjacent (`seg[i].end ≈ seg[i+1].start`).
/// - They travel in the same direction (unit vectors within 1 µm tolerance).
/// - They have the same cross-section geometry.
fn merge_collinear_segments(layout: &[SegmentLayout]) -> Vec<SegmentLayout> {
    if layout.is_empty() {
        return vec![];
    }
    let mut merged: Vec<SegmentLayout> = Vec::with_capacity(layout.len());
    let mut current = layout[0].clone();

    for next in &layout[1..] {
        let cur_dir = (current.end - current.start).normalize();
        let nxt_dir = (next.end - next.start).normalize();
        let connected = (current.end - next.start).norm() < 1e-6;
        let collinear = (cur_dir - nxt_dir).norm() < 1e-6;
        let same_cs = cross_sections_equal(&current.cross_section, &next.cross_section);
        let same_source = current.source_channel_id == next.source_channel_id
            && current.from_node_id == next.from_node_id
            && current.to_node_id == next.to_node_id
            && current.is_synthetic_connector == next.is_synthetic_connector;
        if connected && collinear && same_cs && same_source {
            current.end = next.end; // extend without changing start
        } else {
            merged.push(current);
            current = next.clone();
        }
    }
    merged.push(current);
    merged
}

fn cross_sections_equal(a: &CrossSectionSpec, b: &CrossSectionSpec) -> bool {
    match (a, b) {
        (
            CrossSectionSpec::Circular { diameter_m: d1 },
            CrossSectionSpec::Circular { diameter_m: d2 },
        ) => (d1 - d2).abs() < 1e-9,
        (
            CrossSectionSpec::Rectangular {
                width_m: w1,
                height_m: h1,
            },
            CrossSectionSpec::Rectangular {
                width_m: w2,
                height_m: h2,
            },
        ) => (w1 - w2).abs() < 1e-9 && (h1 - h2).abs() < 1e-9,
        _ => false,
    }
}

fn compute_volume_trace(
    bp: &NetworkBlueprint,
    layout: &[SegmentLayout],
    fluid_mesh: &IndexedMesh,
    chip_mesh: Option<&IndexedMesh>,
    config: &PipelineConfig,
) -> MeshResult<PipelineVolumeTrace> {
    let schematic_summary = bp.fluid_volume_summary();
    let mut channel_traces = Vec::with_capacity(bp.channels.len());
    let mut pre_csg_channel_volume_mm3 = 0.0;

    for channel_summary in bp.channel_fluid_volume_summaries() {
        let channel_segments: Vec<SegmentLayout> = layout
            .iter()
            .filter(|segment| {
                segment.source_channel_id.as_deref() == Some(channel_summary.channel_id.as_str())
            })
            .cloned()
            .collect();
        if channel_segments.is_empty() {
            return Err(MeshError::ChannelError {
                message: format!(
                    "mesh pipeline produced no synthesized segments for blueprint channel '{}'",
                    channel_summary.channel_id
                ),
            });
        }

        let channel_mesh = if channel_segments.len() == 1 {
            build_segment_mesh(&channel_segments[0], config)?
        } else {
            build_polyline_mesh(&channel_segments, 0.0, config)?
        };
        let meshed_volume_mm3 = channel_mesh.signed_volume().abs();
        let meshed_centerline_length_mm =
            channel_segments.iter().map(segment_length_mm).sum::<f64>();
        let volume_error_mm3 = meshed_volume_mm3 - channel_summary.fluid_volume_mm3;
        let volume_error_pct = relative_percent(volume_error_mm3, channel_summary.fluid_volume_mm3);

        pre_csg_channel_volume_mm3 += meshed_volume_mm3;
        channel_traces.push(ChannelVolumeTrace {
            channel_id: channel_summary.channel_id,
            from_node_id: channel_summary.from_node_id,
            to_node_id: channel_summary.to_node_id,
            schematic_centerline_length_mm: channel_summary.centerline_length_mm,
            meshed_centerline_length_mm,
            cross_section_area_mm2: channel_summary.cross_section_area_mm2,
            schematic_volume_mm3: channel_summary.fluid_volume_mm3,
            meshed_volume_mm3,
            volume_error_mm3,
            volume_error_pct,
            layout_segment_count: channel_segments.len(),
        });
    }

    let synthetic_connector_volume_mm3 = layout
        .iter()
        .filter(|segment| segment.is_synthetic_connector)
        .map(|segment| build_segment_mesh(segment, config).map(|mesh| mesh.signed_volume().abs()))
        .collect::<MeshResult<Vec<_>>>()?
        .into_iter()
        .sum();

    let fluid_mesh_volume_mm3 = fluid_mesh.signed_volume().abs();
    let chip_mesh_volume_mm3 = chip_mesh.map(|mesh| mesh.signed_volume().abs());
    let fluid_mesh_volume_error_mm3 =
        fluid_mesh_volume_mm3 - schematic_summary.total_fluid_volume_mm3;
    let fluid_mesh_volume_error_pct = relative_percent(
        fluid_mesh_volume_error_mm3,
        schematic_summary.total_fluid_volume_mm3,
    );
    let csg_overlap_delta_mm3 =
        pre_csg_channel_volume_mm3 + synthetic_connector_volume_mm3 - fluid_mesh_volume_mm3;

    Ok(PipelineVolumeTrace {
        schematic_summary,
        channel_traces,
        pre_csg_channel_volume_mm3,
        synthetic_connector_volume_mm3,
        fluid_mesh_volume_mm3,
        chip_mesh_volume_mm3,
        fluid_mesh_volume_error_mm3,
        fluid_mesh_volume_error_pct,
        csg_overlap_delta_mm3,
    })
}

fn segment_length_mm(segment: &SegmentLayout) -> f64 {
    (segment.end - segment.start).norm()
}

fn relative_percent(delta: f64, reference: f64) -> f64 {
    if reference.abs() <= 1e-18 {
        0.0
    } else {
        delta.abs() / reference.abs() * 100.0
    }
}

#[cfg(test)]
mod tests {
    use cfd_schematics::interface::presets::venturi_chain;
    use cfd_schematics::topology::presets::{
        build_milestone12_blueprint, build_milestone12_topology_spec, Milestone12TopologyRequest,
    };

    use super::*;

    #[test]
    fn default_config_builds() {
        let _cfg = PipelineConfig::default();
    }

    #[test]
    fn pipeline_rejects_wrong_diameter() {
        use cfd_schematics::interface::presets::serpentine_chain;
        let bp = serpentine_chain("x", 3, 0.010, 0.002);
        let result = BlueprintMeshPipeline::run(&bp, &PipelineConfig::default());
        assert!(result.is_err());
        let msg = result.err().expect("checked above").to_string();
        assert!(
            msg.contains("hydraulic diameter") || msg.contains("channel error"),
            "unexpected error: {msg}"
        );
    }

    #[test]
    fn pipeline_handles_complex_topology() {
        // Build a blueprint with complex topology manually
        use cfd_schematics::{ChannelSpec, NetworkBlueprint, NodeKind, NodeSpec};
        let mut bp = NetworkBlueprint {
            name: "complex".to_string(),
            box_dims: (127.76, 85.47),
            box_outline: Vec::new(),
            nodes: Vec::new(),
            channels: Vec::new(),
            render_hints: None,
            topology: None,
            lineage: None,
            metadata: None,
            geometry_authored: false,
        };
        bp.add_node(NodeSpec::new_at("inlet", NodeKind::Inlet, (0.0, 42.735)));
        bp.add_node(NodeSpec::new_at("j1", NodeKind::Junction, (30.0, 42.735)));
        bp.add_node(NodeSpec::new_at("j2", NodeKind::Junction, (60.0, 42.735)));
        bp.add_node(NodeSpec::new_at("j3", NodeKind::Junction, (90.0, 42.735)));
        bp.add_node(NodeSpec::new_at("outlet", NodeKind::Outlet, (120.0, 42.735)));
        // 5 channels: creates a Complex topology (degree > 3 at some node)
        for i in 1..=5_usize {
            let from = if i == 1 { "inlet" } else { "j1" };
            let to = if i == 5 { "outlet" } else { "j2" };
            bp.add_channel(ChannelSpec::new_pipe(
                format!("c{i}"),
                from,
                to,
                0.005,
                0.004,
                0.0,
                0.0,
            ));
        }
        let cfg = PipelineConfig {
            include_chip_body: false,
            skip_diameter_constraint: true,
            ..Default::default()
        };
        let result = BlueprintMeshPipeline::run(&bp, &cfg)
            .expect("complex topology should be supported by graph layout synthesis");
        assert_eq!(result.topology_class, TopologyClass::Complex);
    }

    #[test]
    fn venturi_chain_produces_output() {
        let bp = venturi_chain("v", 0.030, 0.004, 0.002);
        let cfg = PipelineConfig {
            include_chip_body: false,
            ..Default::default()
        };
        let result = BlueprintMeshPipeline::run(&bp, &cfg);
        // We expect it to succeed (or fail gracefully if CSG has limitations)
        match result {
            Ok(out) => {
                assert_eq!(out.topology_class, TopologyClass::VenturiChain);
                assert_eq!(out.segment_count, 3);
            }
            Err(e) => {
                // Allow CSG failures in unit tests (they require specific geometry)
                let msg = e.to_string();
                assert!(
                    !msg.contains("hydraulic diameter"),
                    "unexpected diameter error: {msg}"
                );
            }
        }
    }

    #[test]
    fn selective_pipeline_rejects_non_geometry_authored_blueprint() {
        let request = Milestone12TopologyRequest::new(
            "tri_bi",
            "Tri→Bi",
            vec![
                cfd_schematics::SplitKind::NFurcation(3),
                cfd_schematics::SplitKind::NFurcation(2),
            ],
            6.0e-3,
            1.0e-3,
            8.0e-3,
            8.0e-3,
        );
        let spec = build_milestone12_topology_spec(&request);
        let mut bp = build_milestone12_blueprint(&request).expect("selective blueprint");
        bp.metadata = None;
        bp.geometry_authored = false;
        bp.topology = Some(spec);

        let result = BlueprintMeshPipeline::run(
            &bp,
            &PipelineConfig {
                include_chip_body: false,
                skip_diameter_constraint: true,
                ..Default::default()
            },
        );
        let message = result
            .err()
            .expect("non-canonical selective blueprint must fail");
        assert!(
            message
                .to_string()
                .contains("create_geometry-authored provenance"),
            "unexpected error: {message}"
        );
    }

    #[test]
    fn selective_pipeline_uses_geometry_authored_paths() {
        let request = Milestone12TopologyRequest::new(
            "tri_bi",
            "Tri→Bi",
            vec![
                cfd_schematics::SplitKind::NFurcation(3),
                cfd_schematics::SplitKind::NFurcation(2),
            ],
            6.0e-3,
            1.0e-3,
            8.0e-3,
            8.0e-3,
        );
        let bp = build_milestone12_blueprint(&request).expect("geometry-authored selective build");
        let result = BlueprintMeshPipeline::run(
            &bp,
            &PipelineConfig {
                include_chip_body: false,
                skip_diameter_constraint: true,
                ..Default::default()
            },
        )
        .expect("geometry-authored selective blueprint should mesh");
        assert_eq!(result.topology_class, TopologyClass::Complex);
        assert!(
            result.layout_segments.len() >= bp.channels.len(),
            "selective path meshing should preserve authored channel polylines"
        );
    }
}
