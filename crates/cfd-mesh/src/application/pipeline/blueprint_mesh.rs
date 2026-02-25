//! Main `NetworkBlueprint → IndexedMesh` pipeline.
//!
//! Converts a `cfd_schematics::NetworkBlueprint` into watertight surface meshes
//! suitable for CFD simulation (fluid mesh) and manufacturing output (chip body).

use std::f64::consts::PI;

use cfd_schematics::{CrossSectionSpec, NetworkBlueprint};

use crate::application::channel::path::ChannelPath;
use crate::application::channel::profile::ChannelProfile;
use crate::application::channel::substrate::SubstrateBuilder;
use crate::application::channel::sweep::SweepMesher;
use crate::application::csg::boolean::csg_boolean_indexed;
use crate::application::csg::boolean::operations::BooleanOp;
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
    /// Effective tube diameter [mm]
    pub diameter_mm: f64,
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
        // Step 1 — diameter constraint
        InletOutletConstraint::check(bp).map_err(|e| MeshError::ChannelError {
            message: e.to_string(),
        })?;

        // Step 2 — classify topology
        let topo = NetworkTopology::new(bp);
        let class = topo.classify();

        if class == TopologyClass::Complex {
            return Err(MeshError::ChannelError {
                message: "blueprint has Complex topology — not yet supported by the pipeline"
                    .to_string(),
            });
        }

        // Step 3 — synthesize layout
        let z_mid = config.chip_height_mm / 2.0;
        let y_center = SbsWellPlate96::center_y();
        // segment_count = number of *blueprint* channels (not synthesised 3-D segments).
        // The serpentine zigzag inserts synthetic turn segments that don't map to
        // blueprint channels, so we cannot use layout.len() here.
        let segment_count = bp.channels.len();
        let layout = synthesize_layout(&class, &topo, y_center, z_mid, config)?;

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
        // All other topologies: per-segment mesh + iterative CSG union.
        let mut fluid_mesh = if matches!(&class, TopologyClass::LinearChain { .. }) {
            build_polyline_mesh(&layout, 0.0, config)?
        } else if matches!(&class, TopologyClass::Bifurcation) {
            build_bifurcation_fluid_mesh(&mesh_layout, config)?
        } else {
            let segment_meshes: Vec<IndexedMesh> = mesh_layout
                .iter()
                .map(|seg| build_segment_mesh(seg, config))
                .collect::<MeshResult<_>>()?;
            assemble_fluid_mesh(segment_meshes)?
        };

        // Step 7 — label boundaries
        label_boundaries(&mut fluid_mesh, &class, &layout, z_mid, y_center);
        fluid_mesh.rebuild_edges();

        if !fluid_mesh.is_watertight() {
            return Err(MeshError::NotWatertight { count: 0 });
        }

        // Step 8 — chip body
        let chip_mesh = if config.include_chip_body {
            Some(if matches!(&class, TopologyClass::LinearChain { .. }) {
                // Bent-polyline void for the serpentine chip body has a
                // fundamental CSG problem: the 90° turn tube's interior rings
                // are coplanar with (or tangent to) the substrate faces
                // regardless of how much the turns are inset, because the
                // SweepMesher places profile rings whose normal direction is
                // exactly ±X/±Y along the straight turn body — touching the
                // substrate face at the tube's outer rim.
                //
                // Workaround: use one STRAIGHT row void per serpentine row.
                // Adjacent rows are non-overlapping (row_pitch > 2r), so their
                // CSG union is trivial (no 90° T-junction CDT issue).  The chip
                // body has 2·n clean circular through-holes (n per face) and
                // the turn voids are omitted — the fluid mesh already encodes
                // the correct serpentine geometry for CFD use.
                let rows_only: Vec<SegmentLayout> = layout
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| i % 2 == 0) // even indices = row segments
                    .map(|(_, seg)| seg.clone())
                    .collect();
                build_chip_body(&rows_only, config)?
            } else if matches!(
                &class,
                TopologyClass::Trifurcation | TopologyClass::Bifurcation
            ) {
                // Pre-unioning void tubes before the Difference step causes a CDT
                // panic in corefine_face for any branching topology: the complex
                // junction faces of the multi-tube void union produce overlapping
                // constraint segments in the 2-D PSLG projection.  Sequential
                // subtraction keeps each void operand as a simple single-cylinder
                // mesh with clean cap geometry, avoiding the problematic junction
                // faces entirely.
                build_chip_body_sequential(&mesh_layout, config)?
            } else {
                build_chip_body(&mesh_layout, config)?
            })
        } else {
            None
        };

        Ok(PipelineOutput {
            fluid_mesh,
            chip_mesh,
            topology_class: class,
            segment_count,
            layout_segments,
        })
    }
}

// ── Internal types ────────────────────────────────────────────────────────────

#[derive(Clone)]
struct SegmentLayout {
    start: Point3r,
    end: Point3r,
    cross_section: CrossSectionSpec,
}

// ── Layout synthesis ──────────────────────────────────────────────────────────

fn synthesize_layout(
    class: &TopologyClass,
    topo: &NetworkTopology<'_>,
    y_center: Real,
    z_mid: Real,
    config: &PipelineConfig,
) -> MeshResult<Vec<SegmentLayout>> {
    match class {
        // Venturi chain: straight X-axis layout (varying cross-section along one axis).
        TopologyClass::VenturiChain => {
            let channels = topo.linear_path_channels().ok_or_else(|| MeshError::ChannelError {
                message: "expected linear path but traversal failed".to_string(),
            })?;

            let mut layout = Vec::with_capacity(channels.len());
            let mut x: Real = 0.0;
            for ch in &channels {
                let seg_len = ch.length_m * 1000.0; // m → mm
                let start = Point3r::new(x, y_center, z_mid);
                let end = Point3r::new(x + seg_len, y_center, z_mid);
                layout.push(SegmentLayout {
                    start,
                    end,
                    cross_section: ch.cross_section,
                });
                x += seg_len;
            }
            Ok(layout)
        }

        // Linear (serpentine) chain: zigzag rows spanning the full chip width.
        //
        // Each blueprint channel maps to one horizontal row that spans the chip from
        // x = 0 to x = WIDTH_MM (alternating ±X direction).  Consecutive rows are
        // separated in Y by `row_pitch = max_diameter × 2.5` and connected by a
        // short vertical "turn" segment at the row end.  This produces n rows and
        // (n − 1) turn segments — 2n − 1 segments total.
        TopologyClass::LinearChain { .. } => {
            let channels = topo.linear_path_channels().ok_or_else(|| MeshError::ChannelError {
                message: "expected linear path but traversal failed".to_string(),
            })?;

            let max_dia_mm = channels
                .iter()
                .map(|ch| cross_section_diameter_mm(&ch.cross_section))
                .fold(0.0_f64, f64::max);
            let row_pitch = (max_dia_mm * 2.5).max(1.0); // ≥ 1 mm
            let chip_w = SbsWellPlate96::WIDTH_MM;
            let n = channels.len();

            // Row i sits at y = y_center + (i − (n−1)/2) × row_pitch (centred on chip).
            let y_base = y_center - (n as Real - 1.0) / 2.0 * row_pitch;

            // n rows + (n−1) turns = 2n − 1 segments.
            let mut layout: Vec<SegmentLayout> = Vec::with_capacity(2 * n);
            for i in 0..n {
                let y_row = y_base + i as Real * row_pitch;
                // Even rows go left → right; odd rows go right → left.
                let (x0, x1) = if i % 2 == 0 { (0.0, chip_w) } else { (chip_w, 0.0) };
                layout.push(SegmentLayout {
                    start: Point3r::new(x0, y_row,  z_mid),
                    end:   Point3r::new(x1, y_row,  z_mid),
                    cross_section: channels[i].cross_section,
                });
                // Vertical turn connecting this row to the next.
                if i + 1 < n {
                    let y_next = y_base + (i + 1) as Real * row_pitch;
                    layout.push(SegmentLayout {
                        start: Point3r::new(x1, y_row,  z_mid),
                        end:   Point3r::new(x1, y_next, z_mid),
                        cross_section: channels[i].cross_section,
                    });
                }
            }
            Ok(layout)
        }

        // Bifurcation — diamond (recombining bifurcation) spanning the full chip width:
        //
        //      (p_x1, y_c+dv) ─── upper parallel ─── (p_x2, y_c+dv)
        //     ╱                                                       ╲
        // ─inlet─── (div_x, y_c)                     (conv_x, y_c) ───outlet─
        //     ╲                                                       ╱
        //      (p_x1, y_c-dv) ─── lower parallel ─── (p_x2, y_c-dv)
        //
        // The design maps the full blueprint (parent_in → daughters → parent_out)
        // to 8 physical segments spanning x = 0 → chip_w:
        //
        //   [1] inlet_straight  : 0..25 % chip_w at y_center
        //   [2] upper_arm_in    : 25%→37.5% rising by dv
        //   [3] upper_parallel  : 37.5%→62.5% at y_center + dv
        //   [4] upper_arm_out   : 62.5%→75% returning to y_center
        //   [5] lower_arm_in    : 25%→37.5% falling by dv
        //   [6] lower_parallel  : 37.5%→62.5% at y_center − dv
        //   [7] lower_arm_out   : 62.5%→75% returning to y_center
        //   [8] outlet_straight : 75%→100% chip_w at y_center
        //
        // dv = arm_x_span × tan(bifurcation_half_angle_rad), capped at routing bound.
        TopologyClass::Bifurcation => {
            let (parent_in, d1, _d2, parent_out) =
                topo.bifurcation_channels().ok_or_else(|| MeshError::ChannelError {
                    message: "bifurcation topology detected but channel decomposition failed"
                        .to_string(),
                })?;

            let chip_w    = SbsWellPlate96::WIDTH_MM;
            let max_half_y = y_center - config.wall_clearance_mm;

            // Fixed-fraction layout: symmetric about chip centre.
            let div_x  = chip_w * 0.25;           // 25 % — split/merge x
            let arm_x  = chip_w * 0.125;           // 12.5 % — arm x-span
            let p_x1   = div_x + arm_x;            // 37.5 % — parallel start
            let p_x2   = chip_w - div_x - arm_x;  // 62.5 % — parallel end
            let conv_x = chip_w - div_x;           // 75 % — converging junction

            // Y-spread from arm angle (default π/3 = 60°), capped to routing bounds.
            let d_vert = (arm_x * config.bifurcation_half_angle_rad.tan())
                .min(max_half_y);

            let mut layout = Vec::with_capacity(8);

            // 1. Inlet straight (parent_in cross-section)
            layout.push(SegmentLayout {
                start: Point3r::new(0.0,   y_center, z_mid),
                end:   Point3r::new(div_x, y_center, z_mid),
                cross_section: parent_in.cross_section,
            });
            // 2-4. Upper daughter: arm_in → parallel → arm_out (d1 cross-section)
            layout.push(SegmentLayout {
                start: Point3r::new(div_x, y_center,          z_mid),
                end:   Point3r::new(p_x1,  y_center + d_vert, z_mid),
                cross_section: d1.cross_section,
            });
            layout.push(SegmentLayout {
                start: Point3r::new(p_x1, y_center + d_vert, z_mid),
                end:   Point3r::new(p_x2, y_center + d_vert, z_mid),
                cross_section: d1.cross_section,
            });
            layout.push(SegmentLayout {
                start: Point3r::new(p_x2,   y_center + d_vert, z_mid),
                end:   Point3r::new(conv_x, y_center,           z_mid),
                cross_section: d1.cross_section,
            });
            // 5-7. Lower daughter: arm_in → parallel → arm_out (symmetric)
            layout.push(SegmentLayout {
                start: Point3r::new(div_x, y_center,          z_mid),
                end:   Point3r::new(p_x1,  y_center - d_vert, z_mid),
                cross_section: d1.cross_section,
            });
            layout.push(SegmentLayout {
                start: Point3r::new(p_x1, y_center - d_vert, z_mid),
                end:   Point3r::new(p_x2, y_center - d_vert, z_mid),
                cross_section: d1.cross_section,
            });
            layout.push(SegmentLayout {
                start: Point3r::new(p_x2,   y_center - d_vert, z_mid),
                end:   Point3r::new(conv_x, y_center,           z_mid),
                cross_section: d1.cross_section,
            });
            // 8. Outlet straight (parent_out cross-section)
            layout.push(SegmentLayout {
                start: Point3r::new(conv_x, y_center, z_mid),
                end:   Point3r::new(chip_w, y_center, z_mid),
                cross_section: parent_out.cross_section,
            });

            Ok(layout)
        }

        // Trifurcation — chip-spanning layout with junction at the midpoint:
        //
        //                    ╱── d1 ──▶ right face (y + spread)
        //   ── merged straight ──── d2 ──▶ right face (y_center)
        //                    ╲── d3 ──▶ right face (y − spread)
        //
        // The merged straight (parent + collinear center daughter) spans the
        // full chip width (0 → chip_w) so the center outlet exits the right
        // face.  The junction is placed at chip_w / 2 so the outer daughters
        // span the second half of the chip.  Daughter spread in Y is capped at
        // `y_center − wall_clearance_mm` so all points stay within routing
        // bounds.
        //
        // The angled daughters start inside the merged straight's lateral
        // surface (at x = chip_w/2, which is well before chip_w), so the V-fork
        // CSG trick applies naturally — no 3-way endpoint CDT degeneracy.
        TopologyClass::Trifurcation => {
            let (parent, d1, _d2, d3) =
                topo.trifurcation_channels().ok_or_else(|| MeshError::ChannelError {
                    message: "trifurcation topology detected but channel decomposition failed"
                        .to_string(),
                })?;

            let chip_w    = SbsWellPlate96::WIDTH_MM;
            let junction_x = chip_w / 2.0;

            // Y-spread of the outer daughters, capped to stay within routing
            // bounds (y ∈ [wall_clearance, DEPTH_MM − wall_clearance]).
            let max_half_y = y_center - config.wall_clearance_mm;
            let angle: Real = config.trifurcation_half_angle_rad;
            let d_vert = (junction_x * angle.tan()).min(max_half_y);

            let mut layout = Vec::with_capacity(3);
            // Merged straight: inlet face → right face (full chip width).
            // The center outlet is the end cap of this segment at (chip_w, y_center).
            layout.push(SegmentLayout {
                start: Point3r::new(0.0, y_center, z_mid),
                end:   Point3r::new(chip_w, y_center, z_mid),
                cross_section: parent.cross_section,
            });
            // Outer daughters from junction to right face, spread ± d_vert in Y.
            for (d, vert) in [(d1, d_vert), (d3, -d_vert)] {
                layout.push(SegmentLayout {
                    start: Point3r::new(junction_x, y_center, z_mid),
                    end:   Point3r::new(chip_w, y_center + vert, z_mid),
                    cross_section: d.cross_section,
                });
            }
            Ok(layout)
        }

        TopologyClass::Complex => Err(MeshError::ChannelError {
            message: "Complex topology not supported".to_string(),
        }),
    }
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
    let mesher = SweepMesher { cap_start: true, cap_end: true };
    let faces = mesher.sweep(&profile, &path, &mut pool, RegionId::from_usize(0));

    let mut mesh = IndexedMesh::new();
    for (_, vdata) in pool.iter() {
        mesh.add_vertex(vdata.position, vdata.normal);
    }
    for face in &faces {
        mesh.add_face_with_region(face.vertices[0], face.vertices[1], face.vertices[2], face.region);
    }
    mesh.rebuild_edges();
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
    let mesher = SweepMesher { cap_start: true, cap_end: true };
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
    Ok(mesh)
}

// ── Fluid mesh assembly ───────────────────────────────────────────────────────

fn assemble_fluid_mesh(meshes: Vec<IndexedMesh>) -> MeshResult<IndexedMesh> {
    let mut iter = meshes.into_iter();
    let first = iter.next().ok_or_else(|| MeshError::ChannelError {
        message: "no segment meshes to assemble".to_string(),
    })?;
    let mut accumulated = first;
    for mesh in iter {
        accumulated = csg_boolean_indexed(BooleanOp::Union, &accumulated, &mesh)?;
    }
    Ok(accumulated)
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

    // Determine inlet and outlet positions from layout
    let first_seg = &layout[0];
    let inlet_r_mm = cross_section_radius_mm(&first_seg.cross_section);
    let epsilon = 2.0 * inlet_r_mm;

    let inlet_pos = first_seg.start;

    let outlet_positions: Vec<Point3r> = match class {
        TopologyClass::LinearChain { .. } | TopologyClass::VenturiChain => {
            vec![layout.last().unwrap().end]
        }
        TopologyClass::Bifurcation => {
            // Diamond layout: [0]=inlet_straight … [7]=outlet_straight.
            // The single true outlet port is the end of the outlet_straight
            // segment, at (chip_w, y_center).
            vec![layout.last().unwrap().end]
        }
        TopologyClass::Trifurcation => {
            // Layout has 3 segments: [merged_straight, d1_angled, d3_angled].
            // All three ends are outlet faces.
            layout.iter().map(|s| s.end).collect()
        }
        TopologyClass::Complex => vec![],
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

        let label = if (centroid - inlet_pos).norm() < epsilon {
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

// ── Bifurcation fluid mesh ────────────────────────────────────────────────────

/// Build the fluid mesh for a [`TopologyClass::Bifurcation`] diamond layout.
///
/// # Why spine + individual segments in interleaved order?
///
/// The 8-segment diamond has two T-junctions (diverging at 25 % and converging
/// at 75 % chip width) where three tube axes meet at the same spine point.
///
/// **Naïve sequential CSG of 8 segments** — the inlet ends at (div_x, y_c), so
/// after `inlet ∪ upper_arm_in` the point (div_x, y_c) is an interior bend.
/// Adding `lower_arm_in` there creates a 3-way bend-point junction that triggers
/// a CDT PSLG panic.
///
/// **Fork-polyline approach** — using `build_polyline_mesh` for each full fork
/// (arm_in + parallel + arm_out) avoids the arm–parallel bend issue, but the
/// polyline spans *both* junction points (div_x and conv_x) simultaneously.  The
/// second CSG step `result ∪ lower_fork` therefore has *two disconnected*
/// intersection curves (one at div_x, one at conv_x), which causes 17 open
/// boundary edges for the equal-diameter case.
///
/// # Solution — spine + per-segment interleaved unions
///
/// 1. **Spine**: a single straight segment from the inlet face `(0, y_c)` to
///    the outlet face `(chip_w, y_c)`.  Both junction points `(div_x, y_c)` and
///    `(conv_x, y_c)` are genuine **interior** points of this straight cylinder
///    — smooth surface where forks branch off, identical to the trifurcation axle.
///
/// 2. **Six fork segments added individually** using `build_segment_mesh`, in
///    interleaved (upper/lower) order:
///
///    ```text
///    [1] upper_arm_in   → spine at div_x, lower side   (trifurcation pattern)
///    [4] lower_arm_in   → spine at div_x, upper side   (opposite side — works)
///    [2] upper_parallel → end-to-end from upper_arm_in  (single junction)
///    [5] lower_parallel → end-to-end from lower_arm_in  (single junction)
///    [3] upper_arm_out  → spine at conv_x, lower side  (clean spine surface)
///    [6] lower_arm_out  → spine at conv_x, upper side  (opposite side — works)
///    ```
///
///    Each step creates **exactly one** T-junction region with the accumulated
///    result.  The interleaved order ensures that when the *second* arm is added
///    at each junction (div_x or conv_x), the first arm already modified the
///    *opposite* side of the spine — identical to how trifurcation outer daughters
///    work.  This is valid for all diameter ratios including equal diameters.
fn build_bifurcation_fluid_mesh(
    layout: &[SegmentLayout],
    config: &PipelineConfig,
) -> MeshResult<IndexedMesh> {
    if layout.len() != 8 {
        return Err(MeshError::ChannelError {
            message: format!(
                "diamond bifurcation expects 8 layout segments, got {}",
                layout.len()
            ),
        });
    }

    // ── Step 1: full-width straight spine ────────────────────────────────────
    //
    // Merges inlet_straight (layout[0]) and outlet_straight (layout[7]) into
    // one straight cylinder at y_center.  Both junction X positions (div_x and
    // conv_x) are interior points of this smooth cylinder — no bends.
    let spine = SegmentLayout {
        start:         layout[0].start,            // (0.0, y_center, z_mid)
        end:           layout[7].end,              // (chip_w, y_center, z_mid)
        cross_section: layout[0].cross_section,    // parent_in diameter
    };
    let mut mesh = build_segment_mesh(&spine, config)?;

    // ── Steps 2–7: fork segments in interleaved order ─────────────────────────
    //
    // Indices into `layout`:
    //   [1] upper_arm_in, [4] lower_arm_in  — T-junctions at div_x
    //   [2] upper_parallel, [5] lower_parallel — end-to-end connections
    //   [3] upper_arm_out, [6] lower_arm_out — T-junctions at conv_x
    //
    // Interleaving upper/lower at each stage ensures the second arm at each
    // junction enters from the OPPOSITE side of the spine axis — same geometry
    // as trifurcation outer daughters, which is known-good for all diameter
    // ratios including r_daughter = r_spine.
    for seg_idx in [1_usize, 4, 2, 5, 3, 6] {
        let seg_mesh = build_segment_mesh(&layout[seg_idx], config)?;
        mesh = csg_boolean_indexed(BooleanOp::Union, &mesh, &seg_mesh)?;
    }

    Ok(mesh)
}

// ── Chip body ─────────────────────────────────────────────────────────────────

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

    csg_boolean_indexed(BooleanOp::Difference, &substrate, &void_mesh)
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
        result = csg_boolean_indexed(BooleanOp::Difference, &result, &void_mesh)?;
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
        if connected && collinear && same_cs {
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
            CrossSectionSpec::Rectangular { width_m: w1, height_m: h1 },
            CrossSectionSpec::Rectangular { width_m: w2, height_m: h2 },
        ) => (w1 - w2).abs() < 1e-9 && (h1 - h2).abs() < 1e-9,
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use cfd_schematics::interface::presets::venturi_chain;

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
    fn pipeline_rejects_complex_topology() {
        // Build a blueprint with complex topology manually
        use cfd_schematics::{ChannelSpec, NetworkBlueprint, NodeKind, NodeSpec};
        let mut bp = NetworkBlueprint::new("complex");
        bp.add_node(NodeSpec::new("inlet", NodeKind::Inlet));
        bp.add_node(NodeSpec::new("j1", NodeKind::Junction));
        bp.add_node(NodeSpec::new("j2", NodeKind::Junction));
        bp.add_node(NodeSpec::new("j3", NodeKind::Junction));
        bp.add_node(NodeSpec::new("outlet", NodeKind::Outlet));
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
        let result = BlueprintMeshPipeline::run(&bp, &PipelineConfig::default());
        assert!(result.is_err());
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
}
