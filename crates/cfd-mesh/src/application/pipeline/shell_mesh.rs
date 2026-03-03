//! `InterchangeShellCuboid → IndexedMesh` pipeline.
//!
//! Converts a `cfd_schematics::InterchangeShellCuboid` into watertight surface meshes
//! suitable for CFD simulation (fluid mesh) and manufacturing output (chip body).
//!
//! This pipeline perfectly aligns the TPMS logic with the cavity walls by relying
//! on the watertight SDF-capping logic embedded directly in `build_tpms_box`
//! and `build_tpms_box_graded`. Co-planar CSG unions seamlessly integrate the
//! port connection pipes into the bounded cavity.

use std::f64::consts::PI;

use crate::application::csg::boolean::{csg_boolean_indexed, BooleanOp};
use crate::domain::core::MeshError;
use crate::domain::core::MeshResult;
use crate::domain::core::scalar::Vector3r;
use crate::domain::geometry::primitives::cube::Cube;
use crate::domain::geometry::primitives::PrimitiveMesh;
use crate::domain::geometry::tpms::{
    build_tpms_box, build_tpms_box_graded, FischerKochCY, Frd, Gyroid, Iwp, Lidinoid, Neovius,
    SchwarzD, SchwarzP, SplitP, TpmsBoxParams,
};
use crate::domain::mesh::IndexedMesh;

use cfd_schematics::geometry::types::{InterchangeShellCuboid, TpmsSurfaceKind};

// ── Pipeline configuration ────────────────────────────────────────────────────

/// Configuration for the `ShellMeshPipeline`.
#[derive(Debug, Clone)]
pub struct ShellPipelineConfig {
    /// Desired 3-D height of the fluid cavity (Z-axis) [mm].  Typical: 1–2 mm.
    pub cavity_height_mm: f64,
    /// Desired 3-D height of the physical chip body (Z-axis) [mm].  Must be > cavity_height.
    pub chip_height_mm: f64,
    /// Z-coordinate for the mid-plane of the cavity [mm].
    pub z_mid_mm: f64,
}

impl Default for ShellPipelineConfig {
    fn default() -> Self {
        Self {
            cavity_height_mm: 1.5,
            chip_height_mm: 4.0,
            z_mid_mm: 2.0,
        }
    }
}

// ── Pipeline output ───────────────────────────────────────────────────────────

/// Output of the `ShellMeshPipeline`.
pub struct ShellPipelineOutput {
    /// Mathematical boundary of the fluid domain (cavity + ports).
    pub fluid_mesh: IndexedMesh,
    /// Physical plastic body (substrate - fluid_mesh).
    pub chip_body_mesh: IndexedMesh,
}

// ── Pipeline implementation ───────────────────────────────────────────────────

/// Converts an `InterchangeShellCuboid` into watertight `IndexedMesh` objects.
pub struct ShellMeshPipeline;

impl ShellMeshPipeline {
    /// Run the full pipeline.
    ///
    /// # Steps
    /// 1. Build the outer bounding chip substrate.
    /// 2. Extract bounding box for the inner cavity.
    /// 3. Create the fluid mesh:
    ///    - If TPMS fill requested: build watertight SDF-capped TPMS lattice.
    ///    - If no TPMS fill: build a plain rectangular `Cube`.
    /// 4. Create inlet and outlet port cylinders and CSG-union them.
    /// 5. Create physical chip via CSG-difference: `substrate \ fluid`.
    pub fn run(
        shell: &InterchangeShellCuboid,
        config: &ShellPipelineConfig,
    ) -> MeshResult<ShellPipelineOutput> {
        let (cw, ch) = shell.inner_dims_mm;
        let t = shell.shell_thickness_mm;

        // Origin of cavity is at (t, t) in schematic coords.
        let cvx = t;
        let cvy = t;
        let cvz = config.z_mid_mm - config.cavity_height_mm / 2.0;

        let bounds = [
            cvx,
            cvy,
            cvz,
            cvx + cw,
            cvy + ch,
            cvz + config.cavity_height_mm,
        ];

        // 1. Build the core cavity (either empty box or TPMS)
        let mut fluid_mesh = if let Some(ref fill) = shell.tpms_fill {
            // TPMS branches. The builders ensure perfect `bounds` capping.
            let res = fill.resolution;
            let iso = fill.iso_value;
            match fill.surface {
                TpmsSurfaceKind::Gyroid => Self::build_tpms(&Gyroid, &bounds, fill)?,
                TpmsSurfaceKind::SchwarzP => Self::build_tpms(&SchwarzP, &bounds, fill)?,
                TpmsSurfaceKind::SchwarzD => Self::build_tpms(&SchwarzD, &bounds, fill)?,
                TpmsSurfaceKind::Neovius => Self::build_tpms(&Neovius, &bounds, fill)?,
                TpmsSurfaceKind::Lidinoid => Self::build_tpms(&Lidinoid, &bounds, fill)?,
                TpmsSurfaceKind::Iwp => Self::build_tpms(&Iwp, &bounds, fill)?,
                TpmsSurfaceKind::SplitP => Self::build_tpms(&SplitP, &bounds, fill)?,
                TpmsSurfaceKind::Frd => Self::build_tpms(&Frd, &bounds, fill)?,
                TpmsSurfaceKind::FischerKochCY => Self::build_tpms(&FischerKochCY, &bounds, fill)?,
            }
        } else {
            // Empty cavity
            Cube {
                origin: crate::domain::core::scalar::Point3r::new(cvx, cvy, cvz),
                width: cw,
                height: ch,
                depth: config.cavity_height_mm,
            }
            .build()
            .map_err(|e| MeshError::Other(e.to_string()))?
        };

        // 2. Add port stubs to fluid mesh via Union
        for port in &shell.ports {
            // Port cylinder connects outer wall to cavity wall
            let p1 = Vector3r::new(port.outer_point_mm.0, port.outer_point_mm.1, config.z_mid_mm);
            let p2 = Vector3r::new(port.inner_point_mm.0, port.inner_point_mm.1, config.z_mid_mm);
            
            // To ensure robust overlapping, we extend the inner point slightly into the cavity.
            // Port `end` is on the inner wall. Vector goes from start (outer) to end (inner).
            let mut dir = p2 - p1;
            let len = dir.norm();
            if len > 1e-6 {
                dir.normalize_mut();
            } else {
                continue;
            }
            
            // Extend slightly inwards by 0.5 mm for a clean CSG union interface.
            let p2_extended = p2 + dir * 0.5;
            // Extend slightly outwards by 0.1 mm to ensure it cuts through the substrate perfectly
            let p1_extended = p1 - dir * 0.1;

            // Port diameter defaults to 4mm to interface with 4mm tubing,
            // or 2mm if cavity is smaller.
            let dia = config.cavity_height_mm.min(4.0);
            
            use crate::application::channel::sweep::SweepMesher;
            use crate::application::channel::profile::ChannelProfile;
            use crate::application::channel::path::ChannelPath;
            use crate::infrastructure::storage::vertex_pool::VertexPool;
            use crate::domain::core::index::RegionId;

            let profile = ChannelProfile::Circular {
                radius: dia / 2.0,
                segments: 32,
            };
            let start_ext = crate::domain::core::scalar::Point3r::from(p1_extended);
            let end_ext = crate::domain::core::scalar::Point3r::from(p2_extended);
            let path = ChannelPath::straight(start_ext, end_ext);
            let mut pool = VertexPool::default_millifluidic();
            let mesher = SweepMesher { cap_start: true, cap_end: true };
            let faces = mesher.sweep(&profile, &path, &mut pool, RegionId::new(1));
            
            let mut port_mesh = IndexedMesh::new();
            let mut vmap = std::collections::HashMap::new();
            for (vid, vdata) in pool.iter() {
                let mid = port_mesh.add_vertex(vdata.position, vdata.normal);
                vmap.insert(vid, mid);
            }
            for face in &faces {
                port_mesh.add_face_with_region(
                    *vmap.get(&face.vertices[0]).unwrap(),
                    *vmap.get(&face.vertices[1]).unwrap(),
                    *vmap.get(&face.vertices[2]).unwrap(),
                    face.region
                );
            }
            port_mesh.rebuild_edges();

            fluid_mesh = csg_boolean_indexed(BooleanOp::Union, &fluid_mesh, &port_mesh)?;
        }

        // 3. Build Chip Substrate
        // Uses the Cube primitive to create a perfect outer box.
        let substrate = Cube {
            origin: crate::domain::core::scalar::Point3r::origin(),
            width: shell.outer_dims_mm.0,
            height: shell.outer_dims_mm.1,
            depth: config.chip_height_mm,
        }.build()
        .map_err(|e| MeshError::Other(e.to_string()))?;

        // 4. Physical Chip = Substrate \ Fluid
        let chip_body_mesh = csg_boolean_indexed(BooleanOp::Difference, &substrate, &fluid_mesh)?;

        Ok(ShellPipelineOutput {
            fluid_mesh,
            chip_body_mesh,
        })
    }

    /// Internal helper to dispatch to either uniform or graded TPMS builder.
    fn build_tpms<S: crate::domain::geometry::tpms::Tpms>(
        surface: &S,
        bounds: &[f64; 6],
        fill: &cfd_schematics::geometry::types::TpmsFillSpec,
    ) -> MeshResult<IndexedMesh> {
        if let Some(ref grad) = fill.gradient {
            // Graded TPMS evaluation
            let grad_clone = grad.clone();
            let x0 = bounds[0];
            let y0 = bounds[1];
            let x1 = bounds[3];
            let y1 = bounds[4];
            let w = x1 - x0;
            let h = y1 - y0;

            let period_fn = move |px: f64, py: f64, _pz: f64| {
                // Map world coordinate (px, py) to normalised fractional (x_frac, y_frac)
                let x_frac = ((px - x0) / w).clamp(0.0, 1.0);
                let y_frac = ((py - y0) / h).clamp(0.0, 1.0);
                grad_clone.period_at(x_frac, y_frac)
            };

            let mesh = build_tpms_box_graded(
                surface,
                *bounds,
                fill.resolution,
                fill.iso_value,
                period_fn,
            ).map_err(|e| MeshError::Other(format!("TPMS grading error: {:?}", e)))?;
            Ok(mesh)
        } else {
            // Uniform TPMS
            let params = TpmsBoxParams {
                bounds: *bounds,
                period: fill.period_mm,
                resolution: fill.resolution,
                iso_value: fill.iso_value,
            };
            let mesh = build_tpms_box(surface, &params)
                .map_err(|e| MeshError::Other(format!("TPMS error: {:?}", e)))?;
            Ok(mesh)
        }
    }
}
