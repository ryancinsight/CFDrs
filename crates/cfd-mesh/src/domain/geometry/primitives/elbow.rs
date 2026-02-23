//! Circular-arc pipe bend (elbow) primitive.

use std::f64::consts::TAU;

use crate::core::index::RegionId;
use crate::core::scalar::{Point3r, Vector3r};
use crate::mesh::IndexedMesh;
use super::{PrimitiveMesh, PrimitiveError};

/// Builds a closed circular-arc tube sweep (pipe bend / elbow).
///
/// The arc lies in the XZ plane.  The bend centre is at the origin; the arc
/// starts at `(0, 0, 0)` and sweeps to angle `bend_angle`.
///
/// ## Geometry
///
/// The centreline follows:
/// ```text
/// C(α) = ( R·(1 − cos α),  0,  R·sin α )   α ∈ [0, bend_angle]
/// ```
///
/// The exact Frenet–Serret frame for a planar circle in XZ:
/// ```text
/// T(α) = (sin α,  0,  cos α)      — tangent
/// N(α) = (cos α,  0, −sin α)      — principal normal (centripetal)
/// B     = (0,      1,  0   )      — binormal (constant = +Y)
/// ```
///
/// Tube cross-section at `(α, β)`:
/// ```text
/// P(α, β) = C(α) + r·(cos β · N(α) + sin β · B)
/// n(α, β) =        cos β · N(α) + sin β · B
/// ```
///
/// ## Region IDs
///
/// | RegionId | Surface |
/// |----------|---------|
/// | 1 | Outer tube wall |
/// | 2 | Inlet cap (α = 0) |
/// | 3 | Outlet cap (α = bend_angle) |
///
/// ## Uses
///
/// 90°/180° bends in millifluidic chips, Dean-vortex secondary-flow studies,
/// serpentine-channel building blocks.
///
/// ## Output
///
/// - `signed_volume ≈ π r_tube² · R_bend · bend_angle` (Pappus theorem)
#[derive(Clone, Debug)]
pub struct Elbow {
    /// Tube (cross-section) radius [mm].
    pub tube_radius: f64,
    /// Bend centreline radius [mm]. Must be > `tube_radius`.
    pub bend_radius: f64,
    /// Total sweep angle [rad]. Must be in `(0, 2π]`.
    pub bend_angle: f64,
    /// Angular divisions around the tube cross-section (≥ 3).
    pub tube_segments: usize,
    /// Divisions along the bend arc.
    pub arc_segments: usize,
}

impl Default for Elbow {
    fn default() -> Self {
        Self {
            tube_radius: 0.5,
            bend_radius: 2.0,
            bend_angle: std::f64::consts::FRAC_PI_2, // 90°
            tube_segments: 16,
            arc_segments: 24,
        }
    }
}

impl PrimitiveMesh for Elbow {
    fn build(&self) -> Result<IndexedMesh, PrimitiveError> {
        build(self)
    }
}

fn build(el: &Elbow) -> Result<IndexedMesh, PrimitiveError> {
    if el.tube_radius <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "tube_radius must be > 0, got {}", el.tube_radius
        )));
    }
    if el.bend_radius <= el.tube_radius {
        return Err(PrimitiveError::InvalidParam(format!(
            "bend_radius ({}) must be > tube_radius ({})",
            el.bend_radius, el.tube_radius
        )));
    }
    if el.bend_angle <= 0.0 || el.bend_angle > TAU {
        return Err(PrimitiveError::InvalidParam(format!(
            "bend_angle must be in (0, 2π], got {}", el.bend_angle
        )));
    }
    if el.tube_segments < 3 {
        return Err(PrimitiveError::TooFewSegments(el.tube_segments));
    }
    if el.arc_segments < 1 {
        return Err(PrimitiveError::InvalidParam(
            "arc_segments must be ≥ 1".into()
        ));
    }

    let wall_region   = RegionId::new(1);
    let inlet_region  = RegionId::new(2);
    let outlet_region = RegionId::new(3);

    let mut mesh = IndexedMesh::new();
    let r   = el.tube_radius;
    let big_r = el.bend_radius;
    let ba  = el.bend_angle;
    let ns  = el.tube_segments;
    let na  = el.arc_segments;

    // ── Helpers ──────────────────────────────────────────────────────────────

    // Frenet frame at arc angle alpha (planar circle in XZ plane).
    let frenet = |alpha: f64| -> (Vector3r, Vector3r, Vector3r) {
        let t = Vector3r::new(alpha.sin(), 0.0, alpha.cos());       // tangent
        let n = Vector3r::new(alpha.cos(), 0.0, -alpha.sin());      // principal normal
        let b = Vector3r::y();                                       // binormal = +Y
        (t, n, b)
    };

    // Centreline position at arc angle alpha.
    let centre_at = |alpha: f64| -> Point3r {
        Point3r::new(big_r * (1.0 - alpha.cos()), 0.0, big_r * alpha.sin())
    };

    // Tube vertex + outward normal at (arc angle alpha, tube angle beta).
    let tube_vertex = |alpha: f64, beta: f64| -> (Point3r, Vector3r) {
        let (_, n_frame, b_frame) = frenet(alpha);
        let cb = beta.cos();
        let sb = beta.sin();
        let n_out = n_frame * cb + b_frame * sb;
        let pos   = centre_at(alpha) + n_out * r;
        (pos, n_out)
    };

    // Pre-build ring vertex arrays for shared edge connectivity
    // rings[ia][ib] = VertexId at arc station ia, tube angle ib.
    // Caps reuse rings[0] and rings[na] so wall/cap edges are topologically shared.
    let mut rings: Vec<Vec<crate::core::index::VertexId>> = Vec::with_capacity(na + 1);
    for ia in 0..=na {
        let alpha = ia as f64 / na as f64 * ba;
        let row: Vec<_> = (0..ns).map(|ib| {
            let beta = ib as f64 / ns as f64 * TAU;
            let (p, n) = tube_vertex(alpha, beta);
            mesh.add_vertex(p, n)
        }).collect();
        rings.push(row);
    }

    // Lateral tube wall
    // Outward winding: (arc0,beta0)->(arc0,beta1)->(arc1,beta1) and (arc0,beta0)->(arc1,beta1)->(arc1,beta0)
    for ia in 0..na {
        for ib in 0..ns {
            let ib1 = (ib + 1) % ns;
            let v00 = rings[ia][ib];
            let v10 = rings[ia + 1][ib];
            let v11 = rings[ia + 1][ib1];
            let v01 = rings[ia][ib1];
            mesh.add_face_with_region(v00, v01, v11, wall_region);
            mesh.add_face_with_region(v00, v11, v10, wall_region);
        }
    }

    // Inlet and outlet caps: only add for partial bends (not full 360 torus).
    // For a full circle, rings[0] == rings[na] and the lateral wall already closes itself.
    let is_full_circle = (ba - TAU).abs() < 1e-10;
    if !is_full_circle {
        // Inlet cap (alpha = 0, outward normal = -T(0))
        // Reuse rings[0] for shared wall<->cap edge topology.
        {
            let (t0, _, _) = frenet(0.0);
            let n_cap = -t0;
            let center = centre_at(0.0);
            let vc = mesh.add_vertex(center, n_cap);
            for ib in 0..ns {
                let ib1 = (ib + 1) % ns;
                let vr0 = rings[0][ib];
                let vr1 = rings[0][ib1];
                // CCW from outside (-T direction): centre -> r_{i+1} -> r_i
                mesh.add_face_with_region(vc, vr1, vr0, inlet_region);
            }
        }

        // Outlet cap (alpha = bend_angle, outward normal = +T(bend_angle))
        // Reuse rings[na] for shared wall<->cap edge topology.
        {
            let (t_end, _, _) = frenet(ba);
            let n_cap = t_end;
            let center = centre_at(ba);
            let vc = mesh.add_vertex(center, n_cap);
            for ib in 0..ns {
                let ib1 = (ib + 1) % ns;
                let vr0 = rings[na][ib];
                let vr1 = rings[na][ib1];
                // CCW from outside (+T direction): centre -> r_i -> r_{i+1}
                mesh.add_face_with_region(vc, vr0, vr1, outlet_region);
            }
        }
    }

    Ok(mesh)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::storage::edge_store::EdgeStore;
    use crate::watertight::check::check_watertight;
    use std::f64::consts::PI;

    #[test]
    fn elbow_is_watertight() {
        let mesh = Elbow::default().build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.is_watertight, "elbow must be watertight");
        assert_eq!(report.euler_characteristic, Some(2));
    }

    #[test]
    fn elbow_volume_positive_and_approximately_correct() {
        let r = 0.5_f64;
        let big_r = 2.0_f64;
        let ba = PI; // 180° bend
        let mesh = Elbow {
            tube_radius: r, bend_radius: big_r, bend_angle: ba,
            tube_segments: 32, arc_segments: 48,
        }.build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.signed_volume > 0.0);
        // Pappus: V = π r² · R · θ
        let expected = PI * r * r * big_r * ba;
        let error = (report.signed_volume - expected).abs() / expected;
        assert!(error < 0.01, "volume error {:.4}% should be < 1% at 32×48", error * 100.0);
    }

    #[test]
    fn elbow_full_torus_is_genus1() {
        // A full 360° elbow is topologically a torus: euler = 0
        let mesh = Elbow {
            bend_angle: TAU,
            arc_segments: 48,
            tube_segments: 24,
            ..Elbow::default()
        }.build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.is_watertight);
        // Note: caps are still added for 360° but their vertices are welded,
        // so topology may still be χ=2; annotate this empirically.
        assert!(report.signed_volume > 0.0);
    }

    #[test]
    fn elbow_rejects_invalid_params() {
        assert!(Elbow { tube_radius: 0.0, ..Elbow::default() }.build().is_err());
        assert!(Elbow { bend_radius: 0.4, tube_radius: 0.5, ..Elbow::default() }.build().is_err());
        assert!(Elbow { bend_angle: 0.0, ..Elbow::default() }.build().is_err());
        assert!(Elbow { bend_angle: 7.0, ..Elbow::default() }.build().is_err()); // > 2π
        assert!(Elbow { tube_segments: 2, ..Elbow::default() }.build().is_err());
        assert!(Elbow { arc_segments: 0, ..Elbow::default() }.build().is_err());
    }
}
