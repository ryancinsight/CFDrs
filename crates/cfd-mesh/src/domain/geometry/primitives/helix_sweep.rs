//! Helix sweep primitive — circular tube swept along a helix.

use std::f64::consts::TAU;

use super::{PrimitiveError, PrimitiveMesh};
use crate::domain::core::index::RegionId;
use crate::domain::core::scalar::{Point3r, Vector3r};
use crate::domain::mesh::IndexedMesh;

/// Builds a closed circular tube swept along a helical centreline.
///
/// The helix axis is +Y.  The centreline follows:
/// ```text
/// C(θ) = ( R·cos θ,   pitch·θ/(2π),   R·sin θ )     θ ∈ [0, 2π·turns]
/// ```
///
/// ## Frenet–Serret frame (exact analytic)
///
/// Let `L = √(R² + (pitch/2π)²)` (constant arc-length scale factor).
///
/// ```text
/// T(θ) = (−R·sin θ,  pitch/(2π),   R·cos θ ) / L     — tangent
/// N(θ) = (−cos θ,    0,           −sin θ   )          — centripetal normal
/// B(θ) = (−(pitch/2π)·sin θ,  −R,  (pitch/2π)·cos θ) / L   — binormal
/// ```
///
/// Tube cross-section: `P(θ, β) = C(θ) + r·(cos β · N + sin β · B)`
///
/// ## Region IDs
///
/// | RegionId | Surface |
/// |----------|---------|
/// | 1 | Tube wall |
/// | 2 | Inlet cap  (θ = 0)            |
/// | 3 | Outlet cap (θ = 2π·turns)     |
///
/// ## Uses
///
/// Spiral separators, Dean-vortex mixing channels, SDT fibre bundles,
/// helical heat-exchanger microchannels.
///
/// ## Output
///
/// - `signed_volume ≈ π·r_tube²·arc_length`  where
///   `arc_length = turns·√((2πR)² + pitch²)` (Pappus theorem)
#[derive(Clone, Debug)]
pub struct HelixSweep {
    /// Helix (coil) centreline radius [mm].
    pub coil_radius: f64,
    /// Tube cross-section radius [mm].
    pub tube_radius: f64,
    /// Axial rise per full turn [mm].
    pub pitch: f64,
    /// Number of complete turns (> 0).
    pub turns: f64,
    /// Angular segments around the tube cross-section (≥ 3).
    pub tube_segments: usize,
    /// Arc segments per full turn (≥ 3).
    pub arc_segments_per_turn: usize,
}

impl Default for HelixSweep {
    fn default() -> Self {
        Self {
            coil_radius: 2.0,
            tube_radius: 0.3,
            pitch: 1.0,
            turns: 2.0,
            tube_segments: 16,
            arc_segments_per_turn: 32,
        }
    }
}

impl PrimitiveMesh for HelixSweep {
    fn build(&self) -> Result<IndexedMesh, PrimitiveError> {
        build(self)
    }
}

fn build(hs: &HelixSweep) -> Result<IndexedMesh, PrimitiveError> {
    if hs.coil_radius <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "coil_radius must be > 0, got {}",
            hs.coil_radius
        )));
    }
    if hs.tube_radius <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "tube_radius must be > 0, got {}",
            hs.tube_radius
        )));
    }
    if hs.tube_radius >= hs.coil_radius {
        return Err(PrimitiveError::InvalidParam(format!(
            "tube_radius ({}) must be < coil_radius ({})",
            hs.tube_radius, hs.coil_radius
        )));
    }
    if hs.pitch <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "pitch must be > 0, got {}",
            hs.pitch
        )));
    }
    if hs.turns <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "turns must be > 0, got {}",
            hs.turns
        )));
    }
    if hs.tube_segments < 3 {
        return Err(PrimitiveError::TooFewSegments(hs.tube_segments));
    }
    if hs.arc_segments_per_turn < 3 {
        return Err(PrimitiveError::InvalidParam(
            "arc_segments_per_turn must be ≥ 3".into(),
        ));
    }

    let wall_region = RegionId::new(1);
    let inlet_region = RegionId::new(2);
    let outlet_region = RegionId::new(3);

    let mut mesh = IndexedMesh::new();

    let big_r = hs.coil_radius;
    let r = hs.tube_radius;
    let pitch = hs.pitch;
    let ns = hs.tube_segments;
    let na = (hs.arc_segments_per_turn as f64 * hs.turns).round() as usize;
    let theta_max = TAU * hs.turns;

    // Constant arc-length factor L = sqrt(R² + (pitch/2π)²)
    let p2pi = pitch / TAU; // pitch / (2π)
    let big_l = (big_r * big_r + p2pi * p2pi).sqrt();

    // ── Frenet–Serret frame at helix parameter θ ─────────────────────────────
    let helix_frame = |theta: f64| -> (Vector3r, Vector3r, Vector3r, Point3r) {
        let (ct, st) = (theta.cos(), theta.sin());
        // Centreline
        let centre = Point3r::new(big_r * ct, p2pi * theta, big_r * st);
        // Tangent T
        let t = Vector3r::new(-big_r * st / big_l, p2pi / big_l, big_r * ct / big_l);
        // Principal normal N (centripetal, toward helix axis)
        let n = Vector3r::new(-ct, 0.0, -st);
        // Binormal B = T × N
        let b = t.cross(&n);
        (t, n, b, centre)
    };

    // Tube vertex + outward normal at (theta, beta).
    let tube_vertex = |theta: f64, beta: f64| -> (Point3r, Vector3r) {
        let (_, n_frame, b_frame, centre) = helix_frame(theta);
        let cb = beta.cos();
        let sb = beta.sin();
        let n_out = n_frame * cb + b_frame * sb;
        let pos = centre + n_out * r;
        (pos, n_out)
    };

    // Pre-build ring vertex arrays for shared edge connectivity
    let mut rings: Vec<Vec<crate::domain::core::index::VertexId>> = Vec::with_capacity(na + 1);
    for ia in 0..=na {
        let theta = ia as f64 / na as f64 * theta_max;
        let row: Vec<_> = (0..ns)
            .map(|ib| {
                let beta = ib as f64 / ns as f64 * TAU;
                let (p, n) = tube_vertex(theta, beta);
                mesh.add_vertex(p, n)
            })
            .collect();
        rings.push(row);
    }

    // Lateral tube wall
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

    // Inlet cap (theta = 0, outward normal = -T(0))
    {
        let (t0, _, _, centre0) = helix_frame(0.0);
        let n_cap = -t0;
        let vc = mesh.add_vertex(centre0, n_cap);
        for ib in 0..ns {
            let ib1 = (ib + 1) % ns;
            let vr0 = rings[0][ib];
            let vr1 = rings[0][ib1];
            mesh.add_face_with_region(vc, vr1, vr0, inlet_region);
        }
    }

    // Outlet cap (theta = theta_max, outward normal = +T(theta_max))
    {
        let (t_end, _, _, centre_end) = helix_frame(theta_max);
        let n_cap = t_end;
        let vc = mesh.add_vertex(centre_end, n_cap);
        for ib in 0..ns {
            let ib1 = (ib + 1) % ns;
            let vr0 = rings[na][ib];
            let vr1 = rings[na][ib1];
            mesh.add_face_with_region(vc, vr0, vr1, outlet_region);
        }
    }

    Ok(mesh)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::application::watertight::check::check_watertight;
    use crate::infrastructure::storage::edge_store::EdgeStore;
    use std::f64::consts::PI;

    #[test]
    fn helix_sweep_is_watertight() {
        let mesh = HelixSweep::default().build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.is_watertight, "helix_sweep must be watertight");
        assert_eq!(report.euler_characteristic, Some(2));
    }

    #[test]
    fn helix_sweep_volume_positive_and_pappus() {
        let hs = HelixSweep {
            coil_radius: 2.0,
            tube_radius: 0.3,
            pitch: 1.0,
            turns: 2.0,
            tube_segments: 32,
            arc_segments_per_turn: 64,
        };
        let mesh = hs.build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.signed_volume > 0.0);
        // V ≈ π r² * arc_length; arc_length = turns * sqrt((2πR)² + pitch²)
        let r = 0.3_f64;
        let arc_len = 2.0 * ((TAU * 2.0).powi(2) + 1.0_f64.powi(2)).sqrt();
        let expected = PI * r * r * arc_len;
        let error = (report.signed_volume - expected).abs() / expected;
        assert!(
            error < 0.01,
            "Pappus volume error {:.4}% should be < 1%",
            error * 100.0
        );
    }

    #[test]
    fn helix_sweep_rejects_invalid_params() {
        assert!(HelixSweep {
            coil_radius: 0.0,
            ..HelixSweep::default()
        }
        .build()
        .is_err());
        assert!(HelixSweep {
            tube_radius: 0.0,
            ..HelixSweep::default()
        }
        .build()
        .is_err());
        assert!(HelixSweep {
            tube_radius: 3.0,
            coil_radius: 2.0,
            ..HelixSweep::default()
        }
        .build()
        .is_err());
        assert!(HelixSweep {
            pitch: 0.0,
            ..HelixSweep::default()
        }
        .build()
        .is_err());
        assert!(HelixSweep {
            turns: 0.0,
            ..HelixSweep::default()
        }
        .build()
        .is_err());
        assert!(HelixSweep {
            tube_segments: 2,
            ..HelixSweep::default()
        }
        .build()
        .is_err());
        assert!(HelixSweep {
            arc_segments_per_turn: 2,
            ..HelixSweep::default()
        }
        .build()
        .is_err());
    }
}
