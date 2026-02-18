//! Profile sweep along a path to generate channel geometry.
//!
//! This is the core extrusion engine — the equivalent of blue2mesh's
//! `ExtrusionEngine` but producing indexed mesh faces directly.

use crate::core::index::{VertexId, RegionId};
use crate::storage::face_store::FaceData;
use crate::storage::vertex_pool::VertexPool;
use crate::channel::profile::ChannelProfile;
use crate::channel::path::ChannelPath;

/// Sweep mesher: sweeps a 2D profile along a 3D path.
pub struct SweepMesher {
    /// Whether to cap the start of the sweep.
    pub cap_start: bool,
    /// Whether to cap the end of the sweep.
    pub cap_end: bool,
}

impl SweepMesher {
    /// Create with default settings (both ends capped).
    pub fn new() -> Self {
        Self {
            cap_start: true,
            cap_end: true,
        }
    }

    /// Sweep a profile along a path, producing indexed faces.
    ///
    /// Returns the list of generated faces. New vertices are inserted into
    /// `vertex_pool` via welding.
    pub fn sweep(
        &self,
        profile: &ChannelProfile,
        path: &ChannelPath,
        vertex_pool: &mut VertexPool,
        region: RegionId,
    ) -> Vec<FaceData> {
        let profile_pts = profile.generate_points();
        let frames = path.compute_frames();
        let n_profile = profile_pts.len();
        let n_stations = frames.len();

        // Generate vertex rings at each station
        let mut rings: Vec<Vec<VertexId>> = Vec::with_capacity(n_stations);

        for frame in &frames {
            let mut ring = Vec::with_capacity(n_profile);
            for pt2d in &profile_pts {
                let pos = frame.position
                    + frame.normal * pt2d[0]
                    + frame.binormal * pt2d[1];
                let outward = (pos - frame.position).normalize();
                let vid = vertex_pool.insert_or_weld(pos, outward);
                ring.push(vid);
            }
            rings.push(ring);
        }

        let mut faces = Vec::new();

        // Connect adjacent rings with quad strips (split into triangles)
        for s in 0..(n_stations - 1) {
            let ring_a = &rings[s];
            let ring_b = &rings[s + 1];
            for i in 0..n_profile {
                let j = (i + 1) % n_profile;

                // Quad: ring_a[i], ring_b[i], ring_b[j], ring_a[j]
                // Split into two triangles
                faces.push(FaceData {
                    vertices: [ring_a[i], ring_b[i], ring_b[j]],
                    region,
                });
                faces.push(FaceData {
                    vertices: [ring_a[i], ring_b[j], ring_a[j]],
                    region,
                });
            }
        }

        // Cap start
        if self.cap_start && n_profile >= 3 {
            let center_pos = frames[0].position;
            let center_normal = -frames[0].tangent;
            let center = vertex_pool.insert_or_weld(center_pos, center_normal);
            let ring = &rings[0];
            for i in 0..n_profile {
                let j = (i + 1) % n_profile;
                // Reverse winding for inward-facing cap
                faces.push(FaceData {
                    vertices: [center, ring[j], ring[i]],
                    region,
                });
            }
        }

        // Cap end
        if self.cap_end && n_profile >= 3 {
            let last = n_stations - 1;
            let center_pos = frames[last].position;
            let center_normal = frames[last].tangent;
            let center = vertex_pool.insert_or_weld(center_pos, center_normal);
            let ring = &rings[last];
            for i in 0..n_profile {
                let j = (i + 1) % n_profile;
                faces.push(FaceData {
                    vertices: [center, ring[i], ring[j]],
                    region,
                });
            }
        }

        faces
    }
}

impl Default for SweepMesher {
    fn default() -> Self {
        Self::new()
    }
}
