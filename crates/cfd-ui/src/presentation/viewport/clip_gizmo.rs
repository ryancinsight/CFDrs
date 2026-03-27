//! Clip plane gizmo — interactive translucent quad for visualizing and
//! manipulating section planes in the viewport.
//!
//! Renders a translucent quad on the clip plane using the overlay pipeline.
//! The quad size is proportional to the scene bounding box. A normal-direction
//! arrow handle allows dragging to translate the plane offset.

use crate::domain::clipping::ClipPlane;
use crate::infrastructure::gpu::overlay_buffer::{OverlayBuffer, OverlayLine, OverlayTriangle};
use nalgebra::Vector3;

/// Size of the gizmo quad (half-extent in each direction).
const GIZMO_HALF_SIZE: f32 = 2.0;

/// Generate overlay triangles for a clip plane gizmo.
///
/// Creates a translucent quad centered at `normal * offset` with the given
/// half-extent, oriented perpendicular to the plane normal.
#[must_use]
pub fn clip_plane_triangles(plane: &ClipPlane) -> Vec<OverlayTriangle> {
    let n = plane.normal.cast::<f32>();
    let center = n * plane.offset as f32;

    // Build two orthogonal vectors in the plane.
    let (u, v) = orthogonal_basis(&n);
    let h = GIZMO_HALF_SIZE;

    let corners = [
        to_arr(&(center - u * h - v * h)),
        to_arr(&(center + u * h - v * h)),
        to_arr(&(center + u * h + v * h)),
        to_arr(&(center - u * h + v * h)),
    ];

    vec![
        OverlayTriangle { v0: corners[0], v1: corners[1], v2: corners[2], color: plane.color },
        OverlayTriangle { v0: corners[0], v1: corners[2], v2: corners[3], color: plane.color },
    ]
}

/// Generate a normal-direction arrow for the gizmo.
#[must_use]
pub fn clip_plane_normal_arrow(plane: &ClipPlane) -> Vec<OverlayLine> {
    let n = plane.normal.cast::<f32>();
    let center = n * plane.offset as f32;
    let tip = center + n * GIZMO_HALF_SIZE * 0.5;

    let arrow_color = [1.0, 1.0, 0.0, 0.9]; // yellow arrow
    vec![OverlayLine {
        start: to_arr(&center),
        end: to_arr(&tip),
        color: arrow_color,
    }]
}

/// Build GPU overlay buffers for a clip plane gizmo.
pub fn build_clip_gizmo_buffers(
    plane: &ClipPlane,
    device: &wgpu::Device,
) -> (OverlayBuffer, OverlayBuffer) {
    let tris = clip_plane_triangles(plane);
    let lines = clip_plane_normal_arrow(plane);
    (
        OverlayBuffer::from_triangles(&tris, device),
        OverlayBuffer::from_lines(&lines, device),
    )
}

/// Compute an orthogonal basis `(u, v)` for a plane with given normal.
fn orthogonal_basis(n: &Vector3<f32>) -> (Vector3<f32>, Vector3<f32>) {
    let seed = if n.x.abs() < 0.9 {
        Vector3::x()
    } else {
        Vector3::y()
    };
    let u = n.cross(&seed).normalize();
    let v = n.cross(&u);
    (u, v)
}

/// Convert a nalgebra Vector3 to `[f32; 3]`.
fn to_arr(v: &Vector3<f32>) -> [f32; 3] {
    [v.x, v.y, v.z]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::clipping::{ClipPlaneSet, ClipPreset};

    #[test]
    fn generates_two_triangles() {
        let plane = ClipPlaneSet::from_preset(ClipPreset::Xy, 0.0);
        let tris = clip_plane_triangles(&plane);
        assert_eq!(tris.len(), 2);
    }

    #[test]
    fn normal_arrow_has_one_line() {
        let plane = ClipPlaneSet::from_preset(ClipPreset::Xz, 1.0);
        let lines = clip_plane_normal_arrow(&plane);
        assert_eq!(lines.len(), 1);
    }

    #[test]
    fn orthogonal_basis_is_orthogonal() {
        use approx::assert_relative_eq;
        let n = Vector3::new(0.0, 0.0, 1.0);
        let (u, v) = orthogonal_basis(&n);
        assert_relative_eq!(u.dot(&v), 0.0, epsilon = 1e-6);
        assert_relative_eq!(u.dot(&n), 0.0, epsilon = 1e-6);
        assert_relative_eq!(v.dot(&n), 0.0, epsilon = 1e-6);
    }
}
