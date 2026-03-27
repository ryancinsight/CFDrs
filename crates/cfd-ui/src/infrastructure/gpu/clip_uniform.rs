//! Clip plane GPU uniform — data passed to the fragment shader for clipping.

use crate::domain::clipping::ClipPlaneSet;

/// Clip plane data for the GPU fragment shader.
///
/// Each plane is stored as `vec4(nx, ny, nz, d)` in Hessian normal form.
/// The `enable_mask` bitfield indicates which planes are active: bit `i`
/// enables plane `i`. The fragment shader discards fragments where
/// `dot(normal, world_pos) + d < 0` for any enabled plane.
#[repr(C)]
#[derive(Clone, Copy, Debug, bytemuck::Pod, bytemuck::Zeroable)]
pub struct ClipUniforms {
    /// Up to 6 clip planes, each as `[nx, ny, nz, d]`.
    pub planes: [[f32; 4]; 6],
    /// Bitmask: bit `i` set means plane `i` is active.
    pub enable_mask: u32,
    /// Padding to 16-byte alignment.
    pub _pad: [u32; 3],
}

impl ClipUniforms {
    /// Build GPU uniforms from a `ClipPlaneSet`.
    #[must_use]
    pub fn from_clip_plane_set(set: &ClipPlaneSet) -> Self {
        let mut uniforms = Self::zeroed();
        let mut mask = 0u32;

        for plane in set.active_planes() {
            let i = plane.id.0 as usize;
            if i < 6 {
                uniforms.planes[i] = [
                    plane.normal.x as f32,
                    plane.normal.y as f32,
                    plane.normal.z as f32,
                    plane.offset as f32,
                ];
                mask |= 1 << i;
            }
        }

        uniforms.enable_mask = mask;
        uniforms
    }

    /// A zeroed-out uniform (no clipping).
    #[must_use]
    fn zeroed() -> Self {
        bytemuck::Zeroable::zeroed()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::clipping::{ClipPlaneSet, ClipPreset};

    #[test]
    fn empty_set_produces_zero_mask() {
        let set = ClipPlaneSet::new();
        let u = ClipUniforms::from_clip_plane_set(&set);
        assert_eq!(u.enable_mask, 0);
    }

    #[test]
    fn active_plane_sets_correct_bit() {
        let mut set = ClipPlaneSet::new();
        let plane = ClipPlaneSet::from_preset(ClipPreset::Xy, 1.5);
        let id = set.add(plane).expect("should add");
        let u = ClipUniforms::from_clip_plane_set(&set);
        assert_ne!(u.enable_mask & (1 << id.0), 0);
        assert!((u.planes[id.0 as usize][2] - 1.0).abs() < 1e-6); // Z normal
        assert!((u.planes[id.0 as usize][3] - 1.5).abs() < 1e-6); // offset
    }
}
