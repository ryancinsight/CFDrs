//! Picking application service — dispatches GPU pick readback to domain types.

pub mod highlight;

use crate::domain::picking::{PickQuery, PickResult};
use crate::domain::scene::selection::{SelectionGranularity, SelectionTarget};
use crate::infrastructure::gpu::pick_target::PickTarget;

/// Service that translates screen-space clicks into scene-level pick results.
pub struct PickService;

impl PickService {
    /// Pick the entity at a screen-space location using the GPU pick target.
    pub fn pick_at(
        query: &PickQuery,
        pick_target: &PickTarget,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
    ) -> Option<PickResult> {
        let (node_idx, face_idx) = pick_target.read_pixel(device, queue, query.x, query.y)?;

        let target = match query.granularity {
            SelectionGranularity::Body => SelectionTarget::Body(node_idx as usize),
            SelectionGranularity::Face => SelectionTarget::Face(node_idx as usize, face_idx),
            SelectionGranularity::Edge => SelectionTarget::Edge(node_idx as usize, face_idx),
            SelectionGranularity::Vertex => SelectionTarget::Vertex(node_idx as usize, face_idx),
        };

        Some(PickResult { target })
    }
}
