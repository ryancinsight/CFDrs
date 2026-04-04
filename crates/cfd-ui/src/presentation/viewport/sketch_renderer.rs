//! Sketch renderer — builds overlay geometry from sketch entities,
//! color-coded by constraint status.

use crate::domain::sketch::dof::{DofAnalysis, DofStatus};
use crate::domain::sketch::entity::SketchEntity;
use crate::domain::sketch::sketch::Sketch;
use crate::infrastructure::gpu::overlay_buffer::{OverlayBuffer, OverlayLine};

/// Colors for sketch entities by constraint status.
const FULLY_CONSTRAINED_COLOR: [f32; 4] = [0.2, 0.8, 0.2, 1.0];
const UNDER_CONSTRAINED_COLOR: [f32; 4] = [0.3, 0.5, 1.0, 1.0];
const OVER_CONSTRAINED_COLOR: [f32; 4] = [1.0, 0.2, 0.2, 1.0];
const CONSTRUCTION_COLOR: [f32; 4] = [1.0, 0.6, 0.0, 0.8];

/// Size of point handles in sketch-local units.
const POINT_HANDLE_SIZE: f32 = 0.05;

/// Build overlay lines for all sketch entities.
///
/// Lines are placed in the sketch's local 2D coordinate system (Z = 0).
/// The caller must apply the work plane model matrix to transform to world space.
#[must_use]
pub fn build_sketch_lines(sketch: &Sketch, dof: &DofAnalysis) -> (Vec<OverlayLine>, Vec<OverlayLine>) {
    let color = status_color(dof.status);
    let mut solid_lines = Vec::new();
    let mut dashed_lines = Vec::new();

    for entity in sketch.entities() {
        let target = if entity.is_construction() {
            &mut dashed_lines
        } else {
            &mut solid_lines
        };
        let c = if entity.is_construction() { CONSTRUCTION_COLOR } else { color };

        match entity {
            SketchEntity::Point(p) => {
                // Render as a small cross-hair.
                let h = POINT_HANDLE_SIZE;
                target.push(OverlayLine {
                    start: [p.x as f32 - h, p.y as f32, 0.0],
                    end: [p.x as f32 + h, p.y as f32, 0.0],
                    color: c,
                });
                target.push(OverlayLine {
                    start: [p.x as f32, p.y as f32 - h, 0.0],
                    end: [p.x as f32, p.y as f32 + h, 0.0],
                    color: c,
                });
            }
            SketchEntity::Line(l) => {
                if let (Some(start), Some(end)) = (sketch.point(l.start), sketch.point(l.end)) {
                    target.push(OverlayLine {
                        start: [start.x as f32, start.y as f32, 0.0],
                        end: [end.x as f32, end.y as f32, 0.0],
                        color: c,
                    });
                }
            }
            SketchEntity::Circle(circ) => {
                // Tessellate circle into line segments.
                if let Some(center) = sketch.point(circ.center) {
                    let cx = center.x as f32;
                    let cy = center.y as f32;
                    let r = circ.radius as f32;
                    let segments = 48;
                    for i in 0..segments {
                        let a0 = 2.0 * std::f32::consts::PI * (i as f32) / (segments as f32);
                        let a1 = 2.0 * std::f32::consts::PI * ((i + 1) as f32) / (segments as f32);
                        target.push(OverlayLine {
                            start: [cx + r * a0.cos(), cy + r * a0.sin(), 0.0],
                            end: [cx + r * a1.cos(), cy + r * a1.sin(), 0.0],
                            color: c,
                        });
                    }
                }
            }
            SketchEntity::Arc(_) | SketchEntity::Spline(_) => {
                // Arc and spline rendering deferred to future phases.
            }
        }
    }

    (solid_lines, dashed_lines)
}

/// Build GPU overlay buffers for sketch rendering.
pub fn build_sketch_buffers(
    sketch: &Sketch,
    dof: &DofAnalysis,
    device: &wgpu::Device,
) -> (Option<OverlayBuffer>, Option<OverlayBuffer>) {
    let (solid, dashed) = build_sketch_lines(sketch, dof);
    let solid_buf = if solid.is_empty() {
        None
    } else {
        Some(OverlayBuffer::from_lines(&solid, device))
    };
    let dashed_buf = if dashed.is_empty() {
        None
    } else {
        Some(OverlayBuffer::from_lines(&dashed, device))
    };
    (solid_buf, dashed_buf)
}

/// Pick color based on DOF status.
fn status_color(status: DofStatus) -> [f32; 4] {
    match status {
        DofStatus::FullyConstrained => FULLY_CONSTRAINED_COLOR,
        DofStatus::UnderConstrained => UNDER_CONSTRAINED_COLOR,
        DofStatus::OverConstrained => OVER_CONSTRAINED_COLOR,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::sketch::dof::DofStatus;
    use crate::domain::sketch::entity::SketchPoint;
    use crate::domain::sketch::sketch::{Sketch, SketchId};
    use crate::domain::sketch::work_plane::WorkPlane;

    #[test]
    fn point_generates_cross_hair() {
        let mut sk = Sketch::new(SketchId(0), "test".into(), WorkPlane::xy());
        let id = sk.next_entity_id();
        sk.add_entity(SketchEntity::Point(SketchPoint {
            id, x: 1.0, y: 2.0, construction: false,
        }));
        let dof = DofAnalysis {
            total_dofs: 2, constraint_equations: 0,
            status: DofStatus::UnderConstrained, free_dof_count: 2,
        };
        let (solid, dashed) = build_sketch_lines(&sk, &dof);
        assert_eq!(solid.len(), 2); // horizontal + vertical cross-hair
        assert!(dashed.is_empty());
    }
}
