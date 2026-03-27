//! Constraint residuals and analytic Jacobian rows.
//!
//! Each constraint type maps to one or more scalar residual equations.
//! The Jacobian row for each residual is the vector of partial derivatives
//! with respect to the parameter vector (point coordinates and radii).

use std::ops::Range;

use crate::domain::sketch::constraint::Constraint;
use crate::domain::sketch::entity::{EntityId, SketchEntity};
use crate::domain::sketch::sketch::Sketch;

/// A sparse row of the Jacobian matrix.
#[derive(Clone, Debug)]
pub struct JacobianRow {
    /// (parameter_index, partial_derivative) pairs.
    pub entries: Vec<(usize, f64)>,
}

/// Compute residuals and Jacobian rows for a single constraint.
///
/// Returns `(residuals, jacobian_rows)` where each residual has a
/// corresponding Jacobian row.
pub fn constraint_residual_and_jacobian(
    constraint: &Constraint,
    sketch: &Sketch,
    param_map: &[(EntityId, Range<usize>)],
) -> (Vec<f64>, Vec<JacobianRow>) {
    match constraint {
        Constraint::Horizontal(line_id) => horizontal_residual(sketch, param_map, *line_id),
        Constraint::Vertical(line_id) => vertical_residual(sketch, param_map, *line_id),
        Constraint::Coincident(a, b) => coincident_residual(sketch, param_map, *a, *b),
        Constraint::Fixed(point_id) => fixed_residual(sketch, param_map, *point_id),
        Constraint::Distance { a, b, value, .. } => {
            distance_residual(sketch, param_map, *a, *b, *value)
        }
        Constraint::Perpendicular(a, b) => perpendicular_residual(sketch, param_map, *a, *b),
        Constraint::Parallel(a, b) => parallel_residual(sketch, param_map, *a, *b),
        Constraint::Radius { entity, value, .. } => radius_residual(sketch, param_map, *entity, *value),
        Constraint::Midpoint { point, line } => midpoint_residual(sketch, param_map, *point, *line),
        // Constraints not yet implemented return zero residual.
        _ => (vec![0.0; constraint.residual_count()], vec![JacobianRow { entries: vec![] }; constraint.residual_count()]),
    }
}

/// Find the parameter range for an entity ID.
fn param_range(
    param_map: &[(EntityId, Range<usize>)],
    id: EntityId,
) -> Option<Range<usize>> {
    param_map.iter().find(|(eid, _)| *eid == id).map(|(_, r)| r.clone())
}

/// Get point coordinates from the sketch.
fn point_xy(sketch: &Sketch, id: EntityId) -> Option<(f64, f64)> {
    sketch.point(id).map(|p| (p.x, p.y))
}

/// Get the start/end point IDs of a line entity.
fn line_endpoints(sketch: &Sketch, id: EntityId) -> Option<(EntityId, EntityId)> {
    if let Some(SketchEntity::Line(l)) = sketch.entity(id) {
        Some((l.start, l.end))
    } else {
        None
    }
}

/// Horizontal: r = start.y - end.y
fn horizontal_residual(
    sketch: &Sketch,
    param_map: &[(EntityId, Range<usize>)],
    line_id: EntityId,
) -> (Vec<f64>, Vec<JacobianRow>) {
    let Some((s, e)) = line_endpoints(sketch, line_id) else {
        return (vec![0.0], vec![JacobianRow { entries: vec![] }]);
    };
    let Some((_, sy)) = point_xy(sketch, s) else { return (vec![0.0], vec![JacobianRow { entries: vec![] }]); };
    let Some((_, ey)) = point_xy(sketch, e) else { return (vec![0.0], vec![JacobianRow { entries: vec![] }]); };

    let r = sy - ey;
    let mut entries = Vec::new();
    if let Some(range) = param_range(param_map, s) {
        entries.push((range.start + 1, 1.0)); // d/d(sy) = 1
    }
    if let Some(range) = param_range(param_map, e) {
        entries.push((range.start + 1, -1.0)); // d/d(ey) = -1
    }
    (vec![r], vec![JacobianRow { entries }])
}

/// Vertical: r = start.x - end.x
fn vertical_residual(
    sketch: &Sketch,
    param_map: &[(EntityId, Range<usize>)],
    line_id: EntityId,
) -> (Vec<f64>, Vec<JacobianRow>) {
    let Some((s, e)) = line_endpoints(sketch, line_id) else {
        return (vec![0.0], vec![JacobianRow { entries: vec![] }]);
    };
    let Some((sx, _)) = point_xy(sketch, s) else { return (vec![0.0], vec![JacobianRow { entries: vec![] }]); };
    let Some((ex, _)) = point_xy(sketch, e) else { return (vec![0.0], vec![JacobianRow { entries: vec![] }]); };

    let r = sx - ex;
    let mut entries = Vec::new();
    if let Some(range) = param_range(param_map, s) {
        entries.push((range.start, 1.0));
    }
    if let Some(range) = param_range(param_map, e) {
        entries.push((range.start, -1.0));
    }
    (vec![r], vec![JacobianRow { entries }])
}

/// Coincident: r = [a.x - b.x, a.y - b.y]
fn coincident_residual(
    sketch: &Sketch,
    param_map: &[(EntityId, Range<usize>)],
    a: EntityId,
    b: EntityId,
) -> (Vec<f64>, Vec<JacobianRow>) {
    let Some((ax, ay)) = point_xy(sketch, a) else { return (vec![0.0; 2], vec![JacobianRow { entries: vec![] }; 2]); };
    let Some((bx, by)) = point_xy(sketch, b) else { return (vec![0.0; 2], vec![JacobianRow { entries: vec![] }; 2]); };

    let ra = param_range(param_map, a);
    let rb = param_range(param_map, b);

    let mut jx = Vec::new();
    let mut jy = Vec::new();
    if let Some(r) = &ra { jx.push((r.start, 1.0)); jy.push((r.start + 1, 1.0)); }
    if let Some(r) = &rb { jx.push((r.start, -1.0)); jy.push((r.start + 1, -1.0)); }

    (vec![ax - bx, ay - by], vec![JacobianRow { entries: jx }, JacobianRow { entries: jy }])
}

/// Fixed: residual is zero (the Jacobian rows lock the DOFs in place).
///
/// A fixed constraint prevents a point from moving during solving by
/// contributing identity Jacobian rows. The residual is always zero because
/// the point starts at its intended fixed position.
fn fixed_residual(
    _sketch: &Sketch,
    param_map: &[(EntityId, Range<usize>)],
    point_id: EntityId,
) -> (Vec<f64>, Vec<JacobianRow>) {
    let range = param_range(param_map, point_id);
    let mut jx = Vec::new();
    let mut jy = Vec::new();
    if let Some(r) = &range {
        jx.push((r.start, 1.0));
        jy.push((r.start + 1, 1.0));
    }
    (vec![0.0, 0.0], vec![JacobianRow { entries: jx }, JacobianRow { entries: jy }])
}

/// Distance: r = sqrt(dx^2 + dy^2) - d
fn distance_residual(
    sketch: &Sketch,
    param_map: &[(EntityId, Range<usize>)],
    a: EntityId,
    b: EntityId,
    target: f64,
) -> (Vec<f64>, Vec<JacobianRow>) {
    let Some((ax, ay)) = point_xy(sketch, a) else { return (vec![0.0], vec![JacobianRow { entries: vec![] }]); };
    let Some((bx, by)) = point_xy(sketch, b) else { return (vec![0.0], vec![JacobianRow { entries: vec![] }]); };

    let dx = ax - bx;
    let dy = ay - by;
    let dist = (dx * dx + dy * dy).sqrt();
    let r = dist - target;

    let inv_dist = if dist > 1e-12 { 1.0 / dist } else { 0.0 };
    let mut entries = Vec::new();
    if let Some(range) = param_range(param_map, a) {
        entries.push((range.start, dx * inv_dist));     // d/d(ax)
        entries.push((range.start + 1, dy * inv_dist)); // d/d(ay)
    }
    if let Some(range) = param_range(param_map, b) {
        entries.push((range.start, -dx * inv_dist));
        entries.push((range.start + 1, -dy * inv_dist));
    }
    (vec![r], vec![JacobianRow { entries }])
}

/// Perpendicular: r = dot(dir_a, dir_b) = 0
fn perpendicular_residual(
    sketch: &Sketch,
    _param_map: &[(EntityId, Range<usize>)],
    line_a: EntityId,
    line_b: EntityId,
) -> (Vec<f64>, Vec<JacobianRow>) {
    let dirs = (line_direction(sketch, line_a), line_direction(sketch, line_b));
    let (Some((dax, day)), Some((dbx, dby))) = dirs else {
        return (vec![0.0], vec![JacobianRow { entries: vec![] }]);
    };
    let r = dax * dbx + day * dby;
    // Jacobian is complex (depends on 4 points); provide empty for now.
    // The solver will still converge via numerical approximation fallback.
    (vec![r], vec![JacobianRow { entries: vec![] }])
}

/// Parallel: r = cross(dir_a, dir_b) = dax*dby - day*dbx = 0
fn parallel_residual(
    sketch: &Sketch,
    _param_map: &[(EntityId, Range<usize>)],
    line_a: EntityId,
    line_b: EntityId,
) -> (Vec<f64>, Vec<JacobianRow>) {
    let dirs = (line_direction(sketch, line_a), line_direction(sketch, line_b));
    let (Some((dax, day)), Some((dbx, dby))) = dirs else {
        return (vec![0.0], vec![JacobianRow { entries: vec![] }]);
    };
    let r = dax * dby - day * dbx;
    (vec![r], vec![JacobianRow { entries: vec![] }])
}

/// Radius: r = current_radius - target_value
fn radius_residual(
    sketch: &Sketch,
    param_map: &[(EntityId, Range<usize>)],
    entity: EntityId,
    target: f64,
) -> (Vec<f64>, Vec<JacobianRow>) {
    if let Some(SketchEntity::Circle(c)) = sketch.entity(entity) {
        let r = c.radius - target;
        let mut entries = Vec::new();
        if let Some(range) = param_range(param_map, entity) {
            entries.push((range.start, 1.0));
        }
        (vec![r], vec![JacobianRow { entries }])
    } else {
        (vec![0.0], vec![JacobianRow { entries: vec![] }])
    }
}

/// Midpoint: r = [2*m.x - a.x - b.x, 2*m.y - a.y - b.y]
fn midpoint_residual(
    sketch: &Sketch,
    param_map: &[(EntityId, Range<usize>)],
    point_id: EntityId,
    line_id: EntityId,
) -> (Vec<f64>, Vec<JacobianRow>) {
    let Some((mx, my)) = point_xy(sketch, point_id) else {
        return (vec![0.0; 2], vec![JacobianRow { entries: vec![] }; 2]);
    };
    let Some((s, e)) = line_endpoints(sketch, line_id) else {
        return (vec![0.0; 2], vec![JacobianRow { entries: vec![] }; 2]);
    };
    let Some((sx, sy)) = point_xy(sketch, s) else { return (vec![0.0; 2], vec![JacobianRow { entries: vec![] }; 2]); };
    let Some((ex, ey)) = point_xy(sketch, e) else { return (vec![0.0; 2], vec![JacobianRow { entries: vec![] }; 2]); };

    let rx = 2.0 * mx - sx - ex;
    let ry = 2.0 * my - sy - ey;

    let mut jx = Vec::new();
    let mut jy = Vec::new();
    if let Some(r) = param_range(param_map, point_id) { jx.push((r.start, 2.0)); jy.push((r.start + 1, 2.0)); }
    if let Some(r) = param_range(param_map, s) { jx.push((r.start, -1.0)); jy.push((r.start + 1, -1.0)); }
    if let Some(r) = param_range(param_map, e) { jx.push((r.start, -1.0)); jy.push((r.start + 1, -1.0)); }

    (vec![rx, ry], vec![JacobianRow { entries: jx }, JacobianRow { entries: jy }])
}

/// Direction vector of a line entity (un-normalized).
fn line_direction(sketch: &Sketch, line_id: EntityId) -> Option<(f64, f64)> {
    let (s, e) = line_endpoints(sketch, line_id)?;
    let (sx, sy) = point_xy(sketch, s)?;
    let (ex, ey) = point_xy(sketch, e)?;
    Some((ex - sx, ey - sy))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::sketch::entity::{SketchLine, SketchPoint};
    use crate::domain::sketch::sketch::{Sketch, SketchId};
    use crate::domain::sketch::work_plane::WorkPlane;

    fn two_point_sketch() -> Sketch {
        let mut sk = Sketch::new(SketchId(0), "test".into(), WorkPlane::xy());
        let p0 = sk.next_entity_id();
        let p1 = sk.next_entity_id();
        let ln = sk.next_entity_id();
        sk.add_entity(SketchEntity::Point(SketchPoint { id: p0, x: 0.0, y: 0.0, construction: false }));
        sk.add_entity(SketchEntity::Point(SketchPoint { id: p1, x: 3.0, y: 0.0, construction: false }));
        sk.add_entity(SketchEntity::Line(SketchLine { id: ln, start: p0, end: p1, construction: false }));
        sk
    }

    #[test]
    fn horizontal_residual_is_zero_for_horizontal_line() {
        let sk = two_point_sketch();
        let map = sk.parameter_map();
        let (r, _) = horizontal_residual(&sk, &map, EntityId(2));
        assert!((r[0]).abs() < 1e-12);
    }

    #[test]
    fn distance_residual_correct() {
        let sk = two_point_sketch();
        let map = sk.parameter_map();
        let (r, j) = distance_residual(&sk, &map, EntityId(0), EntityId(1), 3.0);
        assert!((r[0]).abs() < 1e-12); // distance is exactly 3.0
        assert!(!j[0].entries.is_empty());
    }
}
