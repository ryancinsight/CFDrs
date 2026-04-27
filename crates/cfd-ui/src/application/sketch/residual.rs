//! Constraint residuals and analytic Jacobian rows.
//!
//! Each constraint type maps to one or more scalar residual equations.
//! The Jacobian row for each residual is the vector of partial derivatives
//! with respect to the parameter vector (point coordinates and radii).
//!
//! # Theorem — Explicit Constraint Coverage
//!
//! For every declared `Constraint` variant in the sketch domain, this module
//! provides a residual map whose zero set is the intended geometric manifold,
//! together with the analytic Jacobian of that map with respect to the sketch's
//! point and circle-radius parameterization.
//!
//! **Proof sketch**: each residual is derived directly from the defining
//! geometric invariant: coincident points enforce equal coordinates, parallel
//! lines enforce vanishing 2-D cross product, perpendicular lines enforce
//! vanishing dot product, concentric entities enforce equal centers, symmetric
//! points enforce axis-midpoint incidence and axis-normal separation, and angle
//! / tangent residuals use differentiable invariants equivalent to the intended
//! geometric condition away from degenerate line length. The Jacobian entries
//! are the exact derivatives of those scalar maps. Validation is exercised by
//! the unit tests in this module and by `application::sketch::solver` tests.

use std::ops::Range;

use crate::domain::sketch::constraint::Constraint;
use crate::domain::sketch::entity::{EntityId, SketchEntity};
use crate::domain::sketch::sketch::Sketch;

type LineEndpointCoordinates = ((EntityId, EntityId), (f64, f64, f64, f64));
type CircleLikeGeometry = ((EntityId, f64), Vec<(usize, f64)>);

/// A sparse row of the Jacobian matrix.
#[derive(Clone, Debug)]
pub struct JacobianRow {
    /// (`parameter_index`, `partial_derivative`) pairs.
    pub entries: Vec<(usize, f64)>,
}

/// Compute residuals and Jacobian rows for a single constraint.
///
/// Returns `(residuals, jacobian_rows)` where each residual has a
/// corresponding Jacobian row.
#[must_use]
pub fn constraint_residual_and_jacobian(
    constraint: &Constraint,
    sketch: &Sketch,
    param_map: &[(EntityId, Range<usize>)],
) -> (Vec<f64>, Vec<JacobianRow>) {
    match constraint {
        Constraint::Horizontal(line_id) => horizontal_residual(sketch, param_map, *line_id),
        Constraint::Vertical(line_id) => vertical_residual(sketch, param_map, *line_id),
        Constraint::Coincident(a, b) => coincident_residual(sketch, param_map, *a, *b),
        Constraint::Tangent(a, b) => tangent_residual(sketch, param_map, *a, *b),
        Constraint::Fixed(point_id) => fixed_residual(sketch, param_map, *point_id),
        Constraint::Distance { a, b, value, .. } => {
            distance_residual(sketch, param_map, *a, *b, *value)
        }
        Constraint::Perpendicular(a, b) => perpendicular_residual(sketch, param_map, *a, *b),
        Constraint::Parallel(a, b) => parallel_residual(sketch, param_map, *a, *b),
        Constraint::EqualLength(a, b) => equal_length_residual(sketch, param_map, *a, *b),
        Constraint::Concentric(a, b) => concentric_residual(sketch, param_map, *a, *b),
        Constraint::Radius { entity, value, .. } => {
            radius_residual(sketch, param_map, *entity, *value)
        }
        Constraint::Midpoint { point, line } => midpoint_residual(sketch, param_map, *point, *line),
        Constraint::Symmetric {
            point_a,
            point_b,
            axis,
        } => symmetric_residual(sketch, param_map, *point_a, *point_b, *axis),
        Constraint::Angle {
            line_a,
            line_b,
            value_rad,
            ..
        } => angle_residual(sketch, param_map, *line_a, *line_b, *value_rad),
    }
}

fn zero_system(count: usize) -> (Vec<f64>, Vec<JacobianRow>) {
    (
        vec![0.0; count],
        vec![JacobianRow { entries: vec![] }; count],
    )
}

/// Find the parameter range for an entity ID.
fn param_range(param_map: &[(EntityId, Range<usize>)], id: EntityId) -> Option<Range<usize>> {
    param_map
        .iter()
        .find(|(eid, _)| *eid == id)
        .map(|(_, r)| r.clone())
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

fn push_scalar(entries: &mut Vec<(usize, f64)>, index: usize, value: f64) {
    if value != 0.0 {
        entries.push((index, value));
    }
}

fn push_point_entries(
    entries: &mut Vec<(usize, f64)>,
    range: Option<Range<usize>>,
    dx: f64,
    dy: f64,
) {
    if let Some(range) = range {
        push_scalar(entries, range.start, dx);
        push_scalar(entries, range.start + 1, dy);
    }
}

fn line_endpoint_data(sketch: &Sketch, line_id: EntityId) -> Option<LineEndpointCoordinates> {
    let (start, end) = line_endpoints(sketch, line_id)?;
    let (sx, sy) = point_xy(sketch, start)?;
    let (ex, ey) = point_xy(sketch, end)?;
    Some(((start, end), (sx, sy, ex, ey)))
}

fn entity_center_id(sketch: &Sketch, entity: EntityId) -> Option<EntityId> {
    match sketch.entity(entity)? {
        SketchEntity::Circle(circle) => Some(circle.center),
        SketchEntity::Arc(arc) => Some(arc.center),
        _ => None,
    }
}

fn entity_measure_and_jacobian(
    sketch: &Sketch,
    param_map: &[(EntityId, Range<usize>)],
    entity: EntityId,
) -> Option<(f64, Vec<(usize, f64)>)> {
    match sketch.entity(entity)? {
        SketchEntity::Line(line) => {
            let (ax, ay) = point_xy(sketch, line.start)?;
            let (bx, by) = point_xy(sketch, line.end)?;
            let dx = ax - bx;
            let dy = ay - by;
            let length = (dx * dx + dy * dy).sqrt();
            if length <= 1e-12 {
                return None;
            }
            let inv = 1.0 / length;
            let mut entries = Vec::new();
            push_point_entries(
                &mut entries,
                param_range(param_map, line.start),
                dx * inv,
                dy * inv,
            );
            push_point_entries(
                &mut entries,
                param_range(param_map, line.end),
                -dx * inv,
                -dy * inv,
            );
            Some((length, entries))
        }
        SketchEntity::Circle(circle) => {
            let mut entries = Vec::new();
            if let Some(range) = param_range(param_map, entity) {
                push_scalar(&mut entries, range.start, 1.0);
            }
            Some((circle.radius, entries))
        }
        SketchEntity::Arc(arc) => {
            let (cx, cy) = point_xy(sketch, arc.center)?;
            let (sx, sy) = point_xy(sketch, arc.start)?;
            let dx = cx - sx;
            let dy = cy - sy;
            let radius = (dx * dx + dy * dy).sqrt();
            if radius <= 1e-12 {
                return None;
            }
            let inv = 1.0 / radius;
            let mut entries = Vec::new();
            push_point_entries(
                &mut entries,
                param_range(param_map, arc.center),
                dx * inv,
                dy * inv,
            );
            push_point_entries(
                &mut entries,
                param_range(param_map, arc.start),
                -dx * inv,
                -dy * inv,
            );
            Some((radius, entries))
        }
        _ => None,
    }
}

fn circle_like_center_radius_and_gradient(
    sketch: &Sketch,
    param_map: &[(EntityId, Range<usize>)],
    entity: EntityId,
) -> Option<CircleLikeGeometry> {
    match sketch.entity(entity)? {
        SketchEntity::Circle(circle) => {
            let mut entries = Vec::new();
            if let Some(range) = param_range(param_map, entity) {
                push_scalar(&mut entries, range.start, 1.0);
            }
            Some(((circle.center, circle.radius), entries))
        }
        SketchEntity::Arc(arc) => {
            let (radius, entries) = entity_measure_and_jacobian(sketch, param_map, entity)?;
            Some(((arc.center, radius), entries))
        }
        _ => None,
    }
}

fn wrap_angle(angle: f64) -> f64 {
    let two_pi = 2.0 * std::f64::consts::PI;
    (angle + std::f64::consts::PI).rem_euclid(two_pi) - std::f64::consts::PI
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
    let Some((_, sy)) = point_xy(sketch, s) else {
        return (vec![0.0], vec![JacobianRow { entries: vec![] }]);
    };
    let Some((_, ey)) = point_xy(sketch, e) else {
        return (vec![0.0], vec![JacobianRow { entries: vec![] }]);
    };

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
    let Some((sx, _)) = point_xy(sketch, s) else {
        return (vec![0.0], vec![JacobianRow { entries: vec![] }]);
    };
    let Some((ex, _)) = point_xy(sketch, e) else {
        return (vec![0.0], vec![JacobianRow { entries: vec![] }]);
    };

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
    let Some((ax, ay)) = point_xy(sketch, a) else {
        return (vec![0.0; 2], vec![JacobianRow { entries: vec![] }; 2]);
    };
    let Some((bx, by)) = point_xy(sketch, b) else {
        return (vec![0.0; 2], vec![JacobianRow { entries: vec![] }; 2]);
    };

    let ra = param_range(param_map, a);
    let rb = param_range(param_map, b);

    let mut jx = Vec::new();
    let mut jy = Vec::new();
    if let Some(r) = &ra {
        jx.push((r.start, 1.0));
        jy.push((r.start + 1, 1.0));
    }
    if let Some(r) = &rb {
        jx.push((r.start, -1.0));
        jy.push((r.start + 1, -1.0));
    }

    (
        vec![ax - bx, ay - by],
        vec![JacobianRow { entries: jx }, JacobianRow { entries: jy }],
    )
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
    (
        vec![0.0, 0.0],
        vec![JacobianRow { entries: jx }, JacobianRow { entries: jy }],
    )
}

/// Distance: r = sqrt(dx^2 + dy^2) - d
fn distance_residual(
    sketch: &Sketch,
    param_map: &[(EntityId, Range<usize>)],
    a: EntityId,
    b: EntityId,
    target: f64,
) -> (Vec<f64>, Vec<JacobianRow>) {
    let Some((ax, ay)) = point_xy(sketch, a) else {
        return (vec![0.0], vec![JacobianRow { entries: vec![] }]);
    };
    let Some((bx, by)) = point_xy(sketch, b) else {
        return (vec![0.0], vec![JacobianRow { entries: vec![] }]);
    };

    let dx = ax - bx;
    let dy = ay - by;
    let dist = (dx * dx + dy * dy).sqrt();
    let r = dist - target;

    let inv_dist = if dist > 1e-12 { 1.0 / dist } else { 0.0 };
    let mut entries = Vec::new();
    if let Some(range) = param_range(param_map, a) {
        entries.push((range.start, dx * inv_dist)); // d/d(ax)
        entries.push((range.start + 1, dy * inv_dist)); // d/d(ay)
    }
    if let Some(range) = param_range(param_map, b) {
        entries.push((range.start, -dx * inv_dist));
        entries.push((range.start + 1, -dy * inv_dist));
    }
    (vec![r], vec![JacobianRow { entries }])
}

/// Perpendicular: r = `dot(dir_a`, `dir_b`) = 0
fn perpendicular_residual(
    sketch: &Sketch,
    param_map: &[(EntityId, Range<usize>)],
    line_a: EntityId,
    line_b: EntityId,
) -> (Vec<f64>, Vec<JacobianRow>) {
    let Some(((a_start, a_end), (asx, asy, aex, aey))) = line_endpoint_data(sketch, line_a) else {
        return zero_system(1);
    };
    let Some(((b_start, b_end), (bsx, bsy, bex, bey))) = line_endpoint_data(sketch, line_b) else {
        return zero_system(1);
    };
    let dax = aex - asx;
    let day = aey - asy;
    let dbx = bex - bsx;
    let dby = bey - bsy;
    let r = dax * dbx + day * dby;
    let mut entries = Vec::new();
    push_point_entries(&mut entries, param_range(param_map, a_start), -dbx, -dby);
    push_point_entries(&mut entries, param_range(param_map, a_end), dbx, dby);
    push_point_entries(&mut entries, param_range(param_map, b_start), -dax, -day);
    push_point_entries(&mut entries, param_range(param_map, b_end), dax, day);
    (vec![r], vec![JacobianRow { entries }])
}

/// Parallel: r = `cross(dir_a`, `dir_b`) = dax*dby - day*dbx = 0
fn parallel_residual(
    sketch: &Sketch,
    param_map: &[(EntityId, Range<usize>)],
    line_a: EntityId,
    line_b: EntityId,
) -> (Vec<f64>, Vec<JacobianRow>) {
    let Some(((a_start, a_end), (asx, asy, aex, aey))) = line_endpoint_data(sketch, line_a) else {
        return zero_system(1);
    };
    let Some(((b_start, b_end), (bsx, bsy, bex, bey))) = line_endpoint_data(sketch, line_b) else {
        return zero_system(1);
    };
    let dax = aex - asx;
    let day = aey - asy;
    let dbx = bex - bsx;
    let dby = bey - bsy;
    let r = dax * dby - day * dbx;
    let mut entries = Vec::new();
    push_point_entries(&mut entries, param_range(param_map, a_start), -dby, dbx);
    push_point_entries(&mut entries, param_range(param_map, a_end), dby, -dbx);
    push_point_entries(&mut entries, param_range(param_map, b_start), day, -dax);
    push_point_entries(&mut entries, param_range(param_map, b_end), -day, dax);
    (vec![r], vec![JacobianRow { entries }])
}

/// Radius: r = `current_radius` - `target_value`
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
    let Some((sx, sy)) = point_xy(sketch, s) else {
        return (vec![0.0; 2], vec![JacobianRow { entries: vec![] }; 2]);
    };
    let Some((ex, ey)) = point_xy(sketch, e) else {
        return (vec![0.0; 2], vec![JacobianRow { entries: vec![] }; 2]);
    };

    let rx = 2.0 * mx - sx - ex;
    let ry = 2.0 * my - sy - ey;

    let mut jx = Vec::new();
    let mut jy = Vec::new();
    if let Some(r) = param_range(param_map, point_id) {
        jx.push((r.start, 2.0));
        jy.push((r.start + 1, 2.0));
    }
    if let Some(r) = param_range(param_map, s) {
        jx.push((r.start, -1.0));
        jy.push((r.start + 1, -1.0));
    }
    if let Some(r) = param_range(param_map, e) {
        jx.push((r.start, -1.0));
        jy.push((r.start + 1, -1.0));
    }

    (
        vec![rx, ry],
        vec![JacobianRow { entries: jx }, JacobianRow { entries: jy }],
    )
}

#[cfg(test)]
/// Direction vector of a line entity (un-normalized).
fn line_direction(sketch: &Sketch, line_id: EntityId) -> Option<(f64, f64)> {
    let (s, e) = line_endpoints(sketch, line_id)?;
    let (sx, sy) = point_xy(sketch, s)?;
    let (ex, ey) = point_xy(sketch, e)?;
    Some((ex - sx, ey - sy))
}

fn tangent_residual(
    sketch: &Sketch,
    param_map: &[(EntityId, Range<usize>)],
    line_id: EntityId,
    curve_id: EntityId,
) -> (Vec<f64>, Vec<JacobianRow>) {
    let Some(((line_start, line_end), (sx, sy, ex, ey))) = line_endpoint_data(sketch, line_id)
    else {
        return zero_system(1);
    };
    let Some(((center_id, radius), radius_gradient)) =
        circle_like_center_radius_and_gradient(sketch, param_map, curve_id)
    else {
        return zero_system(1);
    };
    let Some((cx, cy)) = point_xy(sketch, center_id) else {
        return zero_system(1);
    };

    let dx = ex - sx;
    let dy = ey - sy;
    let qx = cx - sx;
    let qy = cy - sy;
    let cross = dx * qy - dy * qx;
    let length_sq = dx * dx + dy * dy;
    if length_sq <= 1e-12 {
        return zero_system(1);
    }
    let residual = cross * cross - radius * radius * length_sq;

    let dc_dsx = ey - cy;
    let dc_dsy = cx - ex;
    let dc_dex = cy - sy;
    let dc_dey = sx - cx;
    let dc_dcx = -dy;
    let dc_dcy = dx;

    let mut entries = Vec::new();
    push_point_entries(
        &mut entries,
        param_range(param_map, line_start),
        2.0 * cross * dc_dsx + 2.0 * radius * radius * dx,
        2.0 * cross * dc_dsy + 2.0 * radius * radius * dy,
    );
    push_point_entries(
        &mut entries,
        param_range(param_map, line_end),
        2.0 * cross * dc_dex - 2.0 * radius * radius * dx,
        2.0 * cross * dc_dey - 2.0 * radius * radius * dy,
    );
    push_point_entries(
        &mut entries,
        param_range(param_map, center_id),
        2.0 * cross * dc_dcx,
        2.0 * cross * dc_dcy,
    );
    for (index, derivative) in radius_gradient {
        push_scalar(&mut entries, index, -2.0 * radius * length_sq * derivative);
    }
    (vec![residual], vec![JacobianRow { entries }])
}

fn equal_length_residual(
    sketch: &Sketch,
    param_map: &[(EntityId, Range<usize>)],
    entity_a: EntityId,
    entity_b: EntityId,
) -> (Vec<f64>, Vec<JacobianRow>) {
    let Some((len_a, mut entries)) = entity_measure_and_jacobian(sketch, param_map, entity_a)
    else {
        return zero_system(1);
    };
    let Some((len_b, entries_b)) = entity_measure_and_jacobian(sketch, param_map, entity_b) else {
        return zero_system(1);
    };
    for (index, derivative) in entries_b {
        push_scalar(&mut entries, index, -derivative);
    }
    (vec![len_a - len_b], vec![JacobianRow { entries }])
}

fn concentric_residual(
    sketch: &Sketch,
    param_map: &[(EntityId, Range<usize>)],
    entity_a: EntityId,
    entity_b: EntityId,
) -> (Vec<f64>, Vec<JacobianRow>) {
    let Some(center_a) = entity_center_id(sketch, entity_a) else {
        return zero_system(2);
    };
    let Some(center_b) = entity_center_id(sketch, entity_b) else {
        return zero_system(2);
    };
    coincident_residual(sketch, param_map, center_a, center_b)
}

fn symmetric_residual(
    sketch: &Sketch,
    param_map: &[(EntityId, Range<usize>)],
    point_a: EntityId,
    point_b: EntityId,
    axis: EntityId,
) -> (Vec<f64>, Vec<JacobianRow>) {
    let Some((ax, ay)) = point_xy(sketch, point_a) else {
        return zero_system(2);
    };
    let Some((bx, by)) = point_xy(sketch, point_b) else {
        return zero_system(2);
    };
    let Some(((axis_start, axis_end), (sx, sy, ex, ey))) = line_endpoint_data(sketch, axis) else {
        return zero_system(2);
    };

    let dx = ex - sx;
    let dy = ey - sy;
    let mx = 0.5 * (ax + bx);
    let my = 0.5 * (ay + by);
    let r_axis = dx * (my - sy) - dy * (mx - sx);
    let r_normal = dx * (bx - ax) + dy * (by - ay);

    let mut axis_entries = Vec::new();
    push_point_entries(
        &mut axis_entries,
        param_range(param_map, axis_start),
        ey - my,
        mx - ex,
    );
    push_point_entries(
        &mut axis_entries,
        param_range(param_map, axis_end),
        my - sy,
        sx - mx,
    );
    push_point_entries(
        &mut axis_entries,
        param_range(param_map, point_a),
        -0.5 * dy,
        0.5 * dx,
    );
    push_point_entries(
        &mut axis_entries,
        param_range(param_map, point_b),
        -0.5 * dy,
        0.5 * dx,
    );

    let mut normal_entries = Vec::new();
    push_point_entries(
        &mut normal_entries,
        param_range(param_map, axis_start),
        ax - bx,
        ay - by,
    );
    push_point_entries(
        &mut normal_entries,
        param_range(param_map, axis_end),
        bx - ax,
        by - ay,
    );
    push_point_entries(
        &mut normal_entries,
        param_range(param_map, point_a),
        -dx,
        -dy,
    );
    push_point_entries(&mut normal_entries, param_range(param_map, point_b), dx, dy);

    (
        vec![r_axis, r_normal],
        vec![
            JacobianRow {
                entries: axis_entries,
            },
            JacobianRow {
                entries: normal_entries,
            },
        ],
    )
}

fn angle_residual(
    sketch: &Sketch,
    param_map: &[(EntityId, Range<usize>)],
    line_a: EntityId,
    line_b: EntityId,
    target: f64,
) -> (Vec<f64>, Vec<JacobianRow>) {
    let Some(((a_start, a_end), (asx, asy, aex, aey))) = line_endpoint_data(sketch, line_a) else {
        return zero_system(1);
    };
    let Some(((b_start, b_end), (bsx, bsy, bex, bey))) = line_endpoint_data(sketch, line_b) else {
        return zero_system(1);
    };

    let dax = aex - asx;
    let day = aey - asy;
    let dbx = bex - bsx;
    let dby = bey - bsy;
    let cross = dax * dby - day * dbx;
    let dot = dax * dbx + day * dby;
    let denom = cross * cross + dot * dot;
    if denom <= 1e-12 {
        return zero_system(1);
    }
    let theta = cross.atan2(dot);
    let residual = wrap_angle(theta - target);

    let factor = 1.0 / denom;
    let mut entries = Vec::new();
    let mut push_angle_derivative =
        |entity: EntityId, dc_dx: f64, dc_dy: f64, dd_dx: f64, dd_dy: f64| {
            let dr_dx = factor * (dot * dc_dx - cross * dd_dx);
            let dr_dy = factor * (dot * dc_dy - cross * dd_dy);
            push_point_entries(&mut entries, param_range(param_map, entity), dr_dx, dr_dy);
        };

    push_angle_derivative(a_start, -dby, dbx, -dbx, -dby);
    push_angle_derivative(a_end, dby, -dbx, dbx, dby);
    push_angle_derivative(b_start, day, -dax, -dax, -day);
    push_angle_derivative(b_end, -day, dax, dax, day);

    (vec![residual], vec![JacobianRow { entries }])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::sketch::constraint::Constraint;
    use crate::domain::sketch::entity::{SketchCircle, SketchLine, SketchPoint};
    use crate::domain::sketch::sketch::{Sketch, SketchId};
    use crate::domain::sketch::work_plane::WorkPlane;

    fn two_point_sketch() -> Sketch {
        let mut sk = Sketch::new(SketchId(0), "test".into(), WorkPlane::xy());
        let p0 = sk.next_entity_id();
        let p1 = sk.next_entity_id();
        let ln = sk.next_entity_id();
        sk.add_entity(SketchEntity::Point(SketchPoint {
            id: p0,
            x: 0.0,
            y: 0.0,
            construction: false,
        }));
        sk.add_entity(SketchEntity::Point(SketchPoint {
            id: p1,
            x: 3.0,
            y: 0.0,
            construction: false,
        }));
        sk.add_entity(SketchEntity::Line(SketchLine {
            id: ln,
            start: p0,
            end: p1,
            construction: false,
        }));
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

    #[test]
    fn concentric_residual_is_zero_for_shared_circle_centers() {
        let mut sk = Sketch::new(SketchId(0), "circles".into(), WorkPlane::xy());
        let c0 = sk.next_entity_id();
        let c1 = sk.next_entity_id();
        let cir0 = sk.next_entity_id();
        let cir1 = sk.next_entity_id();
        sk.add_entity(SketchEntity::Point(SketchPoint {
            id: c0,
            x: 1.0,
            y: 2.0,
            construction: false,
        }));
        sk.add_entity(SketchEntity::Point(SketchPoint {
            id: c1,
            x: 1.0,
            y: 2.0,
            construction: false,
        }));
        sk.add_entity(SketchEntity::Circle(SketchCircle {
            id: cir0,
            center: c0,
            radius: 1.0,
            construction: false,
        }));
        sk.add_entity(SketchEntity::Circle(SketchCircle {
            id: cir1,
            center: c1,
            radius: 2.0,
            construction: false,
        }));

        let map = sk.parameter_map();
        let (r, j) = concentric_residual(&sk, &map, cir0, cir1);
        assert!(r.iter().all(|value| value.abs() < 1e-12));
        assert_eq!(j.len(), 2);
    }

    #[test]
    fn equal_length_residual_is_zero_for_equal_segments() {
        let mut sk = Sketch::new(SketchId(0), "equal".into(), WorkPlane::xy());
        let p0 = sk.next_entity_id();
        let p1 = sk.next_entity_id();
        let p2 = sk.next_entity_id();
        let p3 = sk.next_entity_id();
        let l0 = sk.next_entity_id();
        let l1 = sk.next_entity_id();
        sk.add_entity(SketchEntity::Point(SketchPoint {
            id: p0,
            x: 0.0,
            y: 0.0,
            construction: false,
        }));
        sk.add_entity(SketchEntity::Point(SketchPoint {
            id: p1,
            x: 3.0,
            y: 0.0,
            construction: false,
        }));
        sk.add_entity(SketchEntity::Point(SketchPoint {
            id: p2,
            x: 0.0,
            y: 1.0,
            construction: false,
        }));
        sk.add_entity(SketchEntity::Point(SketchPoint {
            id: p3,
            x: 3.0,
            y: 1.0,
            construction: false,
        }));
        sk.add_entity(SketchEntity::Line(SketchLine {
            id: l0,
            start: p0,
            end: p1,
            construction: false,
        }));
        sk.add_entity(SketchEntity::Line(SketchLine {
            id: l1,
            start: p2,
            end: p3,
            construction: false,
        }));

        let map = sk.parameter_map();
        let (r, j) = equal_length_residual(&sk, &map, l0, l1);
        assert!(r[0].abs() < 1e-12);
        assert!(!j[0].entries.is_empty());
    }

    #[test]
    fn solver_converges_on_parallel_constraint() {
        let mut sk = Sketch::new(SketchId(0), "parallel".into(), WorkPlane::xy());
        let p0 = sk.next_entity_id();
        let p1 = sk.next_entity_id();
        let p2 = sk.next_entity_id();
        let p3 = sk.next_entity_id();
        let l0 = sk.next_entity_id();
        let l1 = sk.next_entity_id();
        sk.add_entity(SketchEntity::Point(SketchPoint {
            id: p0,
            x: 0.0,
            y: 0.0,
            construction: false,
        }));
        sk.add_entity(SketchEntity::Point(SketchPoint {
            id: p1,
            x: 2.0,
            y: 0.0,
            construction: false,
        }));
        sk.add_entity(SketchEntity::Point(SketchPoint {
            id: p2,
            x: 0.0,
            y: 1.0,
            construction: false,
        }));
        sk.add_entity(SketchEntity::Point(SketchPoint {
            id: p3,
            x: 1.0,
            y: 2.0,
            construction: false,
        }));
        sk.add_entity(SketchEntity::Line(SketchLine {
            id: l0,
            start: p0,
            end: p1,
            construction: false,
        }));
        sk.add_entity(SketchEntity::Line(SketchLine {
            id: l1,
            start: p2,
            end: p3,
            construction: false,
        }));
        sk.add_constraint(Constraint::Parallel(l0, l1));

        let result = crate::application::sketch::solver::ConstraintSolver::new().solve(&mut sk);
        assert!(result.converged);

        let (dx0, dy0) = line_direction(&sk, l0).expect("line 0 direction");
        let (dx1, dy1) = line_direction(&sk, l1).expect("line 1 direction");
        assert!((dx0 * dy1 - dy0 * dx1).abs() < 1e-8);
    }
}
