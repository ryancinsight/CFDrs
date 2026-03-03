//! Boundary seam repair for CSG boolean operations.
//!
//! After fragment classification, the result face set may contain boundary
//! edges where the two surface CDTs didn't produce matching vertices along
//! the intersection curve.  This module provides two complementary repairs:
//!
//! 1. **T-junction snap-round** ([`snap_round_tjunctions`]) — splits face
//!    edges at boundary vertices that satisfy exact on-segment predicates
//!    with boundary-constrained insertion (and tolerance fallback). This
//!    adds faces (splits) without removing any, avoiding topology cascades.
//!
//! 2. **Boundary loop filling** ([`fill_boundary_loops`]) — fills remaining
//!    closed boundary loops using constrained Delaunay triangulation (CDT)
//!    over an exact-predicate 2D projection, with ear-clipping fallback.
//!
//! ## References
//!
//! - Attene (2010) "A lightweight approach to repairing digitized polygon meshes"
//! - Liepa (2003) "Filling Holes in Meshes"
//! - Held (2001) "FIST: Fast Industrial-Strength Triangulation of Polygons"
//! - Shewchuk (1996) "Triangle: Engineering a 2D Quality Mesh Generator..."

use std::collections::{HashMap, HashSet};

use crate::application::csg::clip::polygon2d::geometry::point_in_polygon;
use crate::application::delaunay::{Cdt, Pslg};
use crate::domain::core::index::VertexId;
use crate::domain::core::scalar::{Point3r, Real};
use crate::domain::geometry::predicates::{orient_2d_arr, Orientation};
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;

/// Maximum boundary loop size to attempt filling.
const MAX_LOOP: usize = 512;

/// Maximum iterations for the fill loop.
const MAX_ITERS: usize = 8;

/// Canonical edge key: (min(a,b), max(a,b)).
fn canon_edge(a: VertexId, b: VertexId) -> (VertexId, VertexId) {
    if a < b {
        (a, b)
    } else {
        (b, a)
    }
}

/// Build canonical edge valence map from a face set.
///
/// Returns a map from canonical edge `(min_id, max_id)` to the number of
/// faces sharing that edge.  Manifold edges have valence 2; boundary edges
/// have valence 1; non-manifold edges have valence > 2.
pub(crate) fn build_canonical_valence(faces: &[FaceData]) -> HashMap<(VertexId, VertexId), u32> {
    let mut valence: HashMap<(VertexId, VertexId), u32> = HashMap::new();
    for face in faces {
        let v = face.vertices;
        for i in 0..3 {
            let j = (i + 1) % 3;
            *valence.entry(canon_edge(v[i], v[j])).or_insert(0) += 1;
        }
    }
    valence
}

/// Check if adding triangle (va, vb, vc) would create a non-manifold edge.
///
/// Returns `true` if any of the triangle's 3 canonical edges already has
/// valence >= 2 in the current mesh (meaning it's already manifold or
/// non-manifold — adding another face would make it worse).
fn would_create_nm(
    va: VertexId,
    vb: VertexId,
    vc: VertexId,
    valence: &HashMap<(VertexId, VertexId), u32>,
) -> bool {
    let edges = [canon_edge(va, vb), canon_edge(vb, vc), canon_edge(vc, va)];
    edges
        .iter()
        .any(|e| valence.get(e).copied().unwrap_or(0) >= 2)
}

/// Update the canonical edge valence map after adding a triangle.
fn record_triangle(
    va: VertexId,
    vb: VertexId,
    vc: VertexId,
    valence: &mut HashMap<(VertexId, VertexId), u32>,
) {
    for &e in &[canon_edge(va, vb), canon_edge(vb, vc), canon_edge(vc, va)] {
        *valence.entry(e).or_insert(0) += 1;
    }
}

/// Fill all boundary loops in the face set by ear-clipping triangulation.
///
/// Iterates until no boundary loops remain or no progress is made.
/// Builds a canonical edge valence map to prevent creating non-manifold
/// edges during fill (triangles whose diagonal edges coincide with
/// existing interior mesh edges are skipped).
pub(crate) fn fill_boundary_loops(faces: &mut Vec<FaceData>, pool: &VertexPool) {
    for _iter in 0..MAX_ITERS {
        let boundary = build_boundary_edges(faces);
        if boundary.is_empty() {
            break;
        }

        let loops = trace_loops(&boundary);
        if loops.is_empty() {
            break;
        }

        // Build canonical edge valence map from current face set.
        let mut valence = build_canonical_valence(faces);

        let mut added = 0_usize;
        for poly in &loops {
            // Modern robust path: CDT first (exact predicates), ear-clip fallback.
            let cdt_added = cdt_fill_loop(poly, pool, faces, &mut valence);
            added += if cdt_added > 0 {
                cdt_added
            } else {
                ear_clip_fill(poly, pool, faces, &mut valence)
            };
        }

        if added == 0 {
            break;
        }

        // Remove duplicate faces after adding fill triangles.
        dedup_faces(faces);

        #[cfg(test)]
        eprintln!(
            "[stitch {}] {} bnd edges, {} loops, {} fill tris",
            _iter,
            boundary.len(),
            loops.len(),
            added,
        );
    }
}

/// Build directed boundary half-edges: edges appearing exactly once with no
/// reverse partner.
fn build_boundary_edges(faces: &[FaceData]) -> Vec<(VertexId, VertexId)> {
    let mut he: HashMap<(VertexId, VertexId), u32> = HashMap::new();
    for face in faces {
        let v = face.vertices;
        for i in 0..3 {
            let j = (i + 1) % 3;
            *he.entry((v[i], v[j])).or_insert(0) += 1;
        }
    }
    he.iter()
        .filter(|&(&(vi, vj), &c)| vi != vj && c == 1 && !he.contains_key(&(vj, vi)))
        .map(|(&e, _)| e)
        .collect()
}

/// Trace closed boundary loops from directed boundary edges using greedy DFS
/// with inner-cycle extraction for figure-8 topologies.
fn trace_loops(boundary: &[(VertexId, VertexId)]) -> Vec<Vec<VertexId>> {
    let mut adj: HashMap<VertexId, Vec<VertexId>> = HashMap::new();
    for &(vi, vj) in boundary {
        adj.entry(vi).or_default().push(vj);
    }
    for v in adj.values_mut() {
        v.sort();
    }

    let mut used: HashSet<(VertexId, VertexId)> = HashSet::new();
    let mut loops: Vec<Vec<VertexId>> = Vec::new();
    let mut starts: Vec<VertexId> = adj.keys().copied().collect();
    starts.sort();

    for start in starts {
        let nexts = match adj.get(&start) {
            Some(s) => s.clone(),
            None => continue,
        };
        for first_next in nexts {
            if used.contains(&(start, first_next)) {
                continue;
            }
            let mut path: Vec<VertexId> = vec![start, first_next];
            used.insert((start, first_next));
            let mut cur = first_next;
            let mut closed = false;

            loop {
                if path.len() > MAX_LOOP {
                    break;
                }
                let succs = match adj.get(&cur) {
                    Some(s) => s,
                    None => break,
                };
                let mut found = false;
                for &n in succs {
                    if used.contains(&(cur, n)) {
                        continue;
                    }
                    used.insert((cur, n));
                    if n == start {
                        closed = true;
                        found = true;
                        break;
                    }
                    // Inner-cycle extraction for figure-8 boundaries.
                    if let Some(pos) = path.iter().position(|&v| v == n) {
                        let inner = path[pos..].to_vec();
                        if inner.len() >= 3 && inner.len() <= MAX_LOOP {
                            loops.push(inner);
                        }
                        path.truncate(pos + 1);
                        cur = n;
                        found = true;
                        break;
                    }
                    path.push(n);
                    cur = n;
                    found = true;
                    break;
                }
                if !found || closed {
                    break;
                }
            }

            if closed && path.len() >= 3 && path.len() <= MAX_LOOP {
                loops.push(path);
            }
        }
    }
    loops
}

/// Fill a single boundary loop using constrained Delaunay triangulation (CDT).
///
/// The 3D boundary polygon is projected to 2D, triangulated with PSLG/CDT,
/// then interior triangles are mapped back to original 3D vertices.
/// This uses exact predicates in the Delaunay pipeline and avoids heuristic
/// ear selection for most loops.
///
/// Returns the number of triangles added, or 0 if CDT is unavailable for
/// this loop (caller should fall back to ear clipping).
pub(crate) fn cdt_fill_loop(
    poly: &[VertexId],
    pool: &VertexPool,
    out: &mut Vec<FaceData>,
    valence: &mut HashMap<(VertexId, VertexId), u32>,
) -> usize {
    let n = poly.len();
    if n < 3 || n > MAX_LOOP {
        return 0;
    }
    if n == 3 {
        let (va, vb, vc) = (poly[0], poly[1], poly[2]);
        let pa = pool.position(va);
        let pb = pool.position(vb);
        let pc = pool.position(vc);
        if (pb - pa).cross(&(pc - pa)).norm_squared() > 1e-30
            && !would_create_nm(va, vb, vc, valence)
        {
            out.push(FaceData::untagged(va, vb, vc));
            record_triangle(va, vb, vc, valence);
            return 1;
        }
        return 0;
    }

    // Compute polygon normal using Newell's method.
    let positions: Vec<_> = poly.iter().map(|&v| *pool.position(v)).collect();
    let (nx, ny, nz) = newell_normal(&positions);
    let nlen_sq = nx * nx + ny * ny + nz * nz;
    if nlen_sq < 1e-30 {
        return 0;
    }
    let inv_nlen = 1.0 / nlen_sq.sqrt();
    let normal = nalgebra::Vector3::new(nx * inv_nlen, ny * inv_nlen, nz * inv_nlen);

    // Build orthogonal 2D basis from normal (Gram-Schmidt).
    let seed = if normal.x.abs() < 0.9 {
        nalgebra::Vector3::new(1.0, 0.0, 0.0)
    } else {
        nalgebra::Vector3::new(0.0, 1.0, 0.0)
    };
    let u_axis = (seed - normal * seed.dot(&normal)).normalize();
    let v_axis = normal.cross(&u_axis);

    // Centroid for projection origin.
    let inv_n = 1.0 / n as Real;
    let cx: Real = positions.iter().map(|p| p.x).sum::<Real>() * inv_n;
    let cy: Real = positions.iter().map(|p| p.y).sum::<Real>() * inv_n;
    let cz: Real = positions.iter().map(|p| p.z).sum::<Real>() * inv_n;
    let centroid = nalgebra::Vector3::new(cx, cy, cz);

    // Project loop vertices to 2D.
    let pts2d: Vec<[Real; 2]> = positions
        .iter()
        .map(|p| {
            let d = p.coords - centroid;
            [d.dot(&u_axis), d.dot(&v_axis)]
        })
        .collect();
    let poly_ccw = match polygon_winding_ccw(&pts2d) {
        Some(w) => w,
        None => return 0,
    };

    // Build PSLG from the projected boundary loop.
    let mut pslg = Pslg::with_capacity(n, n);
    let mut pids = Vec::with_capacity(n);
    for p in &pts2d {
        pids.push(pslg.add_vertex(p[0], p[1]));
    }
    for i in 0..n {
        let j = (i + 1) % n;
        pslg.add_segment(pids[i], pids[j]);
    }

    let cdt = match Cdt::try_from_pslg(&pslg) {
        Ok(cdt) => cdt,
        Err(_) => return 0,
    };

    // Collect candidate triangles from CDT interior cells whose centroid lies
    // inside the boundary polygon (filters convex-hull exterior triangles for
    // concave loops).
    let mut candidates: Vec<[VertexId; 3]> = Vec::new();
    for (_, tri) in cdt.triangulation().interior_triangles() {
        let [a_id, b_id, c_id] = tri.vertices;
        let (ai, bi, ci) = (a_id.idx(), b_id.idx(), c_id.idx());
        if ai >= n || bi >= n || ci >= n {
            // Steiner insertion would require adding new vertices to the pool,
            // which this repair path does not do. Fall back to ear clipping.
            return 0;
        }
        let a2 = pts2d[ai];
        let b2 = pts2d[bi];
        let c2 = pts2d[ci];
        let tcx = (a2[0] + b2[0] + c2[0]) / 3.0;
        let tcy = (a2[1] + b2[1] + c2[1]) / 3.0;
        if !point_in_polygon(tcx, tcy, &pts2d) {
            continue;
        }

        let ori = orient_2d_arr(a2, b2, c2);
        if ori == Orientation::Degenerate {
            continue;
        }
        let tri_ccw = ori == Orientation::Positive;
        let tri_verts = if tri_ccw == poly_ccw {
            [poly[ai], poly[bi], poly[ci]]
        } else {
            [poly[ai], poly[ci], poly[bi]]
        };
        candidates.push(tri_verts);
    }

    if candidates.is_empty() {
        return 0;
    }

    let mut added = 0_usize;
    for [va, vb, vc] in candidates {
        let pa = pool.position(va);
        let pb = pool.position(vb);
        let pc = pool.position(vc);
        if (pb - pa).cross(&(pc - pa)).norm_squared() > 1e-30
            && !would_create_nm(va, vb, vc, valence)
        {
            out.push(FaceData::untagged(va, vb, vc));
            record_triangle(va, vb, vc, valence);
            added += 1;
        }
    }
    added
}

/// Fill a single boundary loop using ear-clipping triangulation.
///
/// Projects the 3D polygon to 2D using the Newell normal, determines winding
/// order, then iteratively clips ears (convex vertices with no interior points)
/// to produce a triangulation that fills the boundary hole.
///
/// The `valence` map tracks canonical edge valence across the entire face set.
/// Before adding each triangle, the function checks that none of its edges
/// would create a non-manifold canonical edge (valence > 2).  Ears whose
/// diagonal edges coincide with existing interior mesh edges are skipped.
///
/// Returns the number of triangles added.
pub(crate) fn ear_clip_fill(
    poly: &[VertexId],
    pool: &VertexPool,
    out: &mut Vec<FaceData>,
    valence: &mut HashMap<(VertexId, VertexId), u32>,
) -> usize {
    let n = poly.len();
    if n < 3 {
        return 0;
    }
    if n == 3 {
        let (va, vb, vc) = (poly[0], poly[1], poly[2]);
        let pa = pool.position(va);
        let pb = pool.position(vb);
        let pc = pool.position(vc);
        if (pb - pa).cross(&(pc - pa)).norm_squared() > 1e-30
            && !would_create_nm(va, vb, vc, valence)
        {
            out.push(FaceData::untagged(va, vb, vc));
            record_triangle(va, vb, vc, valence);
            return 1;
        }
        return 0;
    }

    // Compute polygon normal using Newell's method.
    let positions: Vec<_> = poly.iter().map(|&v| *pool.position(v)).collect();
    let (nx, ny, nz) = newell_normal(&positions);
    let nlen_sq = nx * nx + ny * ny + nz * nz;
    if nlen_sq < 1e-30 {
        return 0; // Degenerate (collinear) polygon — cannot triangulate.
    }
    let inv_nlen = 1.0 / nlen_sq.sqrt();
    let normal = nalgebra::Vector3::new(nx * inv_nlen, ny * inv_nlen, nz * inv_nlen);

    // Build orthogonal 2D basis from normal (Gram-Schmidt).
    let seed = if normal.x.abs() < 0.9 {
        nalgebra::Vector3::new(1.0, 0.0, 0.0)
    } else {
        nalgebra::Vector3::new(0.0, 1.0, 0.0)
    };
    let u_axis = (seed - normal * seed.dot(&normal)).normalize();
    let v_axis = normal.cross(&u_axis);

    // Centroid for projection origin.
    let inv_n = 1.0 / n as Real;
    let cx: Real = positions.iter().map(|p| p.x).sum::<Real>() * inv_n;
    let cy: Real = positions.iter().map(|p| p.y).sum::<Real>() * inv_n;
    let cz: Real = positions.iter().map(|p| p.z).sum::<Real>() * inv_n;
    let centroid = nalgebra::Vector3::new(cx, cy, cz);

    // Project to 2D.
    let pts2d: Vec<[Real; 2]> = positions
        .iter()
        .map(|p| {
            let d = p.coords - centroid;
            [d.dot(&u_axis), d.dot(&v_axis)]
        })
        .collect();

    // Determine polygon orientation with exact predicates.
    let ccw = match polygon_winding_ccw(&pts2d) {
        Some(w) => w,
        None => return 0,
    };

    // Ear-clipping.
    let mut indices: Vec<usize> = (0..n).collect();
    let mut count = 0_usize;
    let mut guard = 0_usize;
    let max_guard = n * n * 2;

    while indices.len() > 3 && guard < max_guard {
        let m = indices.len();
        let mut found_ear = false;

        for i in 0..m {
            let prev = if i == 0 { m - 1 } else { i - 1 };
            let next = (i + 1) % m;

            let a = indices[prev];
            let b = indices[i];
            let c = indices[next];

            // Check convexity using exact orientation predicates.
            let turn = orient_2d_arr(pts2d[a], pts2d[b], pts2d[c]);
            let convex = if ccw {
                turn == Orientation::Positive
            } else {
                turn == Orientation::Negative
            };
            if !convex {
                guard += 1;
                continue;
            }

            // Check no other polygon vertex lies inside this ear triangle.
            let mut ear_ok = true;
            for j in 0..m {
                if j == prev || j == i || j == next {
                    continue;
                }
                if point_in_triangle(&pts2d[indices[j]], &pts2d[a], &pts2d[b], &pts2d[c]) {
                    ear_ok = false;
                    break;
                }
            }

            if ear_ok {
                let (va, vb, vc) = (poly[a], poly[b], poly[c]);
                let pa = pool.position(va);
                let pb = pool.position(vb);
                let pc = pool.position(vc);
                if (pb - pa).cross(&(pc - pa)).norm_squared() > 1e-30
                    && !would_create_nm(va, vb, vc, valence)
                {
                    out.push(FaceData::untagged(va, vb, vc));
                    record_triangle(va, vb, vc, valence);
                    count += 1;
                }
                indices.remove(i);
                found_ear = true;
                break;
            }
            guard += 1;
        }

        if !found_ear {
            // No ear found — polygon may be self-intersecting in projection.
            // Fall back to fan triangulation from vertex 0 for remaining polygon.
            for i in 1..indices.len() - 1 {
                let a = indices[0];
                let b = indices[i];
                let c = indices[i + 1];
                let (va, vb, vc) = (poly[a], poly[b], poly[c]);
                let pa = pool.position(va);
                let pb = pool.position(vb);
                let pc = pool.position(vc);
                if (pb - pa).cross(&(pc - pa)).norm_squared() > 1e-30
                    && !would_create_nm(va, vb, vc, valence)
                {
                    out.push(FaceData::untagged(va, vb, vc));
                    record_triangle(va, vb, vc, valence);
                    count += 1;
                }
            }
            break;
        }
    }

    // Last triangle.
    if indices.len() == 3 {
        let (va, vb, vc) = (poly[indices[0]], poly[indices[1]], poly[indices[2]]);
        let pa = pool.position(va);
        let pb = pool.position(vb);
        let pc = pool.position(vc);
        if (pb - pa).cross(&(pc - pa)).norm_squared() > 1e-30
            && !would_create_nm(va, vb, vc, valence)
        {
            out.push(FaceData::untagged(va, vb, vc));
            record_triangle(va, vb, vc, valence);
            count += 1;
        }
    }

    count
}

/// Compute polygon normal via Newell's method (numerically stable for
/// non-planar polygons).
fn newell_normal(positions: &[nalgebra::Point3<Real>]) -> (Real, Real, Real) {
    let n = positions.len();
    let mut nx: Real = 0.0;
    let mut ny: Real = 0.0;
    let mut nz: Real = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        let pi = &positions[i];
        let pj = &positions[j];
        nx += (pi.y - pj.y) * (pi.z + pj.z);
        ny += (pi.z - pj.z) * (pi.x + pj.x);
        nz += (pi.x - pj.x) * (pi.y + pj.y);
    }
    (nx, ny, nz)
}

/// Signed area of a 2D polygon. Positive = counter-clockwise.
fn polygon_signed_area(pts: &[[Real; 2]]) -> Real {
    let n = pts.len();
    let mut area: Real = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        area += pts[i][0] * pts[j][1] - pts[j][0] * pts[i][1];
    }
    area * 0.5
}

/// Determine polygon winding direction.
///
/// Returns `Some(true)` for CCW and `Some(false)` for CW.
/// Uses exact orientation at an extreme vertex when possible, with
/// shoelace sign as fallback.
fn polygon_winding_ccw(pts: &[[Real; 2]]) -> Option<bool> {
    let n = pts.len();
    if n < 3 {
        return None;
    }

    let mut min_i = 0usize;
    for i in 1..n {
        if pts[i][0] < pts[min_i][0] || (pts[i][0] == pts[min_i][0] && pts[i][1] < pts[min_i][1]) {
            min_i = i;
        }
    }

    for k in 0..n {
        let i = (min_i + k) % n;
        let prev = (i + n - 1) % n;
        let next = (i + 1) % n;
        match orient_2d_arr(pts[prev], pts[i], pts[next]) {
            Orientation::Positive => return Some(true),
            Orientation::Negative => return Some(false),
            Orientation::Degenerate => {}
        }
    }

    let area = polygon_signed_area(pts);
    if area > 0.0 {
        Some(true)
    } else if area < 0.0 {
        Some(false)
    } else {
        None
    }
}

/// Test if point p lies inside (or on edge of) triangle (a, b, c) in 2D.
fn point_in_triangle(p: &[Real; 2], a: &[Real; 2], b: &[Real; 2], c: &[Real; 2]) -> bool {
    let d1 = orient_2d_arr(*a, *b, *p);
    let d2 = orient_2d_arr(*b, *c, *p);
    let d3 = orient_2d_arr(*c, *a, *p);
    let has_neg =
        d1 == Orientation::Negative || d2 == Orientation::Negative || d3 == Orientation::Negative;
    let has_pos =
        d1 == Orientation::Positive || d2 == Orientation::Positive || d3 == Orientation::Positive;
    !(has_neg && has_pos)
}

/// Remove duplicate faces (same 3 vertex IDs in any order).
fn dedup_faces(faces: &mut Vec<FaceData>) {
    let mut seen: HashSet<[VertexId; 3]> = HashSet::with_capacity(faces.len());
    faces.retain(|f| {
        let mut key = f.vertices;
        key.sort();
        seen.insert(key)
    });
}

/// Build undirected boundary vertex adjacency from directed boundary edges.
fn build_boundary_adjacency(boundary: &[(VertexId, VertexId)]) -> HashMap<VertexId, Vec<VertexId>> {
    let mut adj: HashMap<VertexId, Vec<VertexId>> = HashMap::new();
    for &(a, b) in boundary {
        adj.entry(a).or_default().push(b);
        adj.entry(b).or_default().push(a);
    }
    for neigh in adj.values_mut() {
        neigh.sort();
        neigh.dedup();
    }
    adj
}

#[inline]
fn collinear_3d_exact(pa: &Point3r, pb: &Point3r, pv: &Point3r) -> bool {
    orient_2d_arr([pa.x, pa.y], [pb.x, pb.y], [pv.x, pv.y]) == Orientation::Degenerate
        && orient_2d_arr([pa.x, pa.z], [pb.x, pb.z], [pv.x, pv.z]) == Orientation::Degenerate
        && orient_2d_arr([pa.y, pa.z], [pb.y, pb.z], [pv.y, pv.z]) == Orientation::Degenerate
}

/// Exact-predicate interior on-segment check in 3D.
///
/// Uses exact orientation predicates on all three coordinate-plane projections,
/// returning `Some(t)` for interior points.
fn point_on_segment_exact(pa: &Point3r, pb: &Point3r, pv: &Point3r) -> Option<Real> {
    let ab = pb - pa;
    let edge_len_sq = ab.norm_squared();
    if edge_len_sq < 1e-30 {
        return None;
    }
    let av = pv - pa;
    let t = av.dot(&ab) / edge_len_sq;
    if t <= SNAP_EDGE_PARAM_EPS || t >= 1.0 - SNAP_EDGE_PARAM_EPS {
        return None;
    }

    if !collinear_3d_exact(pa, pb, pv) {
        return None;
    }
    Some(t)
}

#[inline]
fn endpoint_constrained(
    v: VertexId,
    a: VertexId,
    b: VertexId,
    bnd_adj: &HashMap<VertexId, Vec<VertexId>>,
) -> bool {
    bnd_adj
        .get(&v)
        .is_some_and(|neigh| neigh.contains(&a) || neigh.contains(&b))
}

// ── T-junction snap-round ────────────────────────────────────────────────────

/// Squared perpendicular distance tolerance for snap-round.
/// A boundary vertex V is considered to lie on face edge [A,B] if:
///   distance(V, line(A,B)) < sqrt(SNAP_TOL_SQ) * |AB|
/// i.e., sin(angle(AV, AB)) < sqrt(SNAP_TOL_SQ) ≈ 0.001
const SNAP_TOL_SQ: Real = 1e-6;
/// Endpoint margin for edge-parameter test `t`.
///
/// Candidate split vertices must be strictly interior to the edge and at least
/// this fraction away from either endpoint.
const SNAP_EDGE_PARAM_EPS: Real = 5e-3;

/// Resolve T-junctions by splitting face edges at boundary vertices.
///
/// A T-junction occurs when boundary vertex V lies on an edge [A,B] of
/// a face [A,B,C] but V is not part of that face.  This creates an
/// unpaired half-edge at V.  Splitting [A,B,C] into [A,V,C] and [V,B,C]
/// inserts V into the mesh topology, pairing the boundary half-edge.
///
/// Candidate splits are constrained to vertices that are boundary-adjacent to
/// one of the edge endpoints.  On-segment detection uses exact orientation
/// predicates first, then falls back to the legacy distance gate.
///
/// Unlike vertex merging, edge splitting:
/// - Only ADDS faces (splits), never removes them
/// - Cannot create degenerate faces (V is verified to be interior to the edge)
/// - Cannot create topology cascades
///
/// Iterates up to `MAX_SNAP_ITERS` times until no more T-junctions are found
/// or boundary edges are eliminated.
pub(crate) fn snap_round_tjunctions(faces: &mut Vec<FaceData>, pool: &VertexPool) {
    const MAX_SNAP_ITERS: usize = 8;

    for _iter in 0..MAX_SNAP_ITERS {
        let boundary = build_boundary_edges(faces);
        if boundary.is_empty() {
            break;
        }

        // Collect unique boundary vertices.
        let mut bnd_verts: Vec<VertexId> = boundary.iter().flat_map(|&(a, b)| [a, b]).collect();
        bnd_verts.sort();
        bnd_verts.dedup();
        let bnd_set: HashSet<VertexId> = bnd_verts.iter().copied().collect();
        let bnd_adj = build_boundary_adjacency(&boundary);
        // For each face, keep the best constrained split candidate:
        // exact-hit first, then shortest residual distance, then centrality.
        let mut best_split_for_face: HashMap<
            usize,
            (bool, Real, Real, VertexId, VertexId, VertexId),
        > = HashMap::new();

        for &v in &bnd_verts {
            let pv = pool.position(v);

            for (fi, face) in faces.iter().enumerate() {
                // Skip faces that already contain this vertex.
                if face.vertices.contains(&v) {
                    continue;
                }

                for edge_idx in 0..3_usize {
                    let a = face.vertices[edge_idx];
                    let b = face.vertices[(edge_idx + 1) % 3];

                    let pa = pool.position(a);
                    let pb = pool.position(b);
                    if !endpoint_constrained(v, a, b, &bnd_adj) {
                        continue;
                    }

                    let mut exact_hit = false;
                    let t;
                    let dist_metric;
                    if let Some(t_exact) = point_on_segment_exact(pa, pb, pv) {
                        exact_hit = true;
                        t = t_exact;
                        dist_metric = 0.0;
                    } else {
                        let ab = pb - pa;
                        let av = pv - pa;
                        let edge_len_sq = ab.norm_squared();
                        if edge_len_sq < 1e-30 {
                            continue;
                        }

                        // Projection parameter: t = dot(AV, AB) / |AB|²
                        t = av.dot(&ab) / edge_len_sq;
                        // V must be strictly interior to edge (not at endpoints).
                        if t <= SNAP_EDGE_PARAM_EPS || t >= 1.0 - SNAP_EDGE_PARAM_EPS {
                            continue;
                        }

                        // Perpendicular distance: |cross(AB, AV)|² / |AB|²
                        let cross = ab.cross(&av);
                        let dist_sq = cross.norm_squared() / edge_len_sq;
                        if dist_sq > SNAP_TOL_SQ * edge_len_sq {
                            continue;
                        }
                        dist_metric = dist_sq / edge_len_sq.max(1e-30);
                    }

                    // Found a T-junction: V lies on edge [A, B] of face fi.
                    // Keep the historical seam-ribbon guard for tolerance-only
                    // hits, but allow exact constrained hits through.
                    if !exact_hit && bnd_set.contains(&a) && bnd_set.contains(&b) {
                        continue;
                    }

                    let center_bias = (t - 0.5).abs();
                    let candidate = (exact_hit, dist_metric, center_bias, a, b, v);
                    match best_split_for_face.get_mut(&fi) {
                        Some(best) => {
                            let best_key = (
                                !best.0,
                                best.1,
                                best.2,
                                best.5.raw(),
                                best.3.raw(),
                                best.4.raw(),
                            );
                            let cand_key = (
                                !candidate.0,
                                candidate.1,
                                candidate.2,
                                candidate.5.raw(),
                                candidate.3.raw(),
                                candidate.4.raw(),
                            );
                            if cand_key < best_key {
                                *best = candidate;
                            }
                        }
                        None => {
                            best_split_for_face.insert(fi, candidate);
                        }
                    }
                }
            }
        }

        let mut splits: Vec<(usize, VertexId, VertexId, VertexId)> = best_split_for_face
            .into_iter()
            .map(|(fi, (_, _, _, a, b, v))| (fi, a, b, v))
            .collect();
        if splits.is_empty() {
            break;
        }

        #[cfg(test)]
        eprintln!(
            "[snap-round {}] {} bnd edges, {} bnd verts, {} splits",
            _iter,
            boundary.len(),
            bnd_verts.len(),
            splits.len(),
        );

        // Apply splits: replace each face with two sub-faces.
        // Process in reverse order of face index to avoid invalidating indices.
        splits.sort_by(|a, b| b.0.cmp(&a.0));
        for (fi, a, b, v) in &splits {
            let face = &faces[*fi];
            // Find the third vertex (not a or b).
            let c = face
                .vertices
                .iter()
                .copied()
                .find(|&vid| vid != *a && vid != *b)
                .expect("face must have 3 distinct vertices");
            let region = face.region;

            // Determine winding order: face has edge [a, b] at some position.
            // The original face is [.., a, b, ..] in some winding.
            // We need [a, v, c] and [v, b, c] with correct winding.
            let verts = face.vertices;
            let (f1, f2) = if (verts[0] == *a && verts[1] == *b)
                || (verts[1] == *a && verts[2] == *b)
                || (verts[2] == *a && verts[0] == *b)
            {
                // Edge a→b is in forward winding: face is [a, b, c]
                // Split into [a, v, c] and [v, b, c]
                (
                    FaceData::new(*a, *v, c, region),
                    FaceData::new(*v, *b, c, region),
                )
            } else {
                // Edge b→a is in forward winding: face is [b, a, c] (effectively)
                // This means the original has a after b in winding.
                // Split into [b, v, c] and [v, a, c]
                (
                    FaceData::new(*b, *v, c, region),
                    FaceData::new(*v, *a, c, region),
                )
            };

            // Replace the original face with f1, append f2.
            faces[*fi] = f1;
            faces.push(f2);
        }

        // Remove any degenerate faces created by numerical edge cases.
        faces.retain(|f| {
            f.vertices[0] != f.vertices[1]
                && f.vertices[1] != f.vertices[2]
                && f.vertices[2] != f.vertices[0]
        });
        dedup_faces(faces);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::core::scalar::{Point3r, Vector3r};

    #[test]
    fn cdt_fill_loop_triangulates_concave_polygon() {
        let mut pool = VertexPool::default_millifluidic();
        let n = Vector3r::new(0.0, 0.0, 1.0);
        let poly = vec![
            pool.insert_or_weld(Point3r::new(0.0, 0.0, 0.0), n),
            pool.insert_or_weld(Point3r::new(2.0, 0.0, 0.0), n),
            pool.insert_or_weld(Point3r::new(2.0, 1.0, 0.0), n),
            pool.insert_or_weld(Point3r::new(1.0, 1.0, 0.0), n),
            pool.insert_or_weld(Point3r::new(1.0, 2.0, 0.0), n),
            pool.insert_or_weld(Point3r::new(0.0, 2.0, 0.0), n),
        ];

        let mut out = Vec::new();
        let mut valence = HashMap::new();
        let added = cdt_fill_loop(&poly, &pool, &mut out, &mut valence);

        assert_eq!(
            added,
            poly.len() - 2,
            "simple concave polygon should triangulate to n-2 triangles"
        );

        let mut area_sum = 0.0_f64;
        for f in &out {
            let a = pool.position(f.vertices[0]);
            let b = pool.position(f.vertices[1]);
            let c = pool.position(f.vertices[2]);
            let area2 = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
            assert!(area2.abs() > 1e-12, "CDT fill produced degenerate triangle");
            area_sum += area2.abs() * 0.5;
        }
        assert!(
            (area_sum - 3.0).abs() < 1e-9,
            "triangulated area must match polygon area (got {area_sum:.12})"
        );
    }

    #[test]
    fn cdt_fill_loop_rejects_collinear_polygon() {
        let mut pool = VertexPool::default_millifluidic();
        let n = Vector3r::new(0.0, 0.0, 1.0);
        let poly = vec![
            pool.insert_or_weld(Point3r::new(0.0, 0.0, 0.0), n),
            pool.insert_or_weld(Point3r::new(1.0, 0.0, 0.0), n),
            pool.insert_or_weld(Point3r::new(2.0, 0.0, 0.0), n),
            pool.insert_or_weld(Point3r::new(3.0, 0.0, 0.0), n),
        ];

        let mut out = Vec::new();
        let mut valence = HashMap::new();
        let added = cdt_fill_loop(&poly, &pool, &mut out, &mut valence);

        assert_eq!(added, 0, "collinear loop must not be triangulated");
        assert!(
            out.is_empty(),
            "no faces should be emitted for collinear loop"
        );
    }

    #[test]
    fn snap_round_splits_exact_tjunction_with_endpoint_constraint() {
        let mut pool = VertexPool::default_millifluidic();
        let n = Vector3r::new(0.0, 0.0, 1.0);
        let a = pool.insert_or_weld(Point3r::new(0.0, 0.0, 0.0), n);
        let b = pool.insert_or_weld(Point3r::new(2.0, 0.0, 0.0), n);
        let c = pool.insert_or_weld(Point3r::new(0.0, 1.0, 0.0), n);
        let m = pool.insert_or_weld(Point3r::new(1.0, 0.0, 0.0), n);
        let d = pool.insert_or_weld(Point3r::new(2.0, 1.0, 0.0), n);

        let mut faces = vec![FaceData::untagged(a, b, c), FaceData::untagged(m, b, d)];
        snap_round_tjunctions(&mut faces, &pool);

        assert_eq!(
            faces.len(),
            3,
            "one constrained exact split should be applied"
        );
        let with_mc = faces
            .iter()
            .filter(|f| f.vertices.contains(&m) && f.vertices.contains(&c))
            .count();
        assert_eq!(
            with_mc, 2,
            "split should replace [a,b,c] with two triangles using split vertex m"
        );
    }

    #[test]
    fn snap_round_rejects_exact_collinear_vertex_without_constraint() {
        let mut pool = VertexPool::default_millifluidic();
        let n = Vector3r::new(0.0, 0.0, 1.0);
        let a = pool.insert_or_weld(Point3r::new(0.0, 0.0, 0.0), n);
        let b = pool.insert_or_weld(Point3r::new(2.0, 0.0, 0.0), n);
        let c = pool.insert_or_weld(Point3r::new(0.0, 1.0, 0.0), n);
        let m = pool.insert_or_weld(Point3r::new(1.0, 0.0, 0.0), n);
        let d = pool.insert_or_weld(Point3r::new(1.0, 1.0, 0.0), n);
        let e = pool.insert_or_weld(Point3r::new(2.0, 1.0, 0.0), n);

        // m lies exactly on [a,b] but is not boundary-adjacent to a or b.
        let mut faces = vec![FaceData::untagged(a, b, c), FaceData::untagged(m, d, e)];
        snap_round_tjunctions(&mut faces, &pool);

        assert_eq!(
            faces.len(),
            2,
            "unconstrained collinear vertices must not trigger edge splits"
        );
    }

    #[test]
    fn point_on_segment_exact_accepts_true_3d_collinear_point() {
        let a = Point3r::new(0.0, 0.0, 0.0);
        let b = Point3r::new(2.0, 4.0, 6.0);
        let p = Point3r::new(1.0, 2.0, 3.0); // exact midpoint
        let t = point_on_segment_exact(&a, &b, &p);
        assert!(t.is_some(), "midpoint on 3D segment must be accepted");
        let tv = t.unwrap();
        assert!(
            (tv - 0.5).abs() < 1e-12,
            "expected midpoint t=0.5, got {tv}"
        );
    }

    #[test]
    fn point_on_segment_exact_rejects_projection_false_positive() {
        let a = Point3r::new(0.0, 0.0, 0.0);
        let b = Point3r::new(2.0, 0.0, 0.0);
        let p = Point3r::new(1.0, 0.0, 1.0); // collinear in XY but not in XZ
        assert!(
            point_on_segment_exact(&a, &b, &p).is_none(),
            "point off the 3D line must be rejected"
        );
    }
}
