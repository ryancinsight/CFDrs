//! Boolean operations: union, difference, intersection.
//!
//! This module exposes two API levels:
//!
//! ## Low-level face-soup API (backward compatible)
//!
//! [`csg_boolean`] operates on `Vec<FaceData>` + a shared `VertexPool`.
//! Internally uses the BSP-tree approach from `csgrs`.
//!
//! ## High-level `IndexedMesh` API (new in Phase 8)
//!
//! [`csg_boolean_indexed`] takes two `IndexedMesh` objects, merges their
//! vertex pools, runs the Boolean pipeline, and reconstructs a fresh mesh.
//!
//! [`CsgNode`] provides a composable CSG tree that evaluates lazily.
//!
//! ## Containment handling
//!
//! The pure BSP approach fails when one mesh is entirely enclosed within the
//! other (no surface intersection exists), because no BSP splitting plane from
//! one mesh ever splits the other's triangles.  The functions here detect this
//! case via a parity ray-cast (`point_in_mesh`) and handle it analytically:
//!
//! | Operation | B ⊂ A        | A ⊂ B        | Disjoint (A outside B) |
//! |-----------|-------------|-------------|------------------------|
//! | A ∪ B     | A           | B           | A + B (concatenate)   |
//! | A ∩ B     | B           | A           | empty (error)         |
//! | A \ B     | A + flip(B) | empty (err) | A (unchanged)         |

use crate::core::error::{MeshError, MeshResult};
use crate::core::index::RegionId;
use crate::core::scalar::{Point3r, Real, Vector3r};
use crate::mesh::IndexedMesh;
use crate::storage::face_store::FaceData;
use crate::storage::vertex_pool::VertexPool;
use crate::csg::bsp::BspNode;
use crate::csg::reconstruct;
use nalgebra::Isometry3;

/// CSG boolean operation type.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BooleanOp {
    /// A ∪ B
    Union,
    /// A ∩ B
    Intersection,
    /// A \ B
    Difference,
}

/// Perform a CSG boolean operation on two sets of faces.
///
/// Both sets share the same `VertexPool`; new intersection vertices are
/// inserted into the pool during splitting.
pub fn csg_boolean(
    op: BooleanOp,
    faces_a: &[FaceData],
    faces_b: &[FaceData],
    pool: &mut VertexPool,
) -> MeshResult<Vec<FaceData>> {
    let result = match op {
        BooleanOp::Union => boolean_union(faces_a, faces_b, pool),
        BooleanOp::Intersection => boolean_intersection(faces_a, faces_b, pool),
        BooleanOp::Difference => boolean_difference(faces_a, faces_b, pool),
    };

    if result.is_empty() {
        return Err(MeshError::EmptyBooleanResult {
            op: format!("{:?}", op),
        });
    }

    Ok(result)
}

// ── Containment detection ─────────────────────────────────────────────────────

/// Spatial relationship between two face soups.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Containment {
    /// The meshes' surfaces actually cross (or cannot be determined otherwise).
    Intersecting,
    /// Mesh B is entirely inside mesh A — no surface crossing.
    BInsideA,
    /// Mesh A is entirely inside mesh B — no surface crossing.
    AInsideB,
    /// The meshes are entirely disjoint — no overlap of any kind.
    Disjoint,
}

/// Test whether `query` is inside a closed triangle mesh using parity ray-
/// casting along +X.
///
/// Fires a ray from `query` in the +X direction and counts triangle crossings.
/// An odd count means the point is inside the solid (Jordan curve theorem).
///
/// ## Edge-hit handling
///
/// When the ray passes exactly through an edge or vertex shared by two
/// triangles, both triangles would be counted, giving an even (wrong) result.
/// This is handled by using the **strict** barycentric bounds `(0, 1)` rather
/// than `[0, 1]`: a hit on the shared boundary `u=0` or `v=0` is only counted
/// once — for the face where that coordinate is strictly positive — so the pair
/// contributes exactly one crossing.
pub(crate) fn point_in_mesh(query: &Point3r, faces: &[FaceData], pool: &VertexPool) -> bool {
    let mut crossings = 0usize;
    let (ox, oy, oz) = (query.x, query.y, query.z);

    for face in faces {
        let a = pool.position(face.vertices[0]);
        let b = pool.position(face.vertices[1]);
        let c = pool.position(face.vertices[2]);

        // Möller–Trumbore with ray direction = +X.
        let edge1 = Vector3r::new(b.x - a.x, b.y - a.y, b.z - a.z);
        let edge2 = Vector3r::new(c.x - a.x, c.y - a.y, c.z - a.z);
        let ray_dir = Vector3r::x();

        let h   = ray_dir.cross(&edge2);
        let det = edge1.dot(&h);
        if det.abs() < 1e-14 { continue; }

        let inv_det = 1.0 / det;
        let s = Vector3r::new(ox - a.x, oy - a.y, oz - a.z);
        let u = inv_det * s.dot(&h);
        // Strict lower bound on u avoids double-counting the u=0 edge boundary
        // when the ray grazes an edge shared by two adjacent triangles.
        // The convention is: the u=0 edge belongs to the *other* triangle.
        if u <= 0.0 || u > 1.0 { continue; }

        let q = s.cross(&edge1);
        let v = inv_det * ray_dir.dot(&q);
        // v >= 0 (inclusive) because u > 0 already handles the vertex ambiguity.
        if v < 0.0 || u + v > 1.0 { continue; }

        let t = inv_det * edge2.dot(&q);
        if t > 1e-10 { crossings += 1; }
    }

    crossings % 2 == 1
}

/// Determine the containment relationship between two face soups in `pool`.
///
/// ## Strategy
///
/// 1. **Disjoint**: AABBs do not overlap → return immediately.
/// 2. **Partial overlap**: neither AABB strictly contains the other → the
///    surfaces must cross → return `Intersecting` immediately (no ray cast
///    needed and avoids sample-point boundary ambiguity).
/// 3. **One AABB inside the other**: do a point-in-mesh ray cast using the
///    centre of the inner AABB as the sample point.  The AABB centre of a
///    convex mesh is strictly interior, so it never lies on a face plane of
///    the enclosing mesh.
fn containment(faces_a: &[FaceData], faces_b: &[FaceData], pool: &VertexPool) -> Containment {
    use crate::geometry::aabb::Aabb;

    let aabb_of = |faces: &[FaceData]| -> Option<Aabb> {
        let mut bb = Aabb::empty();
        let mut any = false;
        for f in faces {
            for &vid in &f.vertices {
                bb.expand(pool.position(vid));
                any = true;
            }
        }
        if any { Some(bb) } else { None }
    };

    let aabb_a = match aabb_of(faces_a) { Some(a) => a, None => return Containment::Disjoint };
    let aabb_b = match aabb_of(faces_b) { Some(b) => b, None => return Containment::Disjoint };

    if !aabb_a.intersects(&aabb_b) {
        return Containment::Disjoint;
    }

    // Check strict AABB containment: if neither AABB is fully inside the other,
    // the two surfaces must geometrically intersect → BSP handles it.
    // `aabb_x.contains_point(min)` and `contains_point(max)` of the other covers all 8 corners.
    let b_in_aabb_a = aabb_a.contains_point(&aabb_b.min) && aabb_a.contains_point(&aabb_b.max);
    let a_in_aabb_b = aabb_b.contains_point(&aabb_a.min) && aabb_b.contains_point(&aabb_a.max);

    if !b_in_aabb_a && !a_in_aabb_b {
        // Partial AABB overlap → surfaces intersect → use BSP.
        return Containment::Intersecting;
    }

    // One AABB is inside the other.  Ray-cast using the inner AABB's centre:
    // it is guaranteed strictly interior to any convex mesh whose AABB contains it,
    // so it will never sit exactly on a face plane of the outer mesh.
    if b_in_aabb_a {
        let b_sample = aabb_b.center();
        if point_in_mesh(&b_sample, faces_a, pool) {
            return Containment::BInsideA;
        }
    }
    if a_in_aabb_b {
        let a_sample = aabb_a.center();
        if point_in_mesh(&a_sample, faces_b, pool) {
            return Containment::AInsideB;
        }
    }

    // AABB containment but ray cast says otherwise (concave mesh, etc.) → BSP.
    Containment::Intersecting
}

/// Reverse the winding of every face (flips the surface normal direction).
fn flip_faces(faces: &[FaceData]) -> Vec<FaceData> {
    faces.iter().map(|f| { let mut ff = *f; ff.flip(); ff }).collect()
}

/// Detect whether a face soup represents a curved (non-flat-faced) mesh.
///
/// Samples up to 32 faces and quantises their outward normals to a coarse
/// grid (0.1 resolution).  A mesh with only 6 distinct quantised normal
/// directions is treated as flat-faced (cube-like); any mesh with more than
/// 8 distinct directions is treated as curved.
///
/// **Complexity:** O(min(n, 32)) — negligible overhead.
///
/// | Mesh type    | Distinct quantised normals |
/// |--------------|---------------------------|
/// | Cube         | 6 (axis-aligned faces)     |
/// | UV sphere    | >> 8 (one per latitude ring) |
/// | Cylinder     | >> 8 (many around the barrel) |
fn is_curved_mesh(faces: &[FaceData], pool: &VertexPool) -> bool {
    use std::collections::HashSet;

    // Cap the sample at 32 faces for O(1) cost.
    let sample = faces.iter().take(32);
    let mut dirs: HashSet<(i32, i32, i32)> = HashSet::new();

    for face in sample {
        let a = pool.position(face.vertices[0]);
        let b = pool.position(face.vertices[1]);
        let c = pool.position(face.vertices[2]);

        // Face normal via cross product.
        let ab = Vector3r::new(b.x - a.x, b.y - a.y, b.z - a.z);
        let ac = Vector3r::new(c.x - a.x, c.y - a.y, c.z - a.z);
        let n = ab.cross(&ac);
        let len = n.norm();
        if len < 1e-20 { continue; }
        let (nx, ny, nz) = (n.x / len, n.y / len, n.z / len);

        // Quantise to nearest 0.1 to group near-identical normals.
        let key = (
            (nx * 10.0).round() as i32,
            (ny * 10.0).round() as i32,
            (nz * 10.0).round() as i32,
        );
        dirs.insert(key);
    }

    // Flat-face meshes (cubes, rectangular prisms) have at most 6 distinct
    // axis directions.  Curved meshes have many more.
    dirs.len() > 8
}

// ── Boolean operations ────────────────────────────────────────────────────────

/// A ∪ B — union of two face soups.
fn boolean_union(faces_a: &[FaceData], faces_b: &[FaceData], pool: &mut VertexPool) -> Vec<FaceData> {
    match containment(faces_a, faces_b, pool) {
        Containment::BInsideA => faces_a.to_vec(),
        Containment::AInsideB => faces_b.to_vec(),
        Containment::Disjoint => {
            let mut r = faces_a.to_vec();
            r.extend_from_slice(faces_b);
            r
        }
        Containment::Intersecting => {
            if is_curved_mesh(faces_a, pool) || is_curved_mesh(faces_b, pool) {
                return crate::csg::arrangement::boolean_intersecting_arrangement(
                    BooleanOp::Union, faces_a, faces_b, pool,
                );
            }
            let mut bsp_a = BspNode::build(faces_a, pool);
            let mut bsp_b = BspNode::build(faces_b, pool);
            bsp_a.clip_to(&bsp_b, pool);
            bsp_b.clip_to(&bsp_a, pool);
            bsp_b.invert();
            bsp_b.clip_to(&bsp_a, pool);
            bsp_b.invert();
            let b_faces = bsp_b.all_faces();
            bsp_a.add_faces(&b_faces, pool);
            bsp_a.all_faces()
        }
    }
}

/// A ∩ B — intersection of two face soups.
fn boolean_intersection(faces_a: &[FaceData], faces_b: &[FaceData], pool: &mut VertexPool) -> Vec<FaceData> {
    match containment(faces_a, faces_b, pool) {
        Containment::BInsideA => faces_b.to_vec(),
        Containment::AInsideB => faces_a.to_vec(),
        Containment::Disjoint => Vec::new(),
        Containment::Intersecting => {
            if is_curved_mesh(faces_a, pool) || is_curved_mesh(faces_b, pool) {
                return crate::csg::arrangement::boolean_intersecting_arrangement(
                    BooleanOp::Intersection, faces_a, faces_b, pool,
                );
            }
            let mut bsp_a = BspNode::build(faces_a, pool);
            let mut bsp_b = BspNode::build(faces_b, pool);
            bsp_a.invert();
            bsp_b.clip_to(&bsp_a, pool);
            bsp_b.invert();
            bsp_a.clip_to(&bsp_b, pool);
            bsp_b.clip_to(&bsp_a, pool);
            let b_faces = bsp_b.all_faces();
            bsp_a.add_faces(&b_faces, pool);
            bsp_a.invert();
            bsp_a.all_faces()
        }
    }
}

/// A \ B — subtract B from A.
fn boolean_difference(faces_a: &[FaceData], faces_b: &[FaceData], pool: &mut VertexPool) -> Vec<FaceData> {
    match containment(faces_a, faces_b, pool) {
        Containment::BInsideA => {
            // B is a cavity inside A: keep A's outer surface + B's surface flipped inward.
            let mut r = faces_a.to_vec();
            r.extend(flip_faces(faces_b));
            r
        }
        Containment::AInsideB => Vec::new(),  // A is swallowed by B → empty
        Containment::Disjoint => faces_a.to_vec(),  // B doesn't touch A → A unchanged
        Containment::Intersecting => {
            if is_curved_mesh(faces_a, pool) || is_curved_mesh(faces_b, pool) {
                return crate::csg::arrangement::boolean_intersecting_arrangement(
                    BooleanOp::Difference, faces_a, faces_b, pool,
                );
            }
            let mut bsp_a = BspNode::build(faces_a, pool);
            let mut bsp_b = BspNode::build(faces_b, pool);
            bsp_a.invert();
            bsp_a.clip_to(&bsp_b, pool);
            bsp_b.clip_to(&bsp_a, pool);
            bsp_b.invert();
            bsp_b.clip_to(&bsp_a, pool);
            bsp_b.invert();
            let b_faces = bsp_b.all_faces();
            bsp_a.add_faces(&b_faces, pool);
            bsp_a.invert();
            bsp_a.all_faces()
        }
    }
}

// ── IndexedMesh API ───────────────────────────────────────────────────────────

/// High-level Boolean operation on two [`IndexedMesh`] objects.
///
/// Merges the vertex pools, runs the Boolean pipeline with containment
/// detection, and reconstructs a fresh deduplicated `IndexedMesh`.
pub fn csg_boolean_indexed(
    op: BooleanOp,
    mesh_a: &IndexedMesh,
    mesh_b: &IndexedMesh,
) -> MeshResult<IndexedMesh> {
    use std::collections::HashMap;
    use crate::core::index::VertexId;

    let mut combined = VertexPool::default_millifluidic();

    let mut remap_a: HashMap<VertexId, VertexId> = HashMap::new();
    for (old_id, _) in mesh_a.vertices.iter() {
        let pos = *mesh_a.vertices.position(old_id);
        let nrm = *mesh_a.vertices.normal(old_id);
        remap_a.insert(old_id, combined.insert_or_weld(pos, nrm));
    }

    let mut remap_b: HashMap<VertexId, VertexId> = HashMap::new();
    for (old_id, _) in mesh_b.vertices.iter() {
        let pos = *mesh_b.vertices.position(old_id);
        let nrm = *mesh_b.vertices.normal(old_id);
        remap_b.insert(old_id, combined.insert_or_weld(pos, nrm));
    }

    let faces_a: Vec<FaceData> = mesh_a.faces.iter().map(|f| {
        FaceData { vertices: f.vertices.map(|vid| remap_a[&vid]), region: f.region }
    }).collect();

    let faces_b: Vec<FaceData> = mesh_b.faces.iter().map(|f| {
        FaceData { vertices: f.vertices.map(|vid| remap_b[&vid]), region: f.region }
    }).collect();

    let result_faces = csg_boolean(op, &faces_a, &faces_b, &mut combined)?;
    Ok(reconstruct::reconstruct_mesh(&result_faces, &combined))
}

// ── CsgNode expression tree ───────────────────────────────────────────────────

/// A composable CSG expression tree over [`IndexedMesh`] operands.
pub enum CsgNode {
    /// A terminal mesh operand.
    Leaf(IndexedMesh),
    /// A ∪ B — mesh union.
    Union {
        /// Left operand.
        left: Box<CsgNode>,
        /// Right operand.
        right: Box<CsgNode>,
    },
    /// A ∩ B — mesh intersection.
    Intersection {
        /// Left operand.
        left: Box<CsgNode>,
        /// Right operand.
        right: Box<CsgNode>,
    },
    /// A \ B — subtract B from A.
    Difference {
        /// Minuend (mesh to subtract from).
        left: Box<CsgNode>,
        /// Subtrahend (mesh to subtract).
        right: Box<CsgNode>,
    },
    /// Apply a rigid-body transform to the sub-tree result.
    Transform {
        /// Sub-tree.
        node: Box<CsgNode>,
        /// Rigid-body isometry.
        iso: Isometry3<Real>,
    },
}

impl CsgNode {
    /// Evaluate the expression tree, consuming `self`.
    pub fn evaluate(self) -> MeshResult<IndexedMesh> {
        match self {
            CsgNode::Leaf(mesh) => Ok(mesh),
            CsgNode::Union { left, right } => {
                csg_boolean_indexed(BooleanOp::Union, &left.evaluate()?, &right.evaluate()?)
            }
            CsgNode::Intersection { left, right } => {
                csg_boolean_indexed(BooleanOp::Intersection, &left.evaluate()?, &right.evaluate()?)
            }
            CsgNode::Difference { left, right } => {
                csg_boolean_indexed(BooleanOp::Difference, &left.evaluate()?, &right.evaluate()?)
            }
            CsgNode::Transform { node, iso } => {
                Ok(transform_mesh(node.evaluate()?, &iso))
            }
        }
    }
}

// ── Transform helper ──────────────────────────────────────────────────────────

/// Apply a rigid-body `Isometry3` transform to all vertices of a mesh.
fn transform_mesh(mesh: IndexedMesh, iso: &Isometry3<Real>) -> IndexedMesh {
    use crate::core::index::VertexId;
    let mut new_mesh = IndexedMesh::new();
    let mut remap: Vec<Option<VertexId>> = vec![None; mesh.vertices.len()];

    for face in mesh.faces.iter() {
        let mut new_verts = [VertexId::default(); 3];
        for (k, &vid) in face.vertices.iter().enumerate() {
            let idx = vid.as_usize();
            let new_id = match remap[idx] {
                Some(id) => id,
                None => {
                    let pos = iso * *mesh.vertices.position(vid);
                    let nrm = iso.rotation * *mesh.vertices.normal(vid);
                    let id = new_mesh.add_vertex(pos, nrm);
                    remap[idx] = Some(id);
                    id
                }
            };
            new_verts[k] = new_id;
        }
        if face.region == RegionId::INVALID {
            new_mesh.add_face(new_verts[0], new_verts[1], new_verts[2]);
        } else {
            new_mesh.add_face_with_region(new_verts[0], new_verts[1], new_verts[2], face.region);
        }
    }
    new_mesh
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a closed 2×2×2 axis-aligned cube centred at (1,1,1) with outward CCW faces.
    fn make_unit_cube(pool: &mut VertexPool) -> Vec<FaceData> {
        let mut p = |x: f64, y: f64, z: f64| {
            pool.insert_or_weld(Point3r::new(x, y, z), Vector3r::zeros())
        };
        let (p000, p100, p110, p010) = (p(0.,0.,0.), p(2.,0.,0.), p(2.,2.,0.), p(0.,2.,0.));
        let (p001, p101, p111, p011) = (p(0.,0.,2.), p(2.,0.,2.), p(2.,2.,2.), p(0.,2.,2.));
        vec![
            // -Z
            FaceData::untagged(p000, p010, p110), FaceData::untagged(p000, p110, p100),
            // +Z
            FaceData::untagged(p001, p101, p111), FaceData::untagged(p001, p111, p011),
            // -Y
            FaceData::untagged(p000, p100, p101), FaceData::untagged(p000, p101, p001),
            // +Y
            FaceData::untagged(p010, p011, p111), FaceData::untagged(p010, p111, p110),
            // -X
            FaceData::untagged(p000, p001, p011), FaceData::untagged(p000, p011, p010),
            // +X
            FaceData::untagged(p100, p110, p111), FaceData::untagged(p100, p111, p101),
        ]
    }

    /// Build a 0.5×0.5×0.5 cube at offset (o,o,o).
    fn make_inner_cube(pool: &mut VertexPool, o: f64) -> Vec<FaceData> {
        let s = 0.5_f64;
        let mut p = |x: f64, y: f64, z: f64| {
            pool.insert_or_weld(Point3r::new(x, y, z), Vector3r::zeros())
        };
        let (p000, p100, p110, p010) = (p(o,o,o), p(o+s,o,o), p(o+s,o+s,o), p(o,o+s,o));
        let (p001, p101, p111, p011) = (p(o,o,o+s), p(o+s,o,o+s), p(o+s,o+s,o+s), p(o,o+s,o+s));
        vec![
            FaceData::untagged(p000, p010, p110), FaceData::untagged(p000, p110, p100),
            FaceData::untagged(p001, p101, p111), FaceData::untagged(p001, p111, p011),
            FaceData::untagged(p000, p100, p101), FaceData::untagged(p000, p101, p001),
            FaceData::untagged(p010, p011, p111), FaceData::untagged(p010, p111, p110),
            FaceData::untagged(p000, p001, p011), FaceData::untagged(p000, p011, p010),
            FaceData::untagged(p100, p110, p111), FaceData::untagged(p100, p111, p101),
        ]
    }

    #[test]
    fn point_in_cube_inside() {
        let mut pool = VertexPool::default_millifluidic();
        let faces = make_unit_cube(&mut pool);
        assert!(point_in_mesh(&Point3r::new(1.0, 1.0, 1.0), &faces, &pool));
    }

    #[test]
    fn point_in_cube_outside() {
        let mut pool = VertexPool::default_millifluidic();
        let faces = make_unit_cube(&mut pool);
        assert!(!point_in_mesh(&Point3r::new(5.0, 5.0, 5.0), &faces, &pool));
    }

    #[test]
    fn containment_b_inside_a_detected() {
        let mut pool = VertexPool::default_millifluidic();
        let faces_a = make_unit_cube(&mut pool);
        let faces_b = make_inner_cube(&mut pool, 0.75); // centred inside big cube
        assert_eq!(containment(&faces_a, &faces_b, &pool), Containment::BInsideA);
    }

    #[test]
    fn containment_disjoint_detected() {
        let mut pool = VertexPool::default_millifluidic();
        let faces_a = make_unit_cube(&mut pool);
        let mut p = |x: f64, y: f64, z: f64| pool.insert_or_weld(Point3r::new(x,y,z), Vector3r::zeros());
        // Tiny cube far away
        let (p0, p1, p2, p3) = (p(10.,10.,10.), p(11.,10.,10.), p(11.,11.,10.), p(10.,11.,10.));
        let (p4, p5, p6, p7) = (p(10.,10.,11.), p(11.,10.,11.), p(11.,11.,11.), p(10.,11.,11.));
        let faces_b = vec![
            FaceData::untagged(p0, p2, p1), FaceData::untagged(p0, p3, p2),
            FaceData::untagged(p4, p5, p6), FaceData::untagged(p4, p6, p7),
            FaceData::untagged(p0, p1, p5), FaceData::untagged(p0, p5, p4),
            FaceData::untagged(p3, p7, p6), FaceData::untagged(p3, p6, p2),
            FaceData::untagged(p0, p4, p7), FaceData::untagged(p0, p7, p3),
            FaceData::untagged(p1, p2, p6), FaceData::untagged(p1, p6, p5),
        ];
        assert_eq!(containment(&faces_a, &faces_b, &pool), Containment::Disjoint);
    }

    #[test]
    fn difference_enclosed_produces_outer_plus_flipped_inner() {
        let mut pool = VertexPool::default_millifluidic();
        let faces_a = make_unit_cube(&mut pool);
        let faces_b = make_inner_cube(&mut pool, 0.75);
        let result = boolean_difference(&faces_a, &faces_b, &mut pool);
        // Outer (12) + inner flipped (12) = 24
        assert_eq!(result.len(), 24);
    }

    #[test]
    fn union_enclosed_returns_outer() {
        let mut pool = VertexPool::default_millifluidic();
        let faces_a = make_unit_cube(&mut pool);
        let faces_b = make_inner_cube(&mut pool, 0.75);
        let result = boolean_union(&faces_a, &faces_b, &mut pool);
        assert_eq!(result.len(), faces_a.len());
    }

    #[test]
    fn intersection_enclosed_returns_inner() {
        let mut pool = VertexPool::default_millifluidic();
        let faces_a = make_unit_cube(&mut pool);
        let faces_b = make_inner_cube(&mut pool, 0.75);
        let result = boolean_intersection(&faces_a, &faces_b, &mut pool);
        assert_eq!(result.len(), faces_b.len());
    }

    #[test]
    fn difference_disjoint_returns_a() {
        let mut pool = VertexPool::default_millifluidic();
        let faces_a = make_unit_cube(&mut pool);
        let mut p = |x: f64, y: f64, z: f64| pool.insert_or_weld(Point3r::new(x,y,z), Vector3r::zeros());
        let (p0, p1, p2, p3) = (p(10.,10.,10.), p(11.,10.,10.), p(11.,11.,10.), p(10.,11.,10.));
        let (p4, p5, p6, p7) = (p(10.,10.,11.), p(11.,10.,11.), p(11.,11.,11.), p(10.,11.,11.));
        let faces_b = vec![
            FaceData::untagged(p0, p2, p1), FaceData::untagged(p0, p3, p2),
            FaceData::untagged(p4, p5, p6), FaceData::untagged(p4, p6, p7),
            FaceData::untagged(p0, p1, p5), FaceData::untagged(p0, p5, p4),
            FaceData::untagged(p3, p7, p6), FaceData::untagged(p3, p6, p2),
            FaceData::untagged(p0, p4, p7), FaceData::untagged(p0, p7, p3),
            FaceData::untagged(p1, p2, p6), FaceData::untagged(p1, p6, p5),
        ];
        let result = boolean_difference(&faces_a, &faces_b, &mut pool);
        assert_eq!(result.len(), faces_a.len());
    }
}
