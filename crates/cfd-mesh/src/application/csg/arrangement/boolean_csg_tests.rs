use super::boolean_csg::{csg_boolean, BooleanOp};
use crate::application::csg::boolean::containment::{containment, Containment};
use crate::domain::core::scalar::{Point3r, Vector3r};
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;
use hashbrown::HashMap;

fn make_cube(pool: &mut VertexPool, offset: Point3r, half_extent: f64) -> Vec<FaceData> {
    let mut vertex =
        |x, y, z| pool.insert_or_weld(offset + Vector3r::new(x, y, z), Vector3r::new(x, y, z));
    let s = half_extent;
    let c000 = vertex(-s, -s, -s);
    let c100 = vertex(s, -s, -s);
    let c010 = vertex(-s, s, -s);
    let c110 = vertex(s, s, -s);
    let c001 = vertex(-s, -s, s);
    let c101 = vertex(s, -s, s);
    let c011 = vertex(-s, s, s);
    let c111 = vertex(s, s, s);

    vec![
        FaceData::untagged(c000, c010, c110),
        FaceData::untagged(c000, c110, c100),
        FaceData::untagged(c001, c101, c111),
        FaceData::untagged(c001, c111, c011),
        FaceData::untagged(c000, c001, c011),
        FaceData::untagged(c000, c011, c010),
        FaceData::untagged(c100, c110, c111),
        FaceData::untagged(c100, c111, c101),
        FaceData::untagged(c010, c011, c111),
        FaceData::untagged(c010, c111, c110),
        FaceData::untagged(c000, c100, c101),
        FaceData::untagged(c000, c101, c001),
    ]
}

fn assert_watertight(faces: &[FaceData]) {
    let mut edges = HashMap::new();
    for face in faces {
        for edge_index in 0..3 {
            let v0 = face.vertices[edge_index];
            let v1 = face.vertices[(edge_index + 1) % 3];
            let key = if v0 < v1 { (v0, v1) } else { (v1, v0) };
            *edges.entry(key).or_insert(0) += 1;
        }
    }

    let non_manifold: Vec<_> = edges.iter().filter(|&(_, &count)| count != 2).collect();
    assert!(
        non_manifold.is_empty(),
        "canonical Boolean arrangement must remain watertight"
    );
}

#[test]
fn adversarial_boolean_union_is_watertight() {
    let mut pool = VertexPool::for_csg();
    let c1 = make_cube(&mut pool, Point3r::new(0.0, 0.0, 0.0), 1.0);
    let c2 = make_cube(&mut pool, Point3r::new(0.5, 0.5, 0.5), 1.0);
    let c3 = make_cube(&mut pool, Point3r::new(-0.25, 0.75, 0.1), 1.0);

    let faces = csg_boolean(BooleanOp::Union, &[c1, c2, c3], &mut pool).expect("union completes");
    assert_watertight(&faces);
}

#[test]
fn adversarial_boolean_intersection_is_watertight() {
    let mut pool = VertexPool::for_csg();
    let c1 = make_cube(&mut pool, Point3r::new(0.0, 0.0, 0.0), 1.0);
    let c2 = make_cube(&mut pool, Point3r::new(0.5, 0.5, 0.5), 1.0);
    let c3 = make_cube(&mut pool, Point3r::new(-0.25, 0.75, 0.1), 1.0);

    let faces = csg_boolean(BooleanOp::Intersection, &[c1, c2, c3], &mut pool)
        .expect("intersection completes");
    assert!(faces.len() > 10, "intersection should exist");
    assert_watertight(&faces);
}

#[test]
fn adversarial_boolean_difference_is_watertight() {
    let mut pool = VertexPool::for_csg();
    let base = make_cube(&mut pool, Point3r::new(0.0, 0.0, 0.0), 1.0);
    let sub1 = make_cube(&mut pool, Point3r::new(1.0, 1.0, 1.0), 1.0);
    let sub2 = make_cube(&mut pool, Point3r::new(0.8, -0.8, 0.8), 1.0);

    let faces = csg_boolean(BooleanOp::Difference, &[base, sub1, sub2], &mut pool)
        .expect("difference completes");
    assert!(faces.len() > 10, "difference should exist");
    assert_watertight(&faces);
}

#[test]
fn binary_disjoint_union_uses_concat_semantics() {
    let mut pool = VertexPool::for_csg();
    let left = make_cube(&mut pool, Point3r::new(-5.0, 0.0, 0.0), 1.0);
    let right = make_cube(&mut pool, Point3r::new(5.0, 0.0, 0.0), 1.0);
    assert_eq!(containment(&left, &right, &pool), Containment::Disjoint);

    let result = csg_boolean(BooleanOp::Union, &[left.clone(), right.clone()], &mut pool)
        .expect("disjoint union should succeed");
    assert_eq!(result.len(), left.len() + right.len());
}

#[test]
fn binary_contained_intersection_returns_inner_faces() {
    let mut pool = VertexPool::for_csg();
    let outer = make_cube(&mut pool, Point3r::new(0.0, 0.0, 0.0), 2.0);
    let inner = make_cube(&mut pool, Point3r::new(0.0, 0.0, 0.0), 0.5);
    assert_eq!(containment(&outer, &inner, &pool), Containment::BInsideA);

    let result = csg_boolean(BooleanOp::Intersection, &[outer, inner.clone()], &mut pool)
        .expect("contained intersection should succeed");
    assert_eq!(result, inner);
}
