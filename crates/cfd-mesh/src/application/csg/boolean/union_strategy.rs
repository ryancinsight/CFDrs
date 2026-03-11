use std::collections::HashMap;

use crate::domain::core::index::{FaceId, VertexId};
use crate::domain::mesh::IndexedMesh;

pub(super) fn spatial_union_order(meshes: &[IndexedMesh]) -> Vec<usize> {
    let n = meshes.len();
    if n <= 1 {
        return (0..n).collect();
    }

    let aabbs: Vec<_> = meshes.iter().map(IndexedMesh::bounding_box).collect();
    let centers: Vec<_> = aabbs.iter().map(|aabb| aabb.center()).collect();

    let seed = (0..n)
        .min_by(|&lhs, &rhs| {
            let a = &aabbs[lhs];
            let b = &aabbs[rhs];
            a.min.x
                .partial_cmp(&b.min.x)
                .unwrap()
                .then(a.min.y.partial_cmp(&b.min.y).unwrap())
                .then(a.min.z.partial_cmp(&b.min.z).unwrap())
                .then(lhs.cmp(&rhs))
        })
        .unwrap();

    let mut order = Vec::with_capacity(n);
    let mut used = vec![false; n];
    let mut accumulated = aabbs[seed];
    let mut accumulated_center = centers[seed];
    used[seed] = true;
    order.push(seed);

    while order.len() < n {
        let next = (0..n)
            .filter(|&idx| !used[idx])
            .min_by(|&lhs, &rhs| {
                let lhs_overlap = overlap_volume(&accumulated, &aabbs[lhs]);
                let rhs_overlap = overlap_volume(&accumulated, &aabbs[rhs]);
                rhs_overlap
                    .partial_cmp(&lhs_overlap)
                    .unwrap()
                    .then_with(|| {
                        center_distance2(&accumulated_center, &centers[lhs])
                            .partial_cmp(&center_distance2(&accumulated_center, &centers[rhs]))
                            .unwrap()
                    })
                    .then(lhs.cmp(&rhs))
            })
            .unwrap();

        used[next] = true;
        order.push(next);
        accumulated = accumulated.union(&aabbs[next]);
        accumulated_center = accumulated.center();
    }

    order
}

pub(super) fn concat_disjoint_meshes(lhs: &IndexedMesh, rhs: &IndexedMesh) -> IndexedMesh {
    let mut merged = lhs.empty_clone();
    append_mesh(&mut merged, lhs);
    append_mesh(&mut merged, rhs);
    merged.rebuild_edges();
    merged
}

fn append_mesh(dst: &mut IndexedMesh, src: &IndexedMesh) {
    let mut vertex_remap: HashMap<VertexId, VertexId> = HashMap::with_capacity(src.vertices.len());
    let mut face_remap: HashMap<FaceId, FaceId> = HashMap::with_capacity(src.faces.len());

    for (old_fid, face) in src.faces.iter_enumerated() {
        let vertices = face.vertices.map(|old_vid| {
            *vertex_remap.entry(old_vid).or_insert_with(|| {
                dst.add_vertex(*src.vertices.position(old_vid), *src.vertices.normal(old_vid))
            })
        });

        let new_fid = dst.add_face_with_region(vertices[0], vertices[1], vertices[2], face.region);
        face_remap.insert(old_fid, new_fid);
    }

    for channel in src.attributes.channel_names() {
        for (&old_fid, &new_fid) in &face_remap {
            if let Some(value) = src.attributes.get(channel, old_fid) {
                dst.attributes.set(channel, new_fid, value);
            }
        }
    }

    for (&old_fid, label) in &src.boundary_labels {
        if let Some(&new_fid) = face_remap.get(&old_fid) {
            dst.boundary_labels.insert(new_fid, label.clone());
        }
    }
}

fn overlap_volume(
    lhs: &crate::domain::geometry::Aabb,
    rhs: &crate::domain::geometry::Aabb,
) -> f64 {
    let dx = (lhs.max.x.min(rhs.max.x) - lhs.min.x.max(rhs.min.x)).max(0.0);
    let dy = (lhs.max.y.min(rhs.max.y) - lhs.min.y.max(rhs.min.y)).max(0.0);
    let dz = (lhs.max.z.min(rhs.max.z) - lhs.min.z.max(rhs.min.z)).max(0.0);
    dx * dy * dz
}

fn center_distance2(
    lhs: &crate::domain::core::scalar::Point3r,
    rhs: &crate::domain::core::scalar::Point3r,
) -> f64 {
    let delta = lhs - rhs;
    delta.norm_squared()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::core::scalar::Point3r;
    use crate::domain::geometry::primitives::{Cube, PrimitiveMesh};

    fn cube(x0: f64, x1: f64) -> IndexedMesh {
        Cube {
            origin: Point3r::new(x0, -1.0, -1.0),
            width: x1 - x0,
            height: 2.0,
            depth: 2.0,
        }
        .build()
        .expect("cube build")
    }

    #[test]
    fn spatial_union_order_prefers_overlapping_cluster_before_far_operand() {
        let meshes = vec![cube(10.0, 12.0), cube(0.0, 2.0), cube(1.0, 3.0)];
        let order = spatial_union_order(&meshes);
        assert_eq!(order.len(), 3);
        assert_eq!(order[0], 1);
        assert_eq!(order[1], 2);
        assert_eq!(order[2], 0);
    }

    #[test]
    fn concat_disjoint_meshes_preserves_watertight_closed_components() {
        let mut merged = concat_disjoint_meshes(&cube(0.0, 2.0), &cube(10.0, 12.0));
        assert!(merged.is_watertight());
        assert!(merged.signed_volume() > 0.0);
    }
}
