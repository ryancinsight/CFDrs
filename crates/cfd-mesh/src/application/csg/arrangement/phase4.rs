//! Phase 4 fragment classification for arrangement CSG.

use super::classify::{
    centroid, classify_fragment_prepared, prepare_classification_faces, tri_normal, FragRecord,
    FragmentClass,
};
use super::dsu::DisjointSet;
use crate::application::csg::boolean::BooleanOp;
use crate::application::csg::predicates3d::triangle_is_degenerate_exact;
use crate::domain::core::index::VertexId;
use crate::domain::core::scalar::{Point3r, Vector3r};
use crate::domain::topology::predicates::{orient3d, Sign};
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;
#[cfg(test)]
use std::collections::BTreeMap;

/// Build connected-component roots for fragment adjacency induced by shared edges.
///
/// Uses a sort-and-scan edge run builder instead of hash-bucket vectors:
/// all canonical undirected edges are collected in a flat vector, sorted, then
/// processed as contiguous runs of equal edge keys.
///
/// # Theorem — Component Equivalence
///
/// For each undirected mesh edge, earlier code connected all incident fragments
/// as a clique (complete graph). Replacing each clique by a star rooted at one
/// arbitrary incident fragment preserves connected components because every
/// clique vertex remains reachable through the hub. Therefore any per-component
/// classification cache keyed by these roots is equivalent to flood-filling the
/// original adjacency graph. ∎
///
/// # Theorem — Sort-Run / Hash-Bucket Equivalence
///
/// Grouping incident fragments by edge key via hash buckets or via sorted equal
/// runs yields the same partition of edge incidences by key. Since each key-group
/// induces the same star-union operations, DSU connectivity is identical. ∎
fn fragment_component_roots(frags: &[FragRecord]) -> Vec<usize> {
    let mut edge_refs: Vec<(VertexId, VertexId, usize, bool)> = Vec::with_capacity(frags.len() * 3);
    for (i, frag) in frags.iter().enumerate() {
        let v = frag.face.vertices;
        for j in 0..3 {
            let a = v[j];
            let b = v[(j + 1) % 3];
            let (mn, mx) = if a < b { (a, b) } else { (b, a) };
            edge_refs.push((mn, mx, i, frag.from_a));
        }
    }
    edge_refs.sort_unstable_by_key(|&(a, b, _, _)| (a.raw(), b.raw()));

    let mut dsu = DisjointSet::new(frags.len());

    let mut run_start = 0usize;
    while run_start < edge_refs.len() {
        let (ea, eb, _, _) = edge_refs[run_start];
        let mut run_end = run_start + 1;
        while run_end < edge_refs.len() && edge_refs[run_end].0 == ea && edge_refs[run_end].1 == eb
        {
            run_end += 1;
        }

        let mut has_a = false;
        let mut has_b = false;
        for &(_, _, _, from_a) in &edge_refs[run_start..run_end] {
            if from_a {
                has_a = true;
            } else {
                has_b = true;
            }
        }
        // Do not connect across A/B seam edges.
        if has_a && has_b {
            run_start = run_end;
            continue;
        }
        if run_end - run_start >= 2 {
            let root = edge_refs[run_start].2;
            for &(_, _, idx, _) in &edge_refs[(run_start + 1)..run_end] {
                dsu.union(root, idx);
            }
        }
        run_start = run_end;
    }

    let mut roots = Vec::with_capacity(frags.len());
    for i in 0..frags.len() {
        roots.push(dsu.find(i));
    }
    roots
}

/// Classify co-refined fragments and return kept (optionally flipped) faces.
pub(crate) fn classify_kept_fragments(
    op: BooleanOp,
    frags: &[FragRecord],
    faces_a: &[FaceData],
    faces_b: &[FaceData],
    pool: &VertexPool,
    coplanar_groups: &[(usize, FaceData)],
) -> Vec<FaceData> {
    // Pre-compute coplanar plane data for coplanar-fragment exclusion below.
    struct CoplanarPlaneInfo {
        a: Point3r,
        b: Point3r,
        c: Point3r,
        valid_plane: bool,
    }
    let coplanar_plane_infos: Vec<CoplanarPlaneInfo> = coplanar_groups
        .iter()
        .map(|(_, rep_face)| {
            let r0 = *pool.position(rep_face.vertices[0]);
            let r1 = *pool.position(rep_face.vertices[1]);
            let r2 = *pool.position(rep_face.vertices[2]);
            CoplanarPlaneInfo {
                a: r0,
                b: r1,
                c: r2,
                valid_plane: !triangle_is_degenerate_exact(&r0, &r1, &r2),
            }
        })
        .collect();

    let component_roots = fragment_component_roots(frags);
    let mut class_cache: Vec<Option<FragmentClass>> = vec![None; frags.len()];
    let prepared_a = prepare_classification_faces(faces_a, pool);
    let prepared_b = prepare_classification_faces(faces_b, pool);
    let mut kept_faces: Vec<FaceData> = Vec::new();

    for (frag_index, frag) in frags.iter().enumerate() {
        let p0 = *pool.position(frag.face.vertices[0]);
        let p1 = *pool.position(frag.face.vertices[1]);
        let p2 = *pool.position(frag.face.vertices[2]);
        let tri = [p0, p1, p2];

        // Exclude fragments fully lying on a coplanar plane; they are handled by
        // Phase 2c coplanar boolean output.
        {
            let mut on_any_coplanar_plane = false;
            for cp in &coplanar_plane_infos {
                if !cp.valid_plane {
                    continue;
                }
                if orient3d(&cp.a, &cp.b, &cp.c, &p0) == Sign::Zero
                    && orient3d(&cp.a, &cp.b, &cp.c, &p1) == Sign::Zero
                    && orient3d(&cp.a, &cp.b, &cp.c, &p2) == Sign::Zero
                {
                    on_any_coplanar_plane = true;
                    break;
                }
            }
            if on_any_coplanar_plane {
                continue;
            }
        }

        let c = centroid(&tri);
        let n = tri_normal(&tri);
        let nlen = n.norm();
        let face_normal = if nlen > 1e-20 {
            n / nlen
        } else {
            Vector3r::zeros()
        };

        // Skip near-degenerate seam slivers.
        {
            let e01_sq = (p1 - p0).norm_squared();
            let e02_sq = (p2 - p0).norm_squared();
            let e12_sq = (p2 - p1).norm_squared();
            let area_sq = n.norm_squared();
            let max_edge_sq = e01_sq.max(e02_sq).max(e12_sq);
            if max_edge_sq > 1e-20 && area_sq < 1e-10 * max_edge_sq {
                continue;
            }
        }

        let mut eval_class = |is_a: bool| -> FragmentClass {
            let comp_root = component_roots[frag_index];
            if let Some(val) = class_cache[comp_root] {
                return val;
            }
            let val = if is_a {
                classify_fragment_prepared(&c, &face_normal, &prepared_b)
            } else {
                classify_fragment_prepared(&c, &face_normal, &prepared_a)
            };
            class_cache[comp_root] = Some(val);
            val
        };

        let (keep, flip) = if frag.from_a {
            let class_b = eval_class(true);
            match op {
                BooleanOp::Union => (
                    class_b == FragmentClass::Outside || class_b == FragmentClass::CoplanarSame,
                    false,
                ),
                BooleanOp::Intersection => (
                    class_b == FragmentClass::Inside || class_b == FragmentClass::CoplanarSame,
                    false,
                ),
                BooleanOp::Difference => (class_b == FragmentClass::Outside, false),
            }
        } else {
            let class_a = eval_class(false);
            match op {
                BooleanOp::Union => (class_a == FragmentClass::Outside, false),
                BooleanOp::Intersection => (class_a == FragmentClass::Inside, false),
                BooleanOp::Difference => (
                    class_a == FragmentClass::Inside || class_a == FragmentClass::CoplanarOpposite,
                    true,
                ),
            }
        };

        if !keep {
            continue;
        }

        let parent_face = if frag.from_a {
            faces_a[frag.parent_idx]
        } else {
            faces_b[frag.parent_idx]
        };

        if flip {
            kept_faces.push(FaceData::new(
                frag.face.vertices[0],
                frag.face.vertices[2],
                frag.face.vertices[1],
                parent_face.region,
            ));
        } else {
            kept_faces.push(FaceData::new(
                frag.face.vertices[0],
                frag.face.vertices[1],
                frag.face.vertices[2],
                parent_face.region,
            ));
        }
    }

    kept_faces
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::infrastructure::storage::face_store::FaceData;

    fn component_signature(frags: &[FragRecord], roots: &[usize]) -> Vec<Vec<[u32; 3]>> {
        let mut groups: BTreeMap<usize, Vec<[u32; 3]>> = BTreeMap::new();
        for (i, frag) in frags.iter().enumerate() {
            let mut key = frag.face.vertices;
            key.sort();
            groups
                .entry(roots[i])
                .or_default()
                .push([key[0].raw(), key[1].raw(), key[2].raw()]);
        }
        let mut out: Vec<Vec<[u32; 3]>> = groups.into_values().collect();
        for g in &mut out {
            g.sort_unstable();
        }
        out.sort_unstable();
        out
    }

    #[test]
    fn fragment_components_respect_same_source_connectivity() {
        let v0 = VertexId::new(0);
        let v1 = VertexId::new(1);
        let v2 = VertexId::new(2);
        let v3 = VertexId::new(3);
        let v4 = VertexId::new(4);
        let v5 = VertexId::new(5);
        let v6 = VertexId::new(6);
        let v7 = VertexId::new(7);
        let v8 = VertexId::new(8);

        // A-chain: f0 --(1,2)-- f1 --(1,3)-- f2
        let f0 = FragRecord {
            face: FaceData::untagged(v0, v1, v2),
            parent_idx: 0,
            from_a: true,
        };
        let f1 = FragRecord {
            face: FaceData::untagged(v2, v1, v3),
            parent_idx: 1,
            from_a: true,
        };
        let f2 = FragRecord {
            face: FaceData::untagged(v1, v3, v4),
            parent_idx: 2,
            from_a: true,
        };

        // B fragments on mixed seam edge (0,1): must not connect to A component.
        let f3 = FragRecord {
            face: FaceData::untagged(v0, v1, v5),
            parent_idx: 3,
            from_a: false,
        };

        // Pure B pair on edge (5,6): should connect.
        let f4 = FragRecord {
            face: FaceData::untagged(v5, v6, v7),
            parent_idx: 4,
            from_a: false,
        };
        let f5 = FragRecord {
            face: FaceData::untagged(v6, v5, v8),
            parent_idx: 5,
            from_a: false,
        };

        let frags = vec![f0, f1, f2, f3, f4, f5];
        let roots = fragment_component_roots(&frags);

        assert_eq!(roots[0], roots[1]);
        assert_eq!(roots[1], roots[2]);
        assert_ne!(roots[3], roots[0]);
        assert_eq!(roots[4], roots[5]);
        assert_ne!(roots[4], roots[0]);
    }

    #[test]
    fn adversarial_permuted_fragment_order_keeps_components() {
        let v0 = VertexId::new(0);
        let v1 = VertexId::new(1);
        let v2 = VertexId::new(2);
        let v3 = VertexId::new(3);
        let v4 = VertexId::new(4);
        let v5 = VertexId::new(5);
        let v6 = VertexId::new(6);
        let v7 = VertexId::new(7);
        let v8 = VertexId::new(8);
        let v9 = VertexId::new(9);

        let frags_a = vec![
            FragRecord {
                face: FaceData::untagged(v0, v1, v2),
                parent_idx: 0,
                from_a: true,
            },
            FragRecord {
                face: FaceData::untagged(v2, v1, v3),
                parent_idx: 1,
                from_a: true,
            },
            FragRecord {
                face: FaceData::untagged(v3, v1, v4),
                parent_idx: 2,
                from_a: true,
            },
            FragRecord {
                face: FaceData::untagged(v5, v6, v7),
                parent_idx: 3,
                from_a: false,
            },
            FragRecord {
                face: FaceData::untagged(v7, v6, v8),
                parent_idx: 4,
                from_a: false,
            },
            // Mixed edge with A side: must remain disconnected.
            FragRecord {
                face: FaceData::untagged(v0, v1, v9),
                parent_idx: 5,
                from_a: false,
            },
        ];
        let mut frags_b = vec![
            FragRecord {
                face: FaceData::untagged(v0, v1, v2),
                parent_idx: 0,
                from_a: true,
            },
            FragRecord {
                face: FaceData::untagged(v2, v1, v3),
                parent_idx: 1,
                from_a: true,
            },
            FragRecord {
                face: FaceData::untagged(v3, v1, v4),
                parent_idx: 2,
                from_a: true,
            },
            FragRecord {
                face: FaceData::untagged(v5, v6, v7),
                parent_idx: 3,
                from_a: false,
            },
            FragRecord {
                face: FaceData::untagged(v7, v6, v8),
                parent_idx: 4,
                from_a: false,
            },
            FragRecord {
                face: FaceData::untagged(v0, v1, v9),
                parent_idx: 5,
                from_a: false,
            },
        ];
        frags_b.reverse();

        let roots_a = fragment_component_roots(&frags_a);
        let roots_b = fragment_component_roots(&frags_b);
        assert_eq!(
            component_signature(&frags_a, &roots_a),
            component_signature(&frags_b, &roots_b),
            "component partition should be invariant to fragment order"
        );
    }
}
