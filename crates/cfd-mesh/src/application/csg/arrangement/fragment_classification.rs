//! Fragment classification for arrangement CSG.

use super::classify::{
    centroid, classify_fragment_prepared, prepare_classification_faces, tri_normal, FragRecord,
    FragmentClass,
};
use super::fragment_analysis::{component_roots_by_source, is_degenerate_sliver_with_normal};
use crate::application::csg::boolean::BooleanOp;
use crate::application::csg::predicates3d::triangle_is_degenerate_exact;
#[cfg(test)]
use crate::domain::core::index::VertexId;
use crate::domain::core::scalar::{Point3r, Vector3r};
use crate::domain::topology::predicates::{orient3d, Sign};
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;
#[cfg(test)]
use std::collections::BTreeMap;

fn fragment_component_roots(frags: &[FragRecord]) -> Vec<usize> {
    component_roots_by_source(frags, |frag| frag.face.vertices, |frag| frag.from_a)
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

        let n = tri_normal(&tri);
        let nlen = n.norm();
        // Scale-relative: compare normal magnitude against edge lengths.
        let e1 = (p1 - p0).norm();
        let e2 = (p2 - p0).norm();
        let edge_product = e1 * e2;
        let face_normal = if nlen > 1e-10 * edge_product {
            n / nlen
        } else {
            Vector3r::zeros()
        };
        let c = centroid(&tri);

        // Skip near-degenerate seam slivers.
        //
        // # Theorem — Threshold Correctness
        //
        // `SLIVER_AREA_RATIO_SQ = 1e-14` skips faces whose altitude-to-longest-edge
        // ratio is below `sqrt(1e-14) = 1e-7`.  The minimum valid millifluidic ratio
        // is ≥ 50 µm / 4 mm × 0.5 correction ≈ 6.25e-3, safely above the threshold.
        // Numerically degenerate slivers from near-parallel face intersections have
        // ratios ~1e-10 to 1e-15, correctly below the threshold. ∎
        {
            if is_degenerate_sliver_with_normal(&tri, &n) {
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
                // Binary intersection can place one source fragment centroid in
                // the coplanar boundary band even when the opposite-source
                // companion fragment does not survive with identical topology.
                // Keeping `CoplanarSame` here prevents seam-adjacent barrel
                // strips from being dropped; exact duplicate faces are removed
                // later by the canonical finalization pass.
                BooleanOp::Intersection => (
                    class_a == FragmentClass::Inside || class_a == FragmentClass::CoplanarSame,
                    false,
                ),
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
    use crate::application::csg::arrangement::fragment_analysis::is_degenerate_sliver;
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

    /// Regression: SLIVER_AREA_RATIO_SQ = 1e-14 must not skip valid millifluidic faces.
    ///
    /// A face with a 4mm longest edge and a 50µm altitude has:
    ///   area_sq = (2 × 0.5 × 4e-3 × 50e-6)² = 4e-14
    ///   max_edge_sq ≈ (4e-3)² = 1.6e-5
    ///   ratio = 4e-14 / 1.6e-5 = 2.5e-9 >> SLIVER_AREA_RATIO_SQ (1e-14)
    ///
    /// The old threshold of 1e-10 would skip faces with altitude ratio < sqrt(1e-10) ≈ 3e-5,
    /// incorrectly eliminating thin-channel seam fragments from near-parallel intersections
    /// (which have altitude ratios of ~1e-6 to 1e-5).  The new threshold 1e-14 is safe. ∎
    #[test]
    fn sliver_area_ratio_keeps_high_aspect_millifluidic_face() {
        let tri = [
            Point3r::new(0.0, 0.0, 0.0),
            Point3r::new(4e-3, 0.0, 0.0),
            Point3r::new(0.0, 50e-6, 0.0),
        ];
        assert!(
            !is_degenerate_sliver(&tri),
            "canonical sliver filter must keep valid 80:1 millifluidic faces"
        );
    }

    /// Regression: genuinely degenerate slivers (altitude ratio < 1e-7) must still be skipped.
    #[test]
    fn sliver_area_ratio_skips_numerically_degenerate_slivers() {
        let tri = [
            Point3r::new(0.0, 0.0, 0.0),
            Point3r::new(1e-3, 0.0, 0.0),
            Point3r::new(0.0, 1e-11, 0.0),
        ];
        assert!(
            is_degenerate_sliver(&tri),
            "canonical sliver filter must skip altitude-ratio 1e-8 artifacts"
        );
    }
}
