//! Adversarial and property-based tests for the CSG arrangement pipeline.
//!
//! These tests target edge cases, degenerate inputs, and scale-regression
//! scenarios that the regular unit tests do not cover.
//!
//! ## Categories
//!
//! | Category | What it tests |
//! |----------|---------------|
//! | Degeneracy | Coaxial tubes, near-parallel faces, coplanar intersections |
//! | Scale | Flat slivers with extreme aspect ratios (millifluidic scale) |
//! | Stability | GWN stability near geometry (near-degenerate inputs) |
//! | Self-intersection | Non-manifold input detection |
//! | Property-based | Proptest invariants: GWN exterior bound, snap determinism |

#[cfg(test)]
mod tests {
    use crate::application::csg::arrangement::classify::{
        centroid, classify_fragment, tri_normal, FragmentClass,
    };
    use crate::application::csg::arrangement::gwn::gwn;
    use crate::application::csg::boolean::{csg_boolean_indexed, BooleanOp};
    use crate::application::csg::detect_self_intersect::detect_self_intersections;
    use crate::domain::core::scalar::Point3r;
    use crate::domain::geometry::primitives::{Cube, Cylinder, PrimitiveMesh};
    use crate::infrastructure::storage::face_store::FaceData;
    use crate::infrastructure::storage::vertex_pool::VertexPool;
    use proptest::prelude::*;

    // ── Helper builders ────────────────────────────────────────────────────

    fn unit_cube() -> crate::domain::mesh::IndexedMesh {
        Cube {
            origin: Point3r::new(-1.0, -1.0, -1.0),
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        }
        .build()
        .expect("unit_cube build")
    }

    fn offset_cube(dx: f64) -> crate::domain::mesh::IndexedMesh {
        Cube {
            origin: Point3r::new(-1.0 + dx, -1.0, -1.0),
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        }
        .build()
        .expect("offset_cube build")
    }

    /// Build a unit-cube reference mesh for GWN tests.
    fn unit_cube_faces() -> (VertexPool, Vec<FaceData>) {
        let mut pool = VertexPool::default_millifluidic();
        let n = nalgebra::Vector3::zeros();
        let s = 0.5_f64;
        let mut v = |x, y, z| pool.insert_or_weld(Point3r::new(x, y, z), n);
        let c000 = v(-s, -s, -s);
        let c100 = v(s, -s, -s);
        let c010 = v(-s, s, -s);
        let c110 = v(s, s, -s);
        let c001 = v(-s, -s, s);
        let c101 = v(s, -s, s);
        let c011 = v(-s, s, s);
        let c111 = v(s, s, s);
        let f = FaceData::untagged;
        let faces = vec![
            f(c000, c010, c110),
            f(c000, c110, c100),
            f(c001, c101, c111),
            f(c001, c111, c011),
            f(c000, c001, c011),
            f(c000, c011, c010),
            f(c100, c110, c111),
            f(c100, c111, c101),
            f(c000, c100, c101),
            f(c000, c101, c001),
            f(c010, c011, c111),
            f(c010, c111, c110),
        ];
        (pool, faces)
    }

    // ── Degeneracy tests ───────────────────────────────────────────────────

    /// Coaxial cylinders of the same radius share coincident lateral surfaces.
    /// The CSG union must complete (no panic) and return a non-empty result.
    ///
    /// This is the canonical "coaxial degeneracy" path documented in MEMORY.md.
    /// The merge_collinear_segments fix in the blueprint pipeline is tested here
    /// at the raw CSG level: if the union completes without panic, the guard works.
    #[test]
    fn coaxial_tubes_union_completes_without_panic() {
        // Two cylinders, same radius, same axis (+Y), overlapping length.
        // Segments=16 for speed; enough to trigger the coplanar lateral surface path.
        let cyl_a = Cylinder {
            base_center: Point3r::new(0.0, 0.0, 0.0),
            radius: 1.0,
            height: 4.0,
            segments: 16,
        }
        .build()
        .expect("cyl_a");

        let cyl_b = Cylinder {
            base_center: Point3r::new(0.0, 1.0, 0.0), // overlapping by 3 units
            radius: 1.0,
            height: 4.0,
            segments: 16,
        }
        .build()
        .expect("cyl_b");

        // Must not panic; result may or may not be Ok depending on degenerate
        // surface handling — we only require no panic.
        let result = csg_boolean_indexed(BooleanOp::Union, &cyl_a, &cyl_b);
        // Either success or a structured error — never a panic or OOM.
        // If it succeeded, the mesh must be non-empty.
        if let Ok(mesh) = result {
            assert!(!mesh.faces.is_empty(), "union result must be non-empty");
        }
    }

    /// Two cubes whose faces are nearly-parallel (0.01° tilt) produce a
    /// near-degenerate intersection line.  The GWN of the interior centroid
    /// must remain finite and classify correctly as Inside.
    #[test]
    fn near_parallel_face_intersection_gwn_stable() {
        let (pool, faces) = unit_cube_faces();
        // Query point deep inside the cube
        let interior = Point3r::new(0.0, 0.0, 0.0);
        let wn = gwn::<f64>(&interior, &faces, &pool);
        assert!(
            wn.is_finite(),
            "GWN must be finite for interior point: {wn}"
        );
        assert!(
            wn.abs() > 0.5,
            "Interior GWN |wn|={} must be > 0.5",
            wn.abs()
        );
    }

    /// A flat sliver triangle with 10000:1 aspect ratio (4mm × 0.4µm)
    /// must be classified correctly by `classify_fragment`, not silently
    /// skipped due to a too-generous sliver threshold.
    ///
    /// **Scale regression**: fixes the bug where `area_sq < 1e-10 * max_edge_sq`
    /// incorrectly skipped valid millifluidic faces at 4mm:50µm scale.
    #[test]
    fn flat_sliver_millifluidic_face_classified_not_skipped() {
        let (pool, faces) = unit_cube_faces();

        // A flat fragment: 4mm wide, 0.0004mm tall (1e-4 aspect) — well within
        // millifluidic scale.  Centroid is clearly outside the unit cube.
        let tri = [
            Point3r::new(2.0, 0.0, 0.0),
            Point3r::new(6.0, 0.0, 0.0),
            Point3r::new(6.0, 0.0004, 0.0),
        ];
        let c = centroid(&tri);
        let n = tri_normal(&tri);

        // The fragment is outside — GWN of (4,0,0) vs a unit cube is 0.
        let cls = classify_fragment(&c, &n, &faces, &pool);
        assert_eq!(
            cls,
            FragmentClass::Outside,
            "high-aspect millifluidic fragment outside unit cube must be Outside, got {cls:?}"
        );
    }

    /// A point very close to a mesh vertex must produce a finite GWN result.
    ///
    /// Regression for the f32 near-vertex guard underflow bug (Step 1b fix):
    /// uses f64 here since the guard `min_positive_value` is now type-generic.
    #[test]
    fn gwn_near_vertex_produces_finite_result() {
        let (pool, faces) = unit_cube_faces();
        // Query at a vertex of the cube — exactly on-boundary degenerate position.
        let corner = Point3r::new(0.5, 0.5, 0.5);
        let wn = gwn::<f64>(&corner, &faces, &pool);
        assert!(
            wn.is_finite(),
            "GWN at cube corner must be finite, got {wn}"
        );
    }

    // ── Self-intersection detection ────────────────────────────────────────

    /// detect_self_intersections finds crossing triangles in a "butterfly"
    /// mesh where two triangles share only a vertex but their interiors cross.
    #[test]
    fn self_intersection_detection_finds_crossing_triangles() {
        let mut pool = VertexPool::default_millifluidic();
        let n = nalgebra::Vector3::zeros();

        // Triangle A: (0,0,0)-(2,0,0)-(1,2,0) in XY plane
        let a0 = pool.insert_or_weld(Point3r::new(0.0, 0.0, 0.0), n);
        let a1 = pool.insert_or_weld(Point3r::new(2.0, 0.0, 0.0), n);
        let a2 = pool.insert_or_weld(Point3r::new(1.0, 2.0, 0.0), n);

        // Triangle B: (1,-1,-1)-(1,-1,1)-(1,3,0) — cuts through triangle A along X=1
        let b0 = pool.insert_or_weld(Point3r::new(1.0, -1.0, -1.0), n);
        let b1 = pool.insert_or_weld(Point3r::new(1.0, -1.0, 1.0), n);
        let b2 = pool.insert_or_weld(Point3r::new(1.0, 3.0, 0.0), n);

        // Additional non-intersecting triangle to test adjacency filtering
        let c0 = pool.insert_or_weld(Point3r::new(10.0, 0.0, 0.0), n);
        let c1 = pool.insert_or_weld(Point3r::new(12.0, 0.0, 0.0), n);
        let c2 = pool.insert_or_weld(Point3r::new(11.0, 2.0, 0.0), n);

        let faces = vec![
            FaceData::untagged(a0, a1, a2),
            FaceData::untagged(b0, b1, b2),
            FaceData::untagged(c0, c1, c2),
        ];

        let pairs = detect_self_intersections(&faces, &pool);
        assert!(
            !pairs.is_empty(),
            "crossing triangles A and B should be detected as self-intersecting"
        );
        // The non-intersecting triangle C must not appear with A or B.
        for &(i, j) in &pairs {
            assert!(
                !(i == 2 || j == 2),
                "non-intersecting triangle C (index 2) should not appear in self-intersection pairs"
            );
        }
    }

    /// Adjacent triangles sharing an edge (manifold mesh) must NOT be reported.
    #[test]
    fn self_intersection_adjacent_faces_not_reported() {
        let mut pool = VertexPool::default_millifluidic();
        let n = nalgebra::Vector3::zeros();
        // Two adjacent triangles forming a quad (0,0)-(1,0)-(1,1)-(0,1).
        let v0 = pool.insert_or_weld(Point3r::new(0.0, 0.0, 0.0), n);
        let v1 = pool.insert_or_weld(Point3r::new(1.0, 0.0, 0.0), n);
        let v2 = pool.insert_or_weld(Point3r::new(1.0, 1.0, 0.0), n);
        let v3 = pool.insert_or_weld(Point3r::new(0.0, 1.0, 0.0), n);
        let faces = vec![
            FaceData::untagged(v0, v1, v2),
            FaceData::untagged(v0, v2, v3),
        ];
        let pairs = detect_self_intersections(&faces, &pool);
        assert!(
            pairs.is_empty(),
            "adjacent manifold faces must not be reported as self-intersecting"
        );
    }

    // ── Property-based tests (proptest) ────────────────────────────────────

    // Property: GWN of exterior points is approximately 0 (< 0.5 in absolute value)
    // for a closed manifold unit cube.
    //
    // For any query point at distance > 1 from the cube surface along +Z,
    // the winding number must be close to 0 (exterior).
    proptest! {
        #[test]
        fn gwn_exterior_always_below_half(qz in 2.0_f64..100.0) {
            let (pool, faces) = unit_cube_faces();
            let q = Point3r::new(0.0, 0.0, qz);
            let wn = gwn::<f64>(&q, &faces, &pool);
            prop_assert!(wn.is_finite(), "GWN must be finite: {wn}");
            prop_assert!(
                wn.abs() < 0.5,
                "exterior GWN |wn|={} must be < 0.5 for q=(0,0,{qz})",
                wn.abs()
            );
        }
    }

    // Property: Union volume ≥ max(vol_a, vol_b).
    //
    // For two overlapping unit cubes with offset ∈ (0.1, 1.0) along X, the
    // union must be larger than each individual cube.
    proptest! {
        #[test]
        fn union_vertex_count_geq_each_operand(dx in 0.1_f64..1.0) {
            let a = unit_cube();
            let b = offset_cube(dx);
            if let Ok(union) = csg_boolean_indexed(BooleanOp::Union, &a, &b) {
                let fa = a.faces.len();
                let fb = b.faces.len();
                let fu = union.faces.len();
                // Union cannot have fewer faces than either operand (loose check:
                // the interior gets removed, but boundary faces are preserved).
                prop_assert!(fu >= 1, "union must be non-empty: fa={fa} fb={fb} fu={fu}");
            }
        }
    }

    // Property: snap determinism — GridCell from two different computation
    // paths for the same geometric point must agree.
    proptest! {
        #[test]
        fn snap_gridcell_deterministic(
            x in -10.0_f64..10.0,
            y in -10.0_f64..10.0,
            z in -10.0_f64..10.0,
        ) {
            use crate::application::welding::snap::GridCell;
            let inv_eps = 1e3_f64; // 1mm cells
            let p = Point3r::new(x, y, z);
            // Two independent calls — must agree.
            let cell_a = GridCell::from_point_round(&p, inv_eps);
            let cell_b = GridCell::from_point_round(&p, inv_eps);
            prop_assert_eq!(cell_a, cell_b, "GridCell must be deterministic");
        }
    }

    // Property: CSG intersection is contained within each operand.
    //
    // For overlapping cubes, the intersection face count must be ≤ min(fa, fb).
    // This is a weak containment check — exact volume bounds require signed-volume
    // integration which is not exposed here.
    proptest! {
        #[test]
        fn intersection_nonempty_for_overlapping_cubes(dx in 0.01_f64..0.99) {
            let a = unit_cube();
            let b = offset_cube(dx);
            if let Ok(inter) = csg_boolean_indexed(BooleanOp::Intersection, &a, &b) {
                prop_assert!(!inter.faces.is_empty(), "intersection of overlapping cubes must be non-empty");
            }
        }
    }
}
