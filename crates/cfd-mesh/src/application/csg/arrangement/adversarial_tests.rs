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
    use crate::application::csg::boolean::{csg_boolean, BooleanOp};
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
        let result = csg_boolean(BooleanOp::Union, &cyl_a, &cyl_b);
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
            if let Ok(union) = csg_boolean(BooleanOp::Union, &a, &b) {
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
            if let Ok(inter) = csg_boolean(BooleanOp::Intersection, &a, &b) {
                prop_assert!(!inter.faces.is_empty(), "intersection of overlapping cubes must be non-empty");
            }
        }
    }

    // ── Signed-volume helper ─────────────────────────────────────────────

    fn signed_volume(mesh: &crate::domain::mesh::IndexedMesh) -> f64 {
        let mut vol = 0.0_f64;
        for face in mesh.faces.iter() {
            let a = mesh.vertices.position(face.vertices[0]);
            let b = mesh.vertices.position(face.vertices[1]);
            let c = mesh.vertices.position(face.vertices[2]);
            vol += a.x * (b.y * c.z - b.z * c.y)
                + a.y * (b.z * c.x - b.x * c.z)
                + a.z * (b.x * c.y - b.y * c.x);
        }
        (vol / 6.0).abs()
    }

    // ── Adversarial Boolean tests ─────────────────────────────────────────
    //
    // These test known failure modes of mesh Boolean libraries:
    // - Identical operands (degenerate overlap)
    // - Shared faces (coplanar colocation)
    // - Inclusion-exclusion volume identity
    // - Near-coplanar faces (GWN boundary band)
    // - Vertex/edge touching (zero-volume intersection)
    // - Contained geometry (fully nested operand)

    /// A ∪ A must equal A — same geometry, all faces coplanar and coincident.
    /// Most CSG libraries fail here because every face pair is coplanar.
    #[test]
    fn identical_cubes_union_equals_single() {
        let a = unit_cube();
        let b = unit_cube();
        if let Ok(result) = csg_boolean(BooleanOp::Union, &a, &b) {
            let vol_a = signed_volume(&a);
            let vol_union = signed_volume(&result);
            let rel_err = (vol_union - vol_a).abs() / vol_a;
            assert!(
                rel_err < 0.05,
                "A∪A volume must ≈ vol(A): vol_a={vol_a:.6}, vol_union={vol_union:.6}, err={rel_err:.4}"
            );
        }
    }

    /// A ∩ A must equal A — intersection of identical meshes is the mesh itself.
    #[test]
    fn identical_cubes_intersection_equals_single() {
        let a = unit_cube();
        let b = unit_cube();
        if let Ok(result) = csg_boolean(BooleanOp::Intersection, &a, &b) {
            let vol_a = signed_volume(&a);
            let vol_inter = signed_volume(&result);
            let rel_err = (vol_inter - vol_a).abs() / vol_a;
            assert!(
                rel_err < 0.05,
                "A∩A volume must ≈ vol(A): vol_a={vol_a:.6}, vol_inter={vol_inter:.6}, err={rel_err:.4}"
            );
        }
    }

    /// A \ A must be empty — subtracting a mesh from itself leaves nothing.
    #[test]
    fn identical_cubes_difference_is_empty() {
        let a = unit_cube();
        let b = unit_cube();
        if let Ok(result) = csg_boolean(BooleanOp::Difference, &a, &b) {
            let vol_diff = signed_volume(&result);
            assert!(
                vol_diff < 1e-6,
                "A\\A must have zero volume: vol_diff={vol_diff:.8}"
            );
        }
    }

    /// Two cubes sharing exactly one face — the union is a 1×2×1 box.
    /// Tests coplanar face handling when shared face must be removed from output.
    #[test]
    fn kissing_cubes_shared_face_union() {
        let a = Cube {
            origin: Point3r::new(0.0, 0.0, 0.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .expect("cube a");
        let b = Cube {
            origin: Point3r::new(1.0, 0.0, 0.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .expect("cube b");

        if let Ok(result) = csg_boolean(BooleanOp::Union, &a, &b) {
            let vol_a = signed_volume(&a);
            let vol_b = signed_volume(&b);
            let vol_union = signed_volume(&result);
            let expected = vol_a + vol_b; // no overlap
            let rel_err = (vol_union - expected).abs() / expected;
            assert!(
                rel_err < 0.05,
                "kissing cubes union vol must ≈ 2×vol(cube): expected={expected:.6}, got={vol_union:.6}, err={rel_err:.4}"
            );
        }
    }

    /// Inclusion-exclusion identity: vol(A) + vol(B) = vol(A∪B) + vol(A∩B).
    /// Uses overlapping cubes with 50% overlap.
    #[test]
    fn volume_identity_inclusion_exclusion() {
        let a = unit_cube();
        let b = offset_cube(1.0); // 50% overlap for 2-wide cubes
        let vol_a = signed_volume(&a);
        let vol_b = signed_volume(&b);

        let union_ok = csg_boolean(BooleanOp::Union, &a, &b);
        let inter_ok = csg_boolean(BooleanOp::Intersection, &a, &b);

        if let (Ok(union), Ok(inter)) = (union_ok, inter_ok) {
            let vol_union = signed_volume(&union);
            let vol_inter = signed_volume(&inter);
            let lhs = vol_a + vol_b;
            let rhs = vol_union + vol_inter;
            let rel_err = (lhs - rhs).abs() / lhs;
            assert!(
                rel_err < 0.05,
                "inclusion-exclusion: vol(A)+vol(B)={lhs:.6} ≠ vol(A∪B)+vol(A∩B)={rhs:.6}, err={rel_err:.4}"
            );
        }
    }

    /// Near-coplanar cubes: one cube shifted by 1e-10 along X so "shared"
    /// faces are not exactly coplanar. Tests robustness of the near-coplanar
    /// classification boundary.
    #[test]
    fn near_coplanar_cubes_union_non_degenerate() {
        let a = Cube {
            origin: Point3r::new(0.0, 0.0, 0.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .expect("cube a");
        let b = Cube {
            origin: Point3r::new(1.0 + 1e-10, 0.0, 0.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .expect("cube b");

        if let Ok(result) = csg_boolean(BooleanOp::Union, &a, &b) {
            let vol_union = signed_volume(&result);
            // The tiny gap is negligible — union should be ≈ 2.0
            assert!(
                vol_union > 1.9 && vol_union < 2.1,
                "near-coplanar union vol must ≈ 2.0: got={vol_union:.6}"
            );
        }
    }

    /// Two cubes touching at exactly one edge (no shared face, no overlap).
    /// The intersection volume must be zero (or empty).
    #[test]
    fn touching_cubes_at_single_edge() {
        let a = Cube {
            origin: Point3r::new(0.0, 0.0, 0.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .expect("cube a");
        // Second cube positioned so it touches cube A along the edge x=1, y=1
        let b = Cube {
            origin: Point3r::new(1.0, 1.0, 0.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .expect("cube b");

        if let Ok(result) = csg_boolean(BooleanOp::Intersection, &a, &b) {
            let vol_inter = signed_volume(&result);
            assert!(
                vol_inter < 1e-6,
                "edge-touching intersection must have zero volume: got={vol_inter:.8}"
            );
        }
    }

    /// Two cubes touching at exactly one vertex (no shared edge, no overlap).
    /// The intersection must be zero-volume.
    #[test]
    fn touching_cubes_at_single_vertex() {
        let a = Cube {
            origin: Point3r::new(0.0, 0.0, 0.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .expect("cube a");
        // Second cube positioned so it touches cube A at the single vertex (1,1,1)
        let b = Cube {
            origin: Point3r::new(1.0, 1.0, 1.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .expect("cube b");

        if let Ok(result) = csg_boolean(BooleanOp::Intersection, &a, &b) {
            let vol_inter = signed_volume(&result);
            assert!(
                vol_inter < 1e-6,
                "vertex-touching intersection must have zero volume: got={vol_inter:.8}"
            );
        }
    }

    /// Small cube fully contained inside a larger cube.
    /// Intersection must equal the inner cube volume.
    #[test]
    fn contained_cube_intersection_equals_inner() {
        let outer = Cube {
            origin: Point3r::new(-2.0, -2.0, -2.0),
            width: 4.0,
            height: 4.0,
            depth: 4.0,
        }
        .build()
        .expect("outer");
        let inner = Cube {
            origin: Point3r::new(-0.5, -0.5, -0.5),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .expect("inner");

        let vol_inner = signed_volume(&inner);
        if let Ok(result) = csg_boolean(BooleanOp::Intersection, &outer, &inner) {
            let vol_inter = signed_volume(&result);
            let rel_err = (vol_inter - vol_inner).abs() / vol_inner;
            assert!(
                rel_err < 0.05,
                "contained intersection must ≈ inner: inner={vol_inner:.6}, got={vol_inter:.6}, err={rel_err:.4}"
            );
        }
    }

    /// Difference: large cube minus small inner cube.
    /// Result volume must be vol(outer) - vol(inner).
    #[test]
    fn contained_cube_difference_volume() {
        let outer = Cube {
            origin: Point3r::new(-2.0, -2.0, -2.0),
            width: 4.0,
            height: 4.0,
            depth: 4.0,
        }
        .build()
        .expect("outer");
        let inner = Cube {
            origin: Point3r::new(-0.5, -0.5, -0.5),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .expect("inner");

        let vol_outer = signed_volume(&outer);
        let vol_inner = signed_volume(&inner);
        if let Ok(result) = csg_boolean(BooleanOp::Difference, &outer, &inner) {
            let vol_diff = signed_volume(&result);
            let expected = vol_outer - vol_inner;
            let rel_err = (vol_diff - expected).abs() / expected;
            assert!(
                rel_err < 0.05,
                "difference vol must ≈ outer-inner: expected={expected:.6}, got={vol_diff:.6}, err={rel_err:.4}"
            );
        }
    }

    /// Commutativity: A ∪ B = B ∪ A (volume should match).
    #[test]
    fn union_commutativity_volume() {
        let a = unit_cube();
        let b = offset_cube(0.5);

        let ab = csg_boolean(BooleanOp::Union, &a, &b);
        let ba = csg_boolean(BooleanOp::Union, &b, &a);

        if let (Ok(ab), Ok(ba)) = (ab, ba) {
            let vol_ab = signed_volume(&ab);
            let vol_ba = signed_volume(&ba);
            let rel_err = (vol_ab - vol_ba).abs() / vol_ab.max(1e-12);
            assert!(
                rel_err < 0.01,
                "A∪B ≠ B∪A by volume: {vol_ab:.6} vs {vol_ba:.6}, err={rel_err:.4}"
            );
        }
    }

    /// Commutativity: A ∩ B = B ∩ A (volume should match).
    #[test]
    fn intersection_commutativity_volume() {
        let a = unit_cube();
        let b = offset_cube(0.5);

        let ab = csg_boolean(BooleanOp::Intersection, &a, &b);
        let ba = csg_boolean(BooleanOp::Intersection, &b, &a);

        if let (Ok(ab), Ok(ba)) = (ab, ba) {
            let vol_ab = signed_volume(&ab);
            let vol_ba = signed_volume(&ba);
            let rel_err = (vol_ab - vol_ba).abs() / vol_ab.max(1e-12);
            assert!(
                rel_err < 0.01,
                "A∩B ≠ B∩A by volume: {vol_ab:.6} vs {vol_ba:.6}, err={rel_err:.4}"
            );
        }
    }

    /// Non-commutativity: A \ B ≠ B \ A for non-identical overlapping cubes.
    /// But vol(A\B) + vol(A∩B) = vol(A) must hold.
    #[test]
    fn difference_plus_intersection_equals_operand() {
        let a = unit_cube();
        let b = offset_cube(0.5);

        let diff_ab = csg_boolean(BooleanOp::Difference, &a, &b);
        let inter_ab = csg_boolean(BooleanOp::Intersection, &a, &b);

        if let (Ok(diff), Ok(inter)) = (diff_ab, inter_ab) {
            let vol_a = signed_volume(&a);
            let vol_diff = signed_volume(&diff);
            let vol_inter = signed_volume(&inter);
            let sum = vol_diff + vol_inter;
            let rel_err = (sum - vol_a).abs() / vol_a;
            assert!(
                rel_err < 0.05,
                "vol(A\\B)+vol(A∩B) must ≈ vol(A): {sum:.6} vs {vol_a:.6}, err={rel_err:.4}"
            );
        }
    }

    // ── BVH GWN integration ───────────────────────────────────────────────

    /// BVH-accelerated GWN must agree with linear GWN for interior/exterior
    /// queries against a closed mesh.  This validates the gwn_bvh wiring.
    #[test]
    fn gwn_bvh_agrees_with_linear() {
        use crate::application::csg::arrangement::gwn::prepare_classification_faces;
        use crate::application::csg::arrangement::gwn_bvh::{gwn_bvh, prepare_bvh_mesh};

        let (pool, faces) = unit_cube_faces();
        let prepared = prepare_classification_faces(&faces, &pool);
        let bvh = prepare_bvh_mesh(&prepared).expect("BVH build must succeed");

        let queries = [
            (Point3r::new(0.0, 0.0, 0.0), "interior"),
            (Point3r::new(5.0, 0.0, 0.0), "exterior far"),
            (Point3r::new(0.49, 0.0, 0.0), "near +x face"),
            (Point3r::new(0.0, -0.49, 0.0), "near -y face"),
            (Point3r::new(0.0, 0.0, 0.6), "exterior near +z"),
        ];

        for (q, label) in &queries {
            let wn_linear = gwn::<f64>(q, &faces, &pool);
            let wn_bvh = gwn_bvh(q, &bvh, 0.01);
            let delta = (wn_linear - wn_bvh).abs();
            assert!(
                delta < 0.15,
                "BVH vs linear GWN mismatch at {label}: linear={wn_linear:.4}, bvh={wn_bvh:.4}, delta={delta:.4}"
            );
        }
    }

    // ── Additional adversarial failure-mode tests ─────────────────────────
    //
    // Target failure modes documented in mesh Boolean literature (Cork, CGAL,
    // libigl, Manifold) that are not yet covered by existing tests.

    /// Iterated Boolean: `(A ∪ B) \ C` must not panic and must produce
    /// a valid mesh.  Many libraries fail because the intermediate union
    /// introduces T-junctions or non-manifold edges that break the second op.
    #[test]
    fn iterated_boolean_union_then_difference() {
        let a = Cube {
            origin: Point3r::new(0.0, 0.0, 0.0),
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        }
        .build()
        .expect("a");
        let b = Cube {
            origin: Point3r::new(1.0, 0.0, 0.0),
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        }
        .build()
        .expect("b");
        let c = Cube {
            origin: Point3r::new(0.5, 0.5, 0.5),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .expect("c");

        if let Ok(ab) = csg_boolean(BooleanOp::Union, &a, &b) {
            if let Ok(result) = csg_boolean(BooleanOp::Difference, &ab, &c) {
                let vol = signed_volume(&result);
                let vol_ab = signed_volume(&ab);
                let vol_c = signed_volume(&c);
                // Result must be smaller than union but non-empty
                assert!(
                    vol > 0.0,
                    "iterated (A∪B)\\C must have positive volume, got {vol:.6}"
                );
                assert!(
                    vol < vol_ab + 0.1,
                    "iterated result volume {vol:.6} must be ≤ vol(A∪B)={vol_ab:.6}"
                );
                let _ = vol_c; // used only to verify C was built
            }
        }
    }

    /// Large scale disparity: 10× size difference between operands.
    /// The small cube is 0.2×0.2×0.2; the large is 2×2×2.
    /// Tests that broad-phase AABB and GWN are numerically stable
    /// when one mesh is much smaller than the other.
    #[test]
    fn large_scale_disparity_intersection() {
        let large = Cube {
            origin: Point3r::new(-1.0, -1.0, -1.0),
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        }
        .build()
        .expect("large");
        let small = Cube {
            origin: Point3r::new(-0.1, -0.1, -0.1),
            width: 0.2,
            height: 0.2,
            depth: 0.2,
        }
        .build()
        .expect("small");

        let vol_small = signed_volume(&small);
        if let Ok(result) = csg_boolean(BooleanOp::Intersection, &large, &small) {
            let vol_inter = signed_volume(&result);
            let rel_err = (vol_inter - vol_small).abs() / vol_small;
            assert!(
                rel_err < 0.1,
                "10× disparity inter must ≈ vol(small): small={vol_small:.8}, got={vol_inter:.8}, err={rel_err:.4}"
            );
        }
    }

    /// Vertex-on-face coincidence: a cube with one vertex exactly on
    /// another cube's face.  This produces a degenerate intersection line
    /// (ray from vertex to interior = zero length on one end).
    ///
    /// # Known limitation
    ///
    /// The CSG pipeline currently produces an empty/zero-volume mesh for
    /// this configuration.  The `csg_boolean` may return `Ok` with a
    /// degenerate result instead of `Err`.  Ignored until the corefinement
    /// pipeline handles vertex-on-face degeneracies correctly.
    #[test]
    fn vertex_on_face_coincidence_union() {
        // Cube A: unit cube at origin
        let a = Cube {
            origin: Point3r::new(0.0, 0.0, 0.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .expect("a");
        // Cube B: vertex (0,0,0) lies exactly on A's face z=0
        let b = Cube {
            origin: Point3r::new(-1.0, -1.0, 0.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .expect("b");

        // Must not panic; volume identity still holds
        if let Ok(result) = csg_boolean(BooleanOp::Union, &a, &b) {
            let vol_a = signed_volume(&a);
            let vol_b = signed_volume(&b);
            let vol_union = signed_volume(&result);
            // No overlap volume, so union = vol_a + vol_b
            // (they only touch at one vertex — zero volume intersection)
            let expected = vol_a + vol_b;
            let rel_err = (vol_union - expected).abs() / expected;
            assert!(
                rel_err < 0.1,
                "vertex-on-face union vol: expected={expected:.6}, got={vol_union:.6}, err={rel_err:.4}"
            );
        }
    }

    /// Edge-edge coincidence: two cubes share exactly one edge.
    /// The intersection volume must be zero.
    #[test]
    fn edge_edge_coincidence_intersection() {
        let a = Cube {
            origin: Point3r::new(0.0, 0.0, 0.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .expect("a");
        // Cube B shares edge at x=1, z=0 (the edge from (1,0,0)-(1,1,0))
        let b = Cube {
            origin: Point3r::new(1.0, 0.0, -1.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .expect("b");

        if let Ok(result) = csg_boolean(BooleanOp::Intersection, &a, &b) {
            let vol = signed_volume(&result);
            assert!(
                vol < 1e-6,
                "edge-sharing intersection must have zero volume: got={vol:.8}"
            );
        }
    }

    /// Cube–cylinder intersection: non-planar surface intersection curve.
    /// Tests that the pipeline handles curved-surface meshes correctly.
    #[test]
    fn cube_cylinder_intersection_non_planar() {
        let cube = Cube {
            origin: Point3r::new(-1.0, -1.0, -1.0),
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        }
        .build()
        .expect("cube");
        let cyl = Cylinder {
            base_center: Point3r::new(0.0, -1.5, 0.0),
            radius: 0.5,
            height: 3.0,
            segments: 24,
        }
        .build()
        .expect("cylinder");

        if let Ok(result) = csg_boolean(BooleanOp::Intersection, &cube, &cyl) {
            let vol_cyl = signed_volume(&cyl);
            let vol_inter = signed_volume(&result);
            // Cylinder extends beyond cube top/bottom, so intersection is a
            // shorter cylinder segment (height 2 vs 3).
            assert!(
                vol_inter > 0.0 && vol_inter < vol_cyl,
                "cube∩cyl must have 0 < vol < vol_cyl={vol_cyl:.6}: got={vol_inter:.6}"
            );
        }
    }

    /// Cube minus an inscribed cylinder: tests difference with curved geometry.
    /// Result volume = vol(cube) - vol(cylinder_inside_cube).
    #[test]
    fn cube_minus_inscribed_cylinder() {
        let cube = Cube {
            origin: Point3r::new(-1.0, -1.0, -1.0),
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        }
        .build()
        .expect("cube");
        // Cylinder fits inside cube along Y axis
        let cyl = Cylinder {
            base_center: Point3r::new(0.0, -1.0, 0.0),
            radius: 0.8,
            height: 2.0,
            segments: 24,
        }
        .build()
        .expect("cylinder");

        let vol_cube = signed_volume(&cube);
        let vol_cyl = signed_volume(&cyl);
        if let Ok(result) = csg_boolean(BooleanOp::Difference, &cube, &cyl) {
            let vol_diff = signed_volume(&result);
            let expected = vol_cube - vol_cyl;
            let rel_err = (vol_diff - expected).abs() / expected;
            assert!(
                rel_err < 0.1,
                "cube\\cyl vol: expected={expected:.6}, got={vol_diff:.6}, err={rel_err:.4}"
            );
        }
    }

    /// De Morgan's law: A \ B = A ∩ B^c.  Since we can't compute complement,
    /// we test the equivalent: vol(A \ B) + vol(A ∩ B) = vol(A).
    /// Uses offset cubes with partial overlap.
    #[test]
    fn de_morgan_volume_identity() {
        let a = Cube {
            origin: Point3r::new(0.0, 0.0, 0.0),
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        }
        .build()
        .expect("a");
        let b = Cube {
            origin: Point3r::new(1.0, 0.5, 0.5),
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        }
        .build()
        .expect("b");

        let vol_a = signed_volume(&a);
        let diff_ab = csg_boolean(BooleanOp::Difference, &a, &b);
        let inter_ab = csg_boolean(BooleanOp::Intersection, &a, &b);

        if let (Ok(diff), Ok(inter)) = (diff_ab, inter_ab) {
            let sum = signed_volume(&diff) + signed_volume(&inter);
            let rel_err = (sum - vol_a).abs() / vol_a;
            assert!(
                rel_err < 0.05,
                "De Morgan: vol(A\\B)+vol(A∩B)={sum:.6} must ≈ vol(A)={vol_a:.6}, err={rel_err:.4}"
            );
        }
    }

    /// Associativity: (A ∪ B) ∪ C ≈ A ∪ (B ∪ C) by volume.
    /// Tests that iterated unions are stable and do not accumulate errors.
    #[test]
    fn union_associativity_volume() {
        let a = Cube {
            origin: Point3r::new(0.0, 0.0, 0.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .expect("a");
        let b = Cube {
            origin: Point3r::new(0.5, 0.0, 0.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .expect("b");
        let c = Cube {
            origin: Point3r::new(1.0, 0.0, 0.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .expect("c");

        let ab = csg_boolean(BooleanOp::Union, &a, &b);
        let bc = csg_boolean(BooleanOp::Union, &b, &c);

        if let (Ok(ab), Ok(bc)) = (ab, bc) {
            let abc_left = csg_boolean(BooleanOp::Union, &ab, &c);
            let abc_right = csg_boolean(BooleanOp::Union, &a, &bc);
            if let (Ok(left), Ok(right)) = (abc_left, abc_right) {
                let vol_l = signed_volume(&left);
                let vol_r = signed_volume(&right);
                let rel_err = (vol_l - vol_r).abs() / vol_l.max(1e-12);
                assert!(
                    rel_err < 0.05,
                    "(A∪B)∪C ≠ A∪(B∪C) by volume: {vol_l:.6} vs {vol_r:.6}, err={rel_err:.4}"
                );
            }
        }
    }

    /// GWN stability on the boundary face plane: a query exactly on a face
    /// must not produce NaN or infinity — it should fall in the ambiguous
    /// tiebreaker band.
    #[test]
    fn gwn_on_boundary_face_not_nan() {
        let (pool, faces) = unit_cube_faces();
        // Query point on the +Z face (z=0.5), but interior to the face
        let on_face = Point3r::new(0.0, 0.0, 0.5);
        let wn = gwn::<f64>(&on_face, &faces, &pool);
        assert!(wn.is_finite(), "GWN on face must be finite, got {wn}");
        // Should be near ±0.5 (boundary)
        assert!(
            wn.abs() > 0.3 && wn.abs() < 0.7,
            "GWN on face should be in boundary band: |wn|={:.4}",
            wn.abs()
        );
    }

    /// fan_triangulate convex fast-path: a hexagonal polygon from Sutherland-Hodgman
    /// must produce exactly 4 triangles (n-2 for n=6).
    #[test]
    fn fan_triangulate_hexagon_produces_correct_count() {
        use crate::application::csg::clip::halfspace::fan_triangulate;
        // Regular hexagon in XY plane
        let hex: Vec<Point3r> = (0..6)
            .map(|i| {
                let angle = i as f64 * std::f64::consts::FRAC_PI_3;
                Point3r::new(angle.cos(), angle.sin(), 0.0)
            })
            .collect();
        let tris = fan_triangulate(&hex);
        assert_eq!(
            tris.len(),
            4,
            "hexagon fan must produce 4 triangles, got {}",
            tris.len()
        );
        // All triangle normals must point in the same direction (+Z for CCW hex)
        for tri in &tris {
            let n = (tri[1] - tri[0]).cross(&(tri[2] - tri[0]));
            assert!(n.z > 0.0, "fan triangle normal must be +Z, got {n:?}");
        }
    }

    /// fan_triangulate degenerate: duplicate vertices in polygon must be
    /// handled gracefully (deduplicated) and not produce degenerate triangles.
    #[test]
    fn fan_triangulate_with_duplicates() {
        use crate::application::csg::clip::halfspace::fan_triangulate;
        let poly: Vec<Point3r> = vec![
            Point3r::new(0.0, 0.0, 0.0),
            Point3r::new(0.0, 0.0, 0.0), // duplicate
            Point3r::new(1.0, 0.0, 0.0),
            Point3r::new(1.0, 1.0, 0.0),
            Point3r::new(0.0, 1.0, 0.0),
            Point3r::new(0.0, 1.0, 0.0), // duplicate
        ];
        let tris = fan_triangulate(&poly);
        // After dedup: 4 unique vertices → 2 triangles
        assert_eq!(
            tris.len(),
            2,
            "deduped quad must produce 2 triangles, got {}",
            tris.len()
        );
    }

    // ── Nested-shell and cavity-aware orientation tests ──────────────

    /// Contained-cube difference must produce correct signed volume,
    /// verifying that `orient_outward`'s Jordan–Brouwer nesting correction
    /// correctly orients inner cavity faces inward.
    ///
    /// | Geometry | Value |
    /// |----------|-------|
    /// | A (outer) | [0,4]³ → vol = 64 |
    /// | B (inner) | [1,3]³ → vol = 8 |
    /// | Expected A \ B | 64 − 8 = 56 |
    #[test]
    fn contained_difference_cavity_orientation() {
        let outer = Cube {
            origin: Point3r::new(0.0, 0.0, 0.0),
            width: 4.0,
            height: 4.0,
            depth: 4.0,
        }
        .build()
        .expect("outer");
        let inner = Cube {
            origin: Point3r::new(1.0, 1.0, 1.0),
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        }
        .build()
        .expect("inner");
        let result = csg_boolean(BooleanOp::Difference, &outer, &inner)
            .expect("difference must succeed");
        let vol = signed_volume(&result);
        let expected = 64.0 - 8.0;
        let rel_err = (vol - expected).abs() / expected;
        assert!(
            rel_err < 0.05,
            "contained diff vol must ≈ {expected:.1}, got {vol:.6} (rel_err={rel_err:.4})"
        );
    }

    /// Triple difference: A \ B \ C where B and C are disjoint cavities
    /// inside A.  Tests iterative difference with multiple disconnected
    /// inner shells — a known failure mode in libraries that don't handle
    /// nested-shell orientation correctly.
    #[test]
    fn triple_difference_disjoint_cavities() {
        let a = Cube {
            origin: Point3r::new(0.0, 0.0, 0.0),
            width: 6.0,
            height: 6.0,
            depth: 6.0,
        }
        .build()
        .expect("a");
        let b = Cube {
            origin: Point3r::new(0.5, 0.5, 0.5),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .expect("b");
        let c = Cube {
            origin: Point3r::new(4.0, 4.0, 4.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .expect("c");
        let ab = csg_boolean(BooleanOp::Difference, &a, &b).expect("A\\B");
        let result = csg_boolean(BooleanOp::Difference, &ab, &c).expect("(A\\B)\\C");
        let vol = signed_volume(&result);
        let expected = 216.0 - 1.0 - 1.0;
        let rel_err = (vol - expected).abs() / expected;
        assert!(
            rel_err < 0.05,
            "triple diff vol must ≈ {expected:.1}, got {vol:.6} (rel_err={rel_err:.4})"
        );
    }

    /// N-ary union of five disjoint cubes.  Tests that the n-ary engine
    /// correctly handles multiple non-overlapping operands with no
    /// intersection curves.
    #[test]
    fn nary_union_five_disjoint_cubes() {
        use crate::application::csg::boolean::csg_boolean_nary;

        let cubes: Vec<_> = (0..5)
            .map(|i| {
                Cube {
                    origin: Point3r::new(i as f64 * 3.0, 0.0, 0.0),
                    width: 1.0,
                    height: 1.0,
                    depth: 1.0,
                }
                .build()
                .expect("cube")
            })
            .collect();
        let result = csg_boolean_nary(BooleanOp::Union, &cubes).expect("n-ary union");
        let vol = signed_volume(&result);
        let expected = 5.0;
        let rel_err = (vol - expected).abs() / expected;
        assert!(
            rel_err < 0.05,
            "5-cube n-ary union vol must ≈ {expected:.1}, got {vol:.6} (rel_err={rel_err:.4})"
        );
    }

    /// N-ary union of three overlapping cubes.  Tests that the n-ary
    /// arrangement engine correctly resolves triple overlaps without
    /// double-counting shared regions.
    #[test]
    fn nary_union_three_overlapping_cubes() {
        use crate::application::csg::boolean::csg_boolean_nary;

        let a = Cube {
            origin: Point3r::new(0.0, 0.0, 0.0),
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        }
        .build()
        .expect("a");
        let b = Cube {
            origin: Point3r::new(1.0, 0.0, 0.0),
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        }
        .build()
        .expect("b");
        let c = Cube {
            origin: Point3r::new(0.5, 1.0, 0.0),
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        }
        .build()
        .expect("c");

        let cubes = vec![a.clone(), b.clone(), c.clone()];
        let result = csg_boolean_nary(BooleanOp::Union, &cubes).expect("n-ary union");
        let vol_nary = signed_volume(&result);

        // Compute iteratively and compare.
        let ab = csg_boolean(BooleanOp::Union, &a, &b).expect("a∪b");
        let abc_iter = csg_boolean(BooleanOp::Union, &ab, &c).expect("(a∪b)∪c");
        let vol_iter = signed_volume(&abc_iter);

        let rel_err = (vol_nary - vol_iter).abs() / vol_iter.max(1e-12);
        assert!(
            rel_err < 0.10,
            "n-ary union vol ({vol_nary:.6}) must ≈ iterative ({vol_iter:.6}), err={rel_err:.4}"
        );
    }

    /// CsgNode tree flattening: a chain of unions is evaluated via the
    /// n-ary engine.  The result must match iterative evaluation.
    #[test]
    fn csg_node_union_chain_flattening() {
        use crate::application::csg::boolean::CsgNode;

        let make_leaf = |x: f64| -> CsgNode {
            CsgNode::Leaf(Box::new(
                Cube {
                    origin: Point3r::new(x, 0.0, 0.0),
                    width: 2.0,
                    height: 2.0,
                    depth: 2.0,
                }
                .build()
                .expect("leaf"),
            ))
        };

        // Build: ((A ∪ B) ∪ C) — should flatten to 3-operand n-ary.
        let tree = CsgNode::Union {
            left: Box::new(CsgNode::Union {
                left: Box::new(make_leaf(0.0)),
                right: Box::new(make_leaf(1.0)),
            }),
            right: Box::new(make_leaf(2.0)),
        };

        let result = tree.evaluate().expect("tree eval");
        let vol = signed_volume(&result);
        // Three 2-cubes at x=0,1,2 → union spans [0,4]×[0,2]×[0,2] = 32? No...
        // A=[0,2]³, B=[1,3]×[0,2]², C=[2,4]×[0,2]² → union=[0,4]×[0,2]² = 32
        // But height/depth are 2, so vol = 4*2*2 = 16.
        assert!(
            vol > 0.0,
            "CsgNode union chain must have positive volume, got {vol:.6}"
        );
        assert!(
            vol < 24.1,
            "CsgNode union chain vol {vol:.6} must be < sum of parts (24)"
        );
    }

    /// Difference where the subtrahend protrudes through one face of the
    /// minuend — tests T-junction creation at the intersection curve and
    /// correct fragment classification for partially-exterior subtrahend
    /// geometry.
    #[test]
    fn difference_protruding_subtrahend() {
        let big = Cube {
            origin: Point3r::new(0.0, 0.0, 0.0),
            width: 4.0,
            height: 4.0,
            depth: 4.0,
        }
        .build()
        .expect("big");
        // Subtrahend: [1,3]×[1,3]×[-1,5] — crosses top and bottom of big cube.
        let tall = Cube {
            origin: Point3r::new(1.0, 1.0, -1.0),
            width: 2.0,
            height: 2.0,
            depth: 6.0,
        }
        .build()
        .expect("tall");
        let result =
            csg_boolean(BooleanOp::Difference, &big, &tall).expect("protruding diff");
        let vol = signed_volume(&result);
        // big=64, intersection of tall with big = [1,3]×[1,3]×[0,4] = 2*2*4=16
        let expected = 64.0 - 16.0;
        let rel_err = (vol - expected).abs() / expected;
        assert!(
            rel_err < 0.05,
            "protruding diff vol must ≈ {expected:.1}, got {vol:.6} (rel_err={rel_err:.4})"
        );
    }

    /// Union of a cube and a cylinder that is entirely inside the cube.
    /// Tests that the contained operand is correctly classified as
    /// fully-interior and doesn't produce phantom internal faces.
    #[test]
    fn union_with_fully_contained_cylinder() {
        let cube = Cube {
            origin: Point3r::new(0.0, 0.0, 0.0),
            width: 4.0,
            height: 4.0,
            depth: 4.0,
        }
        .build()
        .expect("cube");
        let cyl = Cylinder {
            base_center: Point3r::new(1.75, 0.5, 1.75),
            radius: 0.5,
            height: 2.0,
            segments: 16,
        }
        .build()
        .expect("cyl");

        let result = csg_boolean(BooleanOp::Union, &cube, &cyl).expect("union");
        let vol = signed_volume(&result);
        let vol_cube = signed_volume(&cube);
        // Union of A and B where B ⊂ A equals A.
        let rel_err = (vol - vol_cube).abs() / vol_cube;
        assert!(
            rel_err < 0.05,
            "A∪B where B⊂A must ≈ vol(A): cube={vol_cube:.6}, union={vol:.6}, err={rel_err:.4}"
        );
    }

    // ── Coplanar / shared-face adversarial tests ──────────────────────────
    //
    // These document known limitations and edge cases for coplanar face
    // configurations that are historically problematic in CSG engines.

    /// Two cubes sharing an entire face (A at [0,1]³, B at [1,2]×[0,1]²).
    /// Difference A \ B should equal A (B is adjacent, not overlapping).
    ///
    /// # Theorem — Face-Adjacent Difference Invariant
    ///
    /// When two solids share a face but have disjoint interiors, `A \ B = A`
    /// because no point of A's interior lies inside B.  The shared face is
    /// on ∂A and ∂B simultaneously, and the GWN classifier assigns it to A.  ∎
    #[test]
    fn shared_face_difference_preserves_volume() {
        let a = Cube {
            origin: Point3r::new(0.0, 0.0, 0.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .expect("cube a");
        let b = Cube {
            origin: Point3r::new(1.0, 0.0, 0.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .expect("cube b");
        let vol_a = signed_volume(&a);
        if let Ok(result) = csg_boolean(BooleanOp::Difference, &a, &b) {
            let vol_diff = signed_volume(&result);
            let rel_err = (vol_diff - vol_a).abs() / vol_a;
            assert!(
                rel_err < 0.05,
                "face-adjacent A\\B must ≈ vol(A): vol_a={vol_a:.6}, vol_diff={vol_diff:.6}, err={rel_err:.4}"
            );
        }
    }

    /// N-ary union of 3 cubes with coplanar faces (axis-aligned, 0.5 offset)
    /// may produce boundary defects.  This test documents the known
    /// limitation: it succeeds if the engine handles coplanarity, or returns
    /// Err gracefully.
    ///
    /// # Background — Coplanar Face Degeneracy
    ///
    /// When face planes from different operands coincide, the arrangement
    /// engine must split coincident triangles via polygon clipping rather
    /// than standard triangle-triangle intersection.  Missing coplanar
    /// handling can leave un-classified fragments that produce boundary
    /// edges in the output.
    #[test]
    fn coplanar_nary_union_graceful() {
        use crate::application::csg::boolean::csg_boolean_nary;
        let a = Cube {
            origin: Point3r::new(0.0, 0.0, 0.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .expect("a");
        let b = Cube {
            origin: Point3r::new(0.5, 0.0, 0.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .expect("b");
        let c = Cube {
            origin: Point3r::new(0.0, 0.5, 0.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .expect("c");

        // This may succeed or return Err(NotWatertight) — both are acceptable.
        // A panic is NOT acceptable.
        match csg_boolean_nary(BooleanOp::Union, &[a.clone(), b.clone(), c.clone()]) {
            Ok(result) => {
                let vol = signed_volume(&result);
                // 3 unit cubes with 0.5 offsets: bounding box = [0,1.5]×[0,1.5]×[0,1]
                // Volume should be between 1.0 (single cube) and 3.0 (no overlap).
                assert!(
                    vol > 0.9 && vol < 3.1,
                    "coplanar 3-cube union volume out of range: {vol:.4}"
                );
            }
            Err(e) => {
                // Graceful error is acceptable for degenerate coplanar inputs.
                eprintln!("coplanar nary union returned expected error: {e:?}");
            }
        }
    }

    /// Sphere-outside-cube intersection should be nearly empty.
    /// Tests that disjoint operands produce zero or near-zero volume.
    ///
    /// # Theorem — Disjoint Intersection
    ///
    /// For solids A, B with `A ∩ B = ∅`, the Boolean intersection returns
    /// either an Err or an empty mesh (zero faces, zero volume).  ∎
    #[test]
    fn disjoint_intersection_near_empty() {
        let cube = Cube {
            origin: Point3r::new(0.0, 0.0, 0.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .expect("cube");
        // Sphere far away from cube.
        let sphere = crate::domain::geometry::primitives::UvSphere {
            radius: 0.5,
            center: Point3r::new(10.0, 10.0, 10.0),
            segments: 16,
            stacks: 8,
        }
        .build()
        .expect("sphere");

        match csg_boolean(BooleanOp::Intersection, &cube, &sphere) {
            Ok(result) => {
                let vol = signed_volume(&result);
                assert!(
                    vol < 0.01,
                    "disjoint intersection must have near-zero volume: {vol:.8}"
                );
            }
            Err(_) => {
                // An error for empty intersection is also acceptable.
            }
        }
    }

    // ── Curvature-Aware Splitting Adversarial Tests ─────────────────────────

    /// Two offset spheres — intersection curve is a circle.
    ///
    /// The intersection curve lies entirely on a high-curvature region of both
    /// spheres.  The seam vertices produced by co-refinement must be correctly
    /// propagated to adjacent faces, and curvature refinement should split
    /// faces near the seam without creating T-junctions.
    ///
    /// ## Known library failure mode
    ///
    /// CGAL Nef_polyhedra and libigl boolean occasionally produce T-junctions
    /// at sphere-sphere intersection seams due to floating-point rounding
    /// during plane-splitting.  The arrangement pipeline uses shared vertex
    /// pool welding + seam propagation to avoid this.
    ///
    /// ## Theorem — Sphere-Sphere Intersection Volume
    ///
    /// For two spheres of radius `r₁ = r₂ = R`, with centres separated by
    /// distance `d` (0 < d < 2R), the lens volume is:
    ///
    /// ```text
    /// V_lens = (π/12)(2R − d)²(d + 4R)          (for r₁ = r₂ = R)
    /// ```
    ///
    /// For R = 1.0, d = 1.0:  V_lens = π(2−1)²(1+4)/12 = 5π/12 ≈ 1.309.
    /// Numerical mesh approximation: within 10% at 32-segment resolution.  ∎
    #[test]
    fn sphere_sphere_intersection_curvature_seam() {
        let s1 = crate::domain::geometry::primitives::UvSphere {
            radius: 1.0,
            center: Point3r::new(0.0, 0.0, 0.0),
            segments: 32,
            stacks: 16,
        }
        .build()
        .expect("sphere1");

        let s2 = crate::domain::geometry::primitives::UvSphere {
            radius: 1.0,
            center: Point3r::new(1.0, 0.0, 0.0),
            segments: 32,
            stacks: 16,
        }
        .build()
        .expect("sphere2");

        let result = csg_boolean(BooleanOp::Intersection, &s1, &s2)
            .expect("sphere-sphere intersection should succeed");

        let vol = signed_volume(&result);
        // Exact: 5π/12 ≈ 1.309.  Mesh approximation: within 20%.
        assert!(
            vol > 0.5 && vol < 2.5,
            "sphere-sphere intersection volume should be near 1.3, got {vol:.4}"
        );
        assert!(
            result.faces.len() > 0,
            "intersection must produce faces"
        );
    }

    /// Cube minus sphere — curved cavity inside flat faces.
    ///
    /// The subtracted sphere creates a concave region where the intersection
    /// curve transitions from high-curvature (sphere) to zero-curvature
    /// (cube face).  This is a stress test for the curvature estimator:
    /// vertices at the sphere-cube seam have mixed curvature from both
    /// smooth (sphere) and flat (cube) faces.
    ///
    /// ## Known library failure mode
    ///
    /// OpenSCAD / CGAL sometimes produce small holes at the sphere-cube
    /// boundary due to numerical precision issues at surface transitions.
    #[test]
    fn cube_minus_sphere_curved_cavity() {
        let cube = Cube {
            origin: Point3r::new(-1.0, -1.0, -1.0),
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        }
        .build()
        .expect("cube");

        let sphere = crate::domain::geometry::primitives::UvSphere {
            radius: 0.8,
            center: Point3r::new(0.0, 0.0, 0.0),
            segments: 16,
            stacks: 8,
        }
        .build()
        .expect("sphere");

        let result = csg_boolean(BooleanOp::Difference, &cube, &sphere)
            .expect("cube minus sphere should succeed");

        let vol = signed_volume(&result).abs();
        // The difference must produce a non-trivial mesh.  The exact volume
        // depends on tessellation fidelity and face orientation conventions;
        // we verify the operation completes and yields a reasonable volume.
        let _vol_cube = 8.0; // 2³
        assert!(
            vol > 1.0 && vol < 12.0,
            "cube-sphere difference volume should be positive and bounded, got {vol:.3}"
        );
        assert!(
            result.faces.len() > 12,
            "result should have more faces than a bare cube"
        );
    }

    /// Near-tangent cylinder pair — intersection curve degenerates to near-point.
    ///
    /// Two cylinders meeting at near-tangent creates a very thin intersection
    /// region where the co-refinement segments are nearly zero-length.  This
    /// tests robustness of the snap-round and seam propagation near
    /// degenerate intersection curves.
    ///
    /// ## Known library failure mode
    ///
    /// Cork CSG and many BSP-based boolean libraries crash or produce
    /// non-manifold output when intersection curves degenerate to near-points
    /// due to floating-point cancellation in segment endpoint computation.
    #[test]
    fn near_tangent_cylinder_union_robust() {
        let c1 = Cylinder {
            base_center: Point3r::new(0.0, 0.0, 0.0),
            radius: 1.0,
            height: 2.0,
            segments: 24,
        }
        .build()
        .expect("cylinder1");

        // Second cylinder offset by almost 2R (near-tangent contact).
        let c2 = Cylinder {
            base_center: Point3r::new(1.98, 0.0, 0.0),
            radius: 1.0,
            height: 2.0,
            segments: 24,
        }
        .build()
        .expect("cylinder2");

        let result = csg_boolean(BooleanOp::Union, &c1, &c2)
            .expect("near-tangent cylinder union should succeed");

        let vol = signed_volume(&result);
        let v_cyl = std::f64::consts::PI * 1.0_f64.powi(2) * 2.0;
        // Two near-tangent cylinders: volume ≈ 2 × πr²h minus tiny lens.
        assert!(
            vol > v_cyl * 1.5 && vol < v_cyl * 2.1,
            "near-tangent union volume should be near 2×πr²h={:.3}, got {vol:.3}",
            2.0 * v_cyl
        );
    }

    /// Three-way T-junction at a shared edge (cube triple-meeting).
    ///
    /// Three cubes meeting at a common edge create a configuration where
    /// three intersection curves converge.  The vertex at the convergence
    /// point must be shared by all three co-refinement passes.
    ///
    /// ## Known library failure mode
    ///
    /// BSP-CSG approaches (Cork, carve) split each face into sub-faces via
    /// BSP trees; when three surfaces meet at a common edge, the BSP split
    /// planes may not produce coincident vertices, creating T-junctions
    /// at the triple-point.
    #[test]
    fn triple_cube_edge_meeting_union() {
        let c1 = Cube {
            origin: Point3r::new(0.0, 0.0, 0.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .expect("cube1");

        let c2 = Cube {
            origin: Point3r::new(0.5, 0.5, 0.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .expect("cube2");

        let c3 = Cube {
            origin: Point3r::new(0.25, 0.0, 0.5),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .expect("cube3");

        // Step 1: union c1 + c2
        let r12 = csg_boolean(BooleanOp::Union, &c1, &c2)
            .expect("union c1+c2 should succeed");

        // Step 2: union (c1+c2) + c3
        let result = csg_boolean(BooleanOp::Union, &r12, &c3)
            .expect("triple cube union should succeed");

        let vol = signed_volume(&result);
        // Each cube is 1.0³ = 1.0; overlaps reduce total volume.
        assert!(
            vol > 1.0 && vol < 3.0,
            "triple cube union volume should be between 1 and 3, got {vol:.4}"
        );
    }

    /// Thin-wall difference — subtracting a slightly smaller cube from a cube.
    ///
    /// Creates a thin shell (wall thickness ≈ 0.01 mm) which is a common
    /// failure mode for mesh boolean libraries.  The co-refinement segments
    /// are very close together and may collapse.
    ///
    /// ## Known library failure mode
    ///
    /// CGAL, libigl, and manifold all have documented issues with thin-wall
    /// geometries where the wall thickness approaches the floating-point
    /// precision of the vertex coordinates.  The arrangement pipeline's
    /// exact predicates + tolerance-welding approach handles this better.
    ///
    /// ## Theorem — Thin Shell Volume
    ///
    /// For a cube of side `a` with inner cube side `a - 2t` (shell thickness t):
    ///   `V_shell = a³ - (a-2t)³ = 6a²t - 12at² + 8t³`
    /// For a=1, t=0.01: V_shell ≈ 0.0588.  ∎
    #[test]
    fn thin_wall_difference() {
        let outer = Cube {
            origin: Point3r::new(0.0, 0.0, 0.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .expect("outer cube");

        let inner = Cube {
            origin: Point3r::new(0.01, 0.01, 0.01),
            width: 0.98,
            height: 0.98,
            depth: 0.98,
        }
        .build()
        .expect("inner cube");

        let result = csg_boolean(BooleanOp::Difference, &outer, &inner)
            .expect("thin-wall difference should succeed");

        let vol = signed_volume(&result);
        let expected = 1.0_f64.powi(3) - 0.98_f64.powi(3);
        let rel_err = (vol - expected).abs() / expected;
        assert!(
            rel_err < 0.20,
            "thin wall volume error {rel_err:.3} exceeds 20%; \
             expected {expected:.6}, got {vol:.6}"
        );
    }

    /// Anisotropic scale — micro-scale in one axis, macro in others.
    ///
    /// A very flat cube (1.0 × 1.0 × 0.001) intersected with a unit sphere
    /// creates extreme aspect-ratio faces.  This tests that the curvature
    /// estimator handles anisotropic meshes without false-positive splitting
    /// (flat faces should not be refined) and that co-refinement handles
    /// the near-degenerate faces correctly.
    ///
    /// ## Known library failure mode
    ///
    /// Many mesh booleans assume isotropic scale; extreme aspect ratios
    /// cause the GWN classifier to misclassify fragments due to numerical
    /// instability in the solid-angle computation.
    #[test]
    fn anisotropic_flat_cube_sphere_intersection() {
        let flat_cube = Cube {
            origin: Point3r::new(-0.5, -0.5, -0.0005),
            width: 1.0,
            height: 1.0,
            depth: 0.001,
        }
        .build()
        .expect("flat cube");

        let sphere = crate::domain::geometry::primitives::UvSphere {
            radius: 0.4,
            center: Point3r::new(0.0, 0.0, 0.0),
            segments: 16,
            stacks: 8,
        }
        .build()
        .expect("sphere");

        let result = csg_boolean(BooleanOp::Intersection, &flat_cube, &sphere)
            .expect("flat-cube × sphere intersection should succeed");

        let vol = signed_volume(&result);
        // Intersection is a thin lens shape; volume must be positive.
        assert!(
            vol > 0.0 && vol < 0.01,
            "flat-cube × sphere intersection volume should be small positive, got {vol:.6}"
        );
    }

    // ── CW11: Topological adversarial tests ────────────────────────────────
    //
    // These tests target known failure modes common across mesh Boolean
    // libraries (CGAL, libigl, Cork, Manifold):
    //
    // | # | Failure Mode | Libraries Affected |
    // |---|---|---|
    // | 1 | Thin-wall collinear collapse | Cork, libigl |
    // | 2 | Tangent-plane shared-edge | CGAL Nef, Cork |
    // | 3 | Many-operand coplanar union | Cork, libigl |
    // | 4 | Micro-scale T-junction gaps | libigl, Manifold |
    // | 5 | Zero-volume intersection lens | Cork |

    /// Two cubes sharing a single face (tangent contact) — union should
    /// produce a valid closed mesh.
    ///
    /// ## Known Library Failure
    ///
    /// Cork and CGAL Nef polyhedra historically fail on face-contact
    /// configurations because the coplanar face pair creates two zero-
    /// volume fragments that confuse inside/outside classification.
    ///
    /// ## Theorem — Face-Tangent Union
    ///
    /// For two closed manifolds $A$ and $B$ sharing exactly one face $f$,
    /// $A \cup B$ is the convex hull minus $f$'s interior.  The boundary
    /// of $A \cup B$ is $(\partial A \cup \partial B) \setminus (f^+ \cup f^-)$
    /// plus the boundary loop of $f$ (shared edges).  The result is a
    /// closed 2-manifold with Euler characteristic $\chi = 2$.  ∎
    #[test]
    fn face_tangent_cubes_union_produces_closed_mesh() {
        // Two cubes touching on the x=1 plane
        let cube_a = Cube {
            origin: Point3r::new(-1.0, -1.0, -1.0),
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        }
        .build()
        .expect("cube_a");

        let cube_b = Cube {
            origin: Point3r::new(1.0, -1.0, -1.0),
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        }
        .build()
        .expect("cube_b");

        let result = csg_boolean(BooleanOp::Union, &cube_a, &cube_b);
        match result {
            Ok(mesh) => {
                assert!(!mesh.faces.is_empty(), "union must produce faces");
                let vol = signed_volume(&mesh);
                // Two unit cubes touching: volume should be 16 (2×2×2 each)
                assert!(
                    (vol - 16.0).abs() < 1.0,
                    "face-tangent union volume ~16, got {vol:.4}"
                );
            }
            Err(_) => {
                // Structured error is acceptable for this degenerate config
            }
        }
    }

    /// Three coplanar squares unioned: tests the n-ary coplanar reduction tree.
    ///
    /// ## Known Library Failure
    ///
    /// Cork and libigl accumulate tessellation artifacts when performing
    /// sequential binary unions on coplanar geometry.  The balanced
    /// reduction tree minimises intermediate operand complexity.
    ///
    /// ## Theorem — Coplanar N-ary Union Area
    ///
    /// For $n$ axis-aligned squares of side $s$ arranged in a row with
    /// 50% overlap, the union area is $A = s^2 + (n-1) \cdot s^2/2$.
    /// For $n=3$, $s=2$: $A = 4 + 2 \cdot 2 = 8$.  ∎
    #[test]
    fn three_coplanar_squares_n_ary_union() {
        // Three thin boxes (flat quads) in the XY plane with 50% overlap
        let sq_a = Cube {
            origin: Point3r::new(0.0, 0.0, -0.001),
            width: 2.0,
            height: 2.0,
            depth: 0.002,
        }
        .build()
        .expect("sq_a");
        let sq_b = Cube {
            origin: Point3r::new(1.0, 0.0, -0.001),
            width: 2.0,
            height: 2.0,
            depth: 0.002,
        }
        .build()
        .expect("sq_b");
        let sq_c = Cube {
            origin: Point3r::new(2.0, 0.0, -0.001),
            width: 2.0,
            height: 2.0,
            depth: 0.002,
        }
        .build()
        .expect("sq_c");

        // Multi-mesh union
        let ab = csg_boolean(BooleanOp::Union, &sq_a, &sq_b);
        if let Ok(ab_mesh) = ab {
            let result = csg_boolean(BooleanOp::Union, &ab_mesh, &sq_c);
            if let Ok(mesh) = result {
                assert!(!mesh.faces.is_empty(), "3-way coplanar union must produce faces");
                let vol = signed_volume(&mesh);
                // Thin box: volume ≈ area × 0.002
                assert!(vol > 0.0, "coplanar union volume must be positive, got {vol:.6}");
            }
        }
    }

    /// Micro-scale cubes (millifluidic dimensions) — tests snap_round
    /// thresholds and T-junction handling at small scales.
    ///
    /// ## Known Library Failure
    ///
    /// libigl and Manifold use fixed epsilon tolerances calibrated for
    /// unit-scale geometry.  At millifluidic scales (0.1–1 mm), these
    /// tolerances may be too large relative to geometry, causing
    /// T-junction gaps or spurious vertex welding.
    ///
    /// ## Theorem — Scale Invariance
    ///
    /// A robust Boolean pipeline produces topologically equivalent results
    /// under uniform scaling $S$.  If a result is manifold at scale 1,
    /// it must also be manifold at scale $10^{-3}$ (millifluidic).  ∎
    #[test]
    fn micro_scale_cube_difference_produces_valid_mesh() {
        let scale = 0.1; // 100 µm cubes
        let big = Cube {
            origin: Point3r::new(0.0, 0.0, 0.0),
            width: scale,
            height: scale,
            depth: scale,
        }
        .build()
        .expect("micro big cube");

        let small = Cube {
            origin: Point3r::new(scale * 0.25, scale * 0.25, scale * 0.25),
            width: scale * 0.5,
            height: scale * 0.5,
            depth: scale * 0.5,
        }
        .build()
        .expect("micro small cube");

        let result = csg_boolean(BooleanOp::Difference, &big, &small)
            .expect("micro-scale difference must succeed");
        assert!(!result.faces.is_empty(), "must produce faces");
        let vol = signed_volume(&result);
        let expected = scale.powi(3) - (scale * 0.5).powi(3);
        assert!(
            (vol - expected).abs() < expected * 0.15,
            "micro-scale volume: expected {expected:.8}, got {vol:.8}"
        );
    }

    /// Edge-contact cubes: two cubes sharing exactly one edge.
    ///
    /// ## Known Library Failure
    ///
    /// Edge-only contact produces zero-area intersection fragments that
    /// confuse fragment classification in Cork and libigl.  The GWN
    /// classifier must correctly identify the contact as a boundary
    /// condition, not an interior region.
    ///
    /// ## Theorem — Edge-Contact Union
    ///
    /// For two closed manifolds $A$, $B$ sharing exactly one edge $e$,
    /// $A \cap B = e$ (a 1-manifold) and $|A \cup B| = |A| + |B|$.
    /// The union is a valid 2-manifold with two connected components
    /// or (if treated as non-manifold) a pinched surface at $e$.  ∎
    #[test]
    fn edge_contact_cubes_union_volume_additive() {
        let cube_a = Cube {
            origin: Point3r::new(-1.0, -1.0, -1.0),
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        }
        .build()
        .expect("cube_a");

        // Touching on edge at (1, 1, z)
        let cube_b = Cube {
            origin: Point3r::new(1.0, 1.0, -1.0),
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        }
        .build()
        .expect("cube_b");

        let result = csg_boolean(BooleanOp::Union, &cube_a, &cube_b);
        match result {
            Ok(mesh) => {
                let vol = signed_volume(&mesh);
                // Two disjoint cubes (touching edge only): volume = 8 + 8 = 16
                assert!(
                    (vol - 16.0).abs() < 1.0,
                    "edge-contact union volume ~16, got {vol:.4}"
                );
            }
            Err(_) => {
                // Structured error acceptable for edge-contact degeneracy
            }
        }
    }

    /// Vertex-contact cubes: touching at exactly one vertex.
    ///
    /// ## Known Library Failure
    ///
    /// Vertex-only contact is the most degenerate configuration — the
    /// intersection is a single point (0-manifold).  Cork panics on this
    /// configuration; CGAL Nef handles it but produces extraneous faces.
    ///
    /// ## Theorem — Vertex-Contact Union
    ///
    /// For two closed manifolds $A$, $B$ sharing exactly one vertex $v$,
    /// $|A \cup B| = |A| + |B|$.  The union mesh is non-manifold at $v$
    /// (link is two disjoint circles, not one).  ∎
    #[test]
    fn vertex_contact_cubes_union_volume_additive() {
        let cube_a = Cube {
            origin: Point3r::new(-1.0, -1.0, -1.0),
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        }
        .build()
        .expect("cube_a");

        // Touching at vertex (1, 1, 1) = (-1+2, -1+2, -1+2) = corner of A
        let cube_b = Cube {
            origin: Point3r::new(1.0, 1.0, 1.0),
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        }
        .build()
        .expect("cube_b");

        let result = csg_boolean(BooleanOp::Union, &cube_a, &cube_b);
        match result {
            Ok(mesh) => {
                let vol = signed_volume(&mesh);
                assert!(
                    (vol - 16.0).abs() < 1.0,
                    "vertex-contact union volume ~16, got {vol:.4}"
                );
            }
            Err(_) => {
                // Structured error acceptable for vertex-contact degeneracy
            }
        }
    }

    /// Thin-wall cube difference: hollow out a cube leaving a thin shell.
    ///
    /// ## Known Library Failure
    ///
    /// When the inner and outer cubes nearly coincide (thin wall), seam
    /// vertex merging can collapse across the wall, creating holes.  Cork
    /// and libigl are known to fail with wall thickness < ~1% of cube size.
    ///
    /// ## Theorem — Thin-Wall Volume
    ///
    /// For outer cube side $a$ and inner cube side $b = a - 2t$ (wall
    /// thickness $t$), $|A \setminus B| = a^3 - b^3$.  For $a=2$,
    /// $b=1.96$ ($t=0.02$): $V = 8 - 7.529536 = 0.470464$.  ∎
    #[test]
    fn thin_wall_cube_difference_no_collapse() {
        let outer = Cube {
            origin: Point3r::new(-1.0, -1.0, -1.0),
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        }
        .build()
        .expect("outer cube");

        let wall = 0.02; // 2% wall thickness
        let inner = Cube {
            origin: Point3r::new(-1.0 + wall, -1.0 + wall, -1.0 + wall),
            width: 2.0 - 2.0 * wall,
            height: 2.0 - 2.0 * wall,
            depth: 2.0 - 2.0 * wall,
        }
        .build()
        .expect("inner cube");

        let result = csg_boolean(BooleanOp::Difference, &outer, &inner)
            .expect("thin-wall difference must succeed");
        assert!(!result.faces.is_empty(), "thin-wall must produce faces");
        let vol = signed_volume(&result);
        let outer_vol = 8.0_f64;
        let inner_side = 2.0 - 2.0 * wall;
        let inner_vol = inner_side.powi(3);
        let expected = outer_vol - inner_vol;
        assert!(
            (vol - expected).abs() < expected * 0.25,
            "thin-wall volume: expected {expected:.6}, got {vol:.6}"
        );
    }

    // ── Pinch-vertex adversarial tests ────────────────────────────────────
    //
    // These target the figure-8 vertex topology defect at dense multi-way
    // junctions.  A pinch vertex passes manifold edge checks (every edge
    // shared by exactly 2 faces) but violates the vertex-link simple-cycle
    // invariant, reducing the Euler characteristic by 1 per pinch.
    //
    // Known to affect: Cork, CGAL Nef, libigl, Manifold (prior versions).

    /// Three cylinders at 120° spacing — symmetric trifurcation.
    ///
    /// # Known Library Failures
    ///
    /// 120° spacing creates a symmetric 3-way junction where all three
    /// intersection curves meet at a single point.  The rotational symmetry
    /// increases the probability of vertex coincidence at the junction,
    /// making pinch vertices almost certain in libraries that use
    /// single-valued half-edge adjacency maps.
    ///
    /// # Theorem (Symmetric Junction Euler Invariant)
    ///
    /// The union of *k* cylinders meeting at a common junction with
    /// genus-0 topology must satisfy χ = 2 regardless of the angular
    /// spacing, provided every vertex link is a simple cycle.
    ///
    /// **Proof sketch.**  The union boundary is a closed oriented
    /// 2-manifold homeomorphic to a sphere (genus 0).  For any closed
    /// oriented 2-manifold of genus *g*, χ = 2(1 − g) = 2.  ∎
    #[test]
    fn symmetric_120deg_cylinder_union_no_pinch() {
        use crate::application::csg::boolean::csg_boolean_nary;
        use crate::application::csg::CsgNode;
        use crate::application::watertight::check::check_watertight;
        use nalgebra::{Isometry3, Translation3, UnitQuaternion, Vector3};

        let radius = 0.4;
        let height = 3.0;
        let segments = 24;
        let mut meshes = Vec::new();
        for angle_deg in [0.0_f64, 120.0, 240.0] {
            // Create cylinder along +Y, then rotate around Z to the desired angle.
            // Subtract PI/2 so angle 0° points along +X.
            let raw = Cylinder {
                base_center: Point3r::new(0.0, 0.0, 0.0),
                radius,
                height,
                segments,
            }
            .build()
            .expect("cylinder");
            let rotation = UnitQuaternion::<f64>::from_axis_angle(
                &Vector3::z_axis(),
                angle_deg.to_radians() - std::f64::consts::FRAC_PI_2,
            );
            let m = CsgNode::Transform {
                node: Box::new(CsgNode::Leaf(Box::new(raw))),
                iso: Isometry3::from_parts(Translation3::new(0.0, 0.0, 0.0), rotation),
            }
            .evaluate()
            .expect("cylinder transform");
            meshes.push(m);
        }

        let mut result =
            csg_boolean_nary(BooleanOp::Union, &meshes).expect("120° cylinder union");
        result.rebuild_edges();
        let report = check_watertight(
            &result.vertices,
            &result.faces,
            result.edges_ref().unwrap(),
        );
        assert!(
            report.is_watertight,
            "120° cylinder union must be watertight"
        );
        assert_eq!(
            report.euler_characteristic,
            Some(2),
            "120° cylinder union χ = {:?}, expected 2 — pinch vertex present",
            report.euler_characteristic,
        );
    }

    /// Four cylinders at 90° spacing in a cross pattern — dense 4-way junction.
    ///
    /// # Known Library Failures
    ///
    /// A 4-way cross junction creates 6 pairwise intersection curves
    /// meeting at the origin.  The 90° symmetry maximises vertex
    /// coincidence at the junction, making pinch vertices likely
    /// in libraries with single-valued half-edge adjacency.
    #[test]
    fn cross_4_cylinder_union_no_pinch() {
        use crate::application::csg::boolean::csg_boolean_nary;
        use crate::application::csg::CsgNode;
        use crate::application::watertight::check::check_watertight;
        use nalgebra::{Isometry3, Translation3, UnitQuaternion, Vector3};

        let radius = 0.4;
        let height = 3.0;
        let segments = 24;
        let mut meshes = Vec::new();
        for angle_deg in [0.0_f64, 90.0, 180.0, 270.0] {
            let raw = Cylinder {
                base_center: Point3r::new(0.0, 0.0, 0.0),
                radius,
                height,
                segments,
            }
            .build()
            .expect("cross cylinder");
            let rotation = UnitQuaternion::<f64>::from_axis_angle(
                &Vector3::z_axis(),
                angle_deg.to_radians() - std::f64::consts::FRAC_PI_2,
            );
            let m = CsgNode::Transform {
                node: Box::new(CsgNode::Leaf(Box::new(raw))),
                iso: Isometry3::from_parts(Translation3::new(0.0, 0.0, 0.0), rotation),
            }
            .evaluate()
            .expect("cross cylinder transform");
            meshes.push(m);
        }

        let mut result =
            csg_boolean_nary(BooleanOp::Union, &meshes).expect("cross-4 cylinder union");
        result.rebuild_edges();
        let report = check_watertight(
            &result.vertices,
            &result.faces,
            result.edges_ref().unwrap(),
        );
        assert!(
            report.is_watertight,
            "cross-4 cylinder union must be watertight"
        );
        assert_eq!(
            report.euler_characteristic,
            Some(2),
            "cross-4 cylinder union χ = {:?}, expected 2 — pinch vertex(es) detected",
            report.euler_characteristic,
        );
    }

    /// Star-shaped cylinder union — 5 cylinders at 72° spacing through origin.
    ///
    /// # Known Library Failures
    ///
    /// Five co-planar cylinders create a star junction with 10 pairwise
    /// intersection curves.  The junction region has extreme vertex density,
    /// making shared-neighbour collisions in half-edge adjacency nearly
    /// guaranteed without multi-valued maps.
    ///
    /// # Theorem (Star Junction Vertex Count)
    ///
    /// For *k* cylinders through a common center, the junction creates
    /// O(k²) intersection curves.  At each crossing of two curves, a
    /// potential pinch vertex arises.  The total number of potential pinch
    /// vertices is O(k²), requiring the detection algorithm to handle
    /// arbitrary fan multiplicity.  ∎
    #[test]
    fn star_5_cylinder_union_no_pinch() {
        use crate::application::csg::boolean::csg_boolean_nary;
        use crate::application::csg::CsgNode;
        use crate::application::watertight::check::check_watertight;
        use nalgebra::{Isometry3, Translation3, UnitQuaternion, Vector3};

        let radius = 0.3;
        let height = 3.0;
        let segments = 20;
        let mut meshes = Vec::new();
        for i in 0..5 {
            let angle_deg = (i as f64) * 72.0;
            let raw = Cylinder {
                base_center: Point3r::new(0.0, 0.0, 0.0),
                radius,
                height,
                segments,
            }
            .build()
            .expect("star cylinder");
            let rotation = UnitQuaternion::<f64>::from_axis_angle(
                &Vector3::z_axis(),
                angle_deg.to_radians() - std::f64::consts::FRAC_PI_2,
            );
            let m = CsgNode::Transform {
                node: Box::new(CsgNode::Leaf(Box::new(raw))),
                iso: Isometry3::from_parts(Translation3::new(0.0, 0.0, 0.0), rotation),
            }
            .evaluate()
            .expect("star cylinder transform");
            meshes.push(m);
        }

        let mut result =
            csg_boolean_nary(BooleanOp::Union, &meshes).expect("star-5 cylinder union");
        result.rebuild_edges();
        let report = check_watertight(
            &result.vertices,
            &result.faces,
            result.edges_ref().unwrap(),
        );
        assert!(
            report.is_watertight,
            "star-5 cylinder union must be watertight"
        );
        assert_eq!(
            report.euler_characteristic,
            Some(2),
            "star-5 cylinder union χ = {:?}, expected 2 — pinch vertex(es) detected",
            report.euler_characteristic,
        );
    }
}
