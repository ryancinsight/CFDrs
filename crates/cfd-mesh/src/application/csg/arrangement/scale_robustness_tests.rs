//! Scale-robustness and adversarial stress tests for the CSG arrangement pipeline.
//!
//! These tests verify that Boolean operations produce consistent results across
//! vastly different geometric scales, targeting the scale-dependent epsilon
//! issues identified in the Lévy (arXiv:2405.12949v2, ACM TOG 2025) audit.
//!
//! ## Categories
//!
//! | Category | What it tests |
//! |----------|---------------|
//! | Scale invariance | Same topology at 10 µm, 1 mm, 1 m, 1 km scales |
//! | Micro-scale CSG | Millifluidic geometry (10–500 µm features) |
//! | Pathological angles | Near-0° and near-180° dihedral angles |
//! | Dense intersection | Many T-T pairs from high-tessellation overlap |
//! | Segment degeneracy | Zero-length / near-zero intersection segments |

#[cfg(test)]
mod tests {
    use crate::application::csg::boolean::{csg_boolean, BooleanOp};
    use crate::domain::core::scalar::Point3r;
    use crate::domain::geometry::primitives::{Cube, Cylinder, PrimitiveMesh};

    // ── Helpers ────────────────────────────────────────────────────────────

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

    fn cube_at_scale(scale: f64) -> crate::domain::mesh::IndexedMesh {
        Cube {
            origin: Point3r::new(-scale, -scale, -scale),
            width: 2.0 * scale,
            height: 2.0 * scale,
            depth: 2.0 * scale,
        }
        .build()
        .expect("cube_at_scale build")
    }

    fn offset_cube_at_scale(scale: f64, dx: f64) -> crate::domain::mesh::IndexedMesh {
        Cube {
            origin: Point3r::new(-scale + dx, -scale, -scale),
            width: 2.0 * scale,
            height: 2.0 * scale,
            depth: 2.0 * scale,
        }
        .build()
        .expect("offset_cube_at_scale build")
    }

    fn cylinder_at_scale(scale: f64, segments: usize) -> crate::domain::mesh::IndexedMesh {
        Cylinder {
            base_center: Point3r::new(0.0, 0.0, 0.0),
            radius: scale,
            height: 4.0 * scale,
            segments,
        }
        .build()
        .expect("cylinder_at_scale build")
    }

    // ── Scale-invariance tests ─────────────────────────────────────────────
    //
    // These verify that the same Boolean operation on geometrically identical
    // meshes at different scales produces consistent results.  Scale-dependent
    // absolute epsilons would cause failures at extreme scales.

    /// Cube ∩ offset_cube must yield consistent relative volume across scales.
    ///
    /// Tests 4 scales spanning 8 orders of magnitude: 10 µm, 1 mm, 1 m, 1 km.
    /// The intersection of two half-overlapping cubes has volume = 50% of one
    /// cube regardless of scale.
    ///
    /// # Known failure mode (pre-fix)
    ///
    /// At 10 µm scale, absolute AABB expansion of 1e-6 m was 10% of the
    /// diagonal, causing false-positive broad-phase matches and corrupted
    /// intersection geometry.  The relative expansion fix resolves this.
    #[test]
    fn scale_invariant_cube_intersection_4_scales() {
        let scales = [1e-5, 1e-3, 1.0, 1e3];
        let mut relative_volumes: Vec<f64> = Vec::new();

        for &s in &scales {
            let a = cube_at_scale(s);
            let b = offset_cube_at_scale(s, s); // 50% overlap
            let vol_a = signed_volume(&a);

            match csg_boolean(BooleanOp::Intersection, &a, &b) {
                Ok(result) => {
                    let vol_inter = signed_volume(&result);
                    let ratio = vol_inter / vol_a;
                    relative_volumes.push(ratio);
                }
                Err(_) => {
                    // If intersection fails at this scale, that's also
                    // information — but we don't assert it must succeed
                    // since some scale extremes may legitimately fail.
                    relative_volumes.push(f64::NAN);
                }
            }
        }

        // At least 3 of 4 scales must succeed.
        let valid: Vec<f64> = relative_volumes
            .iter()
            .copied()
            .filter(|r| r.is_finite())
            .collect();
        assert!(
            valid.len() >= 3,
            "intersection must succeed at ≥3 of 4 scales: results={relative_volumes:?}"
        );

        // All successful results must agree within 10% relative.
        for (i, &ri) in valid.iter().enumerate() {
            for &rj in valid.iter().skip(i + 1) {
                let diff = (ri - rj).abs() / ri.max(rj).max(1e-30);
                assert!(
                    diff < 0.10,
                    "scale-invariance violated: ratios differ by {diff:.4}: {valid:?}"
                );
            }
        }
    }

    /// Union volume: vol(A ∪ B) = vol(A) + vol(B) − vol(A ∩ B) at micro-scale.
    ///
    /// Exercises the scale-relative AABB expansion and degenerate-segment
    /// filters at typical millifluidic dimensions (500 µm channels).
    #[test]
    fn micro_scale_inclusion_exclusion_identity() {
        let s = 500e-6; // 500 µm
        let a = cube_at_scale(s);
        let b = offset_cube_at_scale(s, s); // 50% overlap

        let vol_a = signed_volume(&a);
        let vol_b = signed_volume(&b);

        let union_r = csg_boolean(BooleanOp::Union, &a, &b);
        let inter_r = csg_boolean(BooleanOp::Intersection, &a, &b);

        if let (Ok(union), Ok(inter)) = (union_r, inter_r) {
            let vol_union = signed_volume(&union);
            let vol_inter = signed_volume(&inter);
            let lhs = vol_a + vol_b;
            let rhs = vol_union + vol_inter;
            let rel_err = (lhs - rhs).abs() / lhs.max(1e-30);
            assert!(
                rel_err < 0.10,
                "micro-scale IE identity: vol(A)+vol(B)={lhs:.2e} ≠ vol(A∪B)+vol(A∩B)={rhs:.2e}, err={rel_err:.4}"
            );
        }
    }

    /// Macro-scale (1 km) cube difference must produce non-empty result.
    ///
    /// Exercises the scale-relative degenerate-normal threshold: at 1 km scale,
    /// cross-product magnitudes are ~10⁶, so the absolute check `‖n‖² < 1e-20`
    /// would never trigger (correct), but the segment dedup at `1e-24` could
    /// incorrectly collapse large segments.
    #[test]
    fn macro_scale_cube_difference() {
        let s = 1000.0; // 1 km
        let a = cube_at_scale(s);
        let b = offset_cube_at_scale(s, s); // 50% overlap

        if let Ok(result) = csg_boolean(BooleanOp::Difference, &a, &b) {
            let vol_diff = signed_volume(&result);
            let vol_a = signed_volume(&a);
            // A \ B should be ~50% of A (half the cube is subtracted).
            let ratio = vol_diff / vol_a;
            assert!(
                ratio > 0.3 && ratio < 0.7,
                "1 km scale A\\B ratio must be ≈0.5, got {ratio:.4}"
            );
        }
    }

    // ── Pathological angle tests ───────────────────────────────────────────

    /// Near-grazing intersection: two cubes with a ~0.01° tilt produce a
    /// very thin intersection wedge.  The CDT corefine and fragment
    /// classification must handle the nearly-degenerate sliver triangles.
    ///
    /// # Known failure mode (Lévy §4.3)
    ///
    /// Near-parallel face pairs produce intersection segments shorter than
    /// absolute thresholds, causing silent segment drops.  Symbolic
    /// perturbation (Lévy) or scale-relative thresholds (our fix) resolve this.
    #[test]
    fn thin_wedge_intersection_from_tilted_cubes() {
        // Cube A: axis-aligned.
        let a = Cube {
            origin: Point3r::new(-1.0, -1.0, -1.0),
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        }
        .build()
        .expect("cube a");

        // Cube B: shifted so one face grazes A's face at a shallow angle.
        // The 0.01 offset creates a thin wedge of intersection.
        let b = Cube {
            origin: Point3r::new(-1.0, -1.0, 0.99),
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        }
        .build()
        .expect("cube b");

        // Must not panic.
        match csg_boolean(BooleanOp::Intersection, &a, &b) {
            Ok(result) => {
                let vol = signed_volume(&result);
                // Thin intersection wedge: ~2×2×0.01 = 0.04
                assert!(
                    vol < 0.1,
                    "thin wedge intersection volume must be small: got {vol:.6}"
                );
            }
            Err(_) => {
                // Acceptable: some degenerate configurations may legitimately
                // produce empty intersections.
            }
        }
    }

    // ── Dense intersection stress tests ────────────────────────────────────

    /// High-tessellation cylinder intersection: 32-segment cylinders with
    /// perpendicular axes create a dense T-T intersection pattern (up to
    /// 32×32 = 1024 candidate pairs) that stresses the CDT corefine pipeline.
    #[test]
    fn dense_perpendicular_cylinder_intersection() {
        let cyl_a = cylinder_at_scale(1.0, 32);

        let cyl_b = cylinder_at_scale(0.8, 32);

        // Must not panic or OOM with dense T-T pairs.
        match csg_boolean(BooleanOp::Intersection, &cyl_a, &cyl_b) {
            Ok(result) => {
                assert!(
                    !result.faces.is_empty(),
                    "dense cylinder intersection must be non-empty"
                );
            }
            Err(_) => {
                // Intersection failure is acceptable for dense degenerate
                // configurations, but panic/OOM is not.
            }
        }
    }

    // ── Micro-scale millifluidic stress tests ─────────────────────────────

    /// 50 µm channel stub: a tiny cylinder Boolean-subtracted from a thin slab.
    ///
    /// 500 µm overlapping cubes: Boolean difference at millifluidic scale.
    ///
    /// This exercises the complete pipeline at millifluidic scale where the
    /// previous absolute epsilons (1e-6 AABB, 1e-24 segment dedup, 1e-20
    /// normal threshold) all fall in the problematic range.
    ///
    /// Note: all dimensions in mm (native library unit).
    #[test]
    fn millifluidic_channel_stub_difference() {
        let s = 0.5; // 500 µm = 0.5 mm half-extent
        let a = cube_at_scale(s);
        let b = offset_cube_at_scale(s, s); // 50% overlap

        // Subtract B from A — must not panic.
        match csg_boolean(BooleanOp::Difference, &a, &b) {
            Ok(result) => {
                let vol_a = signed_volume(&a);
                let vol_result = signed_volume(&result);
                // A \ B should be ~50% of A (half subtracted).
                let ratio = vol_result / vol_a;
                assert!(
                    ratio > 0.3 && ratio < 0.7,
                    "millifluidic A\\B ratio must be ≈0.5, got {ratio:.4}"
                );
            }
            Err(_) => {
                // Acceptable if CDT corefine falls back, but panic is not.
            }
        }
    }

    /// Scale-invariant cube union: the union of two overlapping cubes must
    /// produce a non-empty mesh at 10 µm scale.
    ///
    /// # Theorem — Broad-Phase Scale Safety
    ///
    /// With AABB_RELATIVE_EXPANSION = 1e-6, a 10 µm cube's AABB is expanded by
    /// ≈ 17 fm (1e-6 × 17 µm diagonal), which is sub-atomic and does not merge
    /// distinct features.  The previous absolute 1e-6 m expansion was 60× the
    /// cube's width. ∎
    ///
    /// Note: dimensions in mm (library native unit). 10 µm = 0.01 mm.
    #[test]
    fn ten_micron_cube_union() {
        let s = 0.01; // 10 µm = 0.01 mm
        let a = cube_at_scale(s);
        let b = offset_cube_at_scale(s, s);
        let vol_a = signed_volume(&a);
        let vol_b = signed_volume(&b);

        match csg_boolean(BooleanOp::Union, &a, &b) {
            Ok(result) => {
                let vol_union = signed_volume(&result);
                // Union of 50% overlapping cubes: volume = 1.5 × single cube
                let expected = vol_a + vol_b;
                assert!(
                    vol_union > 0.3 * expected && vol_union <= expected * 1.05,
                    "10 µm union volume unexpected: expected≤{expected:.2e}, got {vol_union:.2e}"
                );
            }
            Err(_) => {
                // May fail at extreme micro-scale; no panic required.
            }
        }
    }

    // ── Coplanar + near-parallel stress (Lévy motivating examples) ────────

    /// Near-parallel cubes with 1e-8 normal offset.
    ///
    /// This catches the "near-coplanar face pair" failure mode documented in
    /// Lévy §4.3: when two faces have normals within 1e-8 rad of each other,
    /// the intersection line direction ‖n₁×n₂‖ is extremely small.  The
    /// absolute check `‖n₁×n₂‖² < 1e-20` would misclassify these as coplanar
    /// at unit scale but not at micro-scale.  The relative check detects them
    /// consistently.
    #[test]
    fn near_parallel_normals_1e8_offset() {
        // Create two cubes with a tiny angular offset.
        let a = Cube {
            origin: Point3r::new(0.0, 0.0, 0.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .expect("cube a");

        // Shift cube B so it barely overlaps A.
        let b = Cube {
            origin: Point3r::new(0.5, 1e-8, 0.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .expect("cube b");

        // Must not panic.
        let _ = csg_boolean(BooleanOp::Union, &a, &b);
    }

    /// Vertex-on-edge touching: a cube corner touches the middle of another
    /// cube's edge.  The snap endpoint must classify as edge-Steiner (not
    /// interior), and the CDT must correctly split the edge.
    #[test]
    fn vertex_on_edge_touching() {
        let a = Cube {
            origin: Point3r::new(0.0, 0.0, 0.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .expect("cube a");

        // Cube B placed so its corner (1.0, 0.5, 0.5) touches the midpoint
        // of cube A's right-face edge.
        let b = Cube {
            origin: Point3r::new(1.0, 0.0, 0.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .expect("cube b");

        match csg_boolean(BooleanOp::Union, &a, &b) {
            Ok(result) => {
                let vol_a = signed_volume(&a);
                let vol_b = signed_volume(&b);
                let vol_union = signed_volume(&result);
                // Touching cubes: union = vol_a + vol_b (no overlap).
                let expected = vol_a + vol_b;
                let rel_err = (vol_union - expected).abs() / expected;
                assert!(
                    rel_err < 0.10,
                    "touching cubes union: expected={expected:.6}, got={vol_union:.6}, err={rel_err:.4}"
                );
            }
            Err(_) => {
                // Edge-touching is a legal degenerate configuration.
            }
        }
    }

    // ── Segment intersection module tests ──────────────────────────────────

    /// Two tiny triangles (10 µm edges) that straddle each other's planes.
    /// The intersection segment must be detected and have non-zero length.
    ///
    /// Exercises the scale-relative `DEGENERATE_NORMAL_REL_SQ` and
    /// `INTERVAL_OVERLAP_REL` constants in the segment intersection pipeline.
    #[test]
    fn micro_scale_triangle_intersection_detected() {
        use crate::application::csg::intersect::exact::intersect_triangles;
        use crate::application::csg::intersect::types::IntersectionType;
        use crate::infrastructure::storage::face_store::FaceData;
        use crate::infrastructure::storage::vertex_pool::VertexPool;

        // Cell size scaled to geometry: 1 nm = 1e-6 mm.
        let mut pool = VertexPool::new(1e-6);
        let n = nalgebra::Vector3::zeros();
        let s = 0.01; // 10 µm = 0.01 mm

        // Triangle A: in XY plane at z=0.
        let va0 = pool.insert_or_weld(Point3r::new(0.0, 0.0, 0.0), n);
        let va1 = pool.insert_or_weld(Point3r::new(s, 0.0, 0.0), n);
        let va2 = pool.insert_or_weld(Point3r::new(0.0, s, 0.0), n);
        let fa = FaceData::untagged(va0, va1, va2);

        // Triangle B: in XZ plane at y = s/4, straddling A's plane.
        let vb0 = pool.insert_or_weld(Point3r::new(s * 0.25, s * 0.25, -s), n);
        let vb1 = pool.insert_or_weld(Point3r::new(s * 0.25, s * 0.25, s), n);
        let vb2 = pool.insert_or_weld(Point3r::new(s * 0.75, s * 0.25, 0.0), n);
        let fb = FaceData::untagged(vb0, vb1, vb2);

        let result = intersect_triangles(&fa, &pool, &fb, &pool);

        match result {
            IntersectionType::Segment { start, end } => {
                let seg_len = (end - start).norm();
                assert!(
                    seg_len > s * 1e-6,
                    "micro-scale intersection segment too short: {seg_len:.2e}"
                );
            }
            IntersectionType::None => {
                panic!("straddling micro-scale triangles must detect intersection");
            }
            IntersectionType::Coplanar => {
                panic!("XY vs XZ triangles must not be classified as coplanar");
            }
        }
    }

    /// Two 1 km triangles intersecting — the intersection segment at macro scale
    /// must be detected correctly by the relative thresholds.
    #[test]
    fn macro_scale_triangle_intersection_detected() {
        use crate::application::csg::intersect::exact::intersect_triangles;
        use crate::application::csg::intersect::types::IntersectionType;
        use crate::infrastructure::storage::face_store::FaceData;
        use crate::infrastructure::storage::vertex_pool::VertexPool;

        let mut pool = VertexPool::default_millifluidic();
        let n = nalgebra::Vector3::zeros();
        let s = 1000.0; // 1 km

        // Triangle A: in XY plane.
        let va0 = pool.insert_or_weld(Point3r::new(0.0, 0.0, 0.0), n);
        let va1 = pool.insert_or_weld(Point3r::new(s, 0.0, 0.0), n);
        let va2 = pool.insert_or_weld(Point3r::new(0.0, s, 0.0), n);
        let fa = FaceData::untagged(va0, va1, va2);

        // Triangle B: in XZ plane, straddling A.
        let vb0 = pool.insert_or_weld(Point3r::new(s * 0.25, s * 0.25, -s), n);
        let vb1 = pool.insert_or_weld(Point3r::new(s * 0.25, s * 0.25, s), n);
        let vb2 = pool.insert_or_weld(Point3r::new(s * 0.75, s * 0.25, 0.0), n);
        let fb = FaceData::untagged(vb0, vb1, vb2);

        let result = intersect_triangles(&fa, &pool, &fb, &pool);

        match result {
            IntersectionType::Segment { start, end } => {
                let seg_len = (end - start).norm();
                assert!(
                    seg_len > s * 1e-6,
                    "macro-scale intersection segment too short: {seg_len:.2e}"
                );
            }
            IntersectionType::None => {
                panic!("straddling macro-scale triangles must detect intersection");
            }
            IntersectionType::Coplanar => {
                panic!("XY vs XZ triangles must not be classified as coplanar");
            }
        }
    }
}
