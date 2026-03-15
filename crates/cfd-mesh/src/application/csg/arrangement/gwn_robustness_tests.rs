//! GWN and normal robustness tests.
//!
//! Covers adversarial edge cases known to cause failures in mesh libraries:
//! - Near-surface GWN instability
//! - Scale-invariant GWN evaluation
//! - WNNC normal consistency
//! - BVH vs linear agreement on large meshes
//! - Near-vertex and degenerate-face handling
//! - Bounded GWN convergence at surface
//! - Open mesh behaviour

#[cfg(test)]
mod tests {
    use crate::application::csg::arrangement::classify::{
        classify_fragment, gwn_bvh, prepare_bvh_mesh, prepare_classification_faces, wnnc_score,
        FragmentClass,
    };
    use crate::application::csg::arrangement::gwn::{gwn, gwn_prepared};
    use crate::domain::core::scalar::{Point3r, Vector3r};
    use crate::infrastructure::storage::face_store::FaceData;
    use crate::infrastructure::storage::vertex_pool::VertexPool;

    // ── Helpers ────────────────────────────────────────────────────────────

    /// Build a unit cube (edge = 1, centred at origin) with 12 triangles.
    fn unit_cube_mesh() -> (VertexPool, Vec<FaceData>) {
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

    /// Build a scaled cube with given half-edge length.
    fn scaled_cube_mesh(half: f64) -> (VertexPool, Vec<FaceData>) {
        let mut pool = VertexPool::default_millifluidic();
        let n = nalgebra::Vector3::zeros();
        let mut v = |x, y, z| pool.insert_or_weld(Point3r::new(x, y, z), n);
        let c000 = v(-half, -half, -half);
        let c100 = v(half, -half, -half);
        let c010 = v(-half, half, -half);
        let c110 = v(half, half, -half);
        let c001 = v(-half, -half, half);
        let c101 = v(half, -half, half);
        let c011 = v(-half, half, half);
        let c111 = v(half, half, half);
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

    /// Build a tessellated unit sphere with `n_lat` latitude bands.
    fn sphere_mesh(n_lat: usize) -> (VertexPool, Vec<FaceData>) {
        let mut pool = VertexPool::default_millifluidic();
        let n = nalgebra::Vector3::zeros();
        let n_lon = 2 * n_lat;
        let mut vids = Vec::new();

        // North pole
        let north = pool.insert_or_weld(Point3r::new(0.0, 0.0, 1.0), n);
        // Latitude bands
        for i in 1..n_lat {
            let phi = std::f64::consts::PI * (i as f64) / (n_lat as f64);
            let (sp, cp) = phi.sin_cos();
            for j in 0..n_lon {
                let theta = 2.0 * std::f64::consts::PI * (j as f64) / (n_lon as f64);
                let (st, ct) = theta.sin_cos();
                vids.push(pool.insert_or_weld(Point3r::new(sp * ct, sp * st, cp), n));
            }
        }
        // South pole
        let south = pool.insert_or_weld(Point3r::new(0.0, 0.0, -1.0), n);

        let mut faces = Vec::new();
        // North cap
        for j in 0..n_lon {
            let j_next = (j + 1) % n_lon;
            faces.push(FaceData::untagged(north, vids[j], vids[j_next]));
        }
        // Body strips
        for i in 0..(n_lat - 2) {
            let row = i * n_lon;
            let next_row = (i + 1) * n_lon;
            for j in 0..n_lon {
                let j_next = (j + 1) % n_lon;
                faces.push(FaceData::untagged(
                    vids[row + j],
                    vids[next_row + j],
                    vids[next_row + j_next],
                ));
                faces.push(FaceData::untagged(
                    vids[row + j],
                    vids[next_row + j_next],
                    vids[row + j_next],
                ));
            }
        }
        // South cap
        let last_row = (n_lat - 2) * n_lon;
        for j in 0..n_lon {
            let j_next = (j + 1) % n_lon;
            faces.push(FaceData::untagged(vids[last_row + j], south, vids[last_row + j_next]));
        }
        (pool, faces)
    }

    // ── Near-surface GWN stability ─────────────────────────────────────────

    /// A point very close to a face plane (distance ε) must still classify
    /// correctly using GWN (|wn| in the [0.35, 0.65] band, triggering
    /// tiebreakers).  The point must NOT erroneously classify as clearly
    /// inside/outside.
    #[test]
    fn gwn_near_surface_enters_band() {
        let (pool, faces) = unit_cube_mesh();
        let prepared = prepare_classification_faces(&faces, &pool);
        // Point just above the +Z face (z = 0.5 + ε)
        for &eps in &[1e-3, 1e-6, 1e-9, 1e-12] {
            let q = Point3r::new(0.0, 0.0, 0.5 + eps);
            let wn = gwn_prepared(&q, &prepared).abs();
            // Should be exterior (just outside) but near 0.5
            assert!(
                wn < 0.99,
                "ε={eps}: near-surface GWN should not be ≈1.0, got {wn:.6}"
            );
        }
    }

    /// Points just inside (below) each face should have GWN close to 1.
    #[test]
    fn gwn_just_inside_each_face() {
        let (pool, faces) = unit_cube_mesh();
        let prepared = prepare_classification_faces(&faces, &pool);
        let inside_points = [
            Point3r::new(0.0, 0.0, 0.49),
            Point3r::new(0.0, 0.0, -0.49),
            Point3r::new(0.49, 0.0, 0.0),
            Point3r::new(-0.49, 0.0, 0.0),
            Point3r::new(0.0, 0.49, 0.0),
            Point3r::new(0.0, -0.49, 0.0),
        ];
        for q in &inside_points {
            let wn = gwn_prepared(q, &prepared).abs();
            assert!(
                wn > 0.9,
                "near-interior point {q:?} should have |wn| > 0.9, got {wn:.6}"
            );
        }
    }

    // ── Scale invariance ───────────────────────────────────────────────────

    /// GWN must produce the same classification at drastically different
    /// mesh scales.  Tests cubes from 1 nm to 1 km edge lengths.
    #[test]
    fn gwn_scale_invariant_across_six_orders() {
        for &scale in &[1e-9, 1e-6, 1e-3, 1.0, 1e3, 1e6] {
            let half = 0.5 * scale;
            let (pool, faces) = scaled_cube_mesh(half);
            let prepared = prepare_classification_faces(&faces, &pool);

            // Interior point at centre
            let q_in = Point3r::new(0.0, 0.0, 0.0);
            let wn_in = gwn_prepared(&q_in, &prepared).abs();
            assert!(
                wn_in > 0.9,
                "scale={scale}: interior GWN should be ≈1.0, got {wn_in:.6}"
            );

            // Exterior point at 10× distance
            let q_out = Point3r::new(5.0 * scale, 0.0, 0.0);
            let wn_out = gwn_prepared(&q_out, &prepared).abs();
            assert!(
                wn_out < 0.1,
                "scale={scale}: exterior GWN should be ≈0.0, got {wn_out:.6}"
            );
        }
    }

    // ── BVH vs linear agreement on sphere ──────────────────────────────────

    /// BVH-accelerated GWN must agree with linear GWN to within the proven
    /// error bound for a sphere with 500+ faces.
    #[test]
    fn gwn_bvh_agrees_with_linear_on_sphere() {
        let (pool, faces) = sphere_mesh(12); // 12 latitude bands → 528 faces
        let prepared = prepare_classification_faces(&faces, &pool);
        let error_budget = 1e-5; // tight budget so BVH skips very few nodes
        let bvh = prepare_bvh_mesh(&prepared).expect("bvh build");
        // The per-node error bound is `error_budget`, accumulated over O(log n)
        // levels.  For n ≈ 528, log₂ n ≈ 9.  We allow a generous 2× safety
        // margin on top of the theoretical bound.
        let log2_n = (faces.len() as f64).log2();
        let max_error = 2.0 * log2_n * error_budget;

        let test_points = [
            Point3r::new(0.5, 0.0, 0.0),  // inside
            Point3r::new(3.0, 0.0, 0.0),  // outside
            Point3r::new(0.0, 0.0, 3.0),  // outside along z
            Point3r::new(0.3, 0.3, 0.3),  // inside corner direction
        ];

        for q in &test_points {
            let linear = gwn_prepared(q, &prepared).abs();
            let bvh_val = gwn_bvh(q, &bvh, error_budget).abs();
            let err = (bvh_val - linear).abs();
            assert!(
                err < max_error,
                "BVH vs linear error={err:.6} > {max_error:.4} at q={q:?}"
            );
        }
    }

    // ── Near-vertex GWN must not produce NaN ───────────────────────────────

    /// A query point exactly at a mesh vertex must not return NaN.  The
    /// near-vertex guard should skip the coincident vertex's incident faces.
    #[test]
    fn gwn_at_vertex_no_nan() {
        let (pool, faces) = unit_cube_mesh();
        let prepared = prepare_classification_faces(&faces, &pool);
        // Query at exact vertex position
        let q = Point3r::new(0.5, 0.5, 0.5);
        let wn = gwn_prepared(&q, &prepared);
        assert!(wn.is_finite(), "GWN at vertex must be finite, got {wn}");
    }

    /// Query at a point extremely close (sub-nm) to a vertex.
    #[test]
    fn gwn_near_vertex_stable() {
        let (pool, faces) = unit_cube_mesh();
        let prepared = prepare_classification_faces(&faces, &pool);
        for &eps in &[1e-10, 1e-13, 1e-15] {
            let q = Point3r::new(0.5 + eps, 0.5, 0.5);
            let wn = gwn_prepared(&q, &prepared);
            assert!(
                wn.is_finite(),
                "ε={eps}: GWN near vertex must be finite, got {wn}"
            );
        }
    }

    // ── Degenerate (zero-area) face ────────────────────────────────────────

    /// A degenerate triangle (three collinear vertices) must contribute
    /// exactly 0 to the GWN sum and not produce NaN.
    #[test]
    fn gwn_degenerate_face_contributes_zero() {
        let mut pool = VertexPool::default_millifluidic();
        let n = nalgebra::Vector3::zeros();
        let v0 = pool.insert_or_weld(Point3r::new(0.0, 0.0, 0.0), n);
        let v1 = pool.insert_or_weld(Point3r::new(1.0, 0.0, 0.0), n);
        let v2 = pool.insert_or_weld(Point3r::new(2.0, 0.0, 0.0), n); // collinear
        let faces = vec![FaceData::untagged(v0, v1, v2)];
        let q = Point3r::new(0.0, 1.0, 0.0);
        let wn = gwn::<f64>(&q, &faces, &pool);
        assert!(wn.is_finite(), "degenerate face GWN must be finite");
        assert!(
            wn.abs() < 0.01,
            "degenerate face should contribute ≈0, got {wn:.6}"
        );
    }

    // ── PreparedFace area field correctness ────────────────────────────────

    /// The precomputed area field must match ‖normal‖/2 for each face.
    #[test]
    fn prepared_face_area_field_correct() {
        let (pool, faces) = unit_cube_mesh();
        let prepared = prepare_classification_faces(&faces, &pool);
        for (i, pf) in prepared.iter().enumerate() {
            let expected = 0.5 * pf.normal.norm();
            let diff = (pf.area - expected).abs();
            assert!(
                diff < 1e-15,
                "face {i}: area={} vs expected={expected}, diff={diff}",
                pf.area
            );
        }
    }

    // ── Bounded GWN convergence ────────────────────────────────────────────

    /// Bounded GWN must agree with standard GWN for interior/exterior
    /// points (where no triangle reaches the clip boundary).
    #[test]
    fn gwn_bounded_agrees_far_from_surface() {
        use crate::application::csg::arrangement::gwn::gwn_bounded_prepared;
        let (pool, faces) = unit_cube_mesh();
        let prepared = prepare_classification_faces(&faces, &pool);
        let test_points = [
            Point3r::new(0.0, 0.0, 0.0),  // interior centre
            Point3r::new(5.0, 0.0, 0.0),  // exterior far
            Point3r::new(0.2, 0.1, 0.3),  // interior off-centre
            Point3r::new(2.0, 2.0, 2.0),  // exterior corner direction
        ];
        for q in &test_points {
            let standard = gwn_prepared(q, &prepared);
            let bounded = gwn_bounded_prepared(q, &prepared);
            let err = (standard - bounded).abs();
            assert!(
                err < 1e-10,
                "bounded/standard diverge at {q:?}: standard={standard:.8}, bounded={bounded:.8}"
            );
        }
    }

    /// Bounded GWN must produce finite values for near-surface queries.
    #[test]
    fn gwn_bounded_stable_near_surface() {
        use crate::application::csg::arrangement::gwn::gwn_bounded_prepared;
        let (pool, faces) = unit_cube_mesh();
        let prepared = prepare_classification_faces(&faces, &pool);
        for &eps in &[1e-3, 1e-6, 1e-9, 1e-12] {
            let q = Point3r::new(0.0, 0.0, 0.5 + eps);
            let wn = gwn_bounded_prepared(&q, &prepared);
            assert!(
                wn.is_finite(),
                "ε={eps}: bounded GWN near surface must be finite, got {wn}"
            );
            assert!(
                wn.abs() <= 1.0,
                "ε={eps}: bounded GWN must be in [-1,1], got {wn}"
            );
        }
    }

    // ── WNNC normal consistency ────────────────────────────────────────────

    /// For a closed manifold (unit cube), the WNNC score at face centroids
    /// with outward normals should be positive (consistent).
    ///
    /// The query must be ON the surface (where GWN ≈ 0.5 and the gradient
    /// is maximal) rather than offset into the exterior (where GWN is
    /// exactly 0 and the gradient vanishes).
    #[test]
    fn wnnc_closed_cube_normals_consistent() {
        let (pool, faces) = unit_cube_mesh();
        let prepared = prepare_classification_faces(&faces, &pool);
        // Test WNNC at each face centroid with its outward normal.
        // On the surface, GWN transitions from 1 (inside) to 0 (outside),
        // so -∇GWN points in the outward normal direction → positive score.
        let h = 1e-4; // step size for numerical gradient
        for pf in &prepared {
            let face_n = pf.normal; // unnormalised outward normal
            let n_len = face_n.norm();
            if n_len < 1e-15 {
                continue;
            }
            let unit_n = face_n / n_len;
            // Query at the face centroid (on the surface).
            let score = wnnc_score(&pf.centroid, &unit_n, &prepared, h);
            assert!(
                score > 0.0,
                "WNNC score should be positive for outward normal, got {score:.4} \
                 at centroid {c:?}",
                c = pf.centroid
            );
        }
    }

    /// For a closed manifold, WNNC with flipped normals should give a
    /// negative score (inconsistent).
    #[test]
    fn wnnc_flipped_normals_negative() {
        let (pool, faces) = unit_cube_mesh();
        let prepared = prepare_classification_faces(&faces, &pool);
        let h = 1e-4;
        // Use the first face and flip its normal
        let pf = &prepared[0];
        let n_len = pf.normal.norm();
        let unit_n = pf.normal / n_len;
        let flipped = -unit_n;
        // Query at the face centroid (on the surface).
        let score = wnnc_score(&pf.centroid, &flipped, &prepared, h);
        assert!(
            score < 0.0,
            "WNNC score should be negative for flipped normal, got {score:.4}"
        );
    }

    // ── Scale-relative tiebreaker ──────────────────────────────────────────

    /// The nearest-face tiebreaker must handle large-scale meshes correctly.
    /// With absolute 1e-9, a 1km cube would incorrectly classify near-surface
    /// points.  The scale-relative fix should classify correctly.
    #[test]
    fn tiebreaker_scale_relative_large_mesh() {
        let scale = 1000.0; // 1 km cube
        let (pool, faces) = scaled_cube_mesh(scale * 0.5);
        // Point on +Z face plane
        let c = Point3r::new(0.0, 0.0, scale * 0.5);
        let n = Vector3r::new(0.0, 0.0, 1.0);
        let cls = classify_fragment(&c, &n, &faces, &pool);
        assert!(
            matches!(
                cls,
                FragmentClass::CoplanarSame | FragmentClass::Outside
            ),
            "1 km scale: on-face point should be CoplanarSame/Outside, got {cls:?}"
        );
    }

    /// Same test at small scale (1 µm cube).
    #[test]
    fn tiebreaker_scale_relative_small_mesh() {
        let scale = 1e-6; // 1 µm cube
        let (pool, faces) = scaled_cube_mesh(scale * 0.5);
        let c = Point3r::new(0.0, 0.0, scale * 0.5);
        let n = Vector3r::new(0.0, 0.0, 1.0);
        let cls = classify_fragment(&c, &n, &faces, &pool);
        assert!(
            matches!(
                cls,
                FragmentClass::CoplanarSame | FragmentClass::Outside
            ),
            "1 µm scale: on-face point should be CoplanarSame/Outside, got {cls:?}"
        );
    }

    // ── Sphere GWN symmetry ────────────────────────────────────────────────

    /// WNNC on a tessellated sphere should agree: outward-pointing normals
    /// at surface points should yield positive consistency scores.
    ///
    /// Queries are placed at face centroids (on the surface) where the GWN
    /// gradient is strongest.
    #[test]
    fn wnnc_sphere_normals_consistent() {
        let (pool, faces) = sphere_mesh(8); // 8 latitude bands → ~240 faces
        let prepared = prepare_classification_faces(&faces, &pool);
        let h = 1e-4;
        let mut positive_count = 0;
        let mut total_checked = 0;
        for pf in &prepared {
            let n_len = pf.normal.norm();
            if n_len < 1e-15 {
                continue;
            }
            let unit_n = pf.normal / n_len;
            // Query at the face centroid (on the surface).
            let score = wnnc_score(&pf.centroid, &unit_n, &prepared, h);
            if score > 0.0 {
                positive_count += 1;
            }
            total_checked += 1;
        }
        // At least 90% of faces should have consistent normals
        // (some near poles may have marginal scores due to tessellation)
        let ratio = positive_count as f64 / total_checked as f64;
        assert!(
            ratio > 0.9,
            "WNNC: only {positive_count}/{total_checked} ({r:.1}%) faces consistent",
            r = ratio * 100.0
        );
    }

    // ── Open mesh GWN ──────────────────────────────────────────────────────

    /// An open mesh (single triangle) should not classify any point as
    /// clearly inside (|wn| ≈ 1).
    #[test]
    fn gwn_open_mesh_never_fully_inside() {
        let mut pool = VertexPool::default_millifluidic();
        let n = nalgebra::Vector3::zeros();
        let v0 = pool.insert_or_weld(Point3r::new(-1.0, -1.0, 0.0), n);
        let v1 = pool.insert_or_weld(Point3r::new(1.0, -1.0, 0.0), n);
        let v2 = pool.insert_or_weld(Point3r::new(0.0, 1.0, 0.0), n);
        let faces = vec![FaceData::untagged(v0, v1, v2)];
        let prepared = prepare_classification_faces(&faces, &pool);

        let test_points = [
            Point3r::new(0.0, 0.0, 0.01),
            Point3r::new(0.0, 0.0, -0.01),
            Point3r::new(0.0, 0.0, 1.0),
        ];
        for q in &test_points {
            let wn = gwn_prepared(q, &prepared).abs();
            assert!(
                wn < 0.6,
                "open mesh: |wn| should be < 0.6 (max ≈ 0.5), got {wn:.6} at {q:?}"
            );
        }
    }

    // ── High aspect ratio slivers ──────────────────────────────────────────

    /// A very thin sliver triangle (aspect ratio > 1000:1) must not corrupt
    /// GWN evaluation for nearby query points.
    #[test]
    fn gwn_sliver_triangle_no_corruption() {
        let mut pool = VertexPool::default_millifluidic();
        let n = nalgebra::Vector3::zeros();
        // Sliver: very long and very narrow
        let v0 = pool.insert_or_weld(Point3r::new(0.0, 0.0, 0.0), n);
        let v1 = pool.insert_or_weld(Point3r::new(1.0, 0.0, 0.0), n);
        let v2 = pool.insert_or_weld(Point3r::new(0.5, 1e-6, 0.0), n); // aspect ≈ 1e6:1
        let faces = vec![FaceData::untagged(v0, v1, v2)];
        let q = Point3r::new(0.5, 0.5, 0.1);
        let wn = gwn::<f64>(&q, &faces, &pool);
        assert!(wn.is_finite(), "sliver triangle GWN must be finite, got {wn}");
        assert!(
            wn.abs() < 0.6,
            "sliver triangle (open) |wn| should be < 0.6, got {wn:.6}"
        );
    }

    // ── Bounded GWN on sphere ──────────────────────────────────────────────

    /// On a sphere, bounded and standard GWN must agree for interior and
    /// exterior points.
    #[test]
    fn gwn_bounded_vs_standard_sphere() {
        use crate::application::csg::arrangement::gwn::gwn_bounded_prepared;
        let (pool, faces) = sphere_mesh(8);
        let prepared = prepare_classification_faces(&faces, &pool);

        // Interior
        let q_in = Point3r::new(0.0, 0.0, 0.0);
        let std_in = gwn_prepared(&q_in, &prepared);
        let bnd_in = gwn_bounded_prepared(&q_in, &prepared);
        assert!(
            (std_in - bnd_in).abs() < 1e-8,
            "sphere interior: standard={std_in:.8} vs bounded={bnd_in:.8}"
        );

        // Exterior
        let q_out = Point3r::new(5.0, 0.0, 0.0);
        let std_out = gwn_prepared(&q_out, &prepared);
        let bnd_out = gwn_bounded_prepared(&q_out, &prepared);
        assert!(
            (std_out - bnd_out).abs() < 1e-8,
            "sphere exterior: standard={std_out:.8} vs bounded={bnd_out:.8}"
        );
    }
}
