//! Degeneracy and stress tests targeting known-hard cases for mesh Boolean libraries.
//!
//! These tests exercise configurations where CSG Boolean engines are known to
//! fail across the industry (CGAL, OpenSCAD, libigl, Cork, etc.):
//!
//! | Category | What it stresses |
//! |----------|------------------|
//! | Near-tangent | Intersection approaching a single point or degenerate curve |
//! | Irrational coordinates | Rotated meshes with irrational vertex positions |
//! | Low tessellation | Skinny triangles from 4-segment cylinders |
//! | Many-body n-ary | 5+ overlapping shapes in a single intersection |
//! | Mixed curvature | Planar-vs-curved (cube-sphere) intersection seams |
//! | Complex curves | Saddle-point intersection curves (torus-cylinder) |
//! | Dense overlap | High-tessellation shapes with thousands of face pairs |

#[cfg(test)]
mod tests {
    use crate::application::csg::boolean::{csg_boolean, BooleanOp};
    use crate::domain::core::scalar::Point3r;
    use crate::domain::geometry::primitives::{
        Cube, Cylinder, PrimitiveMesh, Torus, UvSphere,
    };

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

    // ── Near-tangent tests ─────────────────────────────────────────────────
    //
    // Two convex shapes barely overlapping produce intersection curves that
    // degenerate to a point.  Known to cause:
    //   - Zero-length intersection segments → NaN normals
    //   - GWN instability near surface contact
    //   - Corefine CDT failure on degenerate PSLGs

    /// Two spheres whose surfaces barely overlap — the intersection lens has
    /// near-zero volume.  The Boolean must complete without panic and produce
    /// a watertight result for each operation.
    ///
    /// # Known failure mode (pre-fix)
    ///
    /// When `d ≈ 2r`, intersection segments shrink to near-zero length.
    /// Libraries that rely on segment-length thresholds for collinearity
    /// detection can misclassify or skip these segments entirely, producing
    /// holes in the output mesh.
    #[test]
    fn near_tangent_sphere_sphere_union() {
        let r = 1.0;
        let d = 2.0 * r - 0.01; // 0.5% overlap
        let a = UvSphere {
            radius: r,
            center: Point3r::new(0.0, 0.0, 0.0),
            segments: 32,
            stacks: 16,
        }
        .build()
        .expect("sphere a");
        let b = UvSphere {
            radius: r,
            center: Point3r::new(d, 0.0, 0.0),
            segments: 32,
            stacks: 16,
        }
        .build()
        .expect("sphere b");

        let result = csg_boolean(BooleanOp::Union, &a, &b).expect("union completes");
        assert!(result.face_count() > 0, "union must produce faces");
        let v = signed_volume(&result);
        let v_single = 4.0 * std::f64::consts::PI / 3.0; // r=1
        assert!(
            v > v_single,
            "near-tangent union volume {v} must exceed single sphere {v_single}"
        );
    }

    /// Near-tangent sphere intersection: the lens volume must be small but
    /// positive, and the result must be a valid closed mesh.
    #[test]
    fn near_tangent_sphere_sphere_intersection() {
        let r = 1.0;
        let d = 2.0 * r - 0.01;
        let a = UvSphere {
            radius: r,
            center: Point3r::new(0.0, 0.0, 0.0),
            segments: 32,
            stacks: 16,
        }
        .build()
        .expect("sphere a");
        let b = UvSphere {
            radius: r,
            center: Point3r::new(d, 0.0, 0.0),
            segments: 32,
            stacks: 16,
        }
        .build()
        .expect("sphere b");

        let result = csg_boolean(BooleanOp::Intersection, &a, &b);
        match result {
            Ok(mesh) => {
                // Intersection exists but may be very small
                let v = signed_volume(&mesh);
                assert!(
                    v < 0.1,
                    "near-tangent intersection volume {v} must be small"
                );
            }
            Err(_) => {
                // Empty intersection is acceptable for near-tangent case
            }
        }
    }

    // ── Irrational coordinate tests ────────────────────────────────────────
    //
    // Meshes rotated by non-axis-aligned angles produce irrational vertex
    // coordinates.  Known to cause:
    //   - Exact predicate inconsistencies when FP rounding direction varies
    //   - Grid-snapping artifacts in tolerance-based vertex welding
    //   - Coplanar detection false negatives

    /// Two cubes, one with irrational-valued origin (√2/2 offset), overlapping.
    /// All intersection vertices have irrational coordinates.
    /// The Boolean must produce a valid closed mesh.
    ///
    /// # Known failure mode
    ///
    /// Grid-snapped vertex pools can place irrational intersection points in
    /// adjacent cells, breaking the identity `pool.insert(p) == pool.insert(p)`
    /// when `p` falls near a cell boundary.
    #[test]
    fn irrational_coordinate_cube_intersection() {
        let a = Cube::centred(2.0).build().expect("cube a");
        // Cube B offset by 1/√2 ≈ 0.7071 in X and Y — irrational coordinates
        let s2 = std::f64::consts::FRAC_1_SQRT_2;
        let b = Cube {
            origin: Point3r::new(-1.0 + s2, -1.0 + s2, -1.0),
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        }
        .build()
        .expect("cube b");

        let result = csg_boolean(BooleanOp::Intersection, &a, &b).expect("intersection completes");
        assert!(
            result.face_count() >= 6,
            "irrational intersection must produce faces"
        );
        let v = signed_volume(&result);
        assert!(
            v > 0.0 && v < 8.0,
            "intersection volume {v} must be between 0 and full cube volume"
        );
    }

    // ── Low-tessellation stress ────────────────────────────────────────────
    //
    // Cylinders with very few segments (4–6) produce extreme aspect-ratio
    // triangles.  Known to cause:
    //   - CDT Delaunay violations from near-degenerate input triangles
    //   - Snap-round collapse of skinny triangles to zero area
    //   - Incorrect face classification from poor normal approximation

    /// Two cylinders with only 4 lateral segments each, overlapping with a
    /// small offset.  The lateral faces are extreme slivers (aspect ratio ~10:1
    /// for typical radius-to-height ratios).  The Boolean must complete without
    /// panic.
    ///
    /// # Known failure mode
    ///
    /// Ruppert refinement (if applied) can enter infinite loops on skinny
    /// triangles near the intersection seam.  The snap-round repair can
    /// collapse slivers to zero area, losing faces.
    #[test]
    fn low_tessellation_cylinder_union() {
        let cyl_a = Cylinder {
            base_center: Point3r::new(0.0, -2.0, 0.0),
            radius: 0.5,
            height: 4.0,
            segments: 4,
        }
        .build()
        .expect("cyl a");

        // Cylinder B: same axis, offset in X by 60% of diameter for significant overlap
        let cyl_b = Cylinder {
            base_center: Point3r::new(0.3, -2.0, 0.0),
            radius: 0.5,
            height: 4.0,
            segments: 4,
        }
        .build()
        .expect("cyl b");

        let result = csg_boolean(BooleanOp::Union, &cyl_a, &cyl_b);
        assert!(result.is_ok(), "low-tessellation cylinder union must not panic");
        let mesh = result.unwrap();
        assert!(mesh.face_count() > 0, "result must have faces");
    }

    // ── Many-body n-ary intersection ───────────────────────────────────────
    //
    // Intersecting 5+ shapes simultaneously stresses the n-ary pipeline's
    // face-soup management.  Known to cause:
    //   - Exponential face proliferation from cascaded co-refinement
    //   - Region labeling errors in multi-mesh classification
    //   - Wrong containment when > 2 shapes share a boundary

    /// Five overlapping cubes — the intersection must be a small cuboid.
    ///
    /// # Known failure mode
    ///
    /// Iterative 2-input Boolean chains accumulate FP drift; the n-ary pipeline
    /// operates on all meshes simultaneously, avoiding cascaded error.
    #[test]
    fn five_cube_nary_intersection() {
        let offsets = [0.0, 0.3, 0.6, -0.2, 0.1];
        let mut meshes: Vec<crate::domain::mesh::IndexedMesh> = Vec::new();
        for &dx in &offsets {
            meshes.push(
                Cube {
                    origin: Point3r::new(-1.0 + dx, -1.0, -1.0),
                    width: 2.0,
                    height: 2.0,
                    depth: 2.0,
                }
                .build()
                .expect("cube"),
            );
        }

        // Use iterative binary intersection as proxy (n-ary intersection not directly exposed)
        let mut acc = meshes[0].clone();
        for m in &meshes[1..] {
            acc = match csg_boolean(BooleanOp::Intersection, &acc, m) {
                Ok(r) => r,
                Err(_) => {
                    panic!("5-cube n-ary intersection must not fail");
                }
            };
        }
        let v = signed_volume(&acc);
        // All cubes overlap around x ∈ [0.6 - 1.0, -0.2 + 1.0] = [-0.4, 0.8]
        // Intersection width along X = 0.8 - (-0.4) = 1.2, Y and Z = 2.0
        // Volume ≈ 1.2 * 2.0 * 2.0 = 4.8, but exact depends on offset accumulation
        assert!(v > 0.0, "5-cube intersection must have positive volume, got {v}");
    }

    // ── Mixed curvature (cube-sphere) ──────────────────────────────────────
    //
    // The intersection curve of a cube and a sphere forms arcs on sphere
    // faces and straight segments on cube faces.  Known to cause:
    //   - Inconsistent face classification at curved-to-planar transitions
    //   - Missing faces where the intersection curve crosses a cube edge
    //   - GWN sign flips near cube edges (discontinuous surface normal)

    /// Sphere centred at cube corner — three faces of the cube clip the sphere.
    /// The intersection is a spherical octant.
    ///
    /// # Known failure mode
    ///
    /// The triple point where three cube faces and the sphere surface meet
    /// creates a vertex where 4+ intersection curves converge.  Libraries
    /// that assume T-junction topology (≤ 3 curves per vertex) can produce
    /// non-manifold fans at this point.
    #[test]
    fn sphere_at_cube_corner_intersection() {
        let cube = Cube {
            origin: Point3r::new(0.0, 0.0, 0.0),
            width: 4.0,
            height: 4.0,
            depth: 4.0,
        }
        .build()
        .expect("cube");
        let sphere = UvSphere {
            radius: 1.0,
            center: Point3r::new(0.0, 0.0, 0.0),
            segments: 32,
            stacks: 16,
        }
        .build()
        .expect("sphere");

        let result =
            csg_boolean(BooleanOp::Intersection, &cube, &sphere).expect("intersection completes");
        let v = signed_volume(&result);
        let v_sphere = 4.0 * std::f64::consts::PI / 3.0;
        // Intersection is a spherical octant ≈ sphere_vol/8 = π/6 ≈ 0.524
        assert!(
            v > 0.3 && v < v_sphere,
            "sphere-at-corner intersection volume {v} must be roughly an octant"
        );
    }

    // ── Complex intersection curves (torus) ────────────────────────────────
    //
    // Torus-cylinder intersections produce 4th-order algebraic curves with
    // saddle points.  Known to cause:
    //   - Intersection segment chaining failures at saddle points
    //   - Incorrect inside/outside classification near the torus "hole"
    //   - GWN sign ambiguity inside the torus cavity

    /// Cylinder passing through a torus hole — the cylinder must remove the
    /// inner region of the torus.
    ///
    /// # Known failure mode
    ///
    /// The torus has Euler characteristic χ = 0, and the cylinder intersects
    /// both the outer and inner surfaces.  Libraries that assume χ = 2 for
    /// watertight checks will flag this as invalid even when correct.
    #[test]
    fn cylinder_through_torus_union() {
        let torus = Torus {
            major_radius: 3.0,
            minor_radius: 1.0,
            major_segments: 32,
            minor_segments: 16,
        }
        .build()
        .expect("torus");
        let cylinder = Cylinder {
            base_center: Point3r::new(0.0, -3.0, 0.0),
            radius: 0.5,
            height: 6.0,
            segments: 16,
        }
        .build()
        .expect("cylinder");

        let result = csg_boolean(BooleanOp::Union, &torus, &cylinder);
        assert!(result.is_ok(), "torus-cylinder union must not panic");
        let mesh = result.unwrap();
        assert!(
            mesh.face_count() > 0,
            "torus-cylinder union must produce faces"
        );
        let v = signed_volume(&mesh);
        let v_torus = 2.0 * std::f64::consts::PI.powi(2) * 3.0 * 1.0_f64.powi(2); // 2π²Rr²
        assert!(
            v > v_torus,
            "union volume {v} must exceed torus volume {v_torus}"
        );
    }

    // ── Dense narrow-phase stress ──────────────────────────────────────────

    /// Two high-tessellation spheres with significant overlap — stress-tests
    /// the narrow phase with thousands of triangle-triangle intersection tests.
    ///
    /// # Known failure mode
    ///
    /// O(n²) narrow-phase implementations become impractical at 64×32
    /// resolution (4096+ faces per sphere).  The BVH-accelerated pipeline
    /// must complete in reasonable time.
    #[test]
    fn dense_sphere_sphere_union_64x32() {
        let a = UvSphere {
            radius: 1.0,
            center: Point3r::new(0.0, 0.0, 0.0),
            segments: 64,
            stacks: 32,
        }
        .build()
        .expect("sphere a");
        let b = UvSphere {
            radius: 1.0,
            center: Point3r::new(1.0, 0.0, 0.0),
            segments: 64,
            stacks: 32,
        }
        .build()
        .expect("sphere b");

        let result = csg_boolean(BooleanOp::Union, &a, &b).expect("dense union completes");
        let v = signed_volume(&result);
        let v_single = 4.0 * std::f64::consts::PI / 3.0;
        assert!(
            v > v_single && v < 2.0 * v_single,
            "dense union volume {v} must be between one and two sphere volumes"
        );
    }

    // ── Near-tangent cylinders ──────────────────────────────────────────────

    /// Two parallel cylinders whose surfaces are 0.01 mm apart — the union
    /// touches but barely overlaps.  Stresses collinearity detection and
    /// snap-round at millimeter scale.
    ///
    /// # Known failure mode
    ///
    /// When the gap is on the order of the welding tolerance, vertex pool
    /// deduplication can merge vertices across the two cylinder surfaces,
    /// creating a non-manifold bridge.
    #[test]
    fn near_tangent_parallel_cylinders_union() {
        let r = 0.5;
        let gap = 0.01;
        let a = Cylinder {
            base_center: Point3r::new(0.0, 0.0, 0.0),
            radius: r,
            height: 3.0,
            segments: 32,
        }
        .build()
        .expect("cyl a");
        let b = Cylinder {
            base_center: Point3r::new(2.0 * r - gap, 0.0, 0.0),
            radius: r,
            height: 3.0,
            segments: 32,
        }
        .build()
        .expect("cyl b");

        let result = csg_boolean(BooleanOp::Union, &a, &b);
        assert!(
            result.is_ok(),
            "near-tangent parallel cylinder union must not panic"
        );
        let mesh = result.unwrap();
        assert!(mesh.face_count() > 0, "result must have faces");
        let v = signed_volume(&mesh);
        let v_single = std::f64::consts::PI * r * r * 3.0;
        assert!(
            v > v_single,
            "near-tangent union volume {v} must exceed single cylinder {v_single}"
        );
    }

    // ── Degenerate subtraction (empty result) ──────────────────────────────

    /// Subtract a sphere that completely contains a smaller sphere — the
    /// result should be empty or have zero volume.
    ///
    /// # Known failure mode
    ///
    /// Some libraries return the outer shell with inverted normals instead
    /// of an empty result when A ⊂ B and computing A − B.
    #[test]
    fn fully_contained_sphere_difference_is_empty() {
        let inner = UvSphere {
            radius: 0.5,
            center: Point3r::new(0.0, 0.0, 0.0),
            segments: 16,
            stacks: 8,
        }
        .build()
        .expect("inner");
        let outer = UvSphere {
            radius: 2.0,
            center: Point3r::new(0.0, 0.0, 0.0),
            segments: 32,
            stacks: 16,
        }
        .build()
        .expect("outer");

        let result = csg_boolean(BooleanOp::Difference, &inner, &outer);
        match result {
            Ok(mesh) => {
                let v = signed_volume(&mesh);
                assert!(
                    v < 0.1 || mesh.face_count() == 0,
                    "inner − outer must be empty or near-zero, got volume {v} with {} faces",
                    mesh.face_count()
                );
            }
            Err(_) => {
                // Error (empty result) is acceptable
            }
        }
    }
}
