use super::*;
    use super::*;
    use crate::application::csg::boolean::{csg_boolean, csg_boolean_indexed, BooleanOp};
    use crate::application::watertight::check::{check_watertight, WatertightReport};
    use crate::domain::geometry::primitives::UvSphere;
    use crate::domain::mesh::IndexedMesh;

    /// Convenience wrapper: rebuild edges then return a full watertight report.
    fn watertight_report(mesh: &mut IndexedMesh) -> WatertightReport {
        mesh.rebuild_edges();
        check_watertight(&mesh.vertices, &mesh.faces, mesh.edges_ref().unwrap())
    }

    /// Low-level boolean without the watertight post-check.
    ///
    /// Use this for edge-case regression tests where the mesh is known to be
    /// geometrically correct but may have topological defects (e.g. coplanar
    /// caps that currently produce boundary seams).
    fn boolean_raw(op: BooleanOp, mesh_a: &IndexedMesh, mesh_b: &IndexedMesh) -> IndexedMesh {
        use crate::domain::core::index::VertexId;
        use crate::infrastructure::storage::face_store::FaceData;
        use crate::infrastructure::storage::vertex_pool::VertexPool;
        use std::collections::HashMap;
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
        let faces_a: Vec<FaceData> = mesh_a
            .faces
            .iter()
            .map(|f| FaceData {
                vertices: f.vertices.map(|v| remap_a[&v]),
                region: f.region,
            })
            .collect();
        let faces_b: Vec<FaceData> = mesh_b
            .faces
            .iter()
            .map(|f| FaceData {
                vertices: f.vertices.map(|v| remap_b[&v]),
                region: f.region,
            })
            .collect();
        let result_faces = csg_boolean(op, &faces_a, &faces_b, &mut combined)
            .expect("csg_boolean should not error");
        super::super::reconstruct::reconstruct_mesh(&result_faces, &combined)
    }

    /// Build a UV sphere centred at `(cx, cy, cz)` with radius `r` and the
    /// given latitude/longitude resolution.
    fn make_sphere(
        cx: f64,
        cy: f64,
        cz: f64,
        r: f64,
        stacks: usize,
        segments: usize,
    ) -> IndexedMesh {
        use crate::domain::geometry::primitives::PrimitiveMesh;
        let sphere = UvSphere {
            radius: r,
            center: Point3r::new(cx, cy, cz),
            segments,
            stacks,
        };
        sphere.build().expect("UvSphere::build failed")
    }

    /// Compute the signed volume of an `IndexedMesh` using the divergence theorem.
    ///
    /// `vol = (1/6) * ÃŽÂ£_face  (v0 Ã‚Â· (v1 Ãƒâ€” v2))`
    fn signed_volume(mesh: &IndexedMesh) -> f64 {
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

    /// Sphere-sphere intersection (lens) volume test.
    ///
    /// Two unit spheres (r = 1.0) with centres 1.0 apart.
    /// Analytical lens volume = Ãâ‚¬/12 * (4Ã‚Â·r + d) * (2Ã‚Â·r Ã¢Ë†â€™ d)Ã‚Â²
    ///   where r=1, d=1  Ã¢â€ â€™  Ãâ‚¬/12 * 5 * 1 = 5Ãâ‚¬/12 Ã¢â€°Ë† 1.3090
    ///
    /// Uses 64Ãƒâ€”32 resolution: the lens boundary (intersection seam) is more
    /// finely sampled than the union/difference cases, requiring higher mesh
    /// resolution to reach < 5% discretisation error.
    #[test]
    fn sphere_sphere_intersection_volume() {
        // 64 stacks Ãƒâ€” 32 segments gives dense enough triangulation at the seam
        // to keep discretisation error below 5%.
        let stacks = 64;
        let segments = 32;
        let sphere_a = make_sphere(0.0, 0.0, 0.0, 1.0, stacks, segments);
        let sphere_b = make_sphere(1.0, 0.0, 0.0, 1.0, stacks, segments);

        let result = csg_boolean_indexed(BooleanOp::Intersection, &sphere_a, &sphere_b)
            .expect("sphere-sphere intersection should not fail");

        let vol = signed_volume(&result);
        let expected = 5.0 * std::f64::consts::PI / 12.0; // Ã¢â€°Ë† 1.3090
        let err = (vol - expected).abs() / expected;
        assert!(
            err < 0.05,
            "sphere-sphere intersection volume error {:.1}% > 5% (got {:.4}, expected {:.4})",
            err * 100.0,
            vol,
            expected
        );
    }

    /// Sphere-sphere union volume test.
    ///
    /// vol(A Ã¢Ë†Âª B) = vol(A) + vol(B) Ã¢Ë†â€™ vol(A Ã¢Ë†Â© B)
    ///   = 2 * 4Ãâ‚¬/3 Ã¢Ë†â€™ 5Ãâ‚¬/12 Ã¢â€°Ë† 7.1272
    #[test]
    fn sphere_sphere_union_volume() {
        let stacks = 32;
        let segments = 16;
        let sphere_a = make_sphere(0.0, 0.0, 0.0, 1.0, stacks, segments);
        let sphere_b = make_sphere(1.0, 0.0, 0.0, 1.0, stacks, segments);

        let result = csg_boolean_indexed(BooleanOp::Union, &sphere_a, &sphere_b)
            .expect("sphere-sphere union should not fail");

        let vol = signed_volume(&result);
        let vol_sphere = 4.0 * std::f64::consts::PI / 3.0;
        let vol_lens = 5.0 * std::f64::consts::PI / 12.0;
        let expected = 2.0 * vol_sphere - vol_lens; // Ã¢â€°Ë† 7.1272
        let err = (vol - expected).abs() / expected;
        assert!(
            err < 0.05,
            "sphere-sphere union volume error {:.1}% > 5% (got {:.4}, expected {:.4})",
            err * 100.0,
            vol,
            expected
        );
    }

    /// Sphere-sphere difference volume test.
    ///
    /// vol(A \ B) = vol(A) Ã¢Ë†â€™ vol(A Ã¢Ë†Â© B)
    ///   = 4Ãâ‚¬/3 Ã¢Ë†â€™ 5Ãâ‚¬/12 Ã¢â€°Ë† 2.8798
    #[test]
    fn sphere_sphere_difference_volume() {
        let stacks = 32;
        let segments = 16;
        let sphere_a = make_sphere(0.0, 0.0, 0.0, 1.0, stacks, segments);
        let sphere_b = make_sphere(1.0, 0.0, 0.0, 1.0, stacks, segments);

        let result = csg_boolean_indexed(BooleanOp::Difference, &sphere_a, &sphere_b)
            .expect("sphere-sphere difference should not fail");

        let vol = signed_volume(&result);
        let vol_sphere = 4.0 * std::f64::consts::PI / 3.0;
        let vol_lens = 5.0 * std::f64::consts::PI / 12.0;
        let expected = vol_sphere - vol_lens; // Ã¢â€°Ë† 2.8798
        let err = (vol - expected).abs() / expected;
        assert!(
            err < 0.05,
            "sphere-sphere difference volume error {:.1}% > 5% (got {:.4}, expected {:.4})",
            err * 100.0,
            vol,
            expected
        );
    }

    /// Regression test: equal-height parallel cylinders with exactly coplanar caps.
    ///
    /// Two cylinders of equal radius r = 0.4 mm and equal height h = 3 mm,
    /// with axes offset by d = r = 0.4 mm in X.  Both base caps sit at
    /// exactly y = Ã¢Ë†â€™1.5 and both top caps at y = +1.5.
    ///
    /// Before the fix: the coplanar cap faces were routed through the
    /// Sutherland-Hodgman clip, which duplicated them (orient_3d Degenerate =
    /// "inside"), causing the union volume to exceed V_A + V_B.
    ///
    /// After the fix: coplanar pairs bypass clipping and are classified by
    /// normal alignment.  Co-oriented cap pairs keep exactly one copy, so:
    ///   vol(A Ã¢Ë†Âª B) < V_A + V_B   (overlap correctly subtracted)
    ///   vol(A Ã¢Ë†Â© B) > 0            (lens region positive)
    #[test]
    fn coplanar_caps_no_double_counting() {
        use crate::domain::geometry::primitives::{Cylinder, PrimitiveMesh};
        let r = 0.4_f64;
        let h = 3.0_f64;
        // Both cylinders: Y Ã¢Ë†Ë† [Ã¢Ë†â€™1.5, 1.5] Ã¢â‚¬â€ caps coplanar at y = Ã‚Â±1.5.
        let cyl_a = Cylinder {
            base_center: Point3r::new(-r / 2.0, -h / 2.0, 0.0),
            radius: r,
            height: h,
            segments: 32,
        }
        .build()
        .expect("cyl_a build failed");
        let cyl_b = Cylinder {
            base_center: Point3r::new(r / 2.0, -h / 2.0, 0.0),
            radius: r,
            height: h,
            segments: 32,
        }
        .build()
        .expect("cyl_b build failed");

        let v_each = std::f64::consts::PI * r * r * h; // Ã¢â€°Ë† 1.5080

        // Use the low-level raw path: coplanar caps currently produce boundary
        // seams (known limitation), so the watertight-enforcing high-level API
        // would return an error.  This test validates only the volume invariant.
        let union_mesh = boolean_raw(BooleanOp::Union, &cyl_a, &cyl_b);
        let union_vol = signed_volume(&union_mesh);

        // Union must be strictly less than V_A + V_B (not double-counted).
        assert!(
            union_vol < 2.0 * v_each - 0.01,
            "union vol {:.4} should be < 2Ã‚Â·V_each = {:.4} (coplanar caps double-counted)",
            union_vol,
            2.0 * v_each,
        );
        // And it must be positive.
        assert!(
            union_vol > 0.1,
            "union vol {:.4} should be positive",
            union_vol
        );

        // Intersection must be positive (the lens barrel is non-empty).
        let inter_mesh = boolean_raw(BooleanOp::Intersection, &cyl_a, &cyl_b);
        let inter_vol = signed_volume(&inter_mesh);
        assert!(
            inter_vol > 0.05,
            "intersection vol {:.4} should be positive",
            inter_vol
        );

        // Inclusion-exclusion: vol(A) + vol(B) Ã¢â€°Ë† vol(AÃ¢Ë†ÂªB) + vol(AÃ¢Ë†Â©B)
        // Allow 10% tolerance for triangulation discretisation.
        let ie_lhs = 2.0 * v_each;
        let ie_rhs = union_vol + inter_vol;
        let ie_err = (ie_lhs - ie_rhs).abs() / ie_lhs;
        assert!(
            ie_err < 0.10,
            "inclusion-exclusion error {:.1}% > 10% (lhs={:.4}, rhs={:.4})",
            ie_err * 100.0,
            ie_lhs,
            ie_rhs,
        );
    }

    /// Regression test: ensure faces that are almost but not exactly coplanar
    /// do not get grouped together by a quantization heuristic.
    #[test]
    fn exact_algebraic_coplanarity_no_shatter() {
        let _pool: VertexPool<f64> = VertexPool::default_millifluidic();

        // Let's create `IndexedMesh`es so they run through the full `BooleanOp` pipeline
        // assuring broad-phase overlap tests correctly dispatch to intersection.
        let mut mesh_a = IndexedMesh::new();
        let mut mesh_b = IndexedMesh::new();

        let a1 = mesh_a.add_vertex(Point3r::new(0.0, 0.0, 0.0), Vector3r::z());
        let a2 = mesh_a.add_vertex(Point3r::new(1.0, 0.0, 0.0), Vector3r::z());
        let a3 = mesh_a.add_vertex(Point3r::new(0.0, 1.0, 0.0), Vector3r::z());
        mesh_a.add_face(a1, a2, a3);

        let b1 = mesh_b.add_vertex(Point3r::new(0.0, 0.0, 0.0000001), Vector3r::z());
        let b2 = mesh_b.add_vertex(Point3r::new(1.0, 0.0, 0.0000001), Vector3r::z());
        let b3 = mesh_b.add_vertex(Point3r::new(0.0, 1.0, 0.0), Vector3r::z());
        mesh_b.add_face(b1, b2, b3);

        let result = csg_boolean_indexed(BooleanOp::Intersection, &mesh_a, &mesh_b);

        assert!(result.is_ok(), "Intersection evaluation should not error");
        let result_mesh = result.unwrap();
        // Since they only intersect at exactly one line inside the 2D plane (y=0, z~=0),
        // there is no 3D overlap, and no coplanar 2D area overlap. The intersection mesh correctly
        // produces the exact intersection artifacts (typically 1 or 2 seam fragments) rather than completely
        // vanishing or duplicating due to a coplanar grouping heuristic.
        assert!(
            result_mesh.faces.len() <= 2,
            "Non-coplanar tilted faces should not output massive coplanar intersection"
        );
    }

    /// Diagnostic: perpendicular cylinder union should be watertight.
    ///
    /// Two equal-radius cylinders (r=0.5, h=3) with perpendicular axes crossing
    /// at the centre.  The arrangement pipeline should produce a closed manifold.
    #[test]
    fn perpendicular_cylinder_union_is_watertight() {
        use crate::domain::geometry::primitives::{Cylinder, PrimitiveMesh};
        let r = 0.5_f64;
        let h = 3.0_f64;
        let cyl_a = Cylinder {
            base_center: Point3r::new(0.0, -h / 2.0, 0.0),
            radius: r,
            height: h,
            segments: 16,
        }
        .build()
        .expect("cyl_a build");
        let cyl_b = Cylinder {
            base_center: Point3r::new(-h / 2.0, 0.0, 0.0),
            radius: r,
            height: h,
            segments: 16,
        }
        .build()
        .expect("cyl_b build");

        let mut result =
            csg_boolean_indexed(BooleanOp::Union, &cyl_a, &cyl_b).expect("union should not fail");
        let report = watertight_report(&mut result);
        assert!(
            report.is_watertight,
            "perpendicular cylinder union should be watertight \
             (boundary_edges={}, non_manifold={})",
            report.boundary_edge_count, report.non_manifold_edge_count
        );
    }

    /// Diagnostic: parallel cylinder union (same height Ã¢â€ â€™ coplanar caps) should be watertight.
    #[test]
    fn parallel_cylinder_union_is_watertight() {
        use crate::domain::geometry::primitives::{Cylinder, PrimitiveMesh};
        let r = 0.5_f64;
        let h = 3.0_f64;
        let cyl_a = Cylinder {
            base_center: Point3r::new(0.0, 0.0, 0.0),
            radius: r,
            height: h,
            segments: 16,
        }
        .build()
        .expect("cyl_a build");
        let cyl_b = Cylinder {
            base_center: Point3r::new(r, 0.0, 0.0),
            radius: r,
            height: h,
            segments: 16,
        }
        .build()
        .expect("cyl_b build");

        let mut result =
            csg_boolean_indexed(BooleanOp::Union, &cyl_a, &cyl_b).expect("union should not fail");
        let report = watertight_report(&mut result);

        // Debug: print boundary edge positions if any
        if !report.is_watertight {
            result.rebuild_edges();
            let edges = result.edges_ref().unwrap();
            let mut boundary_positions: Vec<(Point3r, Point3r, f64)> = Vec::new();
            for edge in edges.iter() {
                if edge.valence() == 1 {
                    let pa = *result.vertices.position(edge.vertices.0);
                    let pb = *result.vertices.position(edge.vertices.1);
                    let mid_z = (pa.z + pb.z) / 2.0;
                    boundary_positions.push((pa, pb, mid_z));
                }
            }
            boundary_positions
                .sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal));
            eprintln!("=== Boundary edges ({}) ===", boundary_positions.len());
            for (a, b, _z) in boundary_positions.iter() {
                let mid_y = (a.y + b.y) / 2.0;
                eprintln!(
                    "  yÃ¢â€°Ë†{:.4}  ({:.4},{:.4},{:.4}) Ã¢â€ â€™ ({:.4},{:.4},{:.4})",
                    mid_y, a.x, a.y, a.z, b.x, b.y, b.z
                );
            }
        }

        // Debug: print non-manifold edge positions if any
        if report.non_manifold_edge_count > 0 {
            let edges = result.edges_ref().unwrap();
            for edge in edges.iter() {
                if edge.valence() > 2 {
                    let pa = *result.vertices.position(edge.vertices.0);
                    let pb = *result.vertices.position(edge.vertices.1);
                    eprintln!(
                        "NM_EDGE val={}: ({:.6},{:.6},{:.6})->({:.6},{:.6},{:.6})",
                        edge.valence(),
                        pa.x,
                        pa.y,
                        pa.z,
                        pb.x,
                        pb.y,
                        pb.z
                    );
                }
            }
        }

        assert!(
            report.is_watertight,
            "parallel cylinder union should be watertight \
             (boundary_edges={}, non_manifold={})",
            report.boundary_edge_count, report.non_manifold_edge_count
        );

        // Volume: A u B = 2*V_cyl - V_lens.  d = r => theta = pi/3.
        let theta = (r / (2.0 * r)).acos();
        let a_seg = r * r * (theta - theta.sin() * theta.cos());
        let v_union_expected = 2.0 * std::f64::consts::PI * r * r * h - 2.0 * h * a_seg;
        let vol = result.signed_volume();
        assert!(vol > 0.0, "orientation inverted (vol={vol:.6})");
        assert!(
            (vol - v_union_expected).abs() / v_union_expected < 0.05,
            "volume error >5%: vol={vol:.6} expected={v_union_expected:.6}"
        );
    }

    /// Regression: 64-segment parallel cylinder union volume (r=0.6, h=3.0, d=r).
    ///
    /// Exercises the exact geometry from the `csg_cylinder_cylinder` example.
    /// Verifies cop_faces direct emission produces correct volume for coplanar-cap case.
    #[test]
    fn cylinder_cylinder_union_64seg_volume() {
        use crate::domain::geometry::primitives::{Cylinder, PrimitiveMesh};
        let r = 0.6_f64;
        let h = 3.0_f64;
        let d = r;
        let cyl_a = Cylinder {
            base_center: Point3r::new(-d / 2.0, -h / 2.0, 0.0),
            radius: r,
            height: h,
            segments: 64,
        }
        .build()
        .expect("cyl_a");
        let cyl_b = Cylinder {
            base_center: Point3r::new(d / 2.0, -h / 2.0, 0.0),
            radius: r,
            height: h,
            segments: 64,
        }
        .build()
        .expect("cyl_b");
        let mut result = csg_boolean_indexed(BooleanOp::Union, &cyl_a, &cyl_b).expect("union");
        // Cylindrical lens: d = r => theta = arccos(0.5) = pi/3
        let theta = (d / (2.0 * r)).acos();
        let a_seg = r * r * (theta - theta.sin() * theta.cos());
        let v_intersect = 2.0 * h * a_seg;
        let v_cyl = std::f64::consts::PI * r * r * h;
        let expected = 2.0 * v_cyl - v_intersect;
        let vol = result.signed_volume();
        let report = watertight_report(&mut result);
        assert!(
            report.is_watertight,
            "64-seg union should be watertight (boundary={}, non_manifold={})",
            report.boundary_edge_count, report.non_manifold_edge_count
        );
        assert!(vol > 0.0, "orientation inverted (vol={vol:.6})");
        assert!(
            (vol - expected).abs() / expected < 0.05,
            "volume error >5%: vol={vol:.6} expected={expected:.6}"
        );
    }

    /// Diagnostic: asymmetric cylinder union (different heights Ã¢â€ â€™ non-coplanar caps).
    #[test]
    fn asymmetric_cylinder_union_is_watertight() {
        use crate::domain::geometry::primitives::{Cylinder, PrimitiveMesh};
        let r = 0.6_f64;
        let cyl_a = Cylinder {
            base_center: Point3r::new(-0.3, -1.5, 0.0),
            radius: r,
            height: 3.0,
            segments: 64,
        }
        .build()
        .expect("cyl_a build");
        let cyl_b = Cylinder {
            base_center: Point3r::new(0.3, -2.0, 0.0),
            radius: r,
            height: 4.0,
            segments: 64,
        }
        .build()
        .expect("cyl_b build");

        let mut result =
            csg_boolean_indexed(BooleanOp::Union, &cyl_a, &cyl_b).expect("union should not fail");
        let report = watertight_report(&mut result);

        if !report.is_watertight {
            result.rebuild_edges();
            let edges = result.edges_ref().unwrap();
            let mut boundary_positions: Vec<(Point3r, Point3r)> = Vec::new();
            for edge in edges.iter() {
                if edge.valence() == 1 {
                    let pa = *result.vertices.position(edge.vertices.0);
                    let pb = *result.vertices.position(edge.vertices.1);
                    boundary_positions.push((pa, pb));
                }
            }
            eprintln!(
                "=== Asymmetric boundary edges ({}) ===",
                boundary_positions.len()
            );
            for (a, b) in boundary_positions.iter() {
                eprintln!(
                    "  ({:.6},{:.6},{:.6}) Ã¢â€ â€™ ({:.6},{:.6},{:.6})",
                    a.x, a.y, a.z, b.x, b.y, b.z
                );
            }
        }

        assert!(
            report.is_watertight,
            "asymmetric cylinder union should be watertight \
             (boundary_edges={}, non_manifold={})",
            report.boundary_edge_count, report.non_manifold_edge_count
        );
    }

    /// Regression test: L-shape compound union (stem ∪ elbow ∪ arm) is watertight.
    #[test]
    fn l_shape_compound_union_is_watertight() {
        use crate::application::csg::CsgNode;
        use crate::domain::core::scalar::Real;
        use crate::domain::geometry::primitives::{Cylinder, Elbow, PrimitiveMesh};
        use nalgebra::{Isometry3, Translation3, UnitQuaternion, Vector3};

        let r = 0.5_f64;
        let r_bend = 1.0_f64;
        let h = 3.0_f64;
        let eps = r * 0.05;
        let stem_len = h - r_bend;
        let arm_len = h - r_bend;

        let stem = Cylinder {
            base_center: Point3r::new(0.0, 0.0, 0.0),
            radius: r,
            height: stem_len + eps,
            segments: 32,
        }
        .build()
        .expect("stem build");

        let elbow_raw = Elbow {
            tube_radius: r,
            bend_radius: r_bend,
            bend_angle: std::f64::consts::FRAC_PI_2,
            tube_segments: 32,
            arc_segments: 16,
        }
        .build()
        .expect("elbow build");
        // L-shape example uses -90° about X (not Y): +Z → +Y, +X → +X
        let rot_elbow = UnitQuaternion::<Real>::from_axis_angle(
            &Vector3::x_axis(),
            -std::f64::consts::FRAC_PI_2,
        );
        let elbow = CsgNode::Transform {
            node: Box::new(CsgNode::Leaf(elbow_raw)),
            iso: Isometry3::from_parts(Translation3::new(0.0, stem_len, 0.0), rot_elbow),
        }
        .evaluate()
        .expect("elbow transform");

        let arm_raw = Cylinder {
            base_center: Point3r::new(0.0, 0.0, 0.0),
            radius: r,
            height: arm_len + eps,
            segments: 32,
        }
        .build()
        .expect("arm build");
        let rot_arm = UnitQuaternion::<Real>::from_axis_angle(
            &Vector3::z_axis(),
            -std::f64::consts::FRAC_PI_2,
        );
        let arm_y = h; // arm_y = R_BEND + straight_len = r_bend + stem_len = 1 + 2 = 3
        let arm = CsgNode::Transform {
            node: Box::new(CsgNode::Leaf(arm_raw)),
            iso: Isometry3::from_parts(Translation3::new(r_bend - eps, arm_y, 0.0), rot_arm),
        }
        .evaluate()
        .expect("arm transform");

        let mut stem_elbow = csg_boolean_indexed(BooleanOp::Union, &stem, &elbow)
            .expect("stem ∪ elbow should not fail");
        let rep1 = watertight_report(&mut stem_elbow);
        eprintln!(
            "stem ∪ elbow: boundary={}, non_manifold={}",
            rep1.boundary_edge_count, rep1.non_manifold_edge_count
        );

        let mut result = csg_boolean_indexed(BooleanOp::Union, &stem_elbow, &arm)
            .expect("compound ∪ arm should not fail");
        let report = watertight_report(&mut result);

        if !report.is_watertight {
            result.rebuild_edges();
            let edges = result.edges_ref().unwrap();
            let mut bvec: Vec<(Point3r, Point3r)> = Vec::new();
            for edge in edges.iter() {
                if edge.valence() == 1 {
                    bvec.push((
                        *result.vertices.position(edge.vertices.0),
                        *result.vertices.position(edge.vertices.1),
                    ));
                }
            }
            eprintln!("=== L-shape final boundary edges ({}) ===", bvec.len());
            for (a, b) in bvec.iter() {
                eprintln!(
                    "  ({:.6},{:.6},{:.6}) -> ({:.6},{:.6},{:.6})",
                    a.x, a.y, a.z, b.x, b.y, b.z
                );
            }
        }

        assert!(
            rep1.is_watertight,
            "stem ∪ elbow: boundary={}, non_manifold={}",
            rep1.boundary_edge_count, rep1.non_manifold_edge_count
        );
        assert!(
            report.is_watertight,
            "L-shape compound union: boundary={}, non_manifold={}",
            report.boundary_edge_count, report.non_manifold_edge_count
        );
    }

    /// Regression test: V-shape right_branch (right_elbow Ã¢Ë†Âª right_arm) is watertight.
    ///
    /// Uses the exact same geometry parameters as `cylinder_cylinder_v_shape.rs`
    /// but at reduced resolution (32Ãƒâ€”16) for fast test execution.
    #[test]
    fn v_shape_right_branch_is_watertight() {
        use crate::application::csg::CsgNode;
        use crate::domain::core::scalar::Real;
        use crate::domain::geometry::primitives::{Cylinder, Elbow, PrimitiveMesh};
        use nalgebra::{Isometry3, Translation3, UnitQuaternion, Vector3};

        let r = 0.5_f64;
        let r_bend = 2.0 * r; // = 1.0 mm
        let h = 3.0_f64;
        let theta = std::f64::consts::PI / 6.0; // 30Ã‚Â°
        let (s_th, c_th) = theta.sin_cos();
        let axial_reach = r_bend * s_th;
        let radial_reach = r_bend * (1.0 - c_th);
        let stem_len = h - axial_reach;
        let arm_len = h - radial_reach;
        let eps = r * 0.10;

        // Elbow: tube_segments=32, arc_segments=16 (half of example's 64Ãƒâ€”32)
        let elbow_raw = Elbow {
            tube_radius: r,
            bend_radius: r_bend,
            bend_angle: theta,
            tube_segments: 32,
            arc_segments: 16,
        }
        .build()
        .expect("elbow build");

        // Elbow isometry: -90Ã‚Â° about X + translate to elbow inlet y = -H + stem_len = -axial_reach
        let elbow_inlet_y = -h + stem_len;
        let rot_base = UnitQuaternion::<Real>::from_axis_angle(
            &Vector3::x_axis(),
            -std::f64::consts::FRAC_PI_2,
        );
        let right_elbow = CsgNode::Transform {
            node: Box::new(CsgNode::Leaf(elbow_raw)),
            iso: Isometry3::from_parts(Translation3::new(0.0, elbow_inlet_y, 0.0), rot_base),
        }
        .evaluate()
        .expect("right_elbow transform");

        // Arm: 32 segments, along right branch direction (sinÃŽÂ¸, cosÃŽÂ¸, 0)
        let arm_raw = Cylinder {
            base_center: Point3r::new(0.0, 0.0, 0.0),
            radius: r,
            height: arm_len + eps,
            segments: 32,
        }
        .build()
        .expect("arm build");
        let rot_arm = UnitQuaternion::<Real>::from_axis_angle(&Vector3::z_axis(), -theta);
        let tx = radial_reach - eps * s_th;
        let ty = -eps * c_th;
        let right_arm = CsgNode::Transform {
            node: Box::new(CsgNode::Leaf(arm_raw)),
            iso: Isometry3::from_parts(Translation3::new(tx, ty, 0.0), rot_arm),
        }
        .evaluate()
        .expect("right_arm transform");

        let mut result = csg_boolean_indexed(BooleanOp::Union, &right_elbow, &right_arm)
            .expect("right_elbow Ã¢Ë†Âª right_arm should not fail");
        let report = watertight_report(&mut result);

        if !report.is_watertight {
            result.rebuild_edges();
            let edges = result.edges_ref().unwrap();
            let mut boundary_positions: Vec<(Point3r, Point3r)> = Vec::new();
            for edge in edges.iter() {
                if edge.valence() == 1 {
                    boundary_positions.push((
                        *result.vertices.position(edge.vertices.0),
                        *result.vertices.position(edge.vertices.1),
                    ));
                }
            }
            eprintln!(
                "=== right_branch boundary edges ({}) ===",
                boundary_positions.len()
            );
            for (a, b) in boundary_positions.iter() {
                eprintln!(
                    "  ({:.6},{:.6},{:.6}) Ã¢â€ â€™ ({:.6},{:.6},{:.6})",
                    a.x, a.y, a.z, b.x, b.y, b.z
                );
            }
        }

        assert!(
            report.is_watertight,
            "right_elbow Ã¢Ë†Âª right_arm should be watertight \
             (boundary_edges={}, non_manifold={})",
            report.boundary_edge_count, report.non_manifold_edge_count
        );
    }

    /// Regression test: 90Ã‚Â° elbow union with a straight arm cylinder.
    ///
    /// This tests the elbow (torus-segment) + cylinder union Ã¢â‚¬â€ a more complex
    /// curved mesh operation than cylinder-cylinder.  The arm cylinder cap
    /// penetrates the elbow barrel, requiring `propagate_seam_vertices` to handle
    /// crossings at multiple elbow ring edges.
    #[test]
    fn elbow_cylinder_union_is_watertight() {
        use crate::domain::core::scalar::Real;
        use crate::domain::geometry::primitives::{Cylinder, Elbow, PrimitiveMesh};
        use nalgebra::{Isometry3, Translation3, UnitQuaternion, Vector3};

        let r = 0.5_f64;
        let r_bend = 1.0_f64;
        let h = 2.0_f64;
        let eps = r * 0.10;

        // 90Ã‚Â° elbow: inlet +Z, outlet +X.  Place in canonical position.
        // V-shape parameters: 30Ã‚Â° half-angle
        let theta = std::f64::consts::PI / 6.0; // 30Ã‚Â°
        let (s_th, c_th) = theta.sin_cos();
        let axial_reach = r_bend * s_th;
        let radial_reach = r_bend * (1.0 - c_th);
        let arm_len = h - radial_reach;

        let elbow = Elbow {
            tube_radius: r,
            bend_radius: r_bend,
            bend_angle: theta,
            tube_segments: 64,
            arc_segments: 32,
        }
        .build()
        .expect("elbow build");

        // Apply elbow isometry: -90Ã‚Â° about X, then translate to elbow inlet position.
        let elbow_inlet_y = -(h - axial_reach);
        let rot_base = UnitQuaternion::<Real>::from_axis_angle(
            &Vector3::x_axis(),
            -std::f64::consts::FRAC_PI_2,
        );
        use crate::application::csg::CsgNode;
        let elbow = CsgNode::Transform {
            node: Box::new(CsgNode::Leaf(elbow)),
            iso: Isometry3::from_parts(Translation3::new(0.0, elbow_inlet_y, 0.0), rot_base),
        }
        .evaluate()
        .expect("elbow transform");

        // Arm cylinder: along right branch direction (sinÃŽÂ¸, cosÃŽÂ¸, 0).
        // Base starts eps before elbow outlet.
        let arm_raw = Cylinder {
            base_center: Point3r::new(0.0, 0.0, 0.0),
            radius: r,
            height: arm_len + eps,
            segments: 64,
        }
        .build()
        .expect("arm build");

        // Rotate arm: +Y Ã¢â€ â€™ right branch direction (sinÃŽÂ¸, cosÃŽÂ¸, 0)
        let rot = UnitQuaternion::<Real>::from_axis_angle(&Vector3::z_axis(), -theta);
        let tx = radial_reach - eps * s_th;
        let ty = -eps * c_th;
        let iso = Isometry3::from_parts(Translation3::new(tx, ty, 0.0), rot);
        let arm = CsgNode::Transform {
            node: Box::new(CsgNode::Leaf(arm_raw)),
            iso,
        }
        .evaluate()
        .expect("arm transform");

        let mut result = csg_boolean_indexed(BooleanOp::Union, &elbow, &arm)
            .expect("elbow Ã¢Ë†Âª arm should not fail");
        let report = watertight_report(&mut result);

        if !report.is_watertight {
            result.rebuild_edges();
            let edges = result.edges_ref().unwrap();
            let mut boundary_positions: Vec<(Point3r, Point3r)> = Vec::new();
            for edge in edges.iter() {
                if edge.valence() == 1 {
                    boundary_positions.push((
                        *result.vertices.position(edge.vertices.0),
                        *result.vertices.position(edge.vertices.1),
                    ));
                }
            }
            eprintln!(
                "=== Elbow+Arm boundary edges ({}) ===",
                boundary_positions.len()
            );
            for (a, b) in boundary_positions.iter() {
                eprintln!(
                    "  ({:.6},{:.6},{:.6}) Ã¢â€ â€™ ({:.6},{:.6},{:.6})",
                    a.x, a.y, a.z, b.x, b.y, b.z
                );
            }
        }

        assert!(
            report.is_watertight,
            "elbow ∪ arm should be watertight \
             (boundary_edges={}, non_manifold={})",
            report.boundary_edge_count, report.non_manifold_edge_count
        );
    }

    /// Regression test: V-shape `right_elbow ∪ right_arm` at 64×32 (exact example params).
    ///
    /// Uses the exact geometry from `cylinder_cylinder_v_shape.rs` `run_rounded()`.
    /// R=0.5, H=3.0, THETA=π/6, R_BEND=1.0, tube_segments=64, arc_segments=32.
    #[test]
    fn v_shape_right_branch_64x32_is_watertight() {
        use crate::application::csg::CsgNode;
        use crate::domain::core::scalar::Real;
        use crate::domain::geometry::primitives::{Cylinder, Elbow, PrimitiveMesh};
        use nalgebra::{Isometry3, Translation3, UnitQuaternion, Vector3};

        let r: Real = 0.5;
        let h: Real = 3.0;
        let r_bend: Real = 1.0;
        let theta: Real = std::f64::consts::PI / 6.0;
        let eps: Real = r * 0.10;

        let (s_th, c_th) = theta.sin_cos();
        let axial_reach = r_bend * s_th;
        let radial_reach = r_bend * (1.0 - c_th);
        let arm_len = h - radial_reach;

        // Right elbow: -90° about X, translated to inlet at y = -axial_reach
        let elbow_inlet_y = -h + (h - axial_reach); // = -axial_reach
        let right_elbow_raw = Elbow {
            tube_radius: r,
            bend_radius: r_bend,
            bend_angle: theta,
            tube_segments: 64,
            arc_segments: 32,
        }
        .build()
        .expect("elbow build");
        let rot_base = UnitQuaternion::<Real>::from_axis_angle(
            &Vector3::x_axis(),
            -std::f64::consts::FRAC_PI_2,
        );
        let right_elbow = CsgNode::Transform {
            node: Box::new(CsgNode::Leaf(right_elbow_raw)),
            iso: Isometry3::from_parts(Translation3::new(0.0, elbow_inlet_y, 0.0), rot_base),
        }
        .evaluate()
        .expect("elbow transform");

        // Right arm: 64-segment cylinder, rotated -theta about Z, placed at elbow outlet.
        let arm_raw = Cylinder {
            base_center: Point3r::new(0.0, 0.0, 0.0),
            radius: r,
            height: arm_len + eps,
            segments: 64,
        }
        .build()
        .expect("arm build");
        let rot_arm = UnitQuaternion::<Real>::from_axis_angle(&Vector3::z_axis(), -theta);
        let tx = radial_reach - eps * s_th;
        let ty = -eps * c_th;
        let right_arm = CsgNode::Transform {
            node: Box::new(CsgNode::Leaf(arm_raw)),
            iso: Isometry3::from_parts(Translation3::new(tx, ty, 0.0), rot_arm),
        }
        .evaluate()
        .expect("arm transform");

        let mut result = csg_boolean_indexed(BooleanOp::Union, &right_elbow, &right_arm)
            .expect("right_elbow ∪ right_arm should not fail");
        let report = watertight_report(&mut result);

        if !report.is_watertight {
            result.rebuild_edges();
            let edges = result.edges_ref().unwrap();
            let mut bvec: Vec<(Point3r, Point3r)> = Vec::new();
            let mut nm_vids: HashSet<crate::domain::core::index::VertexId> = HashSet::new();
            // Collect boundary vertex positions
            let mut boundary_vids: HashSet<crate::domain::core::index::VertexId> = HashSet::new();
            for edge in edges.iter() {
                if edge.valence() == 1 {
                    bvec.push((
                        *result.vertices.position(edge.vertices.0),
                        *result.vertices.position(edge.vertices.1),
                    ));
                    boundary_vids.insert(edge.vertices.0);
                    boundary_vids.insert(edge.vertices.1);
                } else if edge.valence() > 2 {
                    nm_vids.insert(edge.vertices.0);
                    nm_vids.insert(edge.vertices.1);
                    let pa = result.vertices.position(edge.vertices.0);
                    let pb = result.vertices.position(edge.vertices.1);
                    eprintln!(
                        "NM_EDGE valence={}: ({:.6},{:.6},{:.6}) -> ({:.6},{:.6},{:.6})",
                        edge.valence(),
                        pa.x,
                        pa.y,
                        pa.z,
                        pb.x,
                        pb.y,
                        pb.z
                    );
                }
            }
            // Print faces containing non-manifold vertices
            for (fi, face) in result.faces.iter_enumerated() {
                let v = face.vertices;
                if v.iter().any(|vid| nm_vids.contains(vid)) {
                    let p0 = result.vertices.position(v[0]);
                    let p1 = result.vertices.position(v[1]);
                    let p2 = result.vertices.position(v[2]);
                    eprintln!(
                        "NM_FACE[{}]: ({:.6},{:.6},{:.6})/({:.6},{:.6},{:.6})/({:.6},{:.6},{:.6})",
                        fi, p0.x, p0.y, p0.z, p1.x, p1.y, p1.z, p2.x, p2.y, p2.z
                    );
                }
            }
            eprintln!(
                "=== V-shape right_branch 64x32 boundary edges ({}) ===",
                bvec.len()
            );
            for (a, b) in bvec.iter() {
                eprintln!(
                    "  ({:.9},{:.9},{:.9}) -> ({:.9},{:.9},{:.9})",
                    a.x, a.y, a.z, b.x, b.y, b.z
                );
            }
        }

        assert!(
            report.is_watertight,
            "V-shape right_elbow ∪ right_arm (64×32) should be watertight \
             (boundary_edges={}, non_manifold={})",
            report.boundary_edge_count, report.non_manifold_edge_count
        );
    }
