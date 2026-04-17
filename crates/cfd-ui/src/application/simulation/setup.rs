//! Simulation setup — boundary conditions and supported reference-flow setup.
//!
//! # Theorem
//! The setup recovered by [`SimulationSetup::from_reference_mesh`] is
//! admissible for an exact Hagen-Poiseuille reference solve if and only if the
//! selected mesh exposes a single inlet patch, a single outlet patch, and a
//! straight circular wall whose radius is constant along the axial direction.
//!
//! **Proof sketch**: a fully developed laminar pipe solution is exact only for
//! straight constant-area circular conduits with incompressible Newtonian flow.
//! The constructor enforces those geometric hypotheses directly from the mesh:
//! boundary labels define inlet/outlet/wall patches, patch circularity rejects
//! non-circular sections, inlet/outlet diameter agreement rejects area changes,
//! and wall-radius spread rejects bent or constricted channels. Any mesh that
//! violates these requirements is rejected before solver execution.

use anyhow::ensure;
use cfd_mesh::domain::core::index::FaceId;
use cfd_mesh::domain::geometry::measure::triangle_area;
use cfd_mesh::IndexedMesh;
use nalgebra::{Point3, Vector3};
use std::collections::HashMap;

/// Canonical Reynolds number used for the UI's deterministic laminar reference
/// solve.
pub const REFERENCE_REYNOLDS_NUMBER: f64 = 100.0;

const MIN_CIRCULARITY: f64 = 0.97;
const MAX_DIAMETER_RELATIVE_MISMATCH: f64 = 0.02;
const MAX_RADIUS_SPREAD_FRACTION: f64 = 0.12;

/// A boundary condition assignment on a set of mesh faces.
#[derive(Clone, Debug)]
pub struct BoundaryAssignment {
    /// The faces this boundary applies to.
    pub faces: Vec<FaceId>,
    /// The label for this boundary (e.g. "inlet", "outlet", "wall").
    pub label: String,
    /// The boundary type.
    pub condition: BoundaryType,
}

/// Supported boundary condition types for simulation.
#[derive(Clone, Debug)]
pub enum BoundaryType {
    /// Velocity inlet (magnitude in m/s).
    VelocityInlet { velocity_m_s: f64 },
    /// Pressure outlet (gauge pressure in Pa).
    PressureOutlet { pressure_pa: f64 },
    /// No-slip wall.
    NoSlipWall,
    /// Slip (free-slip) wall.
    SlipWall,
    /// Symmetry plane.
    Symmetry,
}

/// Complete simulation configuration ready for execution.
#[derive(Clone, Debug)]
pub struct SimulationSetup {
    /// The mesh to simulate on.
    pub mesh_vertex_count: usize,
    /// The mesh face count.
    pub mesh_face_count: usize,
    /// Boundary condition assignments.
    pub boundaries: Vec<BoundaryAssignment>,
    /// Fluid density (kg/m^3).
    pub density_kg_m3: f64,
    /// Fluid dynamic viscosity (Pa.s).
    pub viscosity_pa_s: f64,
    /// Maximum number of solver iterations.
    pub max_iterations: usize,
    /// Convergence tolerance.
    pub convergence_tolerance: f64,
    /// Axial length between inlet and outlet centroids [m].
    pub channel_length_m: f64,
    /// Circular hydraulic diameter [m].
    pub hydraulic_diameter_m: f64,
    /// Area-weighted inlet centroid [m].
    pub inlet_centroid_m: [f64; 3],
    /// Area-weighted outlet centroid [m].
    pub outlet_centroid_m: [f64; 3],
    /// Deterministic inlet velocity used by the reference solve [m/s].
    pub inlet_velocity_m_s: f64,
    /// Deterministic outlet pressure [Pa].
    pub outlet_pressure_pa: f64,
}

#[derive(Clone, Debug)]
struct PatchMetrics {
    face_ids: Vec<FaceId>,
    centroid_m: Point3<f64>,
    circularity: f64,
    equivalent_diameter_m: f64,
}

impl SimulationSetup {
    /// Construct the canonical reference-flow setup from a labeled mesh.
    ///
    /// Failure modes:
    /// - missing or unlabeled inlet/outlet boundary patches
    /// - non-circular inlet or outlet caps
    /// - varying section radius (venturi/constriction) or bent wall geometry
    /// - degenerate geometry or invalid physical parameters
    ///
    /// Validated by:
    /// - `reference_setup_accepts_labeled_cylinder`
    /// - `reference_setup_rejects_variable_radius_channel`
    pub fn from_reference_mesh(mesh: &IndexedMesh<f64>) -> anyhow::Result<Self> {
        let inlet = measure_patch(mesh, "inlet")?;
        let outlet = measure_patch(mesh, "outlet")?;
        let wall_faces = collect_labeled_faces(mesh, "wall");
        ensure!(!wall_faces.is_empty(), "reference flow solve requires labeled wall faces");

        ensure!(
            inlet.circularity >= MIN_CIRCULARITY,
            "reference flow solve requires a circular inlet patch, circularity={:.4}",
            inlet.circularity
        );
        ensure!(
            outlet.circularity >= MIN_CIRCULARITY,
            "reference flow solve requires a circular outlet patch, circularity={:.4}",
            outlet.circularity
        );

        let hydraulic_diameter_m =
            0.5 * (inlet.equivalent_diameter_m + outlet.equivalent_diameter_m);
        ensure!(
            hydraulic_diameter_m.is_finite() && hydraulic_diameter_m > 0.0,
            "reference flow solve requires a positive hydraulic diameter"
        );

        let relative_diameter_mismatch =
            (inlet.equivalent_diameter_m - outlet.equivalent_diameter_m).abs()
                / hydraulic_diameter_m.max(f64::EPSILON);
        ensure!(
            relative_diameter_mismatch <= MAX_DIAMETER_RELATIVE_MISMATCH,
            "reference flow solve requires a constant-area channel; inlet/outlet diameter mismatch={relative_diameter_mismatch:.4}"
        );

        let axis = outlet.centroid_m - inlet.centroid_m;
        let channel_length_m = axis.norm();
        ensure!(
            channel_length_m.is_finite() && channel_length_m > 0.0,
            "reference flow solve requires distinct inlet and outlet centroids"
        );
        let axis_hat = axis / channel_length_m;

        validate_wall_radius_consistency(mesh, &wall_faces, inlet.centroid_m, axis_hat, hydraulic_diameter_m)?;

        let density_kg_m3 = 1000.0;
        let viscosity_pa_s = 1.0e-3;
        let inlet_velocity_m_s =
            REFERENCE_REYNOLDS_NUMBER * viscosity_pa_s / (density_kg_m3 * hydraulic_diameter_m);
        ensure!(
            inlet_velocity_m_s.is_finite() && inlet_velocity_m_s > 0.0,
            "reference inlet velocity must be positive"
        );

        Ok(Self {
            mesh_vertex_count: mesh.vertex_count(),
            mesh_face_count: mesh.face_count(),
            boundaries: vec![
                BoundaryAssignment {
                    faces: inlet.face_ids,
                    label: "inlet".to_owned(),
                    condition: BoundaryType::VelocityInlet {
                        velocity_m_s: inlet_velocity_m_s,
                    },
                },
                BoundaryAssignment {
                    faces: outlet.face_ids,
                    label: "outlet".to_owned(),
                    condition: BoundaryType::PressureOutlet { pressure_pa: 0.0 },
                },
                BoundaryAssignment {
                    faces: wall_faces,
                    label: "wall".to_owned(),
                    condition: BoundaryType::NoSlipWall,
                },
            ],
            density_kg_m3,
            viscosity_pa_s,
            max_iterations: 1,
            convergence_tolerance: 0.0,
            channel_length_m,
            hydraulic_diameter_m,
            inlet_centroid_m: point_to_array(inlet.centroid_m),
            outlet_centroid_m: point_to_array(outlet.centroid_m),
            inlet_velocity_m_s,
            outlet_pressure_pa: 0.0,
        })
    }
}

impl Default for SimulationSetup {
    fn default() -> Self {
        Self {
            mesh_vertex_count: 0,
            mesh_face_count: 0,
            boundaries: Vec::new(),
            density_kg_m3: 1000.0,
            viscosity_pa_s: 1e-3,
            max_iterations: 1000,
            convergence_tolerance: 1e-6,
            channel_length_m: 0.0,
            hydraulic_diameter_m: 0.0,
            inlet_centroid_m: [0.0; 3],
            outlet_centroid_m: [0.0; 3],
            inlet_velocity_m_s: 0.0,
            outlet_pressure_pa: 0.0,
        }
    }
}

fn collect_labeled_faces(mesh: &IndexedMesh<f64>, label: &str) -> Vec<FaceId> {
    mesh.faces
        .iter_enumerated()
        .filter_map(|(face_id, _)| (mesh.boundary_label(face_id) == Some(label)).then_some(face_id))
        .collect()
}

fn measure_patch(mesh: &IndexedMesh<f64>, label: &str) -> anyhow::Result<PatchMetrics> {
    let face_ids = collect_labeled_faces(mesh, label);
    ensure!(
        !face_ids.is_empty(),
        "reference flow solve requires at least one '{label}' boundary face"
    );

    let mut area_m2 = 0.0;
    let mut weighted_centroid = Vector3::zeros();
    let mut edge_counts: HashMap<(u32, u32), usize> = HashMap::new();

    for &face_id in &face_ids {
        let face = mesh.faces.get(face_id);
        let a = mesh.vertices.position(face.vertices[0]);
        let b = mesh.vertices.position(face.vertices[1]);
        let c = mesh.vertices.position(face.vertices[2]);
        let area = triangle_area(a, b, c);
        if !area.is_finite() || area <= 1.0e-18 {
            continue;
        }
        let centroid = (a.coords + b.coords + c.coords) / 3.0;
        area_m2 += area;
        weighted_centroid += centroid * area;

        for (u, v) in face.edges_canonical() {
            let key = (u.0, v.0);
            *edge_counts.entry(key).or_insert(0) += 1;
        }
    }

    ensure!(area_m2.is_finite() && area_m2 > 0.0, "{label} patch must have positive area");
    let centroid_m = Point3::from(weighted_centroid / area_m2);

    let perimeter_m = edge_counts
        .into_iter()
        .filter(|(_, count)| *count == 1)
        .map(|((u, v), _)| {
            let p0 = mesh.vertices.position(cfd_mesh::domain::core::index::VertexId(u));
            let p1 = mesh.vertices.position(cfd_mesh::domain::core::index::VertexId(v));
            (p1 - p0).norm()
        })
        .sum::<f64>();
    ensure!(
        perimeter_m.is_finite() && perimeter_m > 0.0,
        "{label} patch perimeter must be positive"
    );

    let circularity = 4.0 * std::f64::consts::PI * area_m2 / perimeter_m.powi(2);
    let equivalent_diameter_m = 2.0 * (area_m2 / std::f64::consts::PI).sqrt();

    Ok(PatchMetrics {
        face_ids,
        centroid_m,
        circularity,
        equivalent_diameter_m,
    })
}

fn validate_wall_radius_consistency(
    mesh: &IndexedMesh<f64>,
    wall_faces: &[FaceId],
    inlet_centroid_m: Point3<f64>,
    axis_hat: Vector3<f64>,
    hydraulic_diameter_m: f64,
) -> anyhow::Result<()> {
    let reference_radius_m = hydraulic_diameter_m / 2.0;
    let mut min_radius_m = f64::INFINITY;
    let mut max_radius_m = 0.0_f64;

    for &face_id in wall_faces {
        let face = mesh.faces.get(face_id);
        let a = mesh.vertices.position(face.vertices[0]);
        let b = mesh.vertices.position(face.vertices[1]);
        let c = mesh.vertices.position(face.vertices[2]);
        let centroid = Point3::from((a.coords + b.coords + c.coords) / 3.0);
        let offset = centroid - inlet_centroid_m;
        let axial_distance = offset.dot(&axis_hat);
        let radial_vector = offset - axis_hat * axial_distance;
        let radial_distance = radial_vector.norm();
        min_radius_m = min_radius_m.min(radial_distance);
        max_radius_m = max_radius_m.max(radial_distance);
    }

    ensure!(
        min_radius_m.is_finite() && max_radius_m.is_finite(),
        "wall radius statistics must be finite"
    );
    let relative_radius_spread =
        (max_radius_m - min_radius_m).abs() / reference_radius_m.max(f64::EPSILON);
    ensure!(
        relative_radius_spread <= MAX_RADIUS_SPREAD_FRACTION,
        "reference flow solve requires a straight constant-radius wall; relative radius spread={relative_radius_spread:.4}"
    );

    Ok(())
}

fn point_to_array(point: Point3<f64>) -> [f64; 3] {
    [point.x, point.y, point.z]
}

#[cfg(test)]
mod tests {
    use super::*;
    use cfd_mesh::domain::geometry::primitives::{Cylinder, PrimitiveMesh};

    fn label_cylinder_boundaries(mesh: &mut IndexedMesh<f64>, height: f64) {
        let labels: Vec<_> = mesh
            .faces
            .iter_enumerated()
            .map(|(face_id, face)| {
                let a = mesh.vertices.position(face.vertices[0]);
                let b = mesh.vertices.position(face.vertices[1]);
                let c = mesh.vertices.position(face.vertices[2]);
                let centroid_y = (a.y + b.y + c.y) / 3.0;
                let label = if centroid_y < 1.0e-9 {
                    "inlet"
                } else if (centroid_y - height).abs() < 1.0e-9 {
                    "outlet"
                } else {
                    "wall"
                };
                (face_id, label)
            })
            .collect();

        for (face_id, label) in labels {
            mesh.mark_boundary(face_id, label);
        }
    }

    #[test]
    fn reference_setup_accepts_labeled_cylinder() {
        let height = 5.0e-3;
        let mut mesh = Cylinder {
            radius: 5.0e-4,
            height,
            segments: 64,
            ..Cylinder::default()
        }
        .build()
        .expect("cylinder build must succeed");
        label_cylinder_boundaries(&mut mesh, height);

        let setup = SimulationSetup::from_reference_mesh(&mesh)
            .expect("labeled cylinder should be admissible");

        assert_eq!(setup.boundaries.len(), 3);
        assert!(setup.channel_length_m > 0.0);
        assert!(setup.hydraulic_diameter_m > 0.0);
        assert!(setup.inlet_velocity_m_s > 0.0);
    }

    #[test]
    fn reference_setup_rejects_variable_radius_channel() {
        let mesh = crate::application::mesh_ops::channel::build_test_venturi_mesh();
        let error = SimulationSetup::from_reference_mesh(&mesh)
            .expect_err("venturi mesh must be rejected by constant-radius filter");

        assert!(error
            .to_string()
            .contains("constant-radius wall"));
    }
}
