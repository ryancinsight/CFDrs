//! Simulation runner — executes deterministic reference internal-flow analyses.
//!
//! # Theorem
//! For any mesh accepted by [`crate::application::simulation::setup::SimulationSetup::from_reference_mesh`],
//! [`run_reference_internal_flow_analysis`] returns the exact steady
//! Hagen-Poiseuille pressure drop and wall shear stress for the prescribed
//! Reynolds-number reference case.
//!
//! **Proof sketch**: the setup constructor restricts the mesh to straight
//! constant-radius circular channels, so the classical fully developed laminar
//! pipe solution applies exactly. The runner uses the closed-form relations
//! `Δp = 32 μ U L / D²` and `τ_w = 8 μ U / D`, then projects the linear
//! pressure field onto the surface mesh by axial coordinate.

use crate::application::simulation::results::{FieldData, SimulationResults};
use crate::application::simulation::setup::SimulationSetup;
use cfd_mesh::IndexedMesh;
use nalgebra::Point3;
use std::collections::HashMap;
use std::sync::atomic::{AtomicBool, AtomicU32, Ordering};
use std::sync::Arc;

/// State of a simulation run.
#[derive(Clone, Debug)]
pub enum SimulationState {
    /// No simulation running.
    Idle,
    /// Simulation is running.
    Running,
    /// Simulation completed successfully.
    Completed,
    /// Simulation failed with an error message.
    Failed(String),
}

/// Tracks progress and cancellation of a background simulation.
pub struct SimulationRunner {
    /// Current state.
    pub state: SimulationState,
    /// Progress percentage (0..100).
    pub progress: Arc<AtomicU32>,
    /// Cancellation flag.
    pub cancel: Arc<AtomicBool>,
}

/// Summary scalars from a completed reference internal-flow analysis.
#[derive(Clone, Debug)]
pub struct FlowSummary {
    /// Reynolds number of the deterministic reference case.
    pub reynolds_number: f64,
    /// Mean axial velocity [m/s].
    pub mean_velocity_m_s: f64,
    /// Pressure drop [Pa].
    pub pressure_drop_pa: f64,
    /// Wall shear stress [Pa].
    pub wall_shear_pa: f64,
}

impl SimulationRunner {
    /// Create a new idle runner.
    #[must_use]
    pub fn new() -> Self {
        Self {
            state: SimulationState::Idle,
            progress: Arc::new(AtomicU32::new(0)),
            cancel: Arc::new(AtomicBool::new(false)),
        }
    }

    /// Get the current progress percentage.
    #[must_use]
    pub fn progress_percent(&self) -> u32 {
        self.progress.load(Ordering::Relaxed)
    }

    /// Request cancellation of the running simulation.
    pub fn request_cancel(&self) {
        self.cancel.store(true, Ordering::Relaxed);
    }

    /// Check if cancellation was requested.
    #[must_use]
    pub fn is_cancel_requested(&self) -> bool {
        self.cancel.load(Ordering::Relaxed)
    }

    /// Reset the runner to idle state.
    pub fn reset(&mut self) {
        self.state = SimulationState::Idle;
        self.progress.store(0, Ordering::Relaxed);
        self.cancel.store(false, Ordering::Relaxed);
    }
}

impl Default for SimulationRunner {
    fn default() -> Self {
        Self::new()
    }
}

/// Execute the canonical laminar reference-flow analysis on a supported mesh.
///
/// Stores per-face `pressure_pa` and `wall_shear_pa` channels on `mesh`, and
/// returns vertex-averaged fields for UI inspection.
///
/// Validated by:
/// - `reference_analysis_populates_pressure_and_wall_shear_fields`
pub fn run_reference_internal_flow_analysis(
    mesh: &mut IndexedMesh<f64>,
) -> anyhow::Result<(SimulationSetup, SimulationResults, FlowSummary)> {
    let setup = SimulationSetup::from_reference_mesh(mesh)?;
    let axis_start = point_from_array(setup.inlet_centroid_m);
    let axis_end = point_from_array(setup.outlet_centroid_m);
    let axis = axis_end - axis_start;
    let axis_hat = axis / setup.channel_length_m;

    let mean_velocity_m_s = setup.inlet_velocity_m_s;
    let pressure_drop_pa = 32.0 * setup.viscosity_pa_s * mean_velocity_m_s * setup.channel_length_m
        / setup.hydraulic_diameter_m.powi(2);
    let wall_shear_pa = 8.0 * setup.viscosity_pa_s * mean_velocity_m_s / setup.hydraulic_diameter_m;
    let reynolds_number =
        setup.density_kg_m3 * mean_velocity_m_s * setup.hydraulic_diameter_m / setup.viscosity_pa_s;

    let wall_faces: std::collections::HashSet<_> = setup
        .boundaries
        .iter()
        .find(|assignment| assignment.label == "wall")
        .map(|assignment| assignment.faces.iter().copied().collect())
        .unwrap_or_default();

    let mut face_pressures = Vec::with_capacity(mesh.face_count());
    let mut face_wall_shear = Vec::with_capacity(mesh.face_count());

    mesh.attributes.remove_channel("pressure_pa");
    mesh.attributes.remove_channel("wall_shear_pa");

    for (face_id, face) in mesh.faces.iter_enumerated() {
        let a = mesh.vertices.position(face.vertices[0]);
        let b = mesh.vertices.position(face.vertices[1]);
        let c = mesh.vertices.position(face.vertices[2]);
        let centroid = Point3::from((a.coords + b.coords + c.coords) / 3.0);
        let offset = centroid - axis_start;
        let axial_position = (offset.dot(&axis_hat) / setup.channel_length_m).clamp(0.0, 1.0);
        let pressure_pa = setup.outlet_pressure_pa + pressure_drop_pa * (1.0 - axial_position);
        let shear_pa = if wall_faces.contains(&face_id) {
            wall_shear_pa
        } else {
            0.0
        };

        mesh.attributes.set("pressure_pa", face_id, pressure_pa);
        mesh.attributes.set("wall_shear_pa", face_id, shear_pa);
        face_pressures.push((face_id, pressure_pa));
        face_wall_shear.push((face_id, shear_pa));
    }

    let pressure_values = average_face_channel_onto_vertices(mesh, &face_pressures);
    let wall_shear_values = average_face_channel_onto_vertices(mesh, &face_wall_shear);
    let results = SimulationResults {
        fields: vec![
            FieldData::new("Pressure".to_owned(), "Pa".to_owned(), pressure_values),
            FieldData::new(
                "Wall Shear Stress".to_owned(),
                "Pa".to_owned(),
                wall_shear_values,
            ),
        ],
        iterations: 1,
        final_residual: 0.0,
        converged: true,
    };

    let summary = FlowSummary {
        reynolds_number,
        mean_velocity_m_s,
        pressure_drop_pa,
        wall_shear_pa,
    };

    Ok((setup, results, summary))
}

fn average_face_channel_onto_vertices(
    mesh: &IndexedMesh<f64>,
    face_values: &[(cfd_mesh::domain::core::index::FaceId, f64)],
) -> Vec<f64> {
    let mut sums: HashMap<u32, (f64, usize)> = HashMap::new();
    for &(face_id, value) in face_values {
        let face = mesh.faces.get(face_id);
        for vertex_id in face.vertices {
            let entry = sums.entry(vertex_id.0).or_insert((0.0, 0));
            entry.0 += value;
            entry.1 += 1;
        }
    }

    let mut values = Vec::with_capacity(mesh.vertex_count());
    for (vertex_id, _) in mesh.vertices.iter() {
        let (sum, count) = sums.get(&vertex_id.0).copied().unwrap_or((0.0, 0));
        values.push(if count == 0 { 0.0 } else { sum / count as f64 });
    }
    values
}

fn point_from_array(coords: [f64; 3]) -> Point3<f64> {
    Point3::new(coords[0], coords[1], coords[2])
}

#[cfg(test)]
mod tests {
    use super::*;
    use cfd_mesh::domain::geometry::primitives::{Cylinder, PrimitiveMesh};

    fn labeled_cylinder(radius: f64, height: f64) -> IndexedMesh<f64> {
        let mut mesh = Cylinder {
            radius,
            height,
            segments: 64,
            ..Cylinder::default()
        }
        .build()
        .expect("cylinder build must succeed");

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

        mesh
    }

    #[test]
    fn reference_analysis_populates_pressure_and_wall_shear_fields() {
        let mut mesh = labeled_cylinder(5.0e-4, 5.0e-3);

        let (_setup, results, summary) = run_reference_internal_flow_analysis(&mut mesh)
            .expect("labeled cylinder should admit a reference solve");

        assert!(summary.reynolds_number > 0.0);
        assert!(summary.pressure_drop_pa > 0.0);
        assert_eq!(results.fields.len(), 2);
        assert!(mesh.attributes.has_channel("pressure_pa"));
        assert!(mesh.attributes.has_channel("wall_shear_pa"));
        assert!(results.fields[0].max > results.fields[0].min);
        assert!(results.fields[1].max >= results.fields[1].min);
    }
}
