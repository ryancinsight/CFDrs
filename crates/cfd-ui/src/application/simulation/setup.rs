//! Simulation setup — boundary conditions and solver configuration.

use cfd_mesh::domain::core::index::FaceId;

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
}

impl Default for SimulationSetup {
    fn default() -> Self {
        Self {
            mesh_vertex_count: 0,
            mesh_face_count: 0,
            boundaries: Vec::new(),
            density_kg_m3: 1000.0,    // water
            viscosity_pa_s: 1e-3,     // water at 20°C
            max_iterations: 1000,
            convergence_tolerance: 1e-6,
        }
    }
}
