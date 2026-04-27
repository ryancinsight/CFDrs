//! Simulation results — stores computed field data from a CFD run.

/// Scalar field data computed by a simulation.
#[derive(Clone, Debug)]
pub struct FieldData {
    /// Human-readable name (e.g. "Velocity Magnitude", "Pressure").
    pub name: String,
    /// Units string (e.g. "m/s", "Pa").
    pub units: String,
    /// Per-vertex scalar values (indexed by vertex order in the mesh).
    pub values: Vec<f64>,
    /// Minimum value across all vertices.
    pub min: f64,
    /// Maximum value across all vertices.
    pub max: f64,
}

impl FieldData {
    /// Create a new field from per-vertex values.
    #[must_use]
    pub fn new(name: String, units: String, values: Vec<f64>) -> Self {
        let min = values.iter().copied().fold(f64::INFINITY, f64::min);
        let max = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        Self {
            name,
            units,
            values,
            min,
            max,
        }
    }
}

/// Container for all results from a simulation run.
#[derive(Clone, Debug)]
pub struct SimulationResults {
    /// Scalar fields computed by the solver.
    pub fields: Vec<FieldData>,
    /// Total number of iterations performed.
    pub iterations: usize,
    /// Final residual.
    pub final_residual: f64,
    /// Whether the simulation converged.
    pub converged: bool,
}
