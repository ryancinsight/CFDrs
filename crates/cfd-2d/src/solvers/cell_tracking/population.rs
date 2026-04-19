use serde::{Deserialize, Serialize};

/// Cell population type, matching the 1D model's three-population framework.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum CellPopulation {
    /// Circulating tumor cell (diameter 10-15 um, stiff).
    CTC,
    /// White blood cell (diameter 10-12 um, moderate stiffness).
    WBC,
    /// Red blood cell (diameter ~8 um, highly deformable).
    RBC,
}

impl CellPopulation {
    /// Characteristic diameter [m].
    #[must_use]
    pub fn diameter_m(self) -> f64 {
        match self {
            Self::CTC => 12.5e-6,
            Self::WBC => 11.0e-6,
            Self::RBC => 8.0e-6,
        }
    }

    /// Cell density [kg/m3].
    #[must_use]
    pub fn density_kg_m3(self) -> f64 {
        match self {
            Self::CTC => 1068.0,
            Self::WBC => 1070.0,
            Self::RBC => 1100.0,
        }
    }

    /// Di Carlo (2009) inertial lift coefficient.
    /// Stiffer cells experience stronger lift toward equilibrium positions.
    /// These values are for kappa = a/D_h ~ 0.05-0.15 (millifluidic range).
    #[must_use]
    pub fn lift_coefficient(self) -> f64 {
        match self {
            Self::CTC => 0.5,  // stiff, strong focusing
            Self::WBC => 0.35, // moderate deformability
            Self::RBC => 0.15, // high deformability, weaker focusing
        }
    }
}

/// A single tracked cell with position, velocity, and identity.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrackedCell {
    /// The specific cell type (CTC, WBC, RBC).
    pub population: CellPopulation,
    /// X coordinate position [m].
    pub x: f64,
    /// Y coordinate position [m].
    pub y: f64,
    /// Velocity in X direction [m/s].
    pub vx: f64,
    /// Velocity in Y direction [m/s].
    pub vy: f64,
    /// Unique identifier for the cell.
    pub id: usize,
}

/// Trajectory record.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CellTrajectory {
    /// The original cell identifier.
    pub cell_id: usize,
    /// The cell population category.
    pub population: CellPopulation,
    /// Tracked positions over time [x, y, t].
    pub positions: Vec<[f64; 3]>,
    /// Classification of the outlet region through which the cell exited.
    pub exit_outlet: Option<String>,
}

/// Outlet zone definition for classifying cell exit positions.
#[derive(Debug, Clone)]
pub struct OutletZone {
    /// Zone identifier (e.g. "center", "peripheral").
    pub name: String,
    /// x range: cell exits when x >= x_min.
    pub x_min: f64,
    /// Lower y bound for this outlet.
    pub y_lo: f64,
    /// Upper y bound for this outlet.
    pub y_hi: f64,
}

/// Summary of cell routing at outlets.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct CellRoutingSummary {
    /// Total number of CTCs recorded.
    pub ctc_total: usize,
    /// Total number of WBCs recorded.
    pub wbc_total: usize,
    /// Total number of RBCs recorded.
    pub rbc_total: usize,
    /// Total number of CTCs that exited through the center.
    pub ctc_center: usize,
    /// Total number of WBCs that exited through the center.
    pub wbc_center: usize,
    /// Total number of RBCs that exited through the center.
    pub rbc_center: usize,
    /// Fraction of CTCs routed toward the desired outlet.
    pub cancer_center_fraction: f64,
    /// Multiplicative metric representing total separation efficiency.
    pub separation_efficiency: f64,
}
