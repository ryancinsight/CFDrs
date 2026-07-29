use aequitas::systems::si::quantities::{Dimensionless, Length, MassDensity, Time, Velocity};
use serde::{Deserialize, Serialize};

/// Cell population type, matching the 1D model's three-population framework.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum CellPopulation {
    /// Circulating tumor cell.
    CTC,
    /// White blood cell.
    WBC,
    /// Red blood cell.
    RBC,
}

impl CellPopulation {
    /// Return the characteristic cell diameter.
    #[must_use]
    pub fn diameter(self) -> Length {
        Length::from_base(match self {
            Self::CTC => 12.5e-6,
            Self::WBC => 11.0e-6,
            Self::RBC => 8.0e-6,
        })
    }

    /// Return the cell mass density.
    #[must_use]
    pub fn density(self) -> MassDensity {
        MassDensity::from_base(match self {
            Self::CTC => 1068.0,
            Self::WBC => 1070.0,
            Self::RBC => 1100.0,
        })
    }

    /// Return the dimensionless Di Carlo inertial lift coefficient.
    #[must_use]
    pub fn lift_coefficient(self) -> Dimensionless {
        Dimensionless::from_base(match self {
            Self::CTC => 0.5,
            Self::WBC => 0.35,
            Self::RBC => 0.15,
        })
    }
}

/// A single tracked cell with physical position, velocity, and identity.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrackedCell {
    /// The specific cell type.
    pub population: CellPopulation,
    /// X-coordinate position.
    pub x: Length,
    /// Y-coordinate position.
    pub y: Length,
    /// Velocity in the X direction.
    pub vx: Velocity,
    /// Velocity in the Y direction.
    pub vy: Velocity,
    /// Unique identifier for the cell.
    pub id: usize,
}

/// One physical sample in a cell trajectory.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrackedPosition {
    /// X-coordinate position.
    pub x: Length,
    /// Y-coordinate position.
    pub y: Length,
    /// Elapsed trajectory time.
    pub time: Time,
}

/// Trajectory record.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CellTrajectory {
    /// The original cell identifier.
    pub cell_id: usize,
    /// The cell population category.
    pub population: CellPopulation,
    /// Physical positions sampled over the trajectory.
    pub positions: Vec<TrackedPosition>,
    /// Classification of the outlet region through which the cell exited.
    pub exit_outlet: Option<String>,
}

/// Outlet zone definition for classifying cell exit positions.
#[derive(Debug, Clone)]
pub struct OutletZone {
    /// Zone identifier (for example, `center` or `peripheral`).
    pub name: String,
    /// Minimum X-coordinate at which a cell can exit.
    pub x_min: Length,
    /// Lower Y bound for this outlet.
    pub y_lo: Length,
    /// Upper Y bound for this outlet.
    pub y_hi: Length,
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
    pub cancer_center_fraction: Dimensionless,
    /// Dimensionless total separation efficiency.
    pub separation_efficiency: Dimensionless,
}
