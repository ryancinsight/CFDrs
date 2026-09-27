use crate::topology::model::{BranchRole, SplitKind};
use aequitas::systems::si::quantities::Length;

/// Branch specification within a Milestone 12 stage layout.
#[derive(Debug, Clone, PartialEq)]
pub struct Milestone12StageBranchSpec {
    /// Branch label used to derive channel identifiers.
    pub label: String,
    /// Branch role within the treatment network.
    pub role: BranchRole,
    /// Whether this branch carries the treatment path.
    pub treatment_path: bool,
    /// Authored branch width.
    pub width_m: Length<f64>,
}

/// Stage layout within a Milestone 12 topology.
#[derive(Debug, Clone, PartialEq)]
pub struct Milestone12StageLayout {
    /// Split kind of the stage.
    pub split_kind: SplitKind,
    /// Branch specifications for the stage.
    pub branches: Vec<Milestone12StageBranchSpec>,
}
