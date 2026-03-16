use std::path::PathBuf;

use serde::{Deserialize, Serialize};

use crate::analysis::RobustnessReport;
use crate::reporting::{GoalAuditArtifacts, Milestone12ReportDesign, ValidationRow};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Milestone12StageArtifact {
    pub label: String,
    pub path: PathBuf,
}

#[derive(Debug, Clone)]
pub struct Milestone12Option1Run {
    pub total_candidates: usize,
    pub eligible_count: usize,
    pub ranked: Vec<Milestone12ReportDesign>,
    pub audit: GoalAuditArtifacts,
    pub artifacts: Vec<Milestone12StageArtifact>,
}

#[derive(Debug, Clone)]
pub struct Milestone12Option2Run {
    pub total_candidates: usize,
    pub eligible_count: usize,
    pub ranked: Vec<Milestone12ReportDesign>,
    pub robustness: Vec<RobustnessReport>,
    pub audit: GoalAuditArtifacts,
    pub artifacts: Vec<Milestone12StageArtifact>,
}

#[derive(Debug, Clone)]
pub struct Milestone12GaRun {
    pub ranked: Vec<Milestone12ReportDesign>,
    pub best_per_generation: Vec<f64>,
    pub audit: GoalAuditArtifacts,
    pub artifacts: Vec<Milestone12StageArtifact>,
}

#[derive(Debug, Clone)]
pub struct Milestone12ValidationRun {
    pub rows: Vec<ValidationRow>,
    pub artifacts: Vec<Milestone12StageArtifact>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Milestone12RequestedStage {
    Option1,
    Option2,
    Ga,
    Validation,
}

impl Milestone12RequestedStage {
    /// Default report stages. Validation is excluded (run separately via `--stages validation`).
    pub const ALL: [Self; 3] = [Self::Option1, Self::Option2, Self::Ga];

    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::Option1 => "option1",
            Self::Option2 => "option2",
            Self::Ga => "ga",
            Self::Validation => "validation",
        }
    }
}
