use crate::domain::OptimizationGoal;
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct FigureManifestEntry {
    pub figure_id: String,
    pub caption: String,
    pub path: PathBuf,
}

impl FigureManifestEntry {
    #[must_use]
    pub fn exists(&self) -> bool {
        self.path.exists()
    }
}

#[must_use]
pub fn required_figure_ids(goal: OptimizationGoal) -> &'static [&'static str] {
    match goal {
        OptimizationGoal::SelectiveAcousticResidenceSeparation => &[
            "option1_blueprint",
            "option1_residence_map",
            "option1_separation_curve",
        ],
        OptimizationGoal::SelectiveVenturiCavitation => &[
            "option2_blueprint",
            "option2_venturi_screening",
            "option2_cavitation_window",
        ],
        OptimizationGoal::BlueprintGeneticRefinement => {
            &["ga_blueprint", "ga_lineage_delta", "ga_dean_peak_map"]
        }
    }
}

pub fn missing_required_figures(
    goal: OptimizationGoal,
    figures: &[FigureManifestEntry],
) -> Vec<String> {
    required_figure_ids(goal)
        .iter()
        .filter(|required| {
            !figures
                .iter()
                .any(|figure| figure.figure_id == **required && Path::new(&figure.path).exists())
        })
        .map(|required| (*required).to_string())
        .collect()
}
