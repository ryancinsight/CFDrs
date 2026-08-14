#![expect(missing_docs, reason = "ratchet CFDRS-DOCS-1: pre-existing debt")]

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum OptimizationGoal {
    AsymmetricSplitResidenceSeparation,
    AsymmetricSplitVenturiCavitationSelectivity,
    InPlaceDeanSerpentineRefinement,
}
