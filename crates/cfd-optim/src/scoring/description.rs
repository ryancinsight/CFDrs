use super::OptimMode;

// ── Utility ──────────────────────────────────────────────────────────────────

/// Return a human-readable summary line for a scored candidate.
#[must_use]
pub fn score_description(mode: OptimMode) -> &'static str {
    match mode {
        OptimMode::SdtCavitation => "SDT Cavitation",
        OptimMode::UniformExposure => "Uniform Exposure",
        OptimMode::Combined { .. } => "Combined (Cavitation + Exposure)",
        OptimMode::CellSeparation => "Cell Separation + SDT",
        OptimMode::ThreePopSeparation => "Three-Pop Separation (WBC+Cancer→Center, RBC→Wall) + SDT",
        OptimMode::SdtTherapy => "SDT Therapy (Selective Sep + HI + Cav + Dose)",
        OptimMode::PediatricLeukapheresis { .. } => {
            "Paediatric Leukapheresis (WBC Recovery + RBC Removal + Purity + ECV)"
        }
        OptimMode::HydrodynamicCavitationSDT => {
            "Hydrodynamic Cavitation SDT (Cancer-Targeted Cav + Sep3 + RBC Protection + Sono)"
        }
        OptimMode::CombinedSdtLeukapheresis { .. } => {
            "Combined SDT + Leukapheresis (WBC Recovery + RBC Removal + Cancer-Targeted Cavitation)"
        }
        OptimMode::RbcProtectedSdt => {
            "RBC-Protected SDT (Therapeutic Window + Cancer-Cav + Lysis Safety + FDA Compliance)"
        }
    }
}
