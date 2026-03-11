//! Milestone 12 — GA: In-Place Dean–Serpentine Refinement.
//!
//! This track uses a blueprint-native genetic algorithm to refine Option 2
//! venturi designs by inserting **serpentine (curved) segments** on treatment-
//! path channels. The **Dean instability** (De = Re √(D_h / 2R_c), where R_c
//! is the bend radius of curvature) generates counter-rotating secondary
//! vortices that focus larger particles (CTCs) toward the outer wall at bend
//! apices — precisely where centrifugal forces peak.
//!
//! Three classes of architecture-preserving mutations are applied via
//! `BlueprintTopologyMutation`:
//! 1. **Branch width scaling** (treatment ×1.08, bypass ×0.94) — shifts the
//!    Zweifach–Fung flow partition ratio to increase cancer-center fraction.
//! 2. **Serpentine insertion** — adds curved channel segments with specified
//!    bend count and radius, introducing Dean secondary flow. The Dean number
//!    bonus (De_max / 100) in the GA score rewards designs that co-localise
//!    inertial focusing from Dean vortices with hydrodynamic cavitation at
//!    venturi throats positioned at bend apices.
//! 3. **Venturi throat narrowing** (×0.92) — lowers σ at the vena contracta
//!    to increase cavitation intensity per serial stage.
//!
//! GA score is **not comparable** to Option 1 or Option 2 scores — different
//! objective functions weight different physics.

use cfd_optim::{refresh_milestone12_reports, run_milestone12_ga, Milestone12RequestedStage};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    run_milestone12_ga()?;
    let _ = refresh_milestone12_reports(&[Milestone12RequestedStage::Ga])?;
    Ok(())
}
