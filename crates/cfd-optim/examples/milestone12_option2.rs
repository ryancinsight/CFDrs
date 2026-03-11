//! Milestone 12 — Option 2: Selective Venturi Hydrodynamic Cavitation.
//!
//! This track centres on **Rayleigh–Plesset bubble dynamics** in venturi
//! throat constrictions. The Bernoulli cavitation number σ = (p∞ − pᵥ)/(½ρv²)
//! drops below unity when the local throat velocity is sufficient to lower
//! static pressure below the vapour pressure of blood (pᵥ ≈ 6.3 kPa at 37 °C),
//! nucleating vapour/gas microbubbles. Their violent inertial collapse at the
//! divergent exit delivers localised mechanical forces for sonodynamic therapy.
//!
//! Channel sizes are bounded by two physics constraints:
//! - **Lower bound**: throat width ≥ 35 µm to avoid clogging by RBC rouleaux
//!   (∅ ≈ 8 µm) and to remain above the minimum Reynolds number for incipient
//!   cavitation at achievable pump pressures.
//! - **Upper bound**: throat width ≤ 120 µm to maintain sufficient v_throat for
//!   σ < 1 within the ≤ 250 kPa gauge budget of clinical extracorporeal pumps.
//!
//! Asymmetric trifurcation stages upstream of the venturi use Zweifach–Fung
//! flow partitioning to route CTCs preferentially through the treatment lane
//! while diverting WBCs and RBCs to bypass branches, achieving selective
//! cavitation exposure. Scoring couples cavitation intensity I_cav with
//! cancer-center enrichment and RBC shielding in a hybrid additive + synergy
//! formula clamped to \[0.001, 1.0\].

use cfd_optim::{
    refresh_milestone12_reports, run_milestone12_option2, Milestone12RequestedStage,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    run_milestone12_option2()?;
    let _ = refresh_milestone12_reports(&[Milestone12RequestedStage::Option2])?;
    Ok(())
}
