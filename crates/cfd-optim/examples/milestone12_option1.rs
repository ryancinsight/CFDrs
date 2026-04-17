//! Milestone 12 — Option 1: Selective Acoustic Residence/Separation.
//!
//! This track exploits the **Zweifach–Fung effect** at asymmetric bifurcations:
//! when a parent channel (e.g. 6 mm) splits into branches of unequal width
//! (e.g. 1 mm / 2 mm / 3 mm), Fåhræus–Lindqvist plasma skimming and inertial
//! focusing steer larger cells (CTCs ∅ 15–25 µm) toward the wider center
//! branch while smaller WBCs and RBCs preferentially enter narrower peripheral
//! branches. Branch-width ratios are optimised to maximise cancer-center
//! enrichment (cancer_center_fraction) and treatment-lane residence time
//! (τ_res = L/v̄) while minimising healthy-cell treatment exposure.
//!
//! No venturi throats are active — treatment relies on externally applied
//! 412 kHz ultrasound. The acoustic resonance factor (ARF) quantifies how
//! closely each channel's hydraulic diameter D_h matches λ/2 ≈ 1.87 mm
//! (c_blood ≈ 1540 m/s), which creates standing-wave pressure antinodes for
//! preferential bubble trapping and sonosensitiser activation.
//!
//! Scoring uses a hybrid additive + geometric-mean synergy formula clamped
//! to \[0.001, 1.0\], ensuring no valid design collapses to zero.

use cfd_optim::{refresh_milestone12_reports, run_milestone12_option1, Milestone12RequestedStage};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    run_milestone12_option1()?;
    let _ = refresh_milestone12_reports(&[Milestone12RequestedStage::Option1])?;
    Ok(())
}
