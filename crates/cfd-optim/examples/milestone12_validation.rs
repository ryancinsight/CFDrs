//! Milestone 12 multi-fidelity venturi validation.

use cfd_optim::{
    refresh_milestone12_reports, run_milestone12_validation, Milestone12RequestedStage,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    run_milestone12_validation()?;
    let _ = refresh_milestone12_reports(&[Milestone12RequestedStage::Validation])?;
    Ok(())
}
