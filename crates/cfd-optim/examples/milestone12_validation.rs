#![allow(missing_docs)]
//! Milestone 12 multi-fidelity venturi validation.

use cfd_optim::{
    Milestone12RequestedStage, refresh_milestone12_reports, run_milestone12_validation,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    run_milestone12_validation()?;
    let _ = refresh_milestone12_reports(&[Milestone12RequestedStage::Validation])?;
    Ok(())
}
