//! Milestone 12 requirement extraction from contract-facing source documents.
//!
//! Policy: Schedule is authoritative for milestone wording. TDD is used as a
//! cross-check for milestone intent and limits-of-usage traceability.

use std::path::Path;

const SCHEDULE_PATH: &str =
    "report/ATTACHMENT_3_Schedule_of_Milestones_and_Payments_Template_10.03.24.md";
const TDD_PATH: &str = "report/Task_Description_Document_9.27.2024_V5.md";

const M12_DESCRIPTION: &str =
    "Select optimal design(s) that meet desired hydrodynamic and cavitation parameters.";
const M12_EXIT: &str = "Design selected that meets hydrodynamic and cavitation parameters.";
const M12_DELIVERABLE: &str = "Final report detailing the process for selecting optimal designs.";
const TDD_LIMITS_USAGE: &str = "Provide thresholds or limits of usage.";

/// Contract wording used by the narrative report.
#[derive(Debug, Clone)]
pub struct MilestoneContractText {
    pub description: String,
    pub exit_criteria: String,
    pub deliverable: String,
}

/// Load and validate Milestone 12 contract text.
///
/// # Errors
/// Returns an error when required schedule or TDD trace phrases are missing.
pub fn load_m12_contract_text(
    workspace_root: &Path,
) -> Result<MilestoneContractText, Box<dyn std::error::Error>> {
    let schedule = std::fs::read_to_string(workspace_root.join(SCHEDULE_PATH))?;
    ensure_contains(
        &schedule,
        M12_DESCRIPTION,
        "schedule milestone 12 description",
    )?;
    ensure_contains(&schedule, M12_EXIT, "schedule milestone 12 exit criteria")?;
    ensure_contains(
        &schedule,
        M12_DELIVERABLE,
        "schedule milestone 12 deliverable",
    )?;

    let tdd = std::fs::read_to_string(workspace_root.join(TDD_PATH))?;
    ensure_contains(
        &tdd,
        M12_DESCRIPTION,
        "TDD milestone 12 description cross-check",
    )?;
    ensure_contains(&tdd, TDD_LIMITS_USAGE, "TDD limits-of-usage cross-check")?;

    Ok(MilestoneContractText {
        description: M12_DESCRIPTION.to_string(),
        exit_criteria: M12_EXIT.to_string(),
        deliverable: M12_DELIVERABLE.to_string(),
    })
}

fn ensure_contains(
    text: &str,
    needle: &str,
    context: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    if text.contains(needle) {
        Ok(())
    } else {
        Err(format!("missing required phrase for {context}: \"{needle}\"").into())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    #[test]
    fn load_contract_text_from_fixture() {
        let unique_suffix = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system clock must be after the Unix epoch")
            .as_nanos();
        let workspace_root = std::env::temp_dir().join(format!(
            "cfdrs-contract-{}-{unique_suffix}",
            std::process::id()
        ));
        let report = workspace_root.join("report");
        std::fs::create_dir_all(&report).expect("fixture report directory must be created");

        let schedule = format!("{M12_DESCRIPTION}\n{M12_EXIT}\n{M12_DELIVERABLE}\n");
        std::fs::write(
            report.join(SCHEDULE_PATH.strip_prefix("report/").unwrap()),
            schedule,
        )
        .expect("schedule fixture must be written");
        std::fs::write(
            report.join(TDD_PATH.strip_prefix("report/").unwrap()),
            format!("{M12_DESCRIPTION}\n{TDD_LIMITS_USAGE}\n"),
        )
        .expect("TDD fixture must be written");

        let result = load_m12_contract_text(&workspace_root);
        let cleanup_result = std::fs::remove_dir_all(&workspace_root);
        assert!(cleanup_result.is_ok(), "fixture cleanup must succeed");
        let txt = result.expect("contract text must load from the fixture");
        assert_eq!(txt.description, M12_DESCRIPTION);
        assert_eq!(txt.exit_criteria, M12_EXIT);
        assert_eq!(txt.deliverable, M12_DELIVERABLE);
    }
}
