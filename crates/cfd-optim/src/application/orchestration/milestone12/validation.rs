use crate::application::orchestration::{
    ensure_release_reports, init_tracing, resolve_output_directories,
};
use crate::delivery::load_top5_report_json;
use crate::reporting::run_milestone12_validation as run_validation_solver;

use super::types::{Milestone12StageArtifact, Milestone12ValidationRun};

pub fn run_milestone12_validation() -> Result<Milestone12ValidationRun, Box<dyn std::error::Error>>
{
    init_tracing();
    ensure_release_reports()?;
    let (_, out_dir, _) = resolve_output_directories()?;
    let option2_ranked =
        load_top5_report_json(&out_dir.join("two_concept_option2_venturi_top5.json"))?;
    let ga_ranked = load_top5_report_json(&out_dir.join("ga_hydrosdt_top5.json"))?;

    let option2 = option2_ranked
        .first()
        .ok_or("two_concept_option2_venturi_top5.json is empty")?;
    let ga = ga_ranked.first().ok_or("ga_hydrosdt_top5.json is empty")?;
    let rows = run_validation_solver(&out_dir, option2, ga)?;

    Ok(Milestone12ValidationRun {
        rows,
        artifacts: vec![Milestone12StageArtifact {
            label: "validation_rows".to_string(),
            path: out_dir.join("milestone12_validation_rows.json"),
        }],
    })
}
