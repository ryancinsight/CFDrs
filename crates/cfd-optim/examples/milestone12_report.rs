//! Milestone 12 aggregate runner over the shared stage orchestration API.

use cfd_optim::{run_milestone12_report, Milestone12RequestedStage};

#[derive(Default)]
struct RequestedStages {
    stages: Vec<Milestone12RequestedStage>,
}

impl RequestedStages {
    fn from_args() -> Result<Self, Box<dyn std::error::Error>> {
        let mut requested = Self::default();
        let mut args = std::env::args().skip(1);

        while let Some(arg) = args.next() {
            if arg == "--stages" {
                let csv = args
                    .next()
                    .ok_or("--stages requires a comma-separated value")?;
                requested.extend_csv(&csv)?;
            } else {
                return Err(format!("Unsupported argument: {arg}").into());
            }
        }

        if requested.stages.is_empty() {
            requested.stages.extend(Milestone12RequestedStage::ALL);
        }

        Ok(requested)
    }

    fn extend_csv(&mut self, csv: &str) -> Result<(), Box<dyn std::error::Error>> {
        for label in csv
            .split(',')
            .map(str::trim)
            .filter(|label| !label.is_empty())
        {
            let stage = match label.to_ascii_lowercase().as_str() {
                "all" => {
                    self.stages.clear();
                    self.stages.extend(Milestone12RequestedStage::ALL);
                    continue;
                }
                "option1" => Milestone12RequestedStage::Option1,
                "option2" => Milestone12RequestedStage::Option2,
                "ga" => Milestone12RequestedStage::Ga,
                "validation" => Milestone12RequestedStage::Validation,
                "refresh" => Milestone12RequestedStage::Refresh,
                _ => return Err(format!("Unsupported stage: {label}").into()),
            };
            if !self.stages.contains(&stage) {
                self.stages.push(stage);
            }
        }

        Ok(())
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let requested = RequestedStages::from_args()?;
    run_milestone12_report(&requested.stages)?;
    Ok(())
}
