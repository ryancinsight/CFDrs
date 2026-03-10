use crate::BlueprintCandidate;
use cfd_schematics::domain::therapy_metadata::TherapyZone;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Milestone12Stage {
    Option1Base,
    Option2Derived,
    GaRefined,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Milestone12LineageKey {
    stage_sequence_label: String,
    treatment_channel_count: usize,
    terminal_serpentine_segments: usize,
}

#[must_use]
pub fn is_milestone12_lineage_topology(candidate: &BlueprintCandidate) -> bool {
    candidate
        .topology_spec()
        .is_ok_and(|spec| spec.first_stage_is_trifurcation() && spec.is_selective_routing())
}

#[must_use]
pub fn milestone12_lineage_key(candidate: &BlueprintCandidate) -> Option<Milestone12LineageKey> {
    let spec = candidate.topology_spec().ok()?;
    if !spec.first_stage_is_trifurcation() || !spec.is_selective_routing() {
        return None;
    }

    let terminal_segments = spec
        .split_stages
        .iter()
        .rev()
        .flat_map(|stage| stage.branches.iter())
        .find(|branch| branch.treatment_path)
        .and_then(|branch| branch.route.serpentine.as_ref())
        .map_or(0, |serpentine| serpentine.segments);

    Some(Milestone12LineageKey {
        stage_sequence_label: spec.stage_sequence_label(),
        treatment_channel_count: spec.treatment_channel_ids().len(),
        terminal_serpentine_segments: terminal_segments,
    })
}

impl Milestone12LineageKey {
    #[must_use]
    pub fn stage_sequence_label(&self) -> &str {
        &self.stage_sequence_label
    }
}

pub fn validate_milestone12_candidate(
    candidate: &BlueprintCandidate,
    stage: Milestone12Stage,
) -> Result<(), String> {
    let spec = candidate
        .topology_spec()
        .map_err(|error| error.to_string())?;

    if !is_milestone12_lineage_topology(candidate) {
        return Err(format!(
            "{stage:?} must use a tri-first selective-routing scaffold, got {}",
            spec.short_code()
        ));
    }
    if candidate.id.ends_with("-ACS") {
        return Err(format!(
            "{stage:?} candidate {} uses deprecated report-time acoustic clone suffix -ACS",
            candidate.id
        ));
    }

    match stage {
        Milestone12Stage::Option1Base => {
            if !spec.venturi_placements.is_empty() {
                return Err(format!(
                    "Option1Base candidate {} must be ultrasound-only",
                    candidate.id
                ));
            }
        }
        Milestone12Stage::Option2Derived | Milestone12Stage::GaRefined => {
            if spec.venturi_placements.is_empty() {
                return Err(format!(
                    "{stage:?} candidate {} must use venturi treatment",
                    candidate.id
                ));
            }
        }
    }

    candidate.blueprint().validate().map_err(|error| {
        format!(
            "{stage:?} candidate {} produced invalid blueprint: {error}",
            candidate.id
        )
    })?;

    let treatment_lane_count = candidate
        .blueprint()
        .channels
        .iter()
        .filter(|channel| channel.therapy_zone == Some(TherapyZone::CancerTarget))
        .count();
    if treatment_lane_count < spec.treatment_channel_ids().len() {
        return Err(format!(
            "{stage:?} candidate {} collapsed treatment-window lanes: expected at least {}, got {}",
            candidate.id,
            spec.treatment_channel_ids().len(),
            treatment_lane_count
        ));
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::milestone12_lineage_key;
    use crate::domain::fixtures::{
        canonical_option1_candidate, canonical_option2_candidate, operating_point,
    };

    #[test]
    fn milestone12_lineage_key_matches_option1_and_option2_of_same_scaffold() {
        let op = operating_point(2.4e-6, 32_000.0, 0.12);
        let option1 = canonical_option1_candidate("option1", op.clone());
        let option2 = canonical_option2_candidate("option2", op);

        assert_eq!(
            milestone12_lineage_key(&option1),
            milestone12_lineage_key(&option2)
        );
    }
}
