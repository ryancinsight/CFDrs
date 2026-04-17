//! Milestone 12 blueprint-native candidate-space generation.
//!
//! ## Module hierarchy
//!
//! | Module | Responsibility |
//! |--------|----------------|
//! | [`builder`] | Primitive selective blueprint candidate construction |
//! | [`dimensions`] | Milestone 12 selective sweep dimensions |
//! | [`sweep`] | Milestone 12 selective sweep |

mod builder;
mod dimensions;
mod sweep;

use crate::domain::BlueprintCandidate;
pub use sweep::milestone12::CandidateParams;

/// Return lightweight parameter tuples (~100 bytes each) for lazy
/// materialization. This is the memory-efficient alternative to
/// [`build_milestone12_blueprint_candidate_space`].
#[must_use]
pub fn build_milestone12_candidate_params() -> Vec<CandidateParams> {
    sweep::milestone12::build_milestone12_candidate_params_only()
}

/// Build only the canonical primitive selective candidates needed by the
/// Milestone 12 report pipeline.
///
/// This avoids allocating and scanning unrelated topology families when the
/// report only consumes the selective-routing lineage:
/// `Option 1 base -> Option 2 venturi -> GA refinement`.
#[must_use]
pub(crate) fn build_milestone12_candidate_space() -> Vec<BlueprintCandidate> {
    sweep::milestone12::build_milestone12_candidate_space()
}

/// Build the Milestone 12 selective-routing candidate space directly as
/// blueprint-native candidates.
///
/// **Warning:** eagerly materializes ~500K candidates (~7.5 GB). Prefer
/// [`build_milestone12_candidate_params`] for memory-constrained pipelines.
pub fn build_milestone12_blueprint_candidate_space(
) -> Result<Vec<BlueprintCandidate>, crate::error::OptimError> {
    Ok(build_milestone12_candidate_space())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn root_family_label(candidate: &CandidateParams) -> &'static str {
        match candidate.request.split_kinds.first() {
            Some(cfd_schematics::SplitKind::NFurcation(2)) => "Bi",
            Some(cfd_schematics::SplitKind::NFurcation(4)) => "Quad",
            Some(cfd_schematics::SplitKind::NFurcation(5)) => "Penta",
            Some(cfd_schematics::SplitKind::NFurcation(3)) => "Tri",
            Some(cfd_schematics::SplitKind::NFurcation(_)) | None => {
                panic!("unexpected milestone12 root family")
            }
        }
    }

    #[test]
    fn milestone12_candidate_space_is_selective_only() {
        let candidates = build_milestone12_candidate_params();
        assert!(
            !candidates.is_empty(),
            "Milestone 12 candidate space should not be empty"
        );
        assert!(candidates.iter().all(|candidate| {
            candidate
                .request
                .split_kinds
                .first()
                .is_some_and(|split| matches!(split, cfd_schematics::SplitKind::NFurcation(2..=5)))
        }));
        assert!(candidates
            .iter()
            .any(|candidate| root_family_label(candidate) == "Bi"));
        assert!(candidates
            .iter()
            .any(|candidate| root_family_label(candidate) == "Quad"));
        assert!(candidates
            .iter()
            .any(|candidate| root_family_label(candidate) == "Penta"));
    }

    #[test]
    fn milestone12_candidate_space_contains_both_acoustic_and_venturi_modes() {
        let candidates = build_milestone12_candidate_params();
        let has_acoustic = candidates.iter().any(|candidate| !candidate.is_venturi());
        let has_venturi = candidates.iter().any(CandidateParams::is_venturi);
        assert!(
            has_acoustic,
            "Milestone 12 candidate space should include acoustic designs"
        );
        assert!(
            has_venturi,
            "Milestone 12 candidate space should include venturi designs"
        );
    }
}
