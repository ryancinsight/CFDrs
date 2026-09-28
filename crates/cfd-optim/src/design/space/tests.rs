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
    assert!(
        candidates
            .iter()
            .any(|candidate| root_family_label(candidate) == "Bi")
    );
    assert!(
        candidates
            .iter()
            .any(|candidate| root_family_label(candidate) == "Quad")
    );
    assert!(
        candidates
            .iter()
            .any(|candidate| root_family_label(candidate) == "Penta")
    );
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
