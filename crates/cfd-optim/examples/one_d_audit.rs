use cfd_optim::{
    audit_candidate_1d, build_candidate_space, DesignTopology, PrimitiveSplitSequence,
};

fn main() {
    let candidate = build_candidate_space()
        .into_iter()
        .find(|candidate| {
            matches!(
                candidate.topology,
                DesignTopology::PrimitiveSelectiveTree {
                    sequence: PrimitiveSplitSequence::TriTri,
                }
            )
        })
        .expect("primitive Tri->Tri candidate must exist");

    let audit = audit_candidate_1d(&candidate).expect("1D audit must succeed");
    println!(
        "{}",
        serde_json::to_string_pretty(&audit).expect("audit report must serialize")
    );
}
