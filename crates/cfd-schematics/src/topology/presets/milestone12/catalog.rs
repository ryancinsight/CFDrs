use super::super::sequence::MILESTONE12_SWEEP_SEQUENCES;
use super::Milestone12TopologyRequest;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct MirrorVariant {
    suffix: &'static str,
    mirror_x: bool,
    mirror_y: bool,
}

pub(super) const MILESTONE12_MIRROR_VARIANTS: [MirrorVariant; 4] = [
    MirrorVariant {
        suffix: "base",
        mirror_x: false,
        mirror_y: false,
    },
    MirrorVariant {
        suffix: "x",
        mirror_x: true,
        mirror_y: false,
    },
    MirrorVariant {
        suffix: "y",
        mirror_x: false,
        mirror_y: true,
    },
    MirrorVariant {
        suffix: "xy",
        mirror_x: true,
        mirror_y: true,
    },
];

/// Enumerate the canonical Milestone 12 split-tree scaffolds.
///
/// The catalog fixes only the split-stage lineage; callers remain responsible
/// for sweeping operating conditions and Milestone-12-specific geometric
/// parameters such as asymmetric width fractions, venturi throat geometry, and
/// terminal serpentine settings.
#[must_use]
pub fn enumerate_milestone12_topologies() -> Vec<Milestone12TopologyRequest> {
    let mut requests =
        Vec::with_capacity(MILESTONE12_SWEEP_SEQUENCES.len() * MILESTONE12_MIRROR_VARIANTS.len());
    for seq in &MILESTONE12_SWEEP_SEQUENCES {
        let split_kinds = seq.to_split_kinds();
        let family = seq.label().replace('\u{2192}', "");
        for variant in &MILESTONE12_MIRROR_VARIANTS {
            let mut request = Milestone12TopologyRequest::new(
                format!("pst-{}-{}", family.to_ascii_lowercase(), variant.suffix),
                format!("{family}-{}", variant.suffix.to_ascii_uppercase()),
                split_kinds.clone(),
                6.0e-3,
                1.0e-3,
                8.0e-3,
                8.0e-3,
            );
            request.mirror_x = variant.mirror_x;
            request.mirror_y = variant.mirror_y;
            requests.push(request);
        }
    }
    requests
}
