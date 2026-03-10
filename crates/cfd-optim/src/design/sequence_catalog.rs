#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PrimitiveSplitSequence {
    Tri,
    TriBi,
    TriTri,
    TriBiBi,
    TriBiTri,
    TriTriBi,
    TriTriTri,
}

impl PrimitiveSplitSequence {
    #[must_use]
    pub const fn levels(self) -> u8 {
        match self {
            Self::Tri => 1,
            Self::TriBi | Self::TriTri => 2,
            Self::TriBiBi | Self::TriBiTri | Self::TriTriBi | Self::TriTriTri => 3,
        }
    }

    #[must_use]
    pub const fn split_types(self) -> u8 {
        match self {
            Self::Tri => 0b1,
            Self::TriBi => 0b01,
            Self::TriTri => 0b11,
            Self::TriBiBi => 0b001,
            Self::TriBiTri => 0b101,
            Self::TriTriBi => 0b011,
            Self::TriTriTri => 0b111,
        }
    }

    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::Tri => "Tri",
            Self::TriBi => "Tri→Bi",
            Self::TriTri => "Tri→Tri",
            Self::TriBiBi => "Tri→Bi→Bi",
            Self::TriBiTri => "Tri→Bi→Tri",
            Self::TriTriBi => "Tri→Tri→Bi",
            Self::TriTriTri => "Tri→Tri→Tri",
        }
    }
}

/// Restricted tri-first primitive selective lineage used by Milestone 12 flows.
pub(crate) const TRI_FIRST_PRIMITIVE_SELECTIVE_SEQUENCES: [PrimitiveSplitSequence; 7] = [
    PrimitiveSplitSequence::Tri,
    PrimitiveSplitSequence::TriBi,
    PrimitiveSplitSequence::TriTri,
    PrimitiveSplitSequence::TriBiBi,
    PrimitiveSplitSequence::TriBiTri,
    PrimitiveSplitSequence::TriTriBi,
    PrimitiveSplitSequence::TriTriTri,
];

/// Returns `(has_intermediate_tri, has_any_tri, has_any_bi)` for one PST sequence.
#[must_use]
pub(crate) fn primitive_sequence_metadata(sequence: PrimitiveSplitSequence) -> (bool, bool, bool) {
    let bits = sequence.split_types();
    let n_tri = bits.count_ones() as u8;
    let has_any_tri = n_tri > 0;
    let has_any_bi = n_tri < sequence.levels();
    let has_intermediate_tri = n_tri >= 2;
    (has_intermediate_tri, has_any_tri, has_any_bi)
}
